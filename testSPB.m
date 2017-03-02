function models = testSPB(models, varargin)
% TESTSPB Function for comparing performance of various
%   symmetry/ridge/medial axis detection algorithms.

% Default testing options ------------------------------------------------
opts = {'dataset',     'BMAX500',...
        'set',         'val',...   % 'val' or 'test'
        'visualize',   false,...
        'parpoolSize', feature('numcores'),... % set to 0 to run serially
        'nThresh',     30,...      % #thresholds used for computing p-r
        'maxDist',     0.01        % controls max distance of an accurately 
       };                          % detected point from groundtruth.
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths();
if ischar(opts.set) && strcmp(opts.set, 'val')
    opts.imPath = fullfile(paths.bsds500im,'val');
    opts.gtPath = fullfile(paths.symmax500,'val');
    imageList = dir(fullfile(opts.imPath, '*jpg'));
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    opts.imPath = fullfile(paths.bsds500im,'test');
    opts.gtPath = fullfile(paths.symmax500,'test');
    imageList = dir(fullfile(opts.imPath, '*jpg'));
elseif isstruct(opts.set)
    disp('Data provided in struct form')
    imageList = opts.set;
    if strcmp(opts.dataset, 'BMAX500')
        if numel(imageList) == 100, opts.set = 'val'; else opts.set = 'test'; end
    end
else
    error('set can be ''val'', ''test'', or a struct containing test data')
end

% Evaluate models on test images ------------------------------------------
if ~iscell(models), models = {models}; end
for m=1:numel(models)
    models{m} = evaluateModel(models{m},imageList,opts,paths);
end

% Compute dataset-wide stats
for m=1:numel(models)
    [models{m}.stats.odsP,  models{m}.stats.odsR, ...
     models{m}.stats.odsF,  models{m}.stats.odsT, ...
     models{m}.stats.oisP,  models{m}.stats.oisR, ...
     models{m}.stats.oisF,  models{m}.stats.AP] = ...
        computeDatasetStats(models{m}.stats, opts);
    % Create field with dataset-specific stats
    models{m}.(opts.dataset).(opts.set).stats = models{m}.stats;
    models{m}.(opts.dataset).(opts.set).opts = opts;
    models{m} = rmfield(models{m},'stats');
    % And store results
    modelPath = fullfile(paths.spbmil.models, models{m}.name);
    model = models{m}; save(modelPath, 'model')
end

% Plot precision-recall curves --------------------------------------------
plotPrecisionRecall(models,opts.dataset,opts.set) 

% -------------------------------------------------------------------------
function model = evaluateModel(model,imageList,opts,paths)
% -------------------------------------------------------------------------
% Load model
opts.thresh = linspace(1/(opts.nThresh+1),1-1/(opts.nThresh+1),opts.nThresh)';
switch lower(model)
    case 'levinstein'
        opts.thresh = 0.5; opts.nThresh = 1;
        model = loadLevinsteinModel(model,paths);
    case 'amat'
        opts.thresh = 0.5; opts.nThresh = 1;
        model = struct('name',model);
    case 'lindeberg'        
        model = struct('name',model);
    otherwise % load MIL or CNN model
        model = loadModelFromMatFile(model,paths);
end

% Initialize stats
opts.nImages = numel(imageList);
cntP = zeros(opts.nImages, opts.nThresh);
cntR = zeros(opts.nImages, opts.nThresh);
sumP = zeros(opts.nImages, opts.nThresh);
sumR = zeros(opts.nImages, opts.nThresh);
scores = zeros(opts.nImages, 4); % optimal P,R,F,T for each image

modelName = lower(model.name);
ticStart = tic;
% parfor (i=1:opts.nImages, opts.parpoolSize)
for i=1:opts.nImages % keep that just for debugging
    if isfield(imageList(i), 'isdir')
        % Load image and groundtruth data from disk
        [~,iid,~] = fileparts(imageList(i).name);
        gt  = load(fullfile(opts.gtPath,['gt_' iid '.mat' ])); gt = gt.gt;
        img = imread(fullfile(opts.imPath,imageList(i).name));
    else % Read image and groundtruth from struct
        img = imageList(i).img;
        gt  = imageList(i).pts;
        iid = imageList(i).iid;
    end
    
    switch modelName
        case 'levinstein'
            spb = evaluateLevinshtein(model, img);
        case 'lindeberg'
            spb = evaluateLindeberg(img);
        case 'amat'
            spb = evaluateAMAT(img);
        case 'deepskel'
            spb = evaluateDeepSkel(model,img);
        otherwise % MIL
            spb = evaluateModelMIL(model,img,iid,opts);
    end
    [cntP(i,:), sumP(i,:), cntR(i,:), sumR(i,:),scores(i,:)] = ...
        computeImageStats(spb,gt,opts);
    
    msg = sprintf('Testing medial point detection on %s %s set. ', opts.dataset, opts.set);
    progress(msg,i,opts.nImages,ticStart,-1);
end

% Store stats in model struct
model.stats.cntP = cntP;
model.stats.sumP = sumP;
model.stats.cntR = cntR;
model.stats.sumR = sumR;
model.stats.scores = scores;

% -------------------------------------------------------------------------
function spb = evaluateAMAT(img)
% -------------------------------------------------------------------------
[H,W,~] = size(img);
img = imresize(img,0.5,'bilinear');
img = L0Smoothing(img);
mat = amat(img);
spb = any(mat.axis,3); 
spb = imresize(spb,[H,W],'nearest');

% -------------------------------------------------------------------------
function spb = evaluateDeepSkel(model,img)
% -------------------------------------------------------------------------
net = model.net;
img = bsxfun(@minus, single(img), reshape(net.meta.averageImage,1,1,[]));
net.eval({'input',img});
lout= net.getLayer(net.getLayerIndex('softmax'));
spb = net.vars(lout.outputIndexes);
spb = 1-spb.value(:,:,1);

% -------------------------------------------------------------------------
function spb = evaluateModelMIL(model,img,iid,opts)
% -------------------------------------------------------------------------
paths = setPaths();
if strcmp(model.opts.featureSet, 'spectral')
    try
        spectralFeat = load(fullfile(paths.spectral,'spectral_50',...
            opts.set,['spectral_' iid '.mat']));
        spectralFeat = single(spectralFeat.spectral);
    catch
        warning('Was not able to load spectral feature')
    end
else
    spectralFeat = [];
end
spb = spbMIL(img, 'featureSet',model.opts.featureSet, 'fineScale',true,...
             'w',model.w, 'spectralFeature',spectralFeat);
spb = spb.thin;

% -------------------------------------------------------------------------
function ridges = evaluateLindeberg(im)
% -------------------------------------------------------------------------
[~,~,~,contours] = PS00___primal_sketch(double(rgb2gray(im)));
ridges  = PSzz_unzip_contour(contours{1});

% -------------------------------------------------------------------------
function axes = evaluateLevinshtein(model,img)
% -------------------------------------------------------------------------
[parts, part_scales, part_score, image_data] = ...
    DetectSymmetricParts_SingleImage(fullfile(testImDir,testImages(i).name),...
    model.classifier_params,[], [], [], [], true);
%     figure; DisplayParts(im, parts, part_scales, 1:size(parts,3));
[part_labels, part_affinity, part_adjacency] = ...
    GroupSymmetricParts_SingleImage(img, parts, model.part_classifier_params, image_data);
%     figure; DrawPartClustersEllipses(im, parts, part_labels);
selected_parts_indicator = SelectSymmetricPart_SingleImage(parts, ...
    part_labels, part_score, part_affinity, part_adjacency,0.1);
%     figure; DrawPartClustersEllipses(im, parts, part_labels, selected_parts_indicator);
axes = levsym(img,parts, part_labels, selected_parts_indicator);

% -------------------------------------------------------------------------
function [cntP,sumP,cntR,sumR,scores] = computeImageStats(pb,gt,opts)
% -------------------------------------------------------------------------
% For Levinstein's method we do not need to threshold
if islogical(pb), 
    thresh = 0.5; pb = double(pb);
else
    thresh = opts.thresh;
end

% Initialize
cntP = zeros(size(thresh));
sumP = zeros(size(thresh));
cntR = zeros(size(thresh));
sumR = zeros(size(thresh));

% Compute numerator (cntX) and denominator (sumX) for precision and recall.
% Note that because there are multiple groundtruth maps to compare with a
% single machine-generated response, the number of true positives, false
% positives etc. used for computing precision is different than the ones
% that are used for computing recall.
for t = 1:numel(thresh),
    % Threshold probability map and thin to 1-pixel width.
    bmap = (pb >= thresh(t));
%     bmap = bwmorph(bwmorph(bmap,'thin',Inf), 'clean');
    bmap = bwmorph(bmap,'thin',Inf);
    
    % Compute matches between symmetry map and all groundtruth maps
    accP = 0;
    for s=1:size(gt,3)
        [match1,match2] = correspondPixels(double(bmap),double(gt(:,:,s)),opts.maxDist);
        if opts.visualize, plotMatch(1,bmap,gt(:,:,s),match1,match2); end
        % accumulate machine matches
        accP = accP | match1;
        cntR(t) = cntR(t) + nnz(match2>0); % tp (for recall)
    end
    cntP(t) = nnz(accP); % tp (for precision)
    sumP(t) = nnz(bmap); % tp + fp (for precision)
    sumR(t) = nnz(gt);   % tp + fn (for recall)
end

% Compute precision (P), recall (R) and f-measure (F). 
P = cntP ./ max(eps, sumP);
R = cntR ./ max(eps, sumR);

% Use linear interpolation to find best P,R,F combination.
% scores contains the image-specific optimal P,R,F, after computing optimal
% thresholds using linear interpolation.
[bestP,bestR,bestF,bestT] = findBestPRF(P,R,thresh);
scores = [bestP,bestR,bestF,bestT];

% -------------------------------------------------------------------------
function [odsP, odsR, odsF, odsT, oisP, oisR, oisF, AP] = computeDatasetStats(stats,opts)
% Two standard F-based performance metrics are computed:
% i)  ODS: F-measure for a dataset-wide specified optimal threshold.
% ii) OIS: F-measure for an image-specific optimal threshold.
% iii)AP:  Average precision - equivalent to AUC (area under curve).

P = sum(stats.cntP,1) ./ max(eps, sum(stats.sumP,1));
R = sum(stats.cntR,1) ./ max(eps, sum(stats.sumR,1));

if length(P) > 1 % soft probability maps
    % ODS scores (scalars)
    [odsP,odsR,odsF,odsT] = findBestPRF(P,R,opts.thresh);

    % OIS scores (scalars)
    Pi = stats.cntP ./ max(eps, stats.sumP);
    Ri = stats.cntR ./ max(eps, stats.sumR);
    [~,indMaxF] = max(fmeasure(Pi,Ri),[],2);
    indMaxF = sub2ind(size(Pi), (1:size(Pi,1))', indMaxF);
    oisP = sum(stats.cntP(indMaxF)) ./ max(eps, sum(stats.sumP(indMaxF)));
    oisR = sum(stats.cntR(indMaxF)) ./ max(eps, sum(stats.sumR(indMaxF)));
    oisF = fmeasure(oisP,oisR);

    % AP score (scalar)
    AP = interp1(R,P, 0:0.01:1); 
    AP = sum(AP(~isnan(AP)))/100;
else    % binary symmetry maps
    odsP = P; odsR = R; odsF = fmeasure(P,R); odsT = 0.5;
    oisP = odsP; oisR = odsR; oisF = odsF; AP = odsP;
end

% -------------------------------------------------------------------------
function F = fmeasure(P,R), F = 2 .* P .* R ./ max(eps, P+R);

% -------------------------------------------------------------------------
function [bestP,bestR,bestF,bestT] = findBestPRF(P,R,T)
% -------------------------------------------------------------------------
if numel(T) == 1
    bestT = T; bestR = R; bestP = P; bestF = fmeasure(P,R); return
end

bestF = -1;
a = linspace(0,1,100); b = 1-a;
for t = 2:numel(T)
    Rt = a.*R(t) + b.*R(t-1);
    Pt = a.*P(t) + b.*P(t-1);
    Tt = a.*T(t) + b.*T(t-1);
    Ft = fmeasure(Pt,Rt);
    [f,indMax] = max(Ft); 
    if f > bestF
        bestF = f; bestT = Tt(indMax);
        bestP = Pt(indMax); bestR = Rt(indMax); 
    end
end

% -------------------------------------------------------------------------
function model = loadLevinsteinModel(model,paths)
% -------------------------------------------------------------------------
model = struct('name',model);
tmp1 = load('classifier_params_svm.mat');
tmp2 = load('part_classifier_params.mat');
model.classifier_params = tmp1.classifier_params;
model.part_classifier_params = tmp2.part_classifier_params;

% -------------------------------------------------------------------------
function model = loadModelFromMatFile(model,paths)
% -------------------------------------------------------------------------
if exist(fullfile(paths.spbmil.models, model),'file') % MIL model
    tmp = load(fullfile(paths.spbmil.models, model));
elseif exist(model,'file')
    tmp = load(model);
end
if isfield(tmp,'model')     % MIL model
    model = tmp.model;
elseif isfield(tmp,'net')   % DagNN model
    model = struct('trainStats',tmp.stats, 'name','deepskel',...
                   'net',dagnn.DagNN.loadobj(tmp.net));
end

