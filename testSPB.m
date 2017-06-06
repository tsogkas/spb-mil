function models = testSPB(models, varargin)
% TESTSPB Function for comparing performance of various
%   symmetry/ridge/medial axis detection algorithms.

% Default testing options ------------------------------------------------
opts = {'dataset',     'BMAX500',...
        'testSet',     'val',...   % 'val' or 'test'
        'visualize',   false,...
        'parpoolSize', feature('numcores'),... % set to 0 to run serially
        'nThresh',     30,...      % #thresholds used for computing p-r
        'maxDist',     0.01,...    % controls max distance of an accurately detected point from groundtruth.
        'amatws',      1e-4,...    % controls amat coarsenesskai to live
        'ucmthresh',   0.5
       };                          
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths();
switch opts.dataset
    case 'BMAX500'
        disp('Data provided in struct form')
        imageList = opts.testSet;
        if numel(imageList) == 100, opts.testSet = 'val';
        else opts.testSet = 'test';
        end
    case 'SYMMAX500'    % TODO: add support for SYMMAX300
        opts.imPath = fullfile(paths.bsds500, 'images', opts.testSet);
        opts.gtPath = fullfile(paths.symmax500, 'groundtruth', opts.testSet);
        imageList   = dir(fullfile(opts.imPath, '*jpg'));
    case 'SK506'
        opts.testSet = 'test';
        opts.imPath = fullfile(paths.sk506, 'images', 'test');
        opts.gtPath = fullfile(paths.sk506, 'groundTruth', 'test');
        imageList   = dir(fullfile(opts.imPath, '*jpg'));
    case 'WHSYMMAX'
        opts.testSet = 'test';
        opts.imPath = fullfile(paths.whsymmax, 'imgs', 'test');
        opts.gtPath = fullfile(paths.whsymmax, 'gt', 'test');
        imageList   = dir(fullfile(opts.imPath, '*jpg'));        
    otherwise, error('Dataset is not supported')
end

% Evaluate models on test images ------------------------------------------
opts.thresh = linspace(1/(opts.nThresh+1),1-1/(opts.nThresh+1),opts.nThresh)';
if ~iscell(models), models = {models}; end
for m=1:numel(models)
    if ismember(models{m}, {'human','levinshtein','amat','ucm'})
        opts.thresh = 0.5; opts.nThresh = 1;
    end
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
    models{m}.(opts.dataset).(opts.testSet).stats = models{m}.stats;
    models{m}.(opts.dataset).(opts.testSet).opts = opts;
    models{m} = rmfield(models{m},'stats');
    % And store results
    if ~isdir(paths.spbmil.models), mkdir(paths.spbmil.models); end
    modelPath = fullfile(paths.spbmil.models, [models{m}.name '.mat']);
    model = models{m}; save(modelPath, 'model')
end

% Plot precision-recall curves --------------------------------------------
plotPrecisionRecall(models,opts.dataset,opts.testSet) 

% -------------------------------------------------------------------------
function model = evaluateModel(model,imageList,opts,paths)
% -------------------------------------------------------------------------
% Load model
if ischar(model)
    [~,modelName] = fileparts(model);
    if exist(fullfile(paths.spbmil.models, [modelName '.mat']), 'file')
        model = load(fullfile(paths.spbmil.models, [modelName '.mat'])); 
        model = model.model;
    end
else
    switch lower(model)
        case {'human','amat','lindeberg','ucm'}
            model = struct('name',model);
        case 'levinshtein'
            model = loadLevinshteinModel(model,paths);
        otherwise % load MIL or CNN model
            if strfind(model, 'caffemodel')
                model = loadDeepSkelModel(model,paths);
            else
                model = loadModelFromMatFile(model,paths);
            end
    end
end

% Initialize stats
if isfield(model, opts.dataset) && isfield(model.(opts.dataset), opts.testSet)
    model.stats = model.(opts.dataset).(opts.testSet).stats;
end
opts.nImages = numel(imageList);
cntP = zeros(opts.nImages, opts.nThresh);
cntR = zeros(opts.nImages, opts.nThresh);
sumP = zeros(opts.nImages, opts.nThresh);
sumR = zeros(opts.nImages, opts.nThresh);
scores = zeros(opts.nImages, 4); % optimal P,R,F,T for each image

modelName = lower(model.name);
ticStart = tic;
% parfor (i=1:opts.nImages, opts.parpoolSize)
for i=47:opts.nImages % keep that just for debugging
    [img,gt,iid,fgmask] = loadImageGroundtruth(imageList,i,opts);
    switch modelName
        case 'human'
            spb = gt;
        case 'levinshtein'
            spb = evaluateLevinshtein(model, img);
        case 'lindeberg'
            spb = evaluateLindeberg(img);
        case 'amat'
            if strcmp(opts.dataset, 'BMAX500')
                spb = evaluateAMAT(imresize(L0Smoothing(img),0.5),opts.amatws);
                spb = imresize(spb,[size(gt,1),size(gt,2)],'nearest');
            else
                spb = evaluateAMAT(L0Smoothing(img),opts.amatws);
            end
        case {'horse_it14k','sk506_it14k'}
            spb = evaluateDeepSkel(model,img);
        case 'ucm'
            spb = evaluateUCM(iid,opts);
        otherwise % MIL
            spb = evaluateModelMIL(model,img,iid,opts);
    end

    % Optionally filter spb responses
    if ~isempty(fgmask), spb(~fgmask) = 0; end
    
    [cntP(i,:), sumP(i,:), cntR(i,:), sumR(i,:),scores(i,:)] = ...
        computeImageStats(spb,gt,opts);
    
    msg = sprintf('Testing medial point detection on %s %s set. ', opts.dataset, opts.testSet);
    progress(msg,i,opts.nImages,ticStart,-1);
end

% Store stats in model struct
model.stats.cntP = cntP;
model.stats.sumP = sumP;
model.stats.cntR = cntR;
model.stats.sumR = sumR;
model.stats.scores = scores;

% -------------------------------------------------------------------------
function [cntP,sumP,cntR,sumR,scores] = computeImageStats(pb,gt,opts)
% -------------------------------------------------------------------------
if size(pb,3) > 1
    [cntP,sumP,cntR,sumR,scores] = computeImageStatsHuman(gt,opts);
    return
end

% For levinshtein's method we do not need to threshold
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
    bmap = bwmorph(bmap,'thin',Inf);
    
    % Compute matches between symmetry map and all groundtruth maps
    accP = 0;
    for s=1:size(gt,3)
        gt(:,:,s) = bwmorph(gt(:,:,s), 'thin',inf);
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
function [cntP,sumP,cntR,sumR,scores] = computeImageStatsHuman(gt,opts)
% -------------------------------------------------------------------------
% Initialize
thresh  = 0.5;
numSegs = size(gt,3);
cntP = zeros(numSegs,1);
sumP = zeros(numSegs,1); 
cntR = zeros(numSegs,1); 
sumR = zeros(numSegs,1);

% Compute numerator (cntX) and denominator (sumX) for precision and recall.
% To compute human performance we compare one of the annotations at a time
% against all remaining groundtruth maps.
% First make sure all annotations are 1-pixel wide
for t=1:numSegs
    gt(:,:,t) = bwmorph(gt(:,:,t),'thin',inf);
end

% Now compare each gt annotation to all the other gt annotations
indSegs = 1:numSegs;
for t = 1:numSegs
    gtc = gt(:,:,t);
    gtr = gt(:,:,setdiff(indSegs,t));        
    accP = 0;
    for s=1:size(gtr,3)
        [match1,match2] = correspondPixels(double(gtc),double(gtr(:,:,s)),opts.maxDist);
        if opts.visualize, plotMatch(1,gtc,gtr(:,:,s),match1,match2); end
        % accumulate machine matches
        accP = accP | match1;
        cntR(t) = cntR(t) + nnz(match2>0); % tp (for recall)
    end
    cntP(t) = nnz(accP); % tp (for precision)
    sumP(t) = nnz(gtc);  % tp + fp (for precision)
    sumR(t) = nnz(gtr);  % tp + fn (for recall)
end

% Compute precision (P), recall (R) and f-measure (F) for all subjects.
cntP = sum(cntP);
sumP = sum(sumP);
cntR = sum(cntR);
sumR = sum(sumR);
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
function [img,gt,iid,fgmask] = loadImageGroundtruth(imageList,i,opts)
% -------------------------------------------------------------------------
fgmask = []; % optional foreground mask
switch opts.dataset
    case 'BMAX500'  % Read image and groundtruth from struct
        img = imageList(i).img;
        gt  = imageList(i).pts;
        iid = imageList(i).iid;
    case {'SYMMAX500','SYMMAX300'}
        [~,iid] = fileparts(imageList(i).name);
        gt = load(fullfile(opts.gtPath,['gt_' iid '.mat' ]));
        img = imread(fullfile(opts.imPath,imageList(i).name));
    case 'SK506'
        [~,iid] = fileparts(imageList(i).name);
        gt = load(fullfile(opts.gtPath,[iid '.mat' ]));
        img = imread(fullfile(opts.imPath,imageList(i).name));
        % reconstruct foreground mask from edges
        % make sure we only have one connected component
        fgmask = imfill(gt.edge, find(gt.symmetry));
        gt = gt.symmetry;
    case 'WHSYMMAX'
        [~,iid] = fileparts(imageList(i).name);
        gt = load(fullfile(opts.gtPath,[iid '.mat' ]));
        img = imread(fullfile(opts.imPath,imageList(i).name));
        fgmask = gt.segmentation;
        gt = gt.sym > 0;
    otherwise, error('Dataset not supported')
end

% -------------------------------------------------------------------------
function model = loadDeepSkelModel(model,paths)
% -------------------------------------------------------------------------
% MAKE SURE YOU ARE USING THE CORRECT DEPLOY PROTOTXT!
deployFile = '/home/tsogkas/code/DeepSkeleton-master/examples/DeepSkeleton/deploy.prototxt';
[~,modelName] = fileparts(model);
model = struct('name',modelName,'net', caffe.Net(deployFile, model, 'test'));
caffe.set_mode_gpu();
caffe.set_device(0);

% -------------------------------------------------------------------------
function model = loadLevinshteinModel(model,paths)
% -------------------------------------------------------------------------
model = struct('name',model);
tmp1 = load('classifier_params_svm.mat');
tmp2 = load('part_classifier_params.mat');
model.classifier_params = tmp1.classifier_params;
model.part_classifier_params = tmp2.part_classifier_params;

% -------------------------------------------------------------------------
function model = loadModelFromMatFile(model,paths)
% -------------------------------------------------------------------------
% TODO: update this
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

% -------------------------------------------------------------------------
function spb = evaluateAMAT(img,ws)
% -------------------------------------------------------------------------
mat = AMAT(img,'ws',ws).group.simplify;
spb = any(mat.axis,3); 

% -------------------------------------------------------------------------
function spb = evaluateDeepSkel(model,img)
% -------------------------------------------------------------------------
imageMean = [104.00698793,116.66876762,122.67891434]; % BGR
IMAGE_DIM = 500;

% Convert an image returned by Matlab's imread to im_data in caffe's data
% format: W x H x C with BGR channels
img0 = img;
img = img(:, :, [3, 2, 1]);  % permute channels from RGB to BGR
img = permute(img, [2, 1, 3]);  % flip width and height
img = single(img);  % convert from uint8 to single
img = bsxfun(@minus, img, reshape(imageMean, 1,1,[]));
img = imresize(img, [IMAGE_DIM IMAGE_DIM], 'bilinear');  % resize im_data
% Collect symmetry responses and revert RGB
spb = model.net.forward({img});
spb = (1-spb{end}(:,:,1)).'; % symmetry @ all scales = 1-background probability
spb = imresize(spb, [size(img0,1), size(img0,2)], 'bilinear');
% Apply NMS as implemented in Dollar's "Fast Edge Detection Using 
% "Structured Forests".
% [Ox,Oy]=gradient2(convTri(spb,4));
% [Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
% O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
% spb = edgesNmsMex(spb,O,1,5,1.01);
spb = edgeNms(spb,4,1);

% -------------------------------------------------------------------------
function spb = evaluateModelMIL(model,img,iid,opts)
% -------------------------------------------------------------------------
paths = setPaths();
if strcmp(model.opts.featureSet, 'spectral')
    try
        spectralFeat = load(fullfile(paths.spectral,'spectral_50',...
            opts.testSet,['spectral_' iid '.mat']));
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
function ridges = evaluateLindeberg(img)
% -------------------------------------------------------------------------
[~,~,~,contours] = PS00___primal_sketch(double(rgb2gray(img)));
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
function spb = evaluateUCM(iid,opts)
% ------------------------------------------------------------------------
%ucm2 = load(fullfile(paths.bsds500, 'ucm',opts.testSet)); ucm2 = ucm2;
temp = load(fullfile('/home/tsogkas/datasets/BSR/BSDS500/ucm2/',opts.testSet,[iid,'.mat'])); 
ucm2 = temp.ucm2;
seg  = bwlabel(ucm2 <= opts.ucmthresh);
seg  = seg(2:2:end, 2:2:end);
spb  = false(size(seg));
segments = unique(seg);
for i=1:numel(segments)
    skel = skeleton(seg == segments(i));
    spb = spb | (skel >= 50);
end


