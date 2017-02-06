function testSPB(models, varargin)
% TESTSPB Function for comparing performance of various
%   symmetry/ridge/medial axis detection algorithms.
% 


% Default training options ------------------------------------------------
opts = { 'set',       'val',...   % 'val' or 'test'
         'visualize', false,...
         'nThresh',   30,...      % #thresholds used for computing p-r
         'maxDist',   0.01        % controls max distance of an accurately 
       };                         % detected point from groundtruth.
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths();
if ischar(opts.set) && strcmp(opts.set, 'val')
    imageList = dir(fullfile(paths.bsds500imVal, '*jpg'));
    imPath    = paths.bsds500imVal;
    gtPath    = paths.symmax500gtVal;
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    imageList = dir(fullfile(paths.bsds500imTest, '*jpg'));
    imPath    = paths.bsds500imTest;
    gtPath    = paths.symmax500gtTest;
elseif isstruct(opts.set)
    disp('Data provided in struct form')
    imageList = opts.set;
else
    error('set can be ''val'', ''test'', or a struct containing test data')
end

% Load models and initialize stats ----------------------------------------
for m=1:numel(models)
    switch models{m}
        case 'levinstein'
            models{m} = struct('name',models{m});
            tmp1 = load('classifier_params_svm.mat');
            tmp2 = load('part_classifier_params.mat');
            models{m}.classifier_params = tmp1.classifier_params;
            models{m}.part_classifier_params = tmp2.part_classifier_params;
            clear tmp1 tmp2
        case 'lindeberg'
            models{m} = struct('name',models{m});
        otherwise % MIL model
            tmp = load(fullfile(paths.modelsMIL, models{m}.name));
            models{m} = tmp.spbModel;
    end
    models{m}.stats.cntR = zeros(opts.nThresh, opts.nImages);
    models{m}.stats.sumR = zeros(opts.nThresh, opts.nImages);
    models{m}.stats.cntP = zeros(opts.nThresh, opts.nImages);
    models{m}.stats.sumP = zeros(opts.nThresh, opts.nImages);
    models{m}.stats.scores = zeros(opts.nImages, 4); % optimal P,R,F,T for each image
end

% Evaluate models on test images ------------------------------------------
opts.nImages = numel(imageList);
opts.thresh  = linspace(1/(opts.nThresh+1),1-1/(opts.nThresh+1),opts.nThresh)';
for i=1:opts.nImages
    fprintf('Testing on image %d/%d from BSDS500 %s set\n', i, opts.nImages, opts.set);
    if isfield(imageList(i), 'isdir')
        % Load image and groundtruth data from disk
        [~,iid,~] = fileparts(imageList(i).name);
        gt  = load(fullfile(gtPath,['gt_' iid '.mat' ])); gt = gt.gt;
        img = im2double(imread(fullfile(imPath,imageList(i).name)));
    else % Read image and groundtruth from struct
        img = imageList(i).img;
        gt  = imageList(i).seg;
    end
    
    % Compute features and evaluate all models
    features = computeHistogramFeatures(img);
    for m=1:numel(models)
        switch models{m}.name
            case 'levinstein'
                spb = evaluateLevinshtein(models{m}, img);
            case 'lindeberg'
                spb  = evaluateLindeberg(img);
            otherwise % MIL 
                spb = evaluateModelMIL(models{m},img,features,iid,opts);
        end
        [models{m}.stats.cntP(i,:), models{m}.stats.sumP(i,:),...
         models{m}.stats.cntR(i,:), models{m}.stats.sumR(i,:),...
         models{m}.stats.scores(i,:)] = computeImageStats(spb,gt,opts);
    end
end

% Compute dataset-wide stats
for m=1:numel(models)
    [models{m}.stats.oidP,  models{m}.stats.oidR, ...
     models{m}.stats.oidF,  models{m}.stats.oidT, ...
     models{m}.stats.oisP,  models{m}.stats.oisR, ...
     models{m}.stats.oisF,  models{m}.stats.AP] = ...
        computeDatasetStats(models{m}.stats, opts);
end

% Plot precision-recall curves --------------------------------------------
plotPrecisionRecall(0)  % compare methods
plotPrecisionRecall(1)  % compare grouping


% -------------------------------------------------------------------------
function spb = evaluateModelMIL(model,img,histFeatures,iid,opts)
% -------------------------------------------------------------------------
paths = setPaths();
if strcmp(model.featureSet, 'spectral')
    try
        spectralFeat = load(fullfile(paths.spectral,'spectral_50',opts.set,['spectral_' iid '.mat']));
        spectralFeat = single(spectralFeat.spectral);
    catch
        warning('Was not able to load spectral feature')
    end
end
spb = spbMIL(img, 'featureSet',model.featureSet, 'fineScale',true,...
    'w',model.w, 'histFeatures', histFeatures, 'spectralFeat',spectralFeat);
spb = spb.thin;

% -------------------------------------------------------------------------
function ridges = evaluateLindeberg(im)
% -------------------------------------------------------------------------
[~,~,~,contours] = PS00___primal_sketch(double(rgb2gray(im)));
ridges  = PSzz_unzip_contour(contours{1});
%     figure; imshow(ridgesL,[]); title(['Lindeberg iid = ' num2str(iid)]);

% -------------------------------------------------------------------------
function axes = evaluateLevinshtein(model,im)
% -------------------------------------------------------------------------
[parts, part_scales, part_score, image_data] = ...
    DetectSymmetricParts_SingleImage(fullfile(testImDir,testImages(i).name),...
    model.classifier_params,[], [], [], [], true);
%     figure; DisplayParts(im, parts, part_scales, 1:size(parts,3));
[part_labels, part_affinity, part_adjacency] = ...
    GroupSymmetricParts_SingleImage(im, parts, model.part_classifier_params, image_data);
%     figure; DrawPartClustersEllipses(im, parts, part_labels);
selected_parts_indicator = SelectSymmetricPart_SingleImage(parts, ...
    part_labels, part_score, part_affinity, part_adjacency,0.1);
%     figure; DrawPartClustersEllipses(im, parts, part_labels, selected_parts_indicator);
axes = levsym(im,parts, part_labels, selected_parts_indicator);

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
    bmap = bwmorph(bwmorph(bmap,'thin',Inf), 'clean');
    
    % Compute matches between symmetry map and all groundtruth maps
    accP = 0;
    for s=1:size(gt,3)
        [match1,match2] = correspondPixels(double(bmap),double(gt),opts.maxDist);
        if opts.visualize
            plotMatch(1,bmap,gt,match1,match2);
        end
        % accumulate machine matches
        accP = accP | match1;
        cntR(t) = nnz(match2>0); % tp (for recall)
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

% ODS scores
P = sum(stats.cntP,1) ./ max(eps, sum(stats.sumP,1));
R = sum(stats.cntR,1) ./ max(eps, sum(stats.sumR,1));
F = fmeasure(P,R);
[odsP,odsR,odsF,odsT] = findBestPRF(P,R,opts.thresh);

% OIS scores
[~,indMaxF] = max(F,[],2);
oisP = sum(stats.cntP(:,indMaxF)) ./ max(eps, sum(stats.sumP(:,indMaxF)));
oisR = sum(stats.cntR(:,indMaxF)) ./ max(eps, sum(stats.sumR(:,indMaxF)));
oisF = fmeasure(oisP,oisR);

% AP score
AP = interp1(R,P, 0:0.01:1); 
AP = sum(AP(~isnan(AP)))/100;


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



% ridgesCGrouped50  = fpg(ridgesC.thin,50);  % grouping 50 strongest curves
% [statsCGrouped50]  = computeAllScores(ridgesCGrouped50, str2double(iid), gt,...
%     nThresh,[rpbModelColor.name '_grouped50'],statsCGrouped50,i,directories.scores);
% 
% ridgesCGrouped100 = fpg(ridgesC.thin,100); % grouping 100 strongest curves
% [statsCGrouped100] = computeAllScores(ridgesCGrouped100, str2double(iid), gt,...
%     nThresh,[rpbModelColor.name '_grouped100'],statsCGrouped100,i,directories.scores);
% 
% ridgesSGrouped50  = fpg(ridgesS.thin,50);  % grouping 50 strongest curves + spectral
% [statsSGrouped50] = computeAllScores(ridgesSGrouped50, str2double(iid), gt,...
%     nThresh,[rpbModelSpectral.name '_grouped50'],statsSGrouped50,i,directories.scores);
% 
% ridgesSGrouped100  = fpg(ridgesS.thin,100); % grouping 100 strongest curves + spectral
% [statsSGrouped100] = computeAllScores(ridgesSGrouped100, str2double(iid), gt,...
%     nThresh,[rpbModelSpectral.name '_grouped100'],statsSGrouped100,i,directories.scores);
