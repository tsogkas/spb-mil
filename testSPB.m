function testSPB(models, testOpts)

opts.testSet = 'val';
opts.nThresh = 30;
opts.slack   = 0.01;
opts = parseVarargin(opts,testOpts);

% Read test images --------------------------------------------------------
paths = setPaths();
switch opts.testSet
    case 'val'
        imageList = dir(fullfile(paths.bsds500imVal, '*jpg'));
        imPath    = paths.bsds500imVal;
        gtPath    = paths.symmax500gtVal;
    case 'test'
        imageList = dir(fullfile(paths.bsds500imTest, '*jpg'));
        imPath    = paths.bsds500imTest;
        gtPath    = paths.symmax500gtTest;
    otherwise
        error('testSet can be either ''val'' or ''test''')
end

% Load models -------------------------------------------------------------
for i=1:numel(models)
    tmp = load(fullfile(paths.modelsMIL, models{i}.name)); 
    models{i} = tmp.spbModel;
    models{i}.stats = [];
end
statsLind  = [];
statsLevin = [];

% Load parameters used for Levinshtein model
load classifier_params_svm.mat
load part_classifier_params.mat

opts.nImages = numel(imageList);
for i=1:opts.nImages
    fprintf('Testing image %d/%d from BSDS500 %s set\n', i, opts.nImages, opts.testSet);
    [~,iid,~] = fileparts(imageList(i).name);
    gt = load(fullfile(gtPath,['gt_' iid ])); gt = gt.gt; % load image groundtruth
    
    % TODO: add writing results to images
    disp('Running MIL symmetry models...')
    im       = im2double(imread(fullfile(imPath,imageList(i).name)));
    features = computeHistogramFeatures(im);
    for m=1:numel(models)
        spb = evaluateModelMIL(models{m},im,features);
        model{m}.stats = computeStats(spb.thin, gt, model{m}.stats, i, opts);
    end
    disp('Running Lindeberg ridge detector...')
    lind  = evaluateLindeberg(im);
    statsLind = computeStats(lind,gt,statsLind,i,opts);

    disp('Running Levinshtein symmetry axes detector...')
    levin = evaluateLevinshtein(im);
    statsLevin = computeStats(levin,gt,statsLevin,i,opts);
end

disp('Computing precision, recall, f-measure')
for m=1:numel(models)
    model{m}.stats.precision = model{m}.stats.cntP ./ model{m}.stats.sumP;
    model{m}.stats.recall    = model{m}.stats.cntR ./ model{m}.stats.sumR;
    % add ODS, OIS, AP
end
statsLind.precision  = statsLind.cntP  ./statsLind.sumP;
statsLind.recall     = statsLind.cntR  ./statsLind.sumR;
statsLevin.precision = statsLevin.cntP ./statsLevin.sumP;
statsLevin.recall    = statsLevin.cntR ./statsLevin.sumR;

% Plot precision-recall curves --------------------------------------------
plotPrecisionRecall(0)  % compare methods
plotPrecisionRecall(1)  % compare grouping


% -------------------------------------------------------------------------
function spb = evaluateModelMIL(model,im,features)
% -------------------------------------------------------------------------
paths = setPaths();
switch model.featureSet
    case 'spectral'
        try
            spectralFeat = load(fullfile(paths.spectral,'spectral_50',opts.testSet,['spectral_' iid]));
            spectralFeat = single(spectralFeat.spectral);
            spb = spbMIL(im,'spectral',true,model.w,features,spectralFeat);
        catch
            warning('Was not able to load spectral feature')
        end
    case 'color'
        spb = spbMIL(im, 'color',true, model.w, features, []);
    case 'gray'
        spb = spbMIL(im, 'gray', true, model.w, features, []);
    case 'no_texture'
        spb = spbMIL(im,'no_texture',true,model.w,features,[]);
    otherwise
        error('Invalid  feature type')
end

%    ridgesCGrouped50  = fpg(ridgesC.thin,50);  % grouping 50 strongest curves
%    [statsCGrouped50]  = computeAllScores(ridgesCGrouped50, str2double(iid), gt,...
%        nThresh,[rpbModelColor.name '_grouped50'],statsCGrouped50,i,directories.scores);
%
%    ridgesCGrouped100 = fpg(ridgesC.thin,100); % grouping 100 strongest curves
%    [statsCGrouped100] = computeAllScores(ridgesCGrouped100, str2double(iid), gt,...
%        nThresh,[rpbModelColor.name '_grouped100'],statsCGrouped100,i,directories.scores);
%
%     ridgesSGrouped50  = fpg(ridgesS.thin,50);  % grouping 50 strongest curves + spectral
%     [statsSGrouped50] = computeAllScores(ridgesSGrouped50, str2double(iid), gt,...
%         nThresh,[rpbModelSpectral.name '_grouped50'],statsSGrouped50,i,directories.scores);
%
%     ridgesSGrouped100  = fpg(ridgesS.thin,100); % grouping 100 strongest curves + spectral
%     [statsSGrouped100] = computeAllScores(ridgesSGrouped100, str2double(iid), gt,...
%         nThresh,[rpbModelSpectral.name '_grouped100'],statsSGrouped100,i,directories.scores);


% -------------------------------------------------------------------------
function ridges = evaluateLindeberg(im)
% -------------------------------------------------------------------------
[~,~,~,contours] = PS00___primal_sketch(double(rgb2gray(im)));
ridges  = PSzz_unzip_contour(contours{1});
%     figure; imshow(ridgesL,[]); title(['Lindeberg iid = ' num2str(iid)]);

% -------------------------------------------------------------------------
function axes = evaluateLevinshtein(im)
% -------------------------------------------------------------------------
[parts, part_scales, part_score, image_data] = ...
    DetectSymmetricParts_SingleImage(fullfile(testImDir,testImages(i).name),...
    classifier_params,[], [], [], [], true);
%     figure; DisplayParts(im, parts, part_scales, 1:size(parts,3));
[part_labels, part_affinity, part_adjacency] = ...
    GroupSymmetricParts_SingleImage(im, parts, part_classifier_params, image_data);
%     figure; DrawPartClustersEllipses(im, parts, part_labels);
selected_parts_indicator = SelectSymmetricPart_SingleImage(parts, ...
    part_labels, part_score, part_affinity, part_adjacency,0.1);
%     figure; DrawPartClustersEllipses(im, parts, part_labels, selected_parts_indicator);
axes = levsym(im,parts, part_labels, selected_parts_indicator);

% -------------------------------------------------------------------------
function [stats] = computeStats(ridges,gt,stats,iter,opts)
% -------------------------------------------------------------------------
if islogical(ridges), 
    nThresh = 1;
    ridges = double(ridges);
else
    nThresh = max(1,opts.nThresh);
end
thresh = linspace(1/(nThresh+1),1-1/(nThresh+1),nThresh)';

% the first time, initialize by zeroing all counts
if isempty(stats)
    stats.cntR      = zeros(numel(thresh),opts.nImages);
    stats.sumR      = zeros(numel(thresh),opts.nImages);
    stats.cntP      = zeros(numel(thresh),opts.nImages);
    stats.sumP      = zeros(numel(thresh),opts.nImages);
end
for t = 1:nThresh,
    %TODO: add code to handle case where we have multiple gt skeleton maps
    % threshold pb to get binary symmetry map and make 1-pixel thin
    bmap = (ridges >= thresh(t));
    bmap = bwmorph(bwmorph(bmap,'thin',Inf), 'clean');
    [match1,match2] = correspondPixels(double(bmap),double(gt),slack);
    if opts.plot
        plotMatch(1,bmap,gt,match1,match2); 
    end
    stats.sumR(t,iter) = nnz(gt(:));       % tp + fn
    stats.cntR(t,iter) = nnz(match2(:)>0); % tp
    stats.sumP(t,iter) = nnz(bmap(:));     % tp + fp
    stats.cntP(t,iter) = nnz(match1(:)>0); % tp
    assert(stats.cntP(t,iter) == stats.cntR(t,iter))
%     for i = 1:nskels,
%         %fwrite(2,'+');
%         % compute the correspondence
%         [match1,match2] = correspondPixels(bmap,skels(:,:,i),0.01);
%         % accumulate machine matches
%         accP = accP | match1;
%         % compute recall
%         sumR(t) = sumR(t) + sum(sum(skels(:,:,i)));
%         cntR(t) = cntR(t) + sum(match2(:)>0);
%     end
%     % compute precision
%     sumP(t) = sumP(t) + sum(bmap(:));
%     cntP(t) = cntP(t) + sum(accP(:));

end

% -------------------------------------------------------------------------
function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
% -------------------------------------------------------------------------
% Interpolate to find best F and coordinates thereof 
bestT = thresh(1);
bestR = R(1);
bestP = P(1);
bestF = fmeasure(R(1),P(1));
for i = 2:numel(thresh),
  for d = linspace(0,1),
    t = thresh(i)*d + thresh(i-1)*(1-d);
    r = R(i)*d + R(i-1)*(1-d);
    p = P(i)*d + P(i-1)*(1-d);
    f = fmeasure(r,p);
    if f > bestF,
      bestT = t;
      bestR = r;
      bestP = p;
      bestF = f;
    end
  end
end



