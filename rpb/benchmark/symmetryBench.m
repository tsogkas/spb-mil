function symmetryBench(test_iids,groundPath,betaPath,scorePath,nthresh)
% function symmetryBench2(test_iids,groundPath,featPath,betaPath,scorePath,nthresh)
% 
% Evaluate the performance of various ridge/symmetry detectors over the 
% full testing dataset.
% 
% INPUTS
% test_iids:    list with test images.
% groundPath:   path for ground truth.
% betaPath:     path for beta coefficients.
% scorePath:    path for folder where scores are stored.
% nthresh:      length of threshold vector.
% 
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% Setup

if nargin<5, nthresh=30; end

% Load Multiple Instance Learning coefficients
% The names of the .mat files must be set according to the files that
% resulted from training.
betaColor = load(fullfile(betaPath, 'beta_pyramid_balanced_nospb.mat'),'beta');
betaColor = betaColor.beta;
betaGray = load(fullfile(betaPath, 'beta_pyramid_balanced_gray.mat'),'beta');
betaGray = betaGray.beta;
betaNoTexture = load(fullfile(betaPath, 'beta_pyramid_balanced_notexture.mat'),'beta');
betaNoTexture = betaNoTexture.beta;
betaSpectral = load(fullfile(betaPath, 'beta_ncut_balanced.mat'),'beta');
betaSpectral = betaSpectral.beta;

cntR_color = zeros(nthresh,1);
sumR_color = zeros(nthresh,1);
cntP_color = zeros(nthresh,1);
sumP_color = zeros(nthresh,1);
scores_color = zeros(length(test_iids),5);

cntR_gray = zeros(nthresh,1);
sumR_gray = zeros(nthresh,1);
cntP_gray = zeros(nthresh,1);
sumP_gray = zeros(nthresh,1);
scores_gray = zeros(length(test_iids),5);

cntR_notext = zeros(nthresh,1);
sumR_notext = zeros(nthresh,1);
cntP_notext = zeros(nthresh,1);
sumP_notext = zeros(nthresh,1);
scores_notext = zeros(length(test_iids),5);

cntR_spectral = zeros(nthresh,1);
sumR_spectral = zeros(nthresh,1);
cntP_spectral = zeros(nthresh,1);
sumP_spectral = zeros(nthresh,1);
scores_spectral = zeros(length(test_iids),5);

cntR_lind = zeros(nthresh,1);
sumR_lind = zeros(nthresh,1);
cntP_lind = zeros(nthresh,1);
sumP_lind = zeros(nthresh,1);
scores_lind = zeros(length(test_iids),5);


%% Calculate mean precision, recall and f-measure for test data(color, fullcolor and gray)

for i=1:length(test_iids) 
    %% Read image and load groundtruth data
    i
    iid = test_iids(i);
    im = imgRead(iid);
    if exist(fullfile(groundPath,sprintf('groundtruth_%d.mat',iid)),'file')
        load(fullfile(groundPath,sprintf('groundtruth_%d',iid)), 'ridgeUnion')
    
        %% Multiple Instance Learning

        disp('Extracting histogram and spectral features')
        [dlcFeat,drcFeat,dlrFeat] = histGradFeatures(im);
        try
            load(sprintf('/home/tsogkas/symmetry_detector/Features/ncut_50/test/ncut_%d.mat',iid),'sPb')
            sPb = single(sPb);
            % spectral
            disp('Creating non-maximum suppressed image for MIL, with spectral feature')
            [supImSpectral] = rpb(im,betaSpectral,'spectral',dlcFeat,drcFeat,dlrFeat,sPb);
            sflag = 1;
        catch
            sflag = 0;
            disp('Was not able to load spectral feature!!')
        end
        % color
        disp('Creating non-maximum suppressed image for MIL, color features')
        [supImColor] = rpb(im,'color',betaColor,false,dlcFeat,drcFeat,dlrFeat);
        % gray
        disp('Creating non-maximum suppressed image for MIL, gray features')
        [supImGray] = rpb(im,'gray',betaGray,false,dlcFeat,drcFeat,dlrFeat);
        % no texture
        disp('Creating non-maximum suppressed image for MIL, no texture feature')
        [supImNoTexture] = rpb(im,'no_texture',betaNoTexture,false,dlcFeat,drcFeat,dlrFeat);

         %% Lindeberg 
        disp('Lindeberg ridge detector')
        input_image = double(rgb2gray(im));
        [~,~,~,contours] = PS00___primal_sketch(input_image);
        contour_ridge = PSzz_unzip_contour(contours{1}); 


        %% Create precision-recall curves and save scores

        ridgeUnion = ridgeUnion>0; % turn to logical    
        disp('Calculate precision-recall scores')
        % mil color
        [cntR_color,sumR_color,cntP_color,sumP_color,scores_color] = ...
            prScores(supImColor,iid,ridgeUnion,nthresh,'mil_color',...
            cntR_color,sumR_color,cntP_color,sumP_color,scores_color,i,scorePath);
        % mil gray
        [cntR_gray,sumR_gray,cntP_gray,sumP_gray,scores_gray] = ...
            prScores(supImGray,iid,ridgeUnion,nthresh,'boundfeat_noboundval',...
            cntR_gray,sumR_gray,cntP_gray,sumP_gray,scores_gray,i,scorePath);
        % mil no texture
        [cntR_notext,sumR_notext,cntP_notext,sumP_notext,scores_notext] = ...
            prScores(supImNoTexture,iid,ridgeUnion,nthresh,'mil_notext',...
            cntR_notext,sumR_notext,cntP_notext,sumP_notext,scores_notext,i,scorePath);
        if sflag
            % mil spectral
            [cntR_spectral,sumR_spectral,cntP_spectral,sumP_spectral,scores_spectral] = ...
                prScores(supImSpectral,iid,ridgeUnion,nthresh,'mil_spectral',...
                cntR_spectral,sumR_spectral,cntP_spectral,sumP_spectral,scores_spectral,i,scorePath);
        end
        % Lindeberg
        [cntR_lind,sumR_lind,cntP_lind,sumP_lind,scores_lind] = ...
            prScores(contour_ridge,iid,ridgeUnion,nthresh,'lind',...
            cntR_lind,sumR_lind,cntP_lind,sumP_lind,scores_lind,i,scorePath);    
    end
end


%% compute f-measure from recall and precision
function [f] = fmeasure(r,p)
f = 2*p.*r./(p+r+((p+r)==0));


%% interpolate to find best F and coordinates thereof
function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
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

