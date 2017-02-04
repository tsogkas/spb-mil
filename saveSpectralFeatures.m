function saveSpectralFeatures(model, nvec, radius, sigma)
% Compute spectral feature for all images in one of the train, val, test
% sets of the BSDS500 dataset.
% 
%   saveSpectralFeatures(set,nvec,radius,sigma);
% 
%   model:  model used to compute soft symmetry responses that are fed as
%           input to the spectral clustering algorithm.
%   nvec:   number of eigenvectors used {default = 50}.
%   radius: radius for intervening contour cue {default = 5}.
%   sigma:  sigma for intervening contour cue. {default = 0.1}
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: February 2017

if nargin<2, nvec   = 50;      end  % number of eigenvectors computed
if nargin<3, radius = 5;       end  % radius for intervening contour cue
if nargin<4, sigma  = 0.1;     end  % sigma for intervening contour cue

% Load model 
if ischar(model)
    spbModel = load(fullfile(paths.models,model)); spbModel = spbModel.spbModel;
elseif isstruct(model)
    spbModel = model;
else error('Input model must be a path to a mat-file or a struct.')
end

% Create directory where spectral features will be stored
mkdir(fullfile(path.spectral,['eigenvec' num2str(nvec)]));

% Read image paths
imageList = [dir(fullfile(paths.bsds500imTrain,'*jpg'));
             dir(fullfile(paths.bsds500imVal,  '*jpg'));
             dir(fullfile(paths.bsds500imTest, '*jpg'))];
         
% Compute and store spectral features
for i = 1:length(imageList)
    [~,iid,~] = fileparts(imageList(i).name);
    savePath = fullfile(path.spectral, ['eigenvec' num2str(nvec)], ['spectral_' num2str(iid) '.mat']);
    if ~exist(savePath,'file')
        fprintf('Extracting spectral features for image %d/%d (iid = %s)\n',...
            set, i, length(imageList), iid);
        img = im2double(imread(fullfile(imagePath, imageList(i).name)));
        spb = spbMIL(img,spbModel.featureSet,true,spbModel.w);
        [spectral,eigVec] = computeSpectralFeature(spb.thin,nvec,radius,sigma);
        spectral = single(spectral);
%         for iv=1:size(eigVec,3)
%             fig = figure; imagesc(eigVec(:,:,iv)); 
%             print(fig,'-depsc2',sprintf('eigVec%d_%s.eps',iv+1,iid))
%         end
        save(savePath,'spectral')
    end
end

