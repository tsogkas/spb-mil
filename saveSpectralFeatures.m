function saveSpectralFeatures(set, nvec, radius, sigma)
% Compute spectral feature for all images in one of the train, val, test
% sets of the BSDS500 dataset.
% 
%   computeAllSpectralFeatures(set,nvec,radius,sigma);
% 
%   set:    one of {'train', 'val', 'test'}.
%   nvec:   number of eigenvectors used {default = 50}.
%   radius: radius for intervening contour cue {default = 5}.
%   sigma:  sigma for intervening contour cue. {default = 0.1}
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last update: August 2013

if nargin<1, set    = 'train'; end
if nargin<2, nvec   = 50;      end  % number of eigenvectors computed
if nargin<3, radius = 5;       end  % radius for intervening contour cue
if nargin<4, sigma  = 0.1;     end  % sigma for intervening contour cue

% Set paths
directories   = setupDirectories();
switch set
    case 'train'
        imagePath = directories.bsds500imTrain;
    case 'val'
        imagePath = directories.bsds500imVal;
    case 'test'
        imagePath = directories.bsds500imTest;
    otherwise
        error('ERROR:Wrong dataset!')
end
s = load(fullfile(directories.models,'rpbModel_600K_color_loglikeNOR_balanced_trainval.mat'));
w = s.rpbModel.w; 
savePath = fullfile(directories.spectral,['spectral_' num2str(nvec)],set);
imageList = dir(fullfile(imagePath,'*jpg'));
for i = 1:length(imageList)
    [~,iid,~] = fileparts(imageList(i).name);
    if ~exist(fullfile(savePath,['spectral_' num2str(iid) '.mat']),'file')
        fprintf('Extracting spectral features from %s set image %d/%d (iid = %s)\n',...
            set, i, length(imageList), iid);
        im = im2double(imread(fullfile(imagePath, imageList(i).name)));
        ridges = rpb(im,'color',true,w);
        [spectral,eigVec] = computeSpectralFeature(ridges.thin,nvec,radius,sigma);
        spectral = single(spectral);
%             for iv=1:size(eigVec,3)
%                 fig = figure; imagesc(eigVec(:,:,iv)); print(fig,'-depsc2',sprintf('eigVec%d_%s.eps',iv+1,iid))
%             end
        save(fullfile(savePath,['spectral_' num2str(iid)]),'spectral')
    end
end

