function [f,y] = sampleDetectorMIL(iids,groundPath,featPath,numberOfSamples,...
                    samplingMethod,featConfig,boundPerc,buffer,betaRetrain)

% function [f,y] = sampleDetectorMIL(iids,groundPath,featPath,numberOfSamples,...
%                        samplingMethod,featConfig,buffer,boundPerc,betaRetrain)
%
% Sample on and off-symmetry axis pixels from the groundtruth training images,
% returning 0|1 class labels in y with the associated feature 
% vectors in f. 
%
% INPUT
%   iids:           the training images vector
%   groundPath:     the directory path where the groundtruth ridge maps are stored
%   featPath:       the directory path where the features used for training are stored
%	numberOfSamples:Approximate number of samples used. 
%   samplingMethod: one of {'random','balanced','boundary'}
%	featConfig:     one of {'gray','color','color_spectral','color_no_boundary',
%                   'no_texture','no_boundary','gray_no_boundary','fullcolor'}.
%	[buffer=3]:     Buffer zone around ridge pixels where we
%                   don't take off-ridge axis samples.
%   boundPerc:      Percentage of the boundary points that are sampled (n
%                   the case samplingMethod == 'balanced'.
%
% OUTPUT
%	f		Feature vectors (NxKxnInst); N=#samples, K=#features, nInst =
%       	#instances in each bag.
%	y		Vector (Nx1) of 0|1 class labels (1=symmetry axis).
%
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

nPer = ceil(numberOfSamples/numel(iids)); % number of samples per image
fprintf(2,'  %d samples per image...\n', nPer);
y = []; f = [];

for i = 1:numel(iids)
    %% Load pb map and symmetry groundtruth
    tic;  
    iid = iids(i);
    im = imgRead(iid);
    if exist(fullfile(groundPath,sprintf('groundtruth_%d.mat',iid)),'file')
        fprintf(2,'Processing image %d/%d (iid=%d)...\n',i,numel(iids),iid);
        load(fullfile(groundPath,sprintf('groundtruth_%d',iid)), 'ridgeUnion');
        segs = readSegs('color',iid);
        bmap = zeros(size(segs{1}));
        for j = 1:numel(segs),
            bmap = bmap | seg2bmap(segs{j});
        end  
        gt = bwmorph(ridgeUnion,'thin');
        gt(bmap) = 0;  
        dmapSkel = bwdist(gt);
        fprintf(2,'  Sampling...\n');

        %% Extract features
        fprintf(2,'  Extracting features...\n');
        if strcmp(featConfig,'spectral') 
            load(fullfile(featPath,'ncut_50','train',sprintf('ncut_%d.mat',iid)),'sPb')
            [dlcFeat, drcFeat, dlrFeat] = histGradFeatures(im);
            [h,w,norient,nscales,~] = size(dlcFeat);
            sPb = repmat(sPb,[1,1,1,nscales]);
            features = cat(5, ones(h,w,norient,nscales), dlcFeat, drcFeat, dlrFeat, sPb);
        else
            [~,~,~,features] = rpb(im,featConfig,[],true);
        end
            [h,w,norient,nscales,K] = size(features);
            features = reshape(permute(features, [1,2,5,3,4]),[h*w,K,norient*nscales]);

        %% Choose sampling method
        onidx = find(gt);       % skeleton pixels
        offidx = find(dmapSkel>buffer); % non-skeleton/non boundary pixels
        oncnt = numel(onidx);
        offcnt = numel(offidx);
        onused = min(floor(0.5*nPer),oncnt); % 50% positives and 50% negatives
        switch samplingMethod
            case 'random'   % random sampling 
                ind = [ onidx; offidx ];
                cnt = numel(ind);
                idx = randperm(cnt);
                idx = idx(1:min(cnt,nPer));
                ind = ind(idx);
                fprintf(2,'  %d on-pixels and %d off-pixels\n',oncnt,offcnt);
            case 'balanced' % 50% positives, 50% negatives, no boundary samples 
                fprintf(2,'  %d on-pixels and %d off-pixels\n',oncnt,offcnt);
                rest = nPer-onused; % the rest are negatives
                ind = [randsample(onidx,onused); randsample(offidx,rest)];
            case 'boundary' % a percentage of boundary pixels is explicitly used
                if isempty(betaRetrain), boundidx = find(bmap);
                else boundidx = find(bmap & pbth); end
                boundcnt = length(boundidx); % number of boundary pixels
                fprintf(2,'  %d on-pixels %d boundary-pixels and %d off-pixels\n',...
                  oncnt,boundcnt,offcnt);
                boundused = min(floor(boundPerc*(nPer-onused)),boundcnt);
                rest = nPer-onused-boundused;
                ind = [randsample(onidx,onused); randsample(boundidx,boundused);...
                    randsample(offidx,rest)];  % use the strongest boundaries. the rest of the pixels used are off pixels
            otherwise
                error('sampling method must be ''random'', ''balanced'', or ''boundary''')
        end

        pos = sum(gt(ind));
        total = numel(ind);
        y = [ y; gt(ind) ];
        f = cat(1,f,features(ind,:,:));
        fprintf(2,'  %d positives and %d negatives for a total of %d samples.\n',...
            pos,total-pos,total);
    end
toc;
end
