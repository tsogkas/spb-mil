%% Extract Normalized Cuts feature for symmetry

iids = imgList('test'); % BSDS300 image list ('test' or 'train')
%  paths for features, ground-truth and saved ncut features
% featPath = fullfile('/home','tsogkas','symmetry_detector','Features','train_dense_integral');
% groundPath = fullfile('/home','tsogkas','symmetry_detector','Groundtruth','train');
% ncutPath = fullfile('/home','tsogkas','symmetry_detector','Features','ncut_50_pyramid_balanced','train');
% beta = load(fullfile(betaPath,'beta_pyramid_balanced_nospb.mat'),'beta');
% beta = beta.beta;

nvec = 50;       % number of eigenvectors
radius = 5;     % radius for intervening contour cue
sigma = 0.1;    % sigma for intervening contour cue

for iter = 1:length(iids)
    iid = iids(iter);
    fprintf('Extracting features for image %d\n',iid);
    if exist(fullfile(groundPath,sprintf('groundtruth_%d.mat',iid)),'file')
        im = imgRead(iid);
        [sym,pb,orientMap] = rpb(im);
        [sPb eigVec] = ncutsym(sym,nvec,radius,sigma);
        sPb = single(sPb);
        save(fullfile(ncutPath,sprintf('ncut_%d.mat',iid)),'sPb')
    end
end
