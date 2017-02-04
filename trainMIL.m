function spbModel = trainMIL(varargin)
% TRAINMIL Train symmetry detector using Multiple Instance Learning.
% 
%   spbModel = TRAINMIL(trainOpts), trainOpts can either be a cell array,
%              or a struct, containing pairs of arguments/values.
% 
%   INPUTS (default values in brackets)
%   tag:              optional tag for the resulting file and model name {[]}
%   trainSet:         one of ['train','trainval'] {'trainval'}. This can
%                     also be a struct with img/seg fields containing the
%                     images and groundtruth used for training.
%   featureSet:       one of ['spectral','color', 'gray', 'no_texture'] {'color'}. 
%   nSamplesPerImage: total number of samples per image for training {1000}.
%   cost:             cost function to optimize using MIL. One of:
%                     ['nor','max','log'] {'nor'}.
%   nOptIter:         maximum number of optimization iterations {200}.
%   sampling:         one of ['random','balanced','boundary'] {balanced}.
%   slack:            determines the area around groundtruth ridge pixels, 
%                     from where we do not choose training samples.
% 
%   OUTPUTS
%   spbModel: Struct containing learned weights and other training information.
% 
%   ADDITIONAL INFO
%   xstd:     Standard deviation of the vector weights.
%   L:        The final value of the log-likelihood after the optimization.
%   x:        The nSamples x nInstPerBag x nDim array, carrying the features
%             for all bags corresponding to each one of the training samples. 
%             nSamples is the total number of training samples, nDim is the
%             length of the feature vector and nInstPerBag is the number of 
%             instances included in each bag.
%   y:        The 0/1 labels for each sample.
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last upadate: February 2017

% Default training opts ---------------------------------------------------
paths = setPaths;
opts.nSamplesPerImage = 1e3;    % number of samples per image
opts.tag        = [];           % results description
opts.trainSet   = 'trainval';
opts.featureSet = 'color';
opts.cost       = 'nor';        % cost function to optimize
opts.nOptIter   = 200;          % number of optimization iterations
opts.sampling   = 'balanced';
opts.slack      = 3;
opts            = parseVarargin(opts, varargin);

% Setup image files or data used for training. If we provide a struct
% containing image and groundtruth data, then opts.imageList contains the
% data and opts.trainSet is inferred by the size of opts.imageList.
if ischar(opts.trainSet) && strcmp(opts.trainSet, 'train')
    opts.imageList = dir(fullfile(paths.bsds500imTrain, '*jpg'));
elseif ischar(opts.trainSet) && strcmp(opts.trainSet, 'trainval')
    opts.imageList = [dir(fullfile(paths.bsds500imTrain, '*jpg'));...
                      dir(fullfile(paths.bsds500imVal,   '*jpg'))];
elseif isstruct(opts.trainSet) 
    opts.imageList = opts.trainSet;
    if numel(opts.trainSet) == 200, opts.trainSet = 'train';
    elseif numel(opts.trainSet) == 300, opts.trainSet = 'trainval';
    else error(['The training data provided should either be 200 images ',...
            '(train set) or 300 images (trainval set).'])
    end
else error('Invalid training set')
end
    
% Build model name out of parameters --------------------------------------
assert(ismember(opts.cost,{'nor','max','log'}),'Invalid cost function');
spbModel.opts = opts;
spbModel.name = sprintf('spbModel_%dK_%s_%s_%s_%s', ...
    num2str(opts.nSamplesPerImage),opts.featureSet,opts.cost,opts.sampling,opts.trainSet);
if ~isempty(tag), spbModel.name = [spbModel.name '_' tag]; end
spbModelPath  = fullfile(paths.models,spbModel.name);

% Try to load existing model ----------------------------------------------
if exist(spbModelPath,'file')
    warning('You have already trained a model using this configuration. Loading trained model...')
    spbModel = load(spbModelPath,'spbModel'); spbModel = spbModel.spbModel;
    return
else
    disp('Setup successful. Training starting now...')    
end
            
% Sample features and labels (optional retraining) ------------------------
samplesPath = fullfile(paths.models,['samples_' num2str(opts.nSamplesPerImage) 'K_' sampling]);
if ~isempty(tag), samplesPath = [samplesPath '_' tag]; end
try
    tmp = load(samplesPath,'x','y'); x = tmp.x; y = tmp.y; clear tmp;
catch
    [x,y] = getSamples(opts.imageList,opts.featureSet,opts.nSamplesPerImage,opts.sampling,opts.slack);
    save(samplesPath,'x','y','-v7.3');
end
x = x(:,:,getFeatureSubset(featureSet));

% Normalize features to unit variance -------------------------------------
[nSamples,nInstPerBag,nDim] = size(x);
x           = double(x); % double precision is necessary for training     
xReshaped   = reshape(x,[nSamples*nInstPerBag,nDim]);
xstd        = std(xReshaped); xstd = xstd + (xstd==0);  % feature standard deviation
xReshaped   = bsxfun(@rdivide,xReshaped,xstd);
x           = reshape(xReshaped,[nSamples,nInstPerBag,nDim]); 

% Learn weights -----------------------------------------------------------
disp('Fitting model...')
nIterDone = 0;  % sometimes minimize stops prematurely
while nIterDone<20
    w_0 = rand(nDim,1);
    [w,L,nIterDone] = minimize(w_0,'loglikelihood',nOptIter,cost,y,x,xReshaped,true); % weights and log-likelihood vector
end

%     % --- Calculate probabilites for bag instances and bags (Noisy-OR)
%     disp('Calculating posterior probabilities...')
%     s1          = size(xReshaped,1);
%     in_prd      = reshape(xReshaped*w,[s1/nInstPerBag,nInstPerBag]);
%     p_inst      = 1./(1 + exp(-in_prd));
%     % using loglikeMIL or loglikeNOR
%     logp_bags   = sum(log(1-p_inst+eps),2);
%     p_bags      = 1-exp(logp_bags);
%     % % using loglikeMAX
%     [p_bags, instIdx] = max(p_inst,[],2);

% De-normalize beta coefficients ------------------------------------------
xstd            = xstd';
w               = w./xstd;
spbModel.xstd   = xstd;
spbModel.w      = w;
spbModel.w_0    = w_0;
spbModel.loglikelihood = L;
save(spbModelPath,'spbModel')   % save results


% -------------------------------------------------------------------------
function [f,y] = getSamples(trainImages,featureSet,nSamplesPerImage,sampling,slack)
% -------------------------------------------------------------------------
assert(slack>0,'slack must be a positive number!');
assert(nSamplesPerImage>0,'number of samples per image must be positive!');
disp(['Using ' num2str(nSamplesPerImage) ' samples per image'])
nImages = numel(trainImages);
y = []; f = [];

ticStart = tic;
for i = 1:nImages    
    if isfield(trainImages(i), 'isdir') 
        % Read image and groundtruth from disk
        [~,imageName] = fileparts(trainImages(i).name);
        img = im2double(imread(trainImages(i).name));
        sgt = load(['gt_' imageName '.mat']); sgt = sgt.gt; 
        gt  = load([imageName '.mat']); gt = gt.groundTruth;
        sgt = bwmorph(sgt,'thin','inf');
        bgt = false(size(gt{1}.Boundaries));
        for s=1:numel(gt)
            bgt = bgt | gt{i}.Boundaries;
        end
    else 
        % Read image and groundtruth from struct
        img = trainImages(i).img;
        [H,W,~] = size(img);
        sgt = false(H,W);
        bgt = false(H,W);
        for s=1:size(trainImages(i).seg,3)
            sgt = sgt | trainImages(i).seg(:,:,s);
            bgt = bgt | trainImages(i).bnd(:,:,s);
        end
    end
    
    % Compute histogram and spectral features
    if strcmp(featureSet,'gray'), img = rgb2gray(img); end
    fprintf('Sampling image %d/%d (iid = %s)...', i, length(trainImages), imageName);
    histf       = computeHistogramFeatures(img,true);
    [height,width,nOrient,nScales,~] = size(histf.dlc);
    if strcmp(featureSet,'spectral')
        spectralFeat = load(['spectral_' imageName '.mat']);
        spectralFeat = repmat(spectralFeat.spectral,[1 1 1 nScales]);
    else
        spectralFeat = [];
    end
    features = cat(5, ones(height,width,nOrient,nScales,1,'single'),...
                    histf.dlc, histf.drc, histf.dlr,spectralFeat);
    features = reshape(features,height*width,nOrient*nScales,[]);

    % Sample features
    distanceMap = bwdist(sgt);
    ind = getSampleIndexes(sgt,distanceMap,nSamplesPerImage,slack,sampling,boundaries);
    fprintf('Added %d/%d positives and %d/%d negatives for a total of %d samples.\n',...
            nnz(sgt(ind)),nnz(sgt),nnz(~sgt(ind)),nnz(~sgt),numel(ind));

    y = cat(1,y,sgt(ind));
    f = cat(1,f,features(ind,:,:));
    progress('Storing sample features...',i,nImages,ticStart,0);
end


% -------------------------------------------------------------------------
function ind = getSampleIndexes(sgt,distanceMap,nSamplesPerImage,slack,sampling,boundaries)
% -------------------------------------------------------------------------
% Get linear indexes of training instances depending on the sampling scheme
indPos   = find(sgt);    
indNeg   = find(distanceMap>slack);
nPos     = numel(indPos);
nNeg     = numel(indNeg);
switch sampling
    case 'random'   % random sampling 
        ind      = [ indPos; indNeg ];
        ind      = randsample(ind,nSamplesPerImage);
        nPosUsed = nnz(sgt(ind));
        nNegUsed = nSamplesPerImage-nPosUsed;
    case 'balanced' % 50% positives, 50% negatives, no boundary samples 
        nPosUsed = min(floor(0.5*nSamplesPerImage),nPos); % 50% positives and 50% negatives
        nNegUsed = nSamplesPerImage-nPosUsed; % the rest are negatives
        ind      = [randsample(indPos,nPosUsed); randsample(indNeg,nNegUsed)];
    case 'boundary' % a percentage of boundary pixels is explicitly used NOT SURE IF WILL USE THIS
        boundidx = find(boundaries);
        boundcnt = length(boundidx); % number of boundary pixels
        disp([num2str(nPos) ' on-pixels, ' num2str(boundcnt)...
            ' boundary pixels and ' num2str(nNeg) ' off-pixels'])
        boundused = min(floor(boundPerc*(nSamplesPerImage-nPosUsed)),boundcnt);
        rest = nSamplesPerImage-nPosUsed-boundused;
        ind = [randsample(indPos,nPosUsed); randsample(boundidx,boundused);...
            randsample(indNeg,rest)];  % use the strongest boundaries. the rest of the pixels used are off pixels
    otherwise
        error('Sampling must be ''random'', ''balanced'', or ''boundary''')
end

% -------------------------------------------------------------------------
function [X, fX, i, XX] = minimize(X, f, length, P1, P2, P3, P4, P5)
% -------------------------------------------------------------------------
% Written by Carl E. Rasmussen
% 
% Minimize a continuous differentialble multivariate function. Starting point
% is given by "X" (D by 1), and the function named in the string "f", must
% return a function value and a vector of partial derivatives. The Polack-
% Ribiere flavour of conjugate gradients is used to compute search directions,
% and a line search using quadratic and cubic polynomial approximations and the
% Wolfe-Powell stopping criteria is used together with the slope ratio method
% for guessing initial step sizes. Additionally a bunch of checks are made to
% make sure that exploration is taking place and that extrapolation will not
% be unboundedly large. The "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. The
% function returns when either its length is up, or if no further progress can
% be made (ie, we are at a minimum, or so close that due to numerical problems,
% we cannot get any closer). If the function terminates within a few
% iterations, it could be an indication that the function value and derivatives
% are not consistent (ie, there may be a bug in the implementation of your "f"
% function). The function returns the found solution "X", a vector of function
% values "fX" indicating the progress made and "i" the number of iterations
% (line searches or function evaluations, depending on the sign of "length")
% used.
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, P4, P5)
%
% See also: checkgrad
%
% Copyright (C) 2001 by Carl Edward Rasmussen. Date 2001-07-18

RHO = 0.01;                            % a bunch of constants for line searches
SIG = 0.5;       % RHO and SIG are the constants in the Wolfe-Powell conditions
INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                    % extrapolate maximum 3 times the current bracket
MAX = 20;                         % max 20 function evaluations per line search
RATIO = 100;                                      % maximum allowed slope ratio

argstr = [f, '(X'];                      % compose string used to call function
for i = 1:(nargin - 3)
    argstr = [argstr, ',P', int2str(i)];
end
argstr = [argstr, ')'];

if length>0, S=['Linesearch']; else S=['Function evaluation']; end

i = 0;                                            % zero the run length counter
ls_failed = 0;                             % no previous line search has failed
fX = [];
XX = [];
[f1 df1] = eval(argstr);                      % get function value and gradient
XX=vertcat(XX,[1,-f1]);  %by Russ, get first value
i = i + (length<0);                                            % count epochs?!
s = -df1;                                        % search direction is steepest
d1 = -s'*s;                                                 % this is the slope
z1 = 1/(1-d1);                                      % initial step is 1/(|s|+1)

while i < abs(length)                                      % while not finished
    i = i + (length>0);                                      % count iterations?!
    
    X0 = X; f0 = f1; df0 = df1;                   % make a copy of current values
    X = X + z1*s;                                             % begin line search
    [f2 df2] = eval(argstr);
    i = i + (length<0);                                          % count epochs?!
    d2 = df2'*s;
    f3 = f1; d3 = d1; z3 = -z1;             % initialize point 3 equal to point 1
    if length>0, M = MAX; else M = min(MAX, -length-i); end
    success = 0; limit = -1;                     % initialize quanteties
    while 1
        while ((f2 > f1+z1*RHO*d1) || (d2 > -SIG*d1)) && (M > 0)
            limit = z1;                                         % tighten the bracket
            if f2 > f1
                z2 = z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3);                 % quadratic fit
            else
                A = 6*(f2-f3)/z3+3*(d2+d3);                                 % cubic fit
                B = 3*(f3-f2)-z3*(d3+2*d2);
                z2 = (sqrt(B*B-A*d2*z3*z3)-B)/A;       % numerical error possible - ok!
            end
            if isnan(z2) || isinf(z2)
                z2 = z3/2;                  % if we had a numerical problem then bisect
            end
            z2 = max(min(z2, INT*z3),(1-INT)*z3);  % don't accept too close to limits
            z1 = z1 + z2;                                           % update the step
            X = X + z2*s;
            [f2 df2] = eval(argstr);
            M = M - 1; i = i + (length<0);                           % count epochs?!
            d2 = df2'*s;
            z3 = z3-z2;                    % z3 is now relative to the location of z2
        end
        if f2 > f1+z1*RHO*d1 || d2 > -SIG*d1
            break;                                                % this is a failure
        elseif d2 > SIG*d1
            success = 1; break;                                             % success
        elseif M == 0
            break;                                                          % failure
        end
        A = 6*(f2-f3)/z3+3*(d2+d3);                      % make cubic extrapolation
        B = 3*(f3-f2)-z3*(d3+2*d2);
        z2 = -d2*z3*z3/(B+sqrt(B*B-A*d2*z3*z3));        % num. error possible - ok!
        if ~isreal(z2) || isnan(z2) || isinf(z2) || z2 < 0   % num prob or wrong sign?
            if limit < -0.5                               % if we have no upper limit
                z2 = z1 * (EXT-1);                 % the extrapolate the maximum amount
            else
                z2 = (limit-z1)/2;                                   % otherwise bisect
            end
        elseif (limit > -0.5) && (z2+z1 > limit)          % extraplation beyond max?
            z2 = (limit-z1)/2;                                               % bisect
        elseif (limit < -0.5) && (z2+z1 > z1*EXT)       % extrapolation beyond limit
            z2 = z1*(EXT-1.0);                           % set to extrapolation limit
        elseif z2 < -z3*INT
            z2 = -z3*INT;
        elseif (limit > -0.5) && (z2 < (limit-z1)*(1.0-INT))   % too close to limit?
            z2 = (limit-z1)*(1.0-INT);
        end
        f3 = f2; d3 = d2; z3 = -z2;                  % set point 3 equal to point 2
        z1 = z1 + z2; X = X + z2*s;                      % update current estimates
        [f2 df2] = eval(argstr);
        M = M - 1; i = i + (length<0);                             % count epochs?!
        d2 = df2'*s;
    end                                                      % end of line search
    
    if success                                         % if line search succeeded
        f1 = f2; fX = [fX' f1]';
        fprintf('%s %6i;  Value %4.6e\r', S, i, f1);
        XX=vertcat(XX,[i,-f1]); %minus here for our purpuses only
        s = (df2'*df2-df1'*df2)/(df1'*df1)*s - df2;      % Polack-Ribiere direction
        tmp = df1; df1 = df2; df2 = tmp;                         % swap derivatives
        d2 = df1'*s;
        if d2 > 0                                      % new slope must be negative
            s = -df1;                              % otherwise use steepest direction
            d2 = -s'*s;
        end
        z1 = z1 * min(RATIO, d1/(d2-realmin));          % slope ratio but max RATIO
        d1 = d2;
        ls_failed = 0;                              % this line search did not fail
    else
        X = X0; f1 = f0; df1 = df0;  % restore point from before failed line search
        if ls_failed || i > abs(length)          % line search failed twice in a row
            break;                             % or we ran out of time, so we give up
        end
        tmp = df1; df1 = df2; df2 = tmp;                         % swap derivatives
        s = -df1;                                                    % try steepest
        d1 = -s'*s;
        z1 = 1/(1-d1);
        ls_failed = 1;                                    % this line search failed
    end
end
fprintf('\r');

% -------------------------------------------------------------------------
function [L, dL] = loglikelihood(mode,w,y,x,xReshaped,deriv)
% -------------------------------------------------------------------------
% [lli, lliDeriv] = loglikeNOR(mode,w,y,x,xReshaped,deriv)
%
%   INPUTS
%   =======================================================================
%   mode: 'nor','max','log'.
% 	y: 	nSamples x 1 column vector of 0|1 class assignments.
% 	x: 	nSamples x nInstPerBag x nDim array of feature vectors.
%       nInstPerBag = nOrient*nScales is the total number of instances per
%       bag of features. Each bag contains the features for all orientations
%       and scales for the corresponding pixel.
%   w:  classifier coefficients.
%   xReshaped: reshaped (nSamples*nInstPerBag) x nDim version of x to
%              speedup computations.
%   deriv: if true, calculate log-likelihood derivative
%
%   OUTPUTS
%   =======================================================================
%   lli: log likelihood value at point w
%   llgrad: log likelihood gradient with respect to w vector, calculated
%   at point w.
%
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last revision: June 2013

if nargin<5, deriv = true; end

switch mode
    case 'nor'
        %   Log-likelihood and its gradient for the Multiple Instance Learning
        %   paradigm, using noisy-or.
        nInstPerBag = size(x,2);
        s1          = size(xReshaped,1);
        in_prd      = reshape(xReshaped*w,[s1/nInstPerBag,nInstPerBag]);
        p_inst      = 1./(1 + exp(-in_prd)); % bag and bag-instances probabilities (Noisy-OR)
        logp_bags   = sum(log(1-p_inst+eps),2);
        p_bags      = exp(logp_bags); % p_bags = Prod(1-p_inst) so that we can use log1p later
        L           = -y' * log1p(-p_bags) - (1-y)' * logp_bags; % minus log-likelihood
        if deriv
            tc       = (y-(1-p_bags))./(1-p_bags);
            dL = -(tc' * squeeze(sum(bsxfun(@times,x,p_inst), 2)))'; % log-likelihood derivative with respect to w
        end
    case 'max'
        %   Log-likelihood and its gradient for the Multiple Instance Learning
        %   paradigm. Here we don't use the bag approach as in the paper by Viola;
        %   instead we use a simplification, taking into account only the strongest
        %   instance from each bag. (p_i = max(p_ij)).
        [nSamples,nInstPerBag,nDims] = size(x); % get sizes
        s1               = size(xReshaped,1);
        in_prd           = reshape(xReshaped*w,[s1/nInstPerBag,nInstPerBag]);
        p_inst           = 1./(1 + exp(-in_prd));
        [p_bags,indInst] = max(p_inst,[],2);
        linIdx           = [];
        for iDim=1:nDims % get linear indices for instances of strongest instances
            linIdx = [linIdx sub2ind(size(x),1:nSamples,indInst',iDim*ones(1,nSamples))];
        end
        L        = -y' * log(p_bags) - (1-y)' * log(1-p_bags);   % minus log-likelihood
        if deriv
            x       = reshape(x(linIdx'),nSamples,nDims);
            dL = -(x' * (y-p_bags)); % log-likelihood derivative with respect to vector of b coefficients
        end
    case 'log'
        %   Log-likelihood and its gradient for the Multiple Instance Learning
        %   paradigm, using softmax(x,y) = log(exp(x) + exp(y)).
        nInstPerBag = size(x,2);
        s1          = size(xReshaped,1);
        in_prd      = reshape(xReshaped*w,[s1/nInstPerBag,nInstPerBag]);
        p_inst      = 1./(1 + exp(-in_prd)); % bag and bag-instances probabilities (Noisy-OR)
        hardness    = 10;                    % the higher, the more precise the softmax approximation
        expsum      = sum(exp(hardness * p_inst), 2);
        p_bags      = log(expsum) ./ hardness;
        L           = -y' * log(p_bags) - (1-y)' * log(1-p_bags); % minus log-likelihood
        tc          = (y-p_bags) ./ (p_bags .* (1-p_bags) .* expsum);
        dL          = -(tc' * squeeze(sum(bsxfun(@times, x, ...
            exp(hardness*p_inst) .* p_inst .* (1-p_inst)), 2)))'; % log-likelihood derivative with respect w
    otherwise
        error('Mode should be one of ''nor'', ''max'' or ''log''')
end
