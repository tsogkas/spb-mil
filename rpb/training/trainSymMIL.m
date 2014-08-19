function  [beta, xstd, p_bags, lli, y, x] =...
                trainSymMIL(iids,groundPath,featPath,numberOfSamples,...
                    samplingMethod,featConfig,buffer,boundPerc)
% function  [beta, xstd, p_bags, lli, y, x] =...
%                 trainSymMIL(iids,groundPath,featPath,numberOfSamples,...
%                     samplingMethod,featConfig,buffer,boundPerc)
% 
% Function that trains the symmetry classifier using Multiple Instance Learning
% 
% INPUTS
% 
% iids:             the training images iids vector
% groundPath:      the directory path where the groundtruth ridge maps are stored
% featPath:         the directory path where the features used for training are stored
% numberOfSamples:  the total number of samples used in training
% samplingMethod:   sampling method ('random','balanced' or 'boundary').
% featConfig:       determines which features are used for training
% boundPerc:        the percentage of boundary pixels used for training (per image)
% buffer:           determines the area around groundtruth ridge pixels, from where
%                   we do not choose training samples
% 
% OUTPUTS
% 
% beta:     The weights of the detector.
% xstd:     Standard deviation of the vector weights.
% p_bags:   The posterior probabilities calculated as an output by the
%           symmetry detector.
% lli:      The final value of the log-likelihood after the optimization.
% x:        The N x K x nInst array, carrying the features for all bags corresponding
%           to each one of the training samples. N is the number of samples, K is the
%           length of the feature vector and nInst is the number of instances
%           included in each bag.
% y:        The 0/1 labels for each sample.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% Sample features and labels (optional retraining)
fprintf(2,'  Training...\n');
[x,y] = sampleDetectorMIL(iids,groundPath,featPath,numberOfSamples,...
            samplingMethod,featConfig,boundPerc,buffer);

%% Normalize features and calculate inputs for minimization function
[N,K,nInst] = size(x);
% make sure features are double precision (necessary for training)
x = double(x);
% normalize features to unit variance
xReshaped = reshape(permute(x,[1,3,2]),[N*nInst,K]); % x: NxK
xstd = std(xReshaped);
xstd = xstd + (xstd==0);
xReshaped = xReshaped ./ repmat(xstd,size(xReshaped,1),1);
x = permute(reshape(xReshaped,[N,nInst,K]),[1,3,2]); 
x_input{1} = x; % save initial feature matrix for use in minimize function
x_input{2} = xReshaped; % save reshaped and normalized feature matrix

%% Calculate optimal beta coefficients
fprintf(2,'Fitting model...\n');
beta_0 = rand(K,1);
[beta,lliVec] = minimize(beta_0,'loglikeNOR',100,y,x_input);
lli = lliVec(end);

%% Calculate probabilites for bag instances and bags (Noisy-OR)
fprintf(2,'Calculating posterior probabilities...\n');
s1 = size(xReshaped,1);
in_prd =  reshape(xReshaped*beta,[s1/nInst,nInst]);
p_inst = 1./(1 + exp(-in_prd));
logp_bags = log(1-p_inst+eps);
logp_bags = sum(logp_bags,2);
p_bags = 1-exp(logp_bags);

%% De-normalize beta coefficients
xstd = xstd';
beta = beta./xstd;

