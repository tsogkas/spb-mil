function [beta,p_bags,lli] = logist2mil(y,x)
% [beta,p,lli] = logist2mil(y,x) 
%
% 2-class logistic regression with multiple instance learning.  
%
% INPUT
% 	y 	Nx1 colum vector of 0|1 class assignments.
% 	x 	NxKxNINST array of feature vectors. NINST = NORIENT*NSCALES is
% 	the total number of instances per bag of features. Each bag contains
% 	the features for all orientations and scales for the corresponding
% 	pixel.
% 	[w]	Nx1 vector of sample weights.
%
% OUTPUT
% 	beta 	Kx1 column vector of model coefficients
% 	p 	Nx1 column vector of fitted class 1 posteriors
% 	lli 	log likelihood
%
% Class 1 posterior is 1 / (1 + exp(-x*beta))
%


% initial guess for beta: all zeros
beta_0 = zeros(k,1);
%calculate beta coefficients by maximizing log-likelihood
[beta,lliVec] = minimize(beta_0,'loglikeMIL',1000,y,x);
lli = lliVec(end);

% calculate probabilites for bag instances and bags (Noisy-OR)
p_inst = zeros(N,nInst);
p_bags = ones(N,1);
for i=1:nInst
    p_inst(:,i) = 1 ./ (1 + exp(-x(:,:,i)*beta));
    p_bags = p_bags .* (1-p_inst(:,i)');
end
p_bags = 1-p_bags;
