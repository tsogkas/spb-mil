function [lli, lliDeriv] = loglikeNOR(b,y,x_in)
% function [lli, lliDeriv] = loglikeNOR(b,y,x_in)
% 
% Log-likelihood and its gradient for the Multiple Instance Learning 
% Paradigm, using noisy-or.
% 
% INPUTSs
% 	y: 	Nx1 colum vector of 0|1 class assignments.
% 	x: 	NxKxBAGINST array of feature vectors. BAGINST = NORIENT*NSCALES is
% 	the total number of instances per bag of features. Each bag contains
% 	the features for all orientations and scales for the corresponding
% 	pixel.
%   b:  linear classifier coefficients.
% 
% OUTPUT
%   lli: log likelihood value at point b
%   llgrad: log likelihood gradient with respect to b vector, calculated
%   at point b.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

x           = x_in{1};
x_reshaped  = x_in{2};
[N,K,nInst] = size(x); % get sizes
s1 = size(x_reshaped,1);
in_prd =  reshape(x_reshaped*b,[s1/nInst,nInst]);
% Calculate probabilites for bag instances and bags (Noisy-OR)
p_inst = 1./(1 + exp(-in_prd));
logp_bags = log(1-p_inst+eps);
logp_bags = sum(logp_bags,2);
p_bags = exp(logp_bags); % here p_bags = Prod(1-p_inst) for use of log1p later
% minus log-likelihood
lli = -y' * log1p(-p_bags) - (1-y)' * logp_bags;
% coefficients for efficient derivative calculation
tc = (y-(1-p_bags))./(1-p_bags);
ic = reshape(p_inst,[N 1 nInst]);
ic = repmat(ic,[1 K 1]);
% log-likelihood derivative with respect to the vector of b coefficients
lliDeriv = -(tc' * sum(x .* ic, 3))';
