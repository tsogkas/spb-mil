% Script for checking correctness of cost derivative.

nSamples=100000;
nDim=13;
nInstPerBag=96;
% x = single(rand(nSamples,nDim,nInstPerBag));
x = rand(nSamples,nInstPerBag,nDim);
y = zeros(nSamples,1);
idx = randsample(nSamples,0.2*nSamples);
y(idx) = 1;
% normalize features to unit variance
x           = double(x); % double precision is necessary for training     
xReshaped   = reshape(x,[nSamples*nInstPerBag,nDim]);
xstd        = std(xReshaped); xstd = xstd + (xstd==0);  % feature standard deviation
xReshaped   = bsxfun(@rdivide,xReshaped,xstd);
x           = reshape(xReshaped,[nSamples,nInstPerBag,nDim]); 
x_input{1}  = x;         % save initial feature matrix for use in minimize function
x_input{2}  = xReshaped; % save reshaped and normalized feature matrix
beta0 = zeros(nDim,1);
beta1=rand(nDim,1);
beta2=rand(nDim,1);
% profile on
checkgrad('loglikeLOG',beta1,1e-5,y,x,xReshaped)
% profile viewer
