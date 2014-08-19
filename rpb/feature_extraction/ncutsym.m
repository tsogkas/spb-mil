function [sPb,eigVec] = ncutsym(psym,nvec,radius,sigmaW,sigmaOE)
% function [sPb,eigVec] = ncutsym(psym,nvec,radius,sigmaW,sigmaOE)
% 
% Calculate spectral feature. Code based on the code for Berkeley global
% boundary detector (gpb).
% 
% INPUTS
% psym:     probability map for symmetry (extracted from 1st stage detector).
% nvec:     desirable number of eigenvectors.
% radius:   radius defining the area around pixel for affinity calculation.
% sigmaW:   sigma for affinity matrix W.
% sigmaOE:  sigma for oriented smooth Gaussian gradient.
% 
% OUTPUTS
% sPb:      spectral feature at 8 orientations.
% eigVec:   eigenvector images.


%% Calculate affinity matrix

if nargin<5, sigmaOE = 1; end
if nargin<4, sigmaW = 0.1; end
if nargin<3, radius = 5; end
if nargin<2, nvec = 10; end

[h,w,~] = size(psym); 
l{1} = zeros(size(psym, 1) + 1, size(psym, 2));
l{1}(2:end, :) = psym;
l{2} = zeros(size(psym, 1), size(psym, 2) + 1);
l{2}(:, 2:end) = psym;
% build the pairwise affinity matrix for symmetry map
[val,I,J] = buildWexp(l{1},l{2},radius,sigmaW);
Wpsym = sparse(val,I,J);
W = Wpsym;

%% Calculate and visualize eigenvectors

[wx, wy] = size(W);
x = 1 : wx;
S = full(sum(W, 1));
D = sparse(x, x, S, wx, wy);
clear S x;

opts.issym=1;
opts.isreal = 1;
opts.disp=2;
[EigVect, EVal] = eigs(D - W, D, nvec, 'sm',opts);
clear D W opts;

EigVal = diag(EVal);
clear Eval;
EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

disp('This is the Ncut eigenvectors...');
vect = zeros(h, w, nvec);
for iter = 2 : nvec
    vect(:, :, iter) = imresize(reshape(EigVect(:, iter), [w h])',[h,w]);
    vect(:,:,iter)=(vect(:,:,iter)-min(min(vect(:,:,iter))))/(max(max(vect(:,:,iter)))-min(min(vect(:,:,iter))));
end
clear EigVect;

%% Extract spectral symmetry feature

% OE parameters

hil = 0;
deriv = 1;
support = 3;
norient = 8;
dtheta = pi/norient;
% ch_per = [4 3 2 1 8 7 6 5];
ch_per = [2 3 4 5 6 7 8 1];     % change theta sequence to [0,pi) clockwise

eigVec = zeros(h,w,nvec-1);
sPb = zeros(h, w, norient);
for v = 1 : nvec
    if EigVal(v) > 0,
        eigVec = vect(:,:,v)/sqrt(EigVal(v));
        for o = 1 : norient,
            theta = dtheta*o;
            f = oeFilter(sigmaOE, support, theta, deriv, hil);
            sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, eigVec));
        end
    end
end


