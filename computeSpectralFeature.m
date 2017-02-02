function [sPb,eigVec] = computeSpectralFeature(psym,nvec,radius,sigmaW,sigmaOE)
% function [sPb,eigVec] = computeSpectralFeature(psym,nvec,radius,sigmaW,sigmaOE)
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
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last update: August 2013

% Calculate affinity matrix
if nargin<5, sigmaOE = 1;   end
if nargin<4, sigmaW  = 0.1; end
if nargin<3, radius  = 5;   end
if nargin<2, nvec    = 10;  end

[height,width,~] = size(psym); 
l{1}             = zeros(size(psym, 1) + 1, size(psym, 2));
l{2}             = zeros(size(psym, 1), size(psym, 2) + 1);
l{1}(2:end, :)   = psym;
l{2}(:, 2:end)   = psym;
[val,I,J]        = buildWexp(l{1},l{2},radius,sigmaW); % pairwise affinity matrix for symmetry map
Wpsym            = sparse(val,I,J);
W                = Wpsym;

% Calculate (and visualize) eigenvectors
[wx, wy] = size(W);
x        = 1 : wx;
S        = full(sum(W, 1));
D        = sparse(x, x, S, wx, wy);
clear S x;
opts.issym        = 1;
opts.isreal       = 1;
opts.disp         = 2;
[EigVect, EVal]   = eigs(D - W, D, nvec, 'sm',opts); clear D W opts;
EigVal            = diag(EVal); clear Eval;
EigVal(1:end)     = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);
disp('This is the Ncut eigenvectors...');
vect = zeros(height, width, nvec);
for i = 2 : nvec
%     subplot(1,nvec,i-1);    
%     figure; clf;
    vect(:,:,i) = imresize(reshape(EigVect(:, i), [width height])',[height,width]);
    vect(:,:,i) = (vect(:,:,i)-min(min(vect(:,:,i))))/(max(max(vect(:,:,i)))-min(min(vect(:,:,i))));
%     imagesc(vect(:,:,i)); axis image; axis off;
end
clear EigVect;

% Extract spectral symmetry feature
% OE parameters
hil     = 0;
deriv   = 1;
support = 3;
% sigmaOE = 1;
norient = 8;
dtheta  = pi/norient;
% ch_per = [4 3 2 1 8 7 6 5];    % original theta sequence for gpb
ch_per  = [2 3 4 5 6 7 8 1];     % change theta sequence to [0,pi) clockwise

eigVec = zeros(height,width,nvec-1);
sPb    = zeros(height, width, norient);
for v = 1 : nvec
    if EigVal(v) > 0,
        eigVec(:,:,v) = vect(:,:,v)/sqrt(EigVal(v));
%         fig = figure; imagesc(eigVec);
%         print(fig, '-depsc2', sprintf('eigenvec_%d.eps',v))
        for o = 1 : norient,
            theta = dtheta*o;
            f = oeFilter(sigmaOE, support, theta, deriv, hil);
            sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, eigVec(:,:,v)));
        end
    end
end

