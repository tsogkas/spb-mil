function features = computeHistogramFeatures(im)
% Extracts color,brighteness and texture features based on chi-square
% histogram differences.
% 
%   [features] = computeHistogramFeatures(im,fine,b,usechi2,smooth,ratio)
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: February 2017


opts.fine    = true;       % compute features at finer scale
opts.usechi2 = [1 1 1 1];  % use chi-square distance or other
opts.smooth  = 'savgol';   % use savgol filtering 
opts.ratio   = 2;          % ratio between rectangle filter sides
opts         = parseVarargin(opts, varargin);

scales  = [6:2:14, 16:4:28, 32:8:48]; if opts.fine, scales = [4, scales]; end
thetas  = (0:7)*pi/8;
nScales = length(scales);
nOrient = length(thetas);
pad     = 3*7+2;  % maximum padding needed for gaussian pyramid   
[step, pyramidLevel] = assignPyramidLevel(scales);
impyramid = gaussianPyramid(im,pad);

% Create texton map and Lab color map
nTextons = 64;                     % number of texton labels 
nBins = [32*ones(3,1); nTextons];  % number of lab space labels
lmap = cell(4,1);
for isub = 1:4
    cmap = computeColorMap(impyramid{isub},nBins);
    tmap = computeTextonMap(impyramid{isub},nTextons);
    lmap{isub} = uint16(cat(3,cmap,tmap));
end
nChannels = size(lmap{1},3);       % number of channels used
csim = [];

% Collect histogram gradient features
[h,w,~] = size(im); % sizes of the original image
dlc = zeros(h, w, nOrient, nScales, nChannels, 'single');
drc = zeros(h, w, nOrient, nScales, nChannels, 'single');
dlr = zeros(h, w, nOrient, nScales, nChannels, 'single');
for s=1:nScales
    scale   = floor(scales(s)/step(s)); % adjust scale to subsampled size
    hMargin = 3*scale+2;
    wMargin = 3*scale+1;             % adjust image size to speed up convolution
    h0      = floor(h/step(s));    % dimensions for original subsampled image
    w0      = floor(w/step(s));
    for o=1:nOrient
        imrot   = imrotate(lmap{pyramidLevel(s)},rad2deg(thetas(o)));
        imrot   = imrot(pad+1-hMargin:end-pad+hMargin,pad+1-wMargin:end-pad+wMargin,:);
        for c=1:nChannels
            hgrad = computeHistogramGradient(imrot(:,:,c),nBins(c),...
                scale,thetas(o),[h0,w0],opts.usechi2(c),opts.smooth,csim,ratio);            
            if size(hgrad,1)~=h || size(hgrad,2)~=w
                hgrad = imresize(hgrad, [h,w],'bilinear');
            end
            dlc(:,:,o,s,c) = hgrad(:,:,1);
            drc(:,:,o,s,c) = hgrad(:,:,2);
            dlr(:,:,o,s,c) = hgrad(:,:,3);
        end
    end
end
features.scales = scales;
features.thetas = thetas;
features.step   = step;
features.level  = pyramidLevel;
features.dlc    = dlc;
features.drc    = drc;
features.dlr    = dlr;

% --- Compute histogram gradients at a single orientation and scale using
% --- integral image representation for faster computation
function hgrad = computeHistogramGradient(im,nBins,scale,theta,hw,usechi2,smooth,csim,ratio)

% Process options
if nargin<6, usechi2 = true;     end
if nargin<7, smooth  = 'savgol'; end
if nargin<8, csim    = [];       end
if nargin<9, ratio   = 3;        end
scale = floor(max(1,scale));
theta = mod(theta,pi);  % ensure that theta is in the interval [0,pi]

% Calculate histogram differences of rotated integral images
rectArea = (2*scale+1)*(2*ratio*scale+1);
[dlc,drc,dlr] = histgrad(im,nBins,scale,ratio); % use c++ mex file
hgrad      = (0.5/rectArea)*cat(3,dlc,drc,dlr);
% Re-rotate and crop images 
hgrad      = imrotate(hgrad,-rad2deg(theta)); 
% crop images after rotation to initial image size
[uh, uw,~] = size(hgrad); % get size of uncropped image
horMargin  = max(0,floor((uw-hw(2))/2)); 
verMargin  = max(0,floor((uh-hw(1))/2));
hgrad      = hgrad(verMargin+1:end-verMargin,horMargin+1:end-horMargin,:);

% Post-filtering (Savitsky-Golay)
sigma = scale;
switch smooth,
    case 'savgol'
        hgrad(:,:,1) = max(0,fitparab2(hgrad(:,:,1),sigma,sigma,theta));
        hgrad(:,:,2) = max(0,fitparab2(hgrad(:,:,2),sigma,sigma,theta));
        hgrad(:,:,3) = max(0,fitparab2(hgrad(:,:,3),sigma,sigma,theta));
    case 'none'
    otherwise
        error('Filtering method not supported!')
end

% --- Create CIE LAB-space colormap for input image (RGB or grayscale) ----
function cmap = computeColorMap(im,nbins)
% im:       original image
% nbins:    number of bin labels

if nargin<2, nbins = ones(3,1)*32; end
if numel(nbins)==1, nbins = ones(3,1)*nbins; end

if ismatrix(im), % grayscale image
    cmap = max(1,ceil(im*nbins(1)));
else % RGB image
    % convert gamma-corrected image to LAB and scale values into [0,1]
    % min and max values for a,b channels of LAB
    % used to scale values into the unit interval
    abmin = -73;
    abmax = 95;
    lab = applycform(im, makecform('srgb2lab'));
    lab(:,:,1) = lab(:,:,1) ./ 100;
    lab(:,:,2) = (lab(:,:,2) - abmin) ./ (abmax-abmin);
    lab(:,:,3) = (lab(:,:,3) - abmin) ./ (abmax-abmin);
    lab(:,:,2) = max(0,min(1,lab(:,:,2)));
    lab(:,:,3) = max(0,min(1,lab(:,:,3)));
        
    cmap = zeros(size(im,1),size(im,2),3);
    for i=1:3
        cmap(:,:,i) = max(1,ceil(lab(:,:,i)*nbins(i))); 
    end
end

% --- Create texton map of an image (Berkeley Pb code) --------------------
function tmap = computeTextonMap(im,k)
% im:   original image
% k:    number of texton labels

if nargin<2, k = 64; end

no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
fname = sprintf('unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
textonData = load(fname); % defines fb,tex,tsim
if size(im,3)==3, tmapim = rgb2gray(im); else tmapim = im; end
tmap = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);

% --- Get subsampling rates and pyramid level indexes ---------------------
function [step, pyramidLevel] = assignPyramidLevel(scales)
% step: step between rectangle scales in the same level of the pyramid
% pyramidLevel: index of pyramid level according to scale

step = scales;
step(scales<=6)                       = 1;
step(scales>=8  & scales<=14)         = 2; 
step(scales>=16 & scales<=28)         = 4; 
step(scales>=32 & scales<=48)         = 8; 
pyramidLevel(scales<=6)               = 1;
pyramidLevel(scales>=8  & scales<=14) = 2;
pyramidLevel(scales>=16 & scales<=28) = 3; 
pyramidLevel(scales>=32 & scales<=48) = 4; 

% --- Create gaussian pyramid of input image ------------------------------
function impyr = gaussianPyramid(im,pad)
im       = im2double(im);
impyr    = cell(4,1);
impyr{1} = im;
for isub = 2:4
    impyr{isub} = imresize(impyr{isub-1},0.5,'bilinear');
end

% Pad images to avoid wrong symmetry responses at larger scales
for isub = 1:4
    impyr{isub} = padarray(impyr{isub},[pad pad],'replicate');
end
