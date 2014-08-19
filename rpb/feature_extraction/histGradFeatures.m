function [dlc drc dlr scales] = histGradFeatures(im)
% function [ dlc drc dlr scales] = histGradFeatures(im)
% 
% Function that extracts the color/brighteness and texture features based
% on histogram differences.
%
% INPUTS
% im: the original image (RGB or grayscale).
% 
% OUTPUTS
% scales: vector that contains the scales of detection (in pixels)
% symmetryMap: a h x w x nscales array that provides a rough approximation of 
%              the detector output 
% dlc/drc/dlr: Arrays used to store the features of the chi-squared 
% histogram differences between the 3 rectangles we use to extract features
% using integral image representations. dlc: histogram difference between
% left and central rectangle (upper and central), drc: between right and
% central (bottom and central) and dlr between left and right (upper and
% bottom) rectangle areas. These arrays have dimensions
% [h,w,norient,nscales,nchannels], where h and w denote the height and
% width of the input image, norient is typically 8, nscales typically 13
% and nchannels 4 for RGB images and 2 for grayscale.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% SETUP

% Create scale vectors for CG and TG (different scales)
[h,w,~] = size(im); % sizes of the original image
scales = [4:2:14, 16:4:28, 32:8:48];
step = scales;
step(scales<=48) = 8; step(scales<=28) = 4; 
step(scales<=14) = 2; step(scales<=6) = 1;
stepidx = step;
stepidx(stepidx==4) = 3; stepidx(stepidx==8) = 4; 
thetas = (0:7)*pi/8; % orientation vector
nscales = length(scales);
norient = length(thetas);
smooth = 'savgol';
padoff = 3*7+2;  % maximum padding needed for gaussian pyramid   

%% Create gaussian pyramid of images

subsampled = cell(4,1);
subsampled{1} = im;
g = fspecial('gaussian', [3 3], 1); % gaussian filter
for isub = 2:4
    filtered = imfilter(subsampled{isub-1},g,'same','replicate');
    subsampled{isub} = filtered(1:2:end, 1:2:end, :);   % subsample
end
% Pad images to avoid wrong symmetry responces in larger scales
for isub = 1:4
    subsampled{isub} = padarray(subsampled{isub},[padoff padoff],'replicate');
end

%% Create texton map and Lab color map

ntextons = 64;              % number of texton labels 
nbins = [32*ones(3,1); ntextons];  % number of lab space labels
lmap = cell(4,1);
for isub = 1:4
    cmap = createLabColorMap(subsampled{isub},nbins(1));
    tmap = createTextonMap(subsampled{isub},ntextons);
    lmap{isub} = uint16(cat(3,cmap,tmap));
end
nchannels = size(lmap{1},3); % number of channels used

%% Rotate images outside the loops to reduce computational cost

rotsub = cell(norient,4);  % store rotated image versions
for ior = 1:norient 
    for isub = 1:4
        rotsub{ior,isub} = imrotate(lmap{isub},rad2deg(thetas(ior))); % rotated cmap and tmap
    end
end

%% Collect histogram gradient features

% Matrices used to store the brightness/color and texture features
dlc = zeros(h, w, norient, nscales, nchannels, 'single');
drc = zeros(h, w, norient, nscales, nchannels, 'single');
dlr = zeros(h, w, norient, nscales, nchannels, 'single');

fprintf(2, 'Extracting features at scales:\n');
fprintf(2, '[ ');    
for isc=1:nscales
    scale = scales(isc);
    fprintf(2, '%d ',scale);    
    scale = floor(scale/step(isc)); % adjust scale to subsampled size
    hMarg = 3*scale+2; 
    wMarg = 3*scale+1;          % adjust image size to speed up convolution
    ch = floor(h/step(isc));    % dimensions for original subsampled image
    cw = floor(w/step(isc));
    for ich=1:nchannels
        for ior=1:norient
            theta = thetas(ior);
            imrot = rotsub{ior,stepidx(isc)}(padoff+1-hMarg:end-padoff+hMarg,...
                padoff+1-wMarg:end-padoff+wMarg,ich);
            [tg_dlc,tg_drc,tg_dlr] =...
                histGrad(imrot,nbins(ich),scale,theta,[ch,cw],'smooth',smooth);

            % Resize images to original size
            tg_dlc = imresize(tg_dlc, [h,w]);
            tg_drc = imresize(tg_drc, [h,w]);
            tg_dlr = imresize(tg_dlr, [h,w]);

            % store results
            dlc(:,:,ior,isc,ich) = tg_dlc;
            drc(:,:,ior,isc,ich) = tg_drc;
            dlr(:,:,ior,isc,ich) = tg_dlr;
        end
    end
end
fprintf(2, ']\n');    

