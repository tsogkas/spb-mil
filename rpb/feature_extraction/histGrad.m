function [tg_dlc tg_drc tg_dlr] = histGrad(im,nbins,scale,theta,hw,varargin)
% function [tg_dlc tg_drc tg_dlr] = histGrad(im,nbins,scale,theta,hw,varargin)
%
% Compute histogram gradients at a single orientation and scale using the
% integral image representation for faster computation
%
% INPUT
% im:       Pre-rotated image.
% nbins:    Number of bins.
% scale:	Half of rectangle side used to calculate im.
% theta:	Orientation orthogonal to tg.
% 'smooth'	Smoothing method, one of 
%			{'gaussian','savgol','none'}, default 'none'.
% 'ratio'   ratio of rectangle sides.
% 'sigma'	Sigma for smoothing, default to scale.
%
% OUTPUT
%   tg_dlc  Histogram difference between the left and central rectangle.
%   tg_drc  Histogram difference between the right and central rectangle.
%   tg_dlr  Histogram difference between the left and right rectangle.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% Process options
smooth = 'none';
scale = floor(max(1,scale));
theta = mod(theta,pi);  % ensure that theta is in the interval [0,pi]
sigma = scale;
for i = 1:2:numel(varargin),
  opt = varargin{i};
  if ~ischar(opt), error('option names not a string'); end
  if i==numel(varargin), error('option ''%s'' has no value',opt); end
  val = varargin{i+1};
  switch opt,
   case 'smooth',
    switch val,
     case {'none','gaussian','savgol'}, smooth=val;
     otherwise, error('invalid option smooth=''%s''',val);
    end
   case 'sigma', sigma=val;
   otherwise, error('invalid option ''%s''',opt);
  end
end

%% Calculate histogram differences of rotated integral images

[tg_dlc, tg_drc, tg_dlr] = histgrad(im,nbins,scale); % use c++ mex file
tg_d = 0.5*cat(3,tg_dlc,tg_drc,tg_dlr);

%% Re-rotate and crop images 

tg_d = imrotate(tg_d,-rad2deg(theta)); 
% crop images after rotation to initial image size
[uh, uw,~] = size(tg_d); % get size of uncropped image
horMarginRot = max(0,floor((uw-hw(2))/2)); verMarginRot = max(0,floor((uh-hw(1))/2));
tg_d = tg_d(verMarginRot+1:end-verMarginRot,horMarginRot+1:end-horMarginRot,:);

%% Post-filtering (Savitsky-Golay)

switch smooth,
 case 'gaussian',
  f = oeFilter([sigma .5],3,theta+pi/2);
  tg_d(:,:,1) = applyFilter(f,tg_d(:,:,1));
  tg_d(:,:,2) = applyFilter(f,tg_d(:,:,2));
  tg_d(:,:,3) = applyFilter(f,tg_d(:,:,3));
 case 'savgol',
  a = fitparab2(tg_d(:,:,1),sigma,sigma/4,theta);
  tg_d(:,:,1) = max(0,a);
  a = fitparab2(tg_d(:,:,2),sigma,sigma/4,theta);
  tg_d(:,:,2) = max(0,a);
  a = fitparab2(tg_d(:,:,3),sigma,sigma/4,theta);
  tg_d(:,:,3) = max(0,a);
end

tg_dlc = tg_d(:,:,1);
tg_drc = tg_d(:,:,2);
tg_dlr = tg_d(:,:,3);

