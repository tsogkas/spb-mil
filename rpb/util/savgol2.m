function [c] = savgol2(z,d,k,ra,rb,theta)
% function [c] = savgol(z,d,k,ra,rb,theta)
%
% Directional 2D Savitsky-Golay filtering with elliptical support.
% The computation is done with a convolution, so the boundary of the
% output will be biased.  The boundary is of size floor(max(ra,rb)).
%
% INPUT
%	z	Values to fit.
%	d	Degree of fit, usually 2 or 4.
%	k	Coefficient to return in [1,d+1], 1 for smoothing.
%	ra,rb	Radius of elliptical neighborhood, ra=major axis.
%	theta	Orientation of fit (i.e. of minor axis).
%
% OUTPUT
%	c	Coefficient of fit.
%
% Stavros Tsogkas <statsogk@gmail.com>
% December 2011

if d<0, error('d is invalid'); end
if k<1 || k>d+1, error('k is invalid'); end

ra = max(1.5,ra);
rb = max(1.5,rb);
ira2 = 1 / ra^2;
irb2 = 1 / rb^2;
wr = floor(max(ra,rb));
wd = 2*wr+1;
sint = sin(theta);
cost = cos(theta);

% 1. compute linear filters for coefficients
% (a) compute inverse of least-squares problem matrix
[u v] = meshgrid(-wr:wr,-wr:wr);
ai = -u*sint + v*cost; % distance along major axis
bi = u*cost + v*sint; % distance along minor axis
temp = ai(((ai.*ai)*ira2 + (bi.*bi)*irb2)<=1)'; % inside support
temp = [ones(size(temp)); temp(ones(2*d,1),:)];
xx = sum(cumprod(temp,1),2);

A = zeros(d+1,d+1);
for i = 1:d+1, 
  A(:,i) = xx(i:i+d); 
end
A = inv(A);

% (b) solve least-squares problem for delta function at each pixel
temp = ai(:)';
temp = [ones(size(temp)) ; temp(ones(d,1),:)];
yy = cumprod(temp,1);
filt = (A*yy)';
filt(((ai.*ai)*ira2 + (bi.*bi)*irb2)>1,:) = 0; % outside support pixels 
filt = reshape(filt,wd,wd,d+1);

% 2. apply the filter to get the fit coefficient at each pixel
c = conv2(z,filt(:,:,k),'same');

