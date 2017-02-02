function [a,b,c] = fitparab2(z,ra,rb,theta)
% function [a,b,c] = fitparab(z,ra,rb,theta)
%
% Fit cylindrical parabolas to elliptical patches of z at each
% pixel. 
%
% INPUT
%	z	Values to fit.
%	ra,rb	Radius of elliptical neighborhood, ra=major axis.
%	theta	Orientation of fit (i.e. of minor axis).
%
% OUTPUT
%	a,b,c	Coefficients of fit: a + bx + cx^2
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003
% Modified by Stavros Tsogkas November 2011

ra = max(1.5,ra);
rb = max(1.5,rb);

% compute the interior quickly with convolutions
a = savgol2(z,2,1,ra,rb,theta);
if nargout>1, b = savgol2(z,2,2,ra,rb,theta); end
if nargout>2, c = savgol2(z,2,3,ra,rb,theta); end
