% Convert image to CIE LAB color space
% 
%   lab = rgb2lab(rgb)
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last update: January 2013

function lab = rgb2lab(rgb)

lab = applycform(rgb, makecform('srgb2lab'));