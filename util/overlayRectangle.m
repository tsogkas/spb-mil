% Plot rotated rectangle on gray or color image.
% 
%   out = overlayRectangle(im,pos,scale,theta,ratio,color)
% 
%   im: input image.
%   pos: [row, col] coordinates of rectangle center in the original image.
%   scale: a/2 where a is the small side of the rectangle.
%   theta: rotation angle in degrees.
%   ratio: ration between the two sides of the rectangle {b = 3*a}.
%   color: color to use for plotting ('y','r','b','m', etc.)
% 
%   out: original image with drawn rectangle
% 
%   Also check overlay_binary_image.m
% 
% Stavros Tsogkas ,stavros.tsogkas[at]ecp.fr
% Last update: June 2013


function out = overlayRectangle(im,pos,scale,theta,ratio,color)

if nargin<4, theta = 0;   end
if nargin<5, ratio = 3;   end
if nargin<6, color = 'y'; end

[height, width, ~] = size(im);
if isscalar(pos) || iscolumn(pos)
    [row,col] = ind2sub([height,width],pos);
    center    = [row, col];
elseif isrow(pos);
    center = pos;
else
    error('Center position should be either a scalar or a vector')
end


rect = zeros(2*scale+3, 2*ratio*scale+3); 
rect(1:end-1,[1,end-1]) = 1;
rect([1,end-1],1:end-1) = 1;
rotrect     = bwmorph(imrotate(rect,theta),'thin','inf');
rectCenter  = round(size(rotrect)/2);
[rows,cols] = find(rotrect);
distances   = [rows-rectCenter(1), cols-rectCenter(2)]; % find boarder distances from rectangle center
rectim      = zeros(height,width);
for i=1:size(center,1)
    rectim(center(i,1),center(i,2)) = 1;
    rows              = center(i,1) + distances(:,1);
    cols              = center(i,2) + distances(:,2);
    outOfBounds       = rows<1 | rows>height | cols<1 | cols>width;  
    rows(outOfBounds) = [];                 % limit rectangle to image boarder
    cols(outOfBounds) = [];
    ind               = sub2ind([height,width],rows,cols);
    rectim(ind)       = 1;
end
out = overlay_binary_image(im,rectim,color);
