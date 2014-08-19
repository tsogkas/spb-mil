function output = add_segment(input,x1,y1,x2,y2,oper,val)
% function idx = add_segment(x1,y1,x2,y2)
% 
% Draw a line segment that connects points (x1,y1),(x2,y2) in pixel coordinates of
% a 2D image with dimensions h x w. x1 and x2 are the row coordinates and 
% y1,y2 are the column coordinates. After that, a thinning or thickening
% operation can be performed optionally.
% 
% INPUTS
% 
% input: input image in which the linear segment is added
% x1,x2,y1,y2: starting point and end point coordinates
% oper: 'thin' or 'thicken'
% val: times thinning/thickening is performed. Can also be Inf

if nargin<5, error('Insufficient number of coordinates\n');end
if nargin<6, oper = 'thin'; end
if nargin<7 || val<0, val = 0; end
if ~strcmp(oper,'thin')&&~strcmp(oper,'thicken')
    error('oper can either be ''thin'' or ''thicken''!\n')
end 


[h,w] = size(input);
[yy,xx] = meshgrid(1:w,1:h);
x1 = round(x1); x2 = round(x2);
y1 = round(y1); y2 = round(y2);
if x1==x2
    if y2>y1, segment = (yy<=y2)&(yy>=y1)&(xx==x1); 
    else segment = (yy>=y2)&(yy<=y1)&(xx==x1); 
    end
elseif y1==y2
    if x1>x2, segment = (yy==y1)&(xx<=x1)&(xx>=x2);
    else segment = (yy==y1)&(xx>=x1)&(xx<=x2);
    end
else
    a = (y2-y1)/(x2-x1);
    b = y1 - a*x1;
    if x2>x1, segment = (yy<=round(a*(xx+sign(a))+b))&(yy>=round(a*(xx-sign(a))+b))&(xx>=x1)&(xx<=x2);
    else segment = (yy<=round(a*(xx+sign(a))+b))&(yy>=round(a*(xx-sign(a))+b))&(xx>=x2)&(xx<=x1);
    end
end

segment = bwmorph(segment,oper,val);
output = input + segment;

