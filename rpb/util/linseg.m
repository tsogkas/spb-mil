function idx = linseg(size,points)
% function idx = linseg(size,start,end)
% 
% Return linear indices of a line segment connecting two points in a 
% discrete grid of size size(1) by size(2)
% 
% INPUTS
% 
% size:     vector containing sizes of the grid.
% points:   [x1,y1; x2,y2], where (x1,y1) is the start point and (x2,y2) is
%           the endpoint


h = size(1); w = size(2);
[yy,xx] = meshgrid(1:w,1:h);
pstart = points(1,:);
pend = points(2,:);
x1 = pstart(2); y1 = pstart(1);
x2 = pend(2); y2 = pend(1);
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


idx = find(bwmorph(segment,'thin',inf));
