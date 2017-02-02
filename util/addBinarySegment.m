% Draw a line segment that connects points (x1,y1),(x2,y2) in a binary image.
% After that, a thinning or thickening operation can be performed optionally.
% 
%   output = addBinarySegment(input,x1,y1,x2,y2)
% 
%   INPUTS
%   input:       input image in which the linear segment is added
%   x1,x2,y1,y2: starting point and end point coordinates
%   oper:        'thin' or 'thicken'
%   val:         times thinning/thickening is performed. Can also be 'Inf'.
% 
% Stavros Tsogkas, <stavros.tsogkas@ecp.fr>
% Last update: April 2013

function output = addBinarySegment(input,x1,y1,x2,y2,oper,val)

if nargin<5
    error('Insufficient number of coordinates');
end
if nargin<6
    oper = 'thin'; 
    val  = 'inf'; 
end
if nargin<7 || val<0
    val = 1; 
end
if ~any(strcmp(oper,{'thin','thicken'}))
    error('oper argument can either be ''thin'' or ''thicken'' !')
end 

[x,y]        = bresenham(x1,y1,x2,y2);
invalid      = x < 1 | x > size(input,2) | y < 1 | y > size(input,1);
x(invalid)   = [];
y(invalid)   = [];
idx          = sub2ind(size(input),y,x);
segment      = false(size(input));
segment(idx) = true; 
segment      = bwmorph(segment,oper,val);
output       = input | segment;

