function [maxima_points,edges,skeletons,edge_pb] = PS00_get_berkeley_edges(edge_pb,thresh_line_ener);

if nargin==1,
    thresh_line_ener = .1;
end

im_dims = size(edge_pb);
idxs= find(edge_pb>0);
orient = get_orientation_discrete(edge_pb);
%maxima_points.

maxima_points.indexes = idxs;
maxima_points.ener   = edge_pb(idxs);
maxima_points.theta  = orient(idxs)+pi/2;
maxima_points.scl    = 2*ones(size(idxs));
maxima_points.scind  = 2*ones(size(idxs));

[edges, skeletons] = PSL____line_tokens(maxima_points,im_dims,thresh_line_ener);


function res = get_orientation_discrete(edgemap);
[gr_x,gr_y] = meshgrid([-6:6],[-6:6]);
ind=  0;
for orient= 0:pi/8:2*pi,
    ind= ind+1;
    dist_perp  = cos(orient)* gr_x + sin(orient)* gr_y;
    dist_al = -sin(orient)* gr_x + cos(orient)* gr_y;  
    kernel  = exp(-(pow_2(dist_perp)/10 + pow_2(dist_al)/3));
    kernel  = kernel/sum(kernel(:));
    filtered(:,:,ind) = abs(filter2(kernel,double(edgemap)));
end
[mx,idx] = max(filtered,[],3);

or_out = 0;
ind = 0;
for orient = 0:pi/32:2*pi,
    ind = ind+1;
    or_out = or_out + (idx==ind)*orient;
end
res = or_out;