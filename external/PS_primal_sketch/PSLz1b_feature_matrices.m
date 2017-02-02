function [points,skel,ener,scales,orient] = PSLz2_construct_matrices(points,imsize);
% [pointsskel,ener,scales,orient] = PSLz2_construct_matrices(points,imsize)
%
% Some image processing code to refine the set of points considered for
% subsequent processing
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

zs = zeros(imsize);
ind = zs;
ener = zs;
orient = zs;
orient_fine = zs;
scales = zs;
indexes = points.indexes;
% 
ener(indexes) =points.ener;
scales(indexes) = points.scl;
orient(indexes) = points.theta+pi/2;

%% perform nonmaximum suppression on the energy of the line points
%% to obtain a clean feature indicator function
ind = nonmax(ener,orient)>0;

%% some subsequent morphological processing to get a map of 
%% contiguous line segments
boundary_mask = my_patch(peel(ones(size(ener)),1),1,0);
skel = bwmorph(bwmorph(bwmorph(ind,'clean'),'fill'),'thin','inf');
skel = bwmorph(bwmorph(skel.*boundary_mask,'clean'),'bridge');

%% keep only those points that reside on the refined set of locations
points  = keep_points(points,find(skel(points.indexes)));