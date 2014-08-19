function [points,contours]   = PSLz1__preprocess_points(points,imsize);
% [points,contours]   = PSLz1__preprocess_points(points,imsize)
% 
% Takes the set of points delivered by the primal sketch code, and
% (a) removes weaker points lying on the same location with other, stronger points 
% (b) processes  the set of contour locations (morphological thinning) and
%     removes unnecessary points 
% (c) returns the matrix indexes of contour locations.
% 
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

points                           = PSLz1a_isolate_same_location(points,imsize(1));
[points,skel,ener,scales,orient] = PSLz1b_feature_matrices(points,imsize);
contours                         = PSzz_zip_contour(skel,(skel.*orient),(skel.*scales),(skel.*ener));
