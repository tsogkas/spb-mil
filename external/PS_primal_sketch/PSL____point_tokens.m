function  [points] =point_features(points_kept,imsize);
% [points] =point_features(points_kept,imsize)
%
% Converts point structure returned by PS1___scale_space_sketch
% into point structure amenable to subsequent processing, display, etc.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007
 

[c_m,c_n]=  ind2sub(imsize,points_kept.indexes);
reject =find(points_kept.det_hess<0);
curv = points_kept.tr_hess./max(points_kept.det_hess,.000000001);

points.c_m = c_m';
points.c_n = c_n';
points.ratios = ones(size(points.c_m));
points.scales = points_kept.scl';
points.orientations = points_kept.theta';
points.ener = points_kept.ener';
points.curv = curv';
%% will be used later to reject these points 
%% (cuvr>0 is an acceptance condition)
points.curv(reject) = -1;
