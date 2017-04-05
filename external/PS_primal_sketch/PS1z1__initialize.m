function [D_xx,D_yy,D_xy,D_x,D_y,idxs,points,feat_strength_p] = PS1z1__initialize(input_image);
% [D_xx,D_yy,D_xy,D_x,D_y,idxs,points,feat_strength_p] =  PS1z1__initialize(input_image)
%
% Initialization for fields used in primal sketch computation.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

for tp = [1:3],
    D_xx{tp} = [];    D_yy{tp} = [];      D_xy{tp} = [];
    D_y{tp}  = [];    D_x{tp}  = [];      idxs{tp} = [];

    points{tp}.indexes = [];    points{tp}.ener  =[];
    points{tp}.theta = [];      points{tp}.scl = [];
    points{tp}.det_hess = [];   points{tp}.tr_hess  = [];
    points{tp}.scind = [];
    feat_strength_p{tp} = zeros(size(input_image));
end