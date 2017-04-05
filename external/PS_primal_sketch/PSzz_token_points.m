function  [crd_x,crd_y,energy,scale]  = PSzz_token_points(sketch,k);
% [crd_x,crd_y,energy,scale]  = PSzz_token_points(sketch,k);
% 
% Returns the raw sketch information related to each token extracted by
% the PS00___primal_sketch routine
% INPUTS sketch : an edge/ridge structure (returned by PS00___primal_sketch)
%        k      : index into the specific token of interest
% OUTPUTS: crd_x, crd_y: coordinates of points constituting token
%          energy: normalized magnitude of the differential operator
%          scale:  scale at which the feature is detected
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

crd_x = sketch.lines{k}(2,:);
crd_y = sketch.lines{k}(1,:);
energy = sketch.lines{k}(4,:);
scale  = sketch.lines{k}(5,:);
