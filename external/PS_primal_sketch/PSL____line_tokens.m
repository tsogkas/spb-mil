function [line_features,contours,components] = PSL____line_tokens(points,imsize,thresh_line_ener);
% [line_features,contour,components] = PSL____line_tokens(edge_points,imsize,thresh_line_ener)
%
% Converts a list of edge points into a set of straight line tokens.
% INPUT: 
%     points: scale & space maxima locations provided by PS1____scale_space_sketch
%     imsize: image dimensions
%     thresh_line_ener: a threshold to discard very low energy points.  
% OUTPUT:
%     line_features: a structure describing the extracted straight lines.
%     contours: a structure for the set of ridge/edge maxima points and
%               their attributes.
%     components: connected components of edge/ridge curves, before line
%               segmentation is applied.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

%%-------------------------------------------------------------------------
%% Part I: clean up the set of points delivered by the PS1* primal sketch code
%%-------------------------------------------------------------------------
[points,contours]   = PSLz1__preprocess_points(points,imsize);

%%-------------------------------------------------------------------------
%% Part II: find which points are neighbors and track connected components 
%%-------------------------------------------------------------------------

[successors,predecessors,current]           = PSLz2a_get_relatives(points,imsize);
relative_pr                                 = PSLz2b_process_relatives(predecessors,current);
relative_sc                                 = PSLz2b_process_relatives(successors,current);
components                                  = PSLz2c_find_connected_curves(current,relative_sc,relative_pr,imsize);


%%-------------------------------------------------------------------------
%% Part III: break components into straight line segments
%%-------------------------------------------------------------------------
line_parsing        =  PSLz3a_break_into_lines(components,imsize(1));
line_features       =  PSLz3b_embed_lines(line_parsing,thresh_line_ener);




