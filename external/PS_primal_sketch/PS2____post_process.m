function [blob_feats,ridge_feats,edge_feats,contours,conn_components] = PS2___post_process(maxima_points,imsize,settings_tokens);
% [blob_feats,ridge_feats,edge_feats,contours,component_strings] = PS2___post_process(maxima_points,imsize,settings_tokens)
%
% Extracts tokens from the primal sketch maxima by coverting curves into straight
% line segments and then discarding non-salient tokens.
% 
% INPUTS:
%    maxima_points: set of points provided by PS1____scale_space_sketch
%    imsize: size of input image
%    settings_tokens:  setting structure for token extraction (thresholds, etc.)
% OUTPUTS:
%   {ridge,edge,blob}_feats: structures containing coordinates and
%              attributes (scale, orientation, etc) of primal sketch tokens
%   contours: cell array containing locations of ridge/edge maxima points
%   conn_components: connected components of ridge and edge contours
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

structure  = settings_tokens; expand_structure;

%%------------------------------------------------------------------------
%% Get straight edge/ridge line tokens, skeletons, and blobs
%%------------------------------------------------------------------------ 
[ridge_feats,contours{1},conn_components{1}] = PSL____line_tokens(maxima_points{1},imsize,thresh_line_ener);
[edge_feats, contours{2},conn_components{2}] = PSL____line_tokens(maxima_points{2},imsize,thresh_line_ener);
blob_feats  = PSL____point_tokens(maxima_points{3},imsize);

%%------------------------------------------------------------------------
%% Apply thresholds to pick a subset of the tokens
%%------------------------------------------------------------------------ 
ridge_feats = PS2z1__choose_tokens(ridge_feats,{ener_ridge,merit_ridge});
edge_feats  = PS2z1__choose_tokens(edge_feats,{ener_edge,merit_edge});
blob_feats  = PS2z1__choose_tokens(blob_feats,{ener_blob,curv_up,curv_down,scale_cond});
blob_feats  = PS2z2__discard_overlapping(blob_feats,threshold_overlapping);


