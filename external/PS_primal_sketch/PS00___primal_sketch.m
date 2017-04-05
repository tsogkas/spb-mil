function [ridge_feats,edge_feats,blob_feats,contours,conn_components,scale_space_ims,maxima_points] = PS00___primal_sketch(input_image,settings_sketch,settings_tests);
% [ridge_feats,edge_feats,blob_feats,contours,conn_components,scale_space_ims] = PS00___primal_sketch(input_image,settings_sketch,settings_tests)
%
% Gateway routine for primal sketch computation 
% INPUT:
%   input_image:  grayscale image, normalized to lie in [0,1]
%   settings_sketch, settings_tests (optional): user-defined settings,
%               overriding defaults in PS0z{1,2}_settings_{sketch,tokens}
%
% OUTPUT: 
%   {ridge,edge,blob}_feats: structures containing coordinates and
%              attributes (scale, orientation, etc) of primal sketch tokens
%   contours: cell array containing locations of ridge/edge maxima points
%   conn_components: connected components of ridge and edge contours
%   scale_space_ims: cubes containing ridge/edge/blob/intensity in space x scale 
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

%%------------------------------------------------------------------------
%% Get overall settings for primal sketch extraction
%%------------------------------------------------------------------------
if nargin==1, settings_sketch  =[]; settings_tests = []; end
settings_sketch = PS0z1__settings_sketch(input_image,settings_sketch);
settings_tokens = PS0z2__settings_tokens(input_image,settings_tests,settings_sketch);
imsize = size(input_image);

%%------------------------------------------------------------------------
%% Core primal sketch: form scale space, get maxima of 
%% edge, ridge & blob differential operators
%%------------------------------------------------------------------------
[maxima_points,scale_space_ims]                               = PS1____scale_space_sketch(input_image,settings_sketch);
[blob_feats,ridge_feats,edge_feats,contours,conn_components]  = PS2____post_process(maxima_points,imsize,settings_tokens);
