function settings_sketch_out = PS0z1__settings_sketch(input_image,settings_sketch_in);
% settings_sketch_in = PS0z1__settings_sketch(input_image,settings_sketch_in)
%
% Sets default fields for primal sketch computation.
% If settings_sketch_in is not empty, defaults are overriden by
% settings_sketch_in fields.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

%%-------------------------------------------------------------------
%% settings for primal sketch extraction 
%%-------------------------------------------------------------------

%%------------------------------------------------------------------
%% 1. Setup default choices for primal sketch scale space
%%------------------------------------------------------------------

%% starting & ending scales, and scale spacing
min_scale = max(min(size(input_image))/200,.5);
max_scale = min(size(input_image))/4;
nscales   = 50;
linear_spacing = 1;

%% determines whether the dll implementing Deriche's iir gaussian filter &
%% derivatives is used
use_iir  = 1;

%% scale normalization for ridges (check Lindeberg's paper)
norm_r = 'A';
gamma_n_ridge = 3/4;
gamma_n_edge  = 1/2;

%% determines which of the scale space images we will be returning
%% code is : [ridge,edge,blob,image]
keep_scale_space_im = [0,0,0,0];

%%------------------------------------------------------------------
%% 2. override defaults settings 
%%------------------------------------------------------------------
if ~isempty(settings_sketch_in),
    structure = settings_sketch_in; expand_structure;
end

if linear_spacing,
    scales_wt = linspace(min_scale,max_scale,nscales);
else
    scales_wt = logspace(log10(min_scale),log10(max_scale),nscales);
end
    
fields_wt = {'norm_r','min_scale','max_scale','gamma_n_ridge','gamma_n_edge','scales_wt','use_iir','keep_scale_space_im'};
compress_structure;
settings_sketch_out = structure;