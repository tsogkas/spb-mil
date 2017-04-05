function settings_tokens_out = PS0z2__settings_tokens(input_image,settings_tokens_in,settings_sketch);
% settings_tokens_out = PS0z2__settings_tokens(input_image,settings_tokens_in,settings_sketch)
%
% Sets default fields for token extraction & selection.
% If settings_sketch_in is not empty, defaults are overriden by
% the settings_sketch_in fields.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

%%-------------------------------------------------------------------
%% 1. Setup default token thresholds 
%%-------------------------------------------------------------------
[size_m,size_n] = size(input_image);
thresh_line_ener     = .02;   %% a very small threshold on the energy of a ridge/edge point, only to speed up processing

%% thresholds for accepting a ridge token:
%% 1: threshold on averaged energy of differential operator
ener_ridge_threshold  = .13;
%% 2: and threshold on `merit' used in Lowe's '87 cuvre partitioning algorithm
merit_ridge_threshold = 4;

%% same for edges
ener_edge_threshold   = .13;
merit_edge_threshold  = 4;

%% thresholds for accepting a blob token:
%% 1: threshold on energy
ener_blob_threshold   = .12;
%% 2: upper and lower thresholds on the curvature-related 
%% quantity used in Lowe's '04 paper (a high value indicates edges, and should
%% lead to  rejection)
curv_up_threshold     = 10;
curv_down_threshold   = 0;
%% 3: threshold on distance formed by adding blob scale and location 
%% distances used to label a pair of blobs as overlapping
threshold_overlapping = 2;

%%------------------------------------------------------------------
%% 2. override default settings 
%%------------------------------------------------------------------
if ~isempty(settings_tokens_in),
    structure=  settings_tokens_in; expand_structure;
end

%%-------------------------------------------------------------------
%% 3. form string expressions for the feature acceptance conditions
%%-------------------------------------------------------------------
ener_ridge  = sprintf('ener>%.3f',  ener_ridge_threshold);
ener_edge   = sprintf('ener>%.3f',  ener_edge_threshold);
ener_blob   = sprintf('ener>%.3f',  ener_blob_threshold);
merit_ridge = sprintf('merit_geom>%.3f', merit_ridge_threshold);
merit_edge  = sprintf('merit_geom>%.3f', merit_edge_threshold);
curv_up     = sprintf('curv<%.3f',curv_up_threshold);
curv_down   = sprintf('curv>%.3f',curv_down_threshold);

scale_cond  = sprintf('scales>%.3f',settings_sketch.scales_wt(3));

fields_wt= {'ener_ridge','ener_edge','ener_blob','merit_ridge','merit_edge',...
            'curv_up','curv_down','thresh_line_ener','scale_cond','threshold_overlapping',...
            };
compress_structure;
settings_tokens_out = structure;


