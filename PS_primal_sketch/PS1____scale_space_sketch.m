function [points,extras_out_ss] = PS1____scale_space_sketch(input_image,settings_sketch);
% [points,extras_out_ss] = PS1____scale_space_sketch(input_image,settings_sketch)
%
% Receives as input a gray scale image in [0,1] and the settings for primal
% sketch extraction.
% Returns a cell array with maxima points for ridge/edge/blob strengths
% and, if requested, features from the image scale space.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

extras_out_ss = [];
structure = settings_sketch; expand_structure;
[D_xx,D_yy,D_xy,D_x,D_y,idxs,points,feat_strength_p] =  PS1z1__initialize(input_image);


%% loop over scales
for sc_ind = 1:length(scales_wt),
    scale_gauss = scales_wt(sc_ind);
    [gauss_jet] = PS1z2__get_gaussian_jet(input_image,scale_gauss,use_iir);

    for feat_tp = [1:3],                
        [feat_strength_selection,feat_strength_detection,zero_crossing]  = ...
            PS1z3__get_feature_strength(gauss_jet,scale_gauss,feat_tp,settings_sketch);
                           
        %%------------------------------------------------------------------------------
        %% Find feature strength maxima in location and scale. 
        %% NOTE: For memory efficiency, 
        %% we only check whether the maxima at the previous iteration  
        %% (`_p' suffix)  are stronger than the current ones and the ones
        %% before them ('_pp' suffix). That is, there is a 1-iteration delay
        %% in detecting the maxima. 
        %%------------------------------------------------------------------------------

        maxima_indexes =[];
        if feat_tp ==1, sc_min = 3; else sc_min = 2; end  % for ridges it seems best to ignore the finest resolution
        if sc_ind>=sc_min,
            max_str{feat_tp}  = min((feat_strength_p{feat_tp}>feat_strength_selection),(feat_strength_p{feat_tp}>feat_strength_pp{feat_tp}));
            maxima_indexes = find(max_str{feat_tp}.*zc_p{feat_tp});
        end
        if ~isempty(maxima_indexes),
            points{feat_tp} = PS1z4__add_new_points(points{feat_tp},gauss_jet_prev,maxima_indexes,feat_tp,scale_gauss,feat_strength_detection_p,sc_ind);
        end
        
        feat_strength_pp{feat_tp}          = feat_strength_p{feat_tp}; 
        feat_strength_p{feat_tp}           = feat_strength_selection;
        feat_strength_detection_p{feat_tp} = feat_strength_detection;
        zc_p{feat_tp} = zero_crossing;

        %% if the user has asked for it, pack the feature images in a cube.
        if keep_scale_space_im(feat_tp),
            extras_out_ss{feat_tp}(:,:,sc_ind) =  feat_strength_detection;
        end            
    end
    gauss_jet_prev = gauss_jet;    
    if keep_scale_space_im(4),
        extras_out_ss{4}(:,:,sc_ind) =  gauss_jet.im_sc;
    end
end