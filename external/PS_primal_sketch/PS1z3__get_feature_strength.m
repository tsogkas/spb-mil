function [feat_strength_selection,feat_saliency,zc]  = PS1z3__get_feature_strength(gauss_jet,scale_gauss,feat_tp,settings_sketch);
% [feat_strength_selection,feat_saliency,sc_min,zc]  = PS1z3__get_feature_strength(gauss_jet,scale_gauss,feat_tp,settings_sketch)
%
% Given the gaussian jet of an image calculates the strength of different
% features. 
% 
% INPUT:
%   gauss_jet: gaussian jet of the image at scale_gauss
%   feat_tp: feature_type   (1=ridge/2=edge/3=blob) 
%   settings_sketch: settings for sketch (Lindeberg's gamma values etc.)
% OUTPUT:
%   feat_strength_selection:  strength  for scale selection 
%   feat_saliency:            strength  for saliency estimation 
%   zc: zero-crossing / maxima-in-space locations (zc)
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

structure =gauss_jet; expand_structure;
t = scale_gauss.^2;
switch feat_tp,
    case 1,
        %% ridges
        gamma_n_ridge = settings_sketch.gamma_n_ridge;
        %% ridge strength measure used: A (Eq. 51 in Lindeberg)

        differential_operator    =  (pow_2(d_xx - d_yy) + 4*pow_2(d_xy));
        feat_strength_selection  =  t^(2*gamma_n_ridge)*differential_operator;
        feat_saliency            =  2*sqrt(t^2*differential_operator);
    case 2,
        %% edges        
        differential_operator = (pow_2(d_x) + pow_2(d_y));
        feat_strength_selection = t^(settings_sketch.gamma_n_edge).*differential_operator;
        %% as in lindeberg's paper, the saliency measure is obtained at
        %% gamma_n_edge =1 
        feat_saliency = sqrt(t*(2*pi)*differential_operator);        
    case 3,
        %% blob
        gamma_blob  = 1;
        differential_operator = abs(d_xx + d_yy);
        feat_strength_selection = t^(gamma_blob).*differential_operator;
        feat_saliency = feat_strength_selection*3/2;
end

%% get condition for a zero crossing
switch feat_tp,
    case 1, 
        %% ridges
        %% check if image reaches maximum along the ridge's orientation
        [cos_eig,sin_eig,L_1,L_2] =  PS1z3a_get_hessian_orientation(d_xx,d_yy,d_xy);    
        [ismax,ismin]             =  PS1z3b_maxalong_orientation(gauss_jet,cos_eig,sin_eig,L_1,L_2);
        zc                        =  max(ismax,ismin);   
        zc([end-1,end],:) = 0;
    case 2,   
        %% edges
        %% zero crossing of D_vv (Eq. 9 in Lindeberg's paper)
        D_vv  =  pow_2(d_x).*d_xx + 2.*d_x.*d_y.*d_xy + pow_2(d_y).*d_yy;;
        value_right  =  D_vv(:,[2:end,end]);
        value_up     =  D_vv([2:end,2],:);
        zc    = max(((D_vv.*value_right)<0),((D_vv.*value_up)<0));
        zc([end-1,end],:) = 0;
    case {3}   
        %% blobs     
        %% estimate the maximum of the (8+1)-pt neighborhood around each
        %% pixel and check whether it is attained at that pixel
        val_max = max(feat_strength_selection,feat_strength_selection(:,[2:end,end]));
        val_max = max(val_max,feat_strength_selection(:,[1,1:end-1]));
        val_max = max(val_max,feat_strength_selection([1,1:end-1],:));
        val_max = max(val_max,feat_strength_selection([1,1:end-1],[1,1:end-1]));
        val_max = max(val_max,feat_strength_selection([1,1:end-1],[2:end,end]));
        val_max = max(val_max,feat_strength_selection([2:end,end],:));
        val_max = max(val_max,feat_strength_selection([2:end,end],[1,1:end-1]));
        val_max = max(val_max,feat_strength_selection([2:end,end],[2:end,end]));
        zc = (val_max ==feat_strength_selection);
end

