function  points = PS1z4__add_new_points(points,gauss_jet,indexes,feat_tp,scale_gauss,feat_strength_detection_p,sc_ind);
% points = PS1z4__add_new_points(points,gauss_jet,indexes,feat_tp,scale_gauss,feat_strength_detection_p,sc_ind)
%
% Utility function that concatenates previous set of maxima points with the
% ones detected at the current scale.
% 
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

d_xx_wt = gauss_jet.d_xx(indexes);    
d_yy_wt = gauss_jet.d_yy(indexes);    
d_xy_wt = gauss_jet.d_xy(indexes);
d_x_wt =  gauss_jet.d_x(indexes);  
d_y_wt =  gauss_jet.d_y(indexes);

theta =  get_orientation_at_point(d_x_wt,d_y_wt,d_xx_wt,d_yy_wt,d_xy_wt,feat_tp);
trace_hessian_sq     =  pow_2(d_xx_wt + d_yy_wt);
determinant_hessian  = d_xx_wt.*d_yy_wt - pow_2(d_xy_wt);
ener  =  feat_strength_detection_p{feat_tp}(indexes);
scl   =  scale_gauss*ones(length(indexes),1);
scind =  sc_ind*ones(length(indexes),1);


points.theta     = [points.theta;theta];
points.indexes   = [points.indexes;indexes];
points.ener      = [points.ener;ener];
points.scl       = [points.scl;scl];
points.tr_hess   = [points.tr_hess;trace_hessian_sq];
points.det_hess  = [points.det_hess;determinant_hessian];
points.scind     = [points.scind;scind];

%%--------------------------------------------------------------------------------------------
%% Internal functions
%%--------------------------------------------------------------------------------------------
function theta =  get_orientation_at_point(d_x_wt,d_y_wt,d_xx_wt,d_yy_wt,d_xy_wt,feat_tp);
if  feat_tp==1,
    %% for ridges the orientation is determined by the eigenvectors of the
    %% Hessian matrix (local structure is valey/peak)
    common_term = ((d_xx_wt-d_yy_wt).^2 + 4*(d_xy_wt).^2);
    sq_com_term = sqrt(common_term);
    cos_b = sqrt(1/2*( 1 + (d_xx_wt - d_yy_wt)./max(sq_com_term,eps)));
    sin_b = my_sign(d_xy_wt).*sqrt(1/2*( 1 - (d_xx_wt - d_yy_wt)./max(sq_com_term,eps)));
    L_pp_nn = 1/2*(d_xx_wt + d_yy_wt - sq_com_term);
    L_qq_nn = 1/2*(d_xx_wt + d_yy_wt + sq_com_term);
    choose = abs(L_pp_nn)>abs(L_qq_nn);
    theta  = atan2(sin_b,cos_b) + pi/2.*choose;
else
    %% otherwise orientation is determined by the x/y derivatives
    sqren = sqrt(max(pow_2(d_x_wt) + pow_2(d_y_wt),eps));
    cos_b = d_x_wt./sqren;   sin_b = d_y_wt./sqren;
    theta =  atan2(sin_b,cos_b);
end




            