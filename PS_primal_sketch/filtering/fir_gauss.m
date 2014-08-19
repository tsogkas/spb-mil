function [im_sc,d_x,d_xx,d_y,d_yy,d_xy] = fir_gauss(L_in,scale);
%% code for convolution with separable Gaussian and DoG filters.
%% Slow, but controllable alternative to iir filtering.
%% 
%% Iasonas Kokkinos <jkokkin@stat.ucla.edu>

scale_m =  ceil(4*scale);
L_t = wrap_e(L_in,scale_m);
domain = [-scale_m:scale_m];

gauss       = 1/(sqrt(2*pi)*scale)*(exp( - domain.^2/(2*scale^2)));
gauss_d     = -(gauss.*(-domain/(scale^2)));
gauss_dd    = ((-1/scale^2).*gauss + (domain.^2/scale^4).*gauss);

filt_v_smooth  =  filter1(gauss,L_t,1,scale_m);
filt_h_smooth  =  filter1(gauss,L_t,2,scale_m);

im_sc   =  peel(filter1(gauss,filt_h_smooth,1,scale_m),scale_m);
d_x     =  peel(filter1(gauss_d,filt_v_smooth,2,scale_m),scale_m);
d_y     =  peel(filter1(gauss_d,filt_h_smooth,1,scale_m),scale_m);
d_xx    =  peel(filter1(gauss_dd  ,filt_v_smooth,2,scale_m),scale_m);
d_yy    =  peel(filter1(gauss_dd  ,filt_h_smooth,1,scale_m),scale_m);
d_xy    =  peel(filter1(gauss_d,filter1(gauss_d,L_t,1,scale_m),2,scale_m),scale_m);

