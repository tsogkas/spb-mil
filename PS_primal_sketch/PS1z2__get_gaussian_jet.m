function [gauss_jet] = PS1z2__get_gaussian_jet(input_image,scale_gauss,use_iir);
% [gauss_jet] = PS1z2__get_gaussian_jet(input_image,scale_gauss,use_iir)
%
% Gets up to second-order derivatives of gaussian at scale_gauss.
% The results may vary depending on whether the iir filter is
% used or not, due to the different behavior of the fir/iir filters on 
% the boundaries. 
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

if use_iir,
    [im_sc,d_x,d_xx,d_y,d_yy,d_xy] = iir_gauss(input_image,scale_gauss);
else
    [im_sc,d_x,d_xx,d_y,d_yy,d_xy] = fir_gauss(input_image,scale_gauss);
end
fields_wt = {'im_sc','d_x','d_y','d_xx','d_xy','d_yy'};
compress_structure;
gauss_jet = structure;
    