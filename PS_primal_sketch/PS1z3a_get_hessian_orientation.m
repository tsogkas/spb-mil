function [cos_hessian,sin_hessian,L_pp,L_qq]   = PS1z3a_get_hessian_orientation(d_xx,d_yy,d_xy);
% [cos_hessian,sin_hessian,L_pp,L_qq]   = PS1z3a_get_hessian_orientation(d_xx,d_yy,d_xy)
%
% Calculate eigenvalues and eigenvectors of the image hessian for orientation estimation.
% The eigenvector corresponding to the maximal eigenvalue gives the orientation 
% along which dark features (valleys) lie. Ridges lie along the perpendicular
% orientation.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

sq_common_term = sqrt(pow_2(d_xx - d_yy) + 4*pow_2(d_xy));
%% eigenvalues 
L_pp = 1/2*(d_xx + d_yy - sq_common_term);
L_qq = 1/2*(d_xx + d_yy + sq_common_term);

%% First eigenvector (corresponding to L_pp):
coeff_1 =  -(1/2*d_yy-1/2*d_xx) + 1/2*sq_common_term;           
coeff_2 = d_xy;
ener =  max(sqrt(pow_2(coeff_1) + pow_2(coeff_2)),eps);
nener = ener==eps; 

%% for the special case where `ener' is too low, instead of 
%% normalizing the eignevector elements, set them to  [1,0]

cos_hessian = (coeff_1./ener).*(1-nener) + nener.*1;
sin_hessian = (coeff_2./ener).*(1-nener) + nener.*0;

function submolic_verify();
%% symbolic toolbox code for a quick  verification
a_11= sym('d_xx'); a_12  = sym('d_xy'); a_21= a_12; a_22 = sym('d_yy');
[vec,val] = eig([a_11,a_12;a_21,a_22]);
vec, val