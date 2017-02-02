function   [ismax,ismin] = PS1z3b_maxalong_orientation(derivstr,cos_eig,sin_eig,L_1,L_2);
% [ismax,ismin] = PS1z3b_maxalong_orientation(derivstr,cos_eig,sin_eig,L_1,L_2)
%
% Finds ridges and edges of image.
% INPUT: 
%   derivstr: structure containing image derivatives up to 2nd order
%   cos_eig,sin_eig: elements of the largest eigenvector of the Hessian
%    L_1, L_2: eigenvalues of the Hessian. 
% OUTPUT: 
%   ismax: indicator function of valleys
%   ismin: indicator function of ridges
% 
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

d_x = derivstr.d_x; d_y = derivstr.d_y;
for it =[0,1],
    im_an = sin_eig;
    if it ==0,
        %% it =0:  bright features (ridges)-> find zero-crossings of derivative
        %% along the orientation perpendicular to [cos_hessian,sin_hessian],
        %% and require that second derivative is negative (maximum condition)
        cos_ang = - sin_eig;
        sin_ang =   cos_eig;
        cond = (abs(L_1)>abs(L_2)).*(L_1<=0) + (abs(L_1)<abs(L_2)).*(L_2<=0); %*(L_1<L_2);
    else
        %% it =1:  dark features (valleys)-> take derivative along
        %% along  [cos_hessian,sin_hessian],
        %% and require that second derivative is positive (minimum condition)
        cos_ang = cos_eig;
        sin_ang = sin_eig;
        cond =  (abs(L_1)<abs(L_2)).*(L_2>=0) + (abs(L_1)>abs(L_2)).*(L_1>=0); %(max(L_1,L_2)>=0); %.*(L_1>L_2);
    end
    
    %% locate the image locations where sin_ang changes sign 
    %% due to the orientation ambiguity (look below for details)
    [zc_hs,zc_vs] = find_orientation_ambiguity(im_an);

    %% form directional derivative 
    im_dr = (cos_ang.*d_x + sin_ang.*d_y);
    
    %% find vertical/horizontal zero crossings of im_dr
    %% by fetching the values of neighboring points 
    [shift_dr_up,shift_dr_left] =  shift_matrix_up_left(im_dr);

    %% for locations with {-1/1} sin_ang changes,
    %% invert derivative sing for neighboring points 
    shift_dr_up   = shift_dr_up.*(1- 2*zc_vs);
    shift_dr_left = shift_dr_left.*(1- 2*zc_hs);
    
    %% find locations where the derivative changes sign
    zc_v  = shift_dr_up.*im_dr<0;
    zc_h  = shift_dr_left.*im_dr<0;

    zc_locations  = (zc_h|zc_v);
    if it==0,
        ismax = (zc_locations).*cond;
    else
        ismin = (zc_locations).*cond;
    end
end

function   [zc_hs,zc_vs] = find_orientation_ambiguity(im_an);
%% a small patch that fixes the problem of orientation ambiguity:
%% Observe that around -1/1 sin_eig may change abruptly, even though the orientation
%% does not really change 
%% ( [-1,0] and [1,0] encode the same orientation, but with different direction)
%% For this we directly find locations where this happens, by identifying
%% zero crossings of sin_eig, and keeping those where sin_eig has large
%% magnitude. The rest (sin_eig~=0) are natural.
%% zc_vs, zc_hs thus indicate whether a vertical/horizontal zero crossing
%% of sin_eig is an artifact of this ambiguity.

[shift_an_up,shift_an_left] =  shift_matrix_up_left(im_an);
zc_vs = ((shift_an_up.*im_an)<0);
zc_vs = zc_vs&(abs(im_an)>.8);
zc_hs = ((shift_an_left.*im_an)<0);
zc_hs = zc_hs&(abs(im_an)>.8);

function [shift_up,shift_left] = shift_matrix_up_left(inp);
[sz_m,sz_n]  = size(inp);
shift_up    =  [[inp(2:end,1:end,:)];zeros(1,sz_n)];
shift_left  =  [[inp(1:end,2:end,:),zeros(sz_m,1)]];

