function [res,ind]  = pick_ps_point(input_image,ps_points);
if ~isempty(input_image),
end
[x,y,z] = impixel;
dists  = (x - ps_points.c_n).^2  +(y-ps_points.c_m).^2;
[cl,ind]  = min(dists);

res = keep_points(ps_points,ind);
show_ellipses_on_image([],res,[],[],'color_ell','k');