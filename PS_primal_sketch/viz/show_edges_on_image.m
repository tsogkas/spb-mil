function show_ridges_on_image(input_image,keyps,sc_k,selection,varargin)
color_ell = [1,0,0];

color_arr = 'g';
lw = 3;
if ~isempty(input_image),
    imshow(input_image)
end
if ~isempty(varargin),
    expand_varargin;
end
structure = keyps; expand_structure;

sz = 10;
pt_x_along = [-sz:sz]/(sz);
pt_y_along_1 = zeros(size(pt_x_along));
pt_y_along_2 = ones(size(pt_x_along));
pt_y_along_3 = -ones(size(pt_x_along));

pt_y_across = [-sz:sz]/(sz);
pt_x_across_1 = ones(size(pt_y_across));
pt_x_across_2 = -ones(size(pt_y_across));

if (nargin <4)|isempty(selection)
    selection = 1:length(c_m);
end
try sc_k = sc_k; catch sc_k  = 1; end
c_m  =round(c_m); c_n = round(c_n);
for k_ind= 1:length(selection),
    k  = selection(k_ind);
    orn = orientations(k);
    sc =  scales(k);
    try,
        ratio = ratios(k);
    catch
        ratio =1;
    end
    orn  = orn;
    pts_x_along = (pt_x_along)*sc;
    pts_y_along_1 = (pt_y_along_1)*ratio*sc;
    pts_y_along_2 = (pt_y_along_2)*ratio*sc;
    pts_y_along_3 = (pt_y_along_3)*ratio*sc;
    pts_y_across = (pt_y_across)*ratio*sc;
    pts_x_across_1 = (pt_x_across_1)*sc;
    pts_x_across_2 = (pt_x_across_2)*sc;
    
    rotation = [cos(orn),-sin(orn);sin(orn),cos(orn)];
    center = [c_n(k);c_m(k)];
    hold on;
    pts_xs_along_1 =  center(1) + rotation(1,1)*pts_x_along +  rotation(1,2)*pts_y_along_1;
    pts_xs_along_2 =  center(1) + rotation(1,1)*pts_x_along +  rotation(1,2)*pts_y_along_2;
    pts_xs_along_3 =  center(1) + rotation(1,1)*pts_x_along +  rotation(1,2)*pts_y_along_3;
    
    pts_ys_along_1 =  center(2) + rotation(2,1)*pts_x_along +  rotation(2,2)*pts_y_along_1;
    pts_ys_along_2 =  center(2) + rotation(2,1)*pts_x_along +  rotation(2,2)*pts_y_along_2;
    pts_ys_along_3 =  center(2) + rotation(2,1)*pts_x_along +  rotation(2,2)*pts_y_along_3;

   
    pts_xs_across_1 =  center(1) + rotation(1,1)*pts_x_across_1 +  rotation(1,2)*pts_y_across;
    pts_xs_across_2 =  center(1) + rotation(1,1)*pts_x_across_2 +  rotation(1,2)*pts_y_across;
    pts_ys_across_1 =  center(2) + rotation(2,1)*pts_x_across_1 +  rotation(2,2)*pts_y_across;
    pts_ys_across_2 =  center(2) + rotation(2,1)*pts_x_across_2 +  rotation(2,2)*pts_y_across;
    
    h_i1=  line(pts_xs_along_1,pts_ys_along_1);

%    h_o1=  line(pts_xs_across_1,pts_ys_across_1);
%    h_o2=  line(pts_xs_across_2,pts_ys_across_2);
    
    
    set(h_i1,'color',color_arr,'linewidth',lw);
    %set(h_i2,'color',color_ell,'linewidth',lw);
    %set(h_i3,'color',color_ell,'linewidth',lw);
%    set(h_o1,'color',color_ell,'linewidth',lw);
%    set(h_o2,'color',color_ell,'linewidth',lw);   
    hold on;
%    h = text(c_n(k),c_m(k),sprintf('%i',selection(k_ind)));
    %set(h,'linewidth',3);
%    h = text(c_n(k),c_m(1,k),sprintf('%i',idxs(k)));
end
