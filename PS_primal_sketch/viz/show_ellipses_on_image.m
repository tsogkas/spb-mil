function show_ellipses_on_image(input_image,keyps,sc_k,selection,varargin)
color_ell = 'k';
color_arr = [1,1/3,0];
lw = 2.5;
selection_given = exist('selection');
if ~isempty(input_image),
    imshow(input_image)
end
if ~isempty(varargin),
    expand_varargin;
end
no_arrow = 0;
try,
    if ~all(color_arr==[1,1,1]),
        error,
    end
    no_arrow=  1;
end
structure = keyps; expand_structure;

if ~exist('selection')|isempty(selection)
    selection = 1:length(c_m);
end

if ~exist('sc_k')|(isempty(sc_k))
    sc_k  = 1;
end
c_m  =round(c_m); c_n = round(c_n);
for k_ind= 1:length(selection),
    k   = selection(k_ind);
    orn = orientations(k);
    sc  = scales(k);

    try,
        ratio = ratios(k);
    catch
        ratio =1;
    end
   
    rotation = [cos(orn),-sin(orn);sin(orn),cos(orn)];
    eigs = sc*[1,0;0,ratio];
    center = [c_n(k);c_m(k)];
    hold on;
    h = my_draw_ellipse(center,rotation,diag(eigs));
    set(h,'color',color_ell);
    set(h,'linewidth',lw);
    hold on;
    if ~no_arrow,
        sc_u =sc_k*sc;
        h =draw_arrow([c_n(k);c_m(k)],sc_u,orn,color_arr);
        set(h,'linewidth',lw);
    end
end
