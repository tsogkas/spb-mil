function show_ridges_on_image(input_image,keyps,sc_k,selection,varargin)


if exist('sc_k'),
    try,
        show_ellipses_on_image(input_image,keyps,sc_k,selection,'color_ell',varargin{2});
    catch
        show_ellipses_on_image(input_image,keyps,sc_k,selection);
    end
else
    show_ellipses_on_image(input_image,keyps,[],[],'color_ell','b')
end
