function show_intermediate_curves_on_image(input_image,tokens,color)
if nargin==2,
    color = 'r';
end
if ~isempty(input_image),
    imshow(input_image)
end
is_component = isfield(tokens,'imsize');
if is_component,
    tokens= keep_points(tokens,find(tokens.lst>5));    
    nlines = length(tokens.lst);
else
    nlines = length(tokens.lines);
end

for k_ind = 1:nlines,
    hold on,
    if is_component,
        [crd_x,crd_y]  = PSzz_conn_comp_points(tokens,k_ind);
    else
        [crd_x,crd_y]  = PSzz_token_points(tokens,k_ind);
    end
    plot(crd_x,crd_y,color,'linewidth',2);
    hold on,
    scatter(crd_x(1),crd_y(1),16,'k','filled');
    hold on,
    scatter(crd_x(end),crd_y(end),16,'k','filled');
end
