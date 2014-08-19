function colored_sketch  =  colored_contours(input_image,ridge_skel,edge_skel);

ridge_skel = clip(1.5*ridge_skel,0,1);
edge_skel = clip(1.5*edge_skel,0,1);
sz=  3;
ridge_skel_u  = ordfilt2(ridge_skel,sz^2,ones(sz))/max(max(abs(ridge_skel(:))),eps);
edge_skel_u   = ordfilt2(edge_skel,sz^2,ones(sz))/max(max(abs(edge_skel(:))),eps);
fct =1.5;

sum_other_colors = edge_skel_u + ridge_skel_u;
colored_sketch(:,:,3) = max(input_image + ridge_skel_u - edge_skel_u,0);
colored_sketch(:,:,2) = max(input_image + edge_skel_u - ridge_skel_u,0);
colored_sketch(:,:,1) = max(input_image - sum_other_colors,0);