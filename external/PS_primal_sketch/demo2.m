%% demo.2
%% Shows how to modify the standard code for primal sketch computation
%% by overriding default choices made inside the code.
%% Defaults are in PS0z1_settings_sketch and PS0z2_settings_tokens
%% 
%% Iasonas Kokkinos <jkokkin@stat.ucla.edu>

input_image = double(imread('horse010.jpg'))/256;

%% override defaults for primal sketch extraction 
override_sketch.nscales             = 30;  
override_sketch.linear_spacing      = 0;     % require geometrical spacing
override_sketch.min_scale           = .3;

%% Request the cubes formed by the intensity and feature-strength images 
%% as the scale increases. Binary entries, corresponding to (ridge, edge, blob, intensity)
override_sketch.keep_scale_space_im = [1,0,0,1];

%% override defaults for criteria used to select token subset
override_tests.ener_blob_threshold = .2;
override_tests.ener_edge_threshold = .1;

input_image = double(imread('horse010.jpg'))/256;
[ridges,edges,blobs,contours,conn_components,scale_space_ims] = PS00___primal_sketch(input_image,override_sketch,override_tests);
[ridge_markers] = PSzz_components_to_markers(conn_components{1},input_image);

contour_ridge       = PSzz_unzip_contour(contours{1}); 
contour_edge        = PSzz_unzip_contour(contours{2});
colored_sketch      =  colored_contours(input_image,contour_ridge,contour_edge);
%% show sketch tokens
figure,
subplot(2,2,1); show_ridges_on_image(input_image,ridges); title('Ridges')
subplot(2,2,2); show_edges_on_image(input_image,edges);   title('Edges');
subplot(2,2,3); show_ellipses_on_image(input_image,blobs);title('Blobs')
subplot(2,2,4); imshow(colored_sketch); title('Ridge and Edge contours');