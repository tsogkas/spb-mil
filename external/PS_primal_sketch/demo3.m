%% demo3.m
%% A demonstration of the intermediate
%% results used for primal sketch computation
%% Shows: scale space construction, maxima locations, curve segments,
%% marker image generation.
%%
%% Iasonas Kokkinos <jkokkin@stat.ucla.edu>

input_image = double(imread('horse010.jpg'))/256;

%% Request the cubes formed by the intensity and feature-strength images
%% as the scale increases. Binary entries, corresponding to (ridge, edge, blob, intensity)
override_sketch.keep_scale_space_im = [1,1,0,1];
override_tests  = [];

[ridges,edges,blobs,contours,connected_components,scale_space_ims] =...
    PS00___primal_sketch(input_image,override_sketch,override_tests);

for it =[1,2],  % iterate (ridges /edges)
    if it==1,
        tokens = ridges;
        location_maxima     = PSzz_unzip_contour(contours{1});        color = 'b';
        colored_sketch      =  colored_contours(input_image,location_maxima,zeros(size(input_image)));
    else
        tokens = edges;
        location_maxima     = PSzz_unzip_contour(contours{2});          color = 'g';
        colored_sketch      =  colored_contours(input_image,zeros(size(input_image)),location_maxima);
    end

    %%  Visualize feature strength-versus-scale as a cube.
    figure,show_3d_volume(scale_space_ims{it});
    
    %%  Show the 3 processing stages involved in token computation

    figure,
    subplot(2,2,1);
    imshow(colored_sketch);
    title('Maxima Locations.');

    subplot(2,2,2);
    show_intermediate_curves_on_image(input_image,connected_components{it},color);
    title('Larger than 5 pixel-long Connected Components.');

    subplot(2,2,3);
    show_intermediate_curves_on_image(input_image,tokens,color);
    title('Curve Segments, after thresholding.');

    subplot(2,2,4);
    if it==1,
        show_ridges_on_image(input_image,tokens);
    else
        show_edges_on_image(input_image,tokens);
    end
    title('Straight tokens.');
end

%% turn the connected components of the ridge 
%% contours into markers
[ridge_markers] = PSzz_components_to_markers(connected_components{1},input_image);
figure,imshow(min(input_image,1-ridge_markers)); title('Ridge markers superimposed on image');