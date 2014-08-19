%% demo1.m
%% basic primal sketch functionality:
%% extract edge/ridge/blob tokens from a grayscale image
%% 
%% Iasonas Kokkinos <jkokkin@stat.ucla.edu>

%% load image and convert to double in [0,1]
% input_image = double(imread('horse010.jpg'))/256;
input_image = double(rgb2gray(imread('E:\Thesis\images\BSDS300\images\train\46076.jpg')))/256;
[ridges,edges,blobs,contours] = PS00___primal_sketch(input_image);

contour_ridge   = PSzz_unzip_contour(contours{1}); 
contour_edge    = PSzz_unzip_contour(contours{2});
colored_sketch  = colored_contours(input_image,contour_ridge,contour_edge);

%% show sketch tokens
figure,
subplot(2,2,1); show_ridges_on_image(input_image,ridges); title('Ridges')
subplot(2,2,2); show_edges_on_image(input_image,edges);   title('Edges');
subplot(2,2,3); show_ellipses_on_image(input_image,blobs);title('Blobs')
subplot(2,2,4); imshow(colored_sketch); title('Ridge and Edge contours');

figure, imshow(contour_ridge)