function [coord_m,coord_n] = make_coords(size_n,size_m);
[coord_m,coord_n] =  meshgrid([1:size_m],[1:size_n]);