function [ridge_markers] = PSzz_components_to_markers(ridge_components,input_image);
% [ridge_markers] = PSzz_components_to_markers(ridge_components,input_image)
% 
% Turns connected ridges into a marker image. 
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

ridge_markers = zeros(size(input_image));

%% apply some  thresholds to rule out weak ridges 
ridge_components =  keep_points(ridge_components,find(ridge_components.lst>4));
for k=1:length(ridge_components.lstrings),
    mean_ener(k) = mean(ridge_components.attribs{k}.ener);
    mean_scale(k)= mean(ridge_components.attribs{k}.scl);
end
ridge_components = keep_points(ridge_components,find((mean_ener>.1)&(mean_scale>2)));

%% and then put a single label for each connected component
for k= 1:length(ridge_components.lstrings),
     ridge_markers(ridge_components.lstrings{k}(3:end-2)) = k;
end
