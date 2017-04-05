function a_s = my_patch(a,band_width,new_values);
[size_x,size_y]=size(a);
if nargin ==2,
    a_s=[a(:,band_width:-1:1),a,a(:,(size_y):-1:(size_y-(band_width-1)))];
    a_s=[a_s(band_width:-1:1,:);a_s;a_s((size_x):-1:(size_x-(band_width-1)),:)];
else
    a_s=[new_values*ones(size(a,1),band_width),a,new_values*ones(size(a,1),band_width)];
    a_s=[new_values*ones(band_width,size(a_s,2));a_s;new_values*ones(band_width,size(a_s,2))];
end

