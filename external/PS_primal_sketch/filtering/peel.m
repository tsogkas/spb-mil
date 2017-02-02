function peeled  = peel(input,patch_size,p_s_2,p_s_3,p_s_4,inverse);

[size_x,size_y,nd] = size(input);
if nargin==2,
    peeled = input(patch_size+1:size_x-(patch_size),patch_size+1:size_y-(patch_size),:);
elseif nargin == 5
    peeled = input(patch_size+1:(size_x - p_s_2), (p_s_3+1):(size_y-p_s_4),:);
else
    peeled =  -1*ones(size_x + p_s_2+patch_size,size_y+p_s_3+p_s_4,nd);
    [size_x,size_y,nd]  =size(peeled);
    peeled(patch_size+1:(size_x - p_s_2), (p_s_3+1):(size_y-p_s_4),:) =  input;
end
