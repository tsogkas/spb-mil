function res= filter1(filt,input,dim,patch_size);

if nargin<4,
    patch_size = 0;
else
    [size_m,size_n]  =size(input);
    patch_size_m = patch_size*(dim==1);
    patch_size_n = patch_size*(dim==2);

    dim_wt_m = [(patch_size_m + 1):(patch_size_m+size_m)];
    dim_wt_n = [(patch_size_n + 1):(patch_size_n+size_n)];
    
    input = [zeros(patch_size_m,size_n);input;zeros(patch_size_m,size_n)];
    input = [zeros(size_m+2*patch_size_m,patch_size_n),input,...
             zeros(size_m+2*patch_size_m,patch_size_n)];
end

[size_m,size_n] = size(input);
filt = reshape(filt,[1, length(filt)]);
l_f = ceil((length(filt)+1)/2);

input = [[input;zeros(l_f*(dim==1),size_n)],zeros(size_m + l_f*(dim==1),l_f*(dim==2))];
res = filter(fliplr(filt),1,input,[],dim);
res = res((1:size_m)+(dim==1)*(l_f-1),(1:size_n)  + (dim==2)*(l_f-1));
if nargin==4,
    res = res(dim_wt_m,dim_wt_n);
end