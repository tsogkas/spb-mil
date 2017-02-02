function  [crd_x,crd_y]  = PSzz_token_raw_information(component,k_ind);
imsize  = component.imsize;
indexes = component.lstrings{k_ind};
[crd_y,crd_x] = ind2sub(imsize,indexes);
