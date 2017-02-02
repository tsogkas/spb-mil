function sp = permuteThetas(sPb)

temp = sPb;
sp = sPb;
% ch_theta = [1,8,7,6,5,4,3,2];
ch_theta = [5,4,3,2,1,8,7,6];
for iter=1:size(sPb,3)
    sp(:,:,iter,:) = temp(:,:,ch_theta(iter),:);
end