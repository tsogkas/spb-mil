function  ewrap = wrap_e(e_mn,maxwrap);
[N_m,N_n,nf] = size(e_mn);
%ewrap = zeros(N_m+2*maxwrap,
ewrap(maxwrap+(1:N_m),maxwrap+(1:N_n),:)=e_mn;
ewrap(1:maxwrap,:,:) = ewrap(maxwrap+N_m +[-maxwrap + 1:0],:,:);
ewrap(:,1:maxwrap,:) = ewrap(:,maxwrap+N_n +[-maxwrap + 1:0],:);
ewrap(maxwrap+N_m + (1:maxwrap),:,:)=ewrap(maxwrap + (1:maxwrap),:,:);
ewrap(:,maxwrap+N_n + (1:maxwrap),:)=ewrap(:,maxwrap + (1:maxwrap),:);
