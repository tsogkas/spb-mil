function res = PSzz_zip_skeleton(skeleton,orien,scl,ener);
% res = PSzz_zip_skeleton(skeleton,orien,scl,ener)
%
% Form strcuture containing information about the skeleton
% only on its sparse non-zero points
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

wt = find(skeleton);
res.idx = wt';
res.val = skeleton(wt)';
res.orn = orien(wt)';
res.scl = scl(wt)';
res.ener =ener(wt)';
res.sz =  size(skeleton);
