function [res,orn,scl] = PSzz_unzip_skeleton(skeleton);
% [res,orn] = PSzz_unzip_skeleton(skeleton)
%
% Turn zipped skeleton into image.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

res= zeros(skeleton.sz);
try,
    res(skeleton.idx) = skeleton.ener;
catch
    res(skeleton.idx) = skeleton.val;
end
if nargout>=2,
    orn = zeros(skeleton.sz);
    orn(skeleton.idx) = skeleton.orn;
end

if nargout>=3,
    scl = zeros(skeleton.sz);
    scl(skeleton.idx) = skeleton.scl;
end