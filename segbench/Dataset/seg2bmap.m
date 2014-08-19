
function [bmap] = seg2bmap(seg)
% [bmap] = seg2bmap(seg)
%
% INPUTS
%	seg	Segments labeled from 1..k.
%
% OUTPUTS
%	bmap	Binary boundary map.
%

% check segmentation
if min(seg(:))<1, error('bug'); end
if length(unique(seg(:)))~=max(seg(:)), error('bug'); end

[h,w] = size(seg);
k = max(seg(:));
bmap = zeros(size(seg));

lut = makelut(@f,3);
for i = 1:k,
  bmap = bmap | applylut((seg==i),lut);
end
bmap(1,:)=0; bmap(h,:)=0; % clear image boundary
bmap(:,1)=0; bmap(:,w)=0;

function y = f(x)
y = (sum(x(:)) < 9) & (sum(x(:)) > 0);
