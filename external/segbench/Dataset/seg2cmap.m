function [cmap,cid2sid] = seg2cmap(seg,bmap)
% [cmap,cid2sid] = seg2cmap(seg,bmap)
%
% INPUTS
%	seg	Segments labeled from 1..k.
%	bmap	Binary boundary map.
%
% OUTPUTS
%	cmap	Contours labeled from 1..n.
%		Junctions labeled -1.
%		Other pixels labeled 0.
%	cid2sid 
%		2xn array mapping contour ID to segment IDs.
%		e.g. cid2sid(:,i) gives segment IDs for contour i.
%
% David Martin <dmartin@eecs.berkeley.edu>
% November, 2002

% check segmentation
if min(seg(:))<1, error('bug'); end
if length(unique(seg(:)))~=max(seg(:)), error('bug'); end

[h,w] = size(seg);
k = max(seg(:));

ids = zeros(k,k); % map seg ID pairs to contour ID
cid2sid = zeros(2,0); % map contour ID to seg ID pairs
cmap = zeros(size(seg));

nextID = 1;

for x = 1:w,
  for y = 1:h,
    if bmap(y,x)==0, continue; end
    x1 = max(1,x-1); x2 = min(w,x+1);
    y1 = max(1,y-1); y2 = min(h,y+1);
    s = seg(y,x); % segment of current pixel
    other = 0; % neighboring segment
    morethanone = 0; % if we found more than 1 neighboring segment
    for i = x1:x2,
      for j = y1:y2,
        s2 = seg(j,i);
        if s2 ~= s, % different segment
          if other==0 | other==s2, % just one other so far
            other = s2;
          else % more than one other segment
            morethanone = 1;
          end
        end
      end
    end
    if morethanone, % it's a junction
      cmap(y,x) = -1; % mark it
    else % not a junction
      a = min(s,other);
      b = max(s,other);
      if ids(a,b)==0,
        ids(a,b) = nextID;
        cid2sid(:,nextID) = [a;b];
        nextID = nextID + 1;
      end
      cmap(y,x) = ids(a,b);
    end
  end
end

% break up disconnected contours in cmap
ncontours = nextID-1;
for i = 1:ncontours,
  m = bwlabel(cmap==i);
  count = max(m(:));
  for j = 2:count,
    cmap(find(m==j)) = nextID;
    cid2sid(:,nextID) = cid2sid(:,i);
    nextID = nextID + 1;
  end
end

% make sure cmap is non-zero wherever bmap is non-zero
if (cmap~=0) ~= (bmap~=0),
  error('bug');
end

if size(cid2sid,1) ~= 2, error('bug'); end
if size(cid2sid,2) ~= max(cmap(:)), error('bug'); end
if min(cmap(:)) < -1, error('bug'); end
