function points_kept = PSLz1a_isolate_same_location(points,size_m);
% points_kept = PSLz1a_isolate_same_location(points,size_m)
%
% Removes edge/ridge points lying on identical image locations 
% (happens due to different structures residing at different locations).
% Simple remedy: keep only the strongest edge point.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

ener = points.ener;
indexes = points.indexes;
%%-------------------------------------------------------------------------
%% identify edge points at identical locations
%%-------------------------------------------------------------------------
[srt,ind] = sort(indexes); 
% identical points should form blocks of indexes with identical value 
% start of identical: previous is not the same  but next is identical. 
wt_start =  find((srt(2:end-1) ==srt(3:end)).*(srt(1:end-2)~=srt(2:end-1)))+1;
% end of identical: previous is the same but next is not the same
wt_end   = find((srt(2:end-1) ~=srt(3:end)).*(srt(1:end-2)==srt(2:end-1)))+1;

if ~isempty(wt_start)  
    if wt_end(1)<wt_start(1), 
        wt_end = wt_end(2:end); 
    end    
    if wt_start(end)>wt_end(end),
        wt_start = wt_start(1:end-1);
    end
    
    %% for each cluster of identical points, pick the one with
    %% maximum energy and put the others in the remove list
    
    remove = [];
    for k=1:length(wt_start)
        indexes_same  = ind(wt_start(k):wt_end(k));
        ener_same     = ener(indexes_same);
        [maxen,indh] =max(ener_same);    
        indexes_same(indh) =[];
        remove =[remove,indexes_same']; 
    end
    keep = setdiff(1:length(indexes),remove);
    points_kept = keep_points(points,keep);
else
    points_kept = points;
end