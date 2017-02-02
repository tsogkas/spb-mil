function  res=  PSPz2_discard_overlapping(input_points,threshold);
% res=  PSPz3_discard_overlapping(input_points,threshold)
% 
% Simple heuristic code to discard overlapping blobs.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

c_m = input_points.c_m;
c_n = input_points.c_n;
coords = [c_m;c_n];
scales = input_points.scales;
eners = input_points.ener;
scale_kron =kron(scales,scales');
dists  = sqrt(dist2(coords',coords'));
dists_scales = sqrt(dist2(log(scales)',log(scales)'));
dists = dists./sqrt(scale_kron);
dists = dists + dists_scales*3;

together = sum(dists<threshold,2);
covered = find(together>1);
keep_ind = setdiff([1:length(c_m)],covered);
cnt = 0;
while ~isempty(covered)
    cnt  = cnt+1;
    ener_covered =  eners(covered);
    [max_ener,ind]  = max(ener_covered);
    ind = covered(ind);
    keep_ind  = [keep_ind,ind];
    neighbors = find(dists(ind,:)<threshold);
    covered = setdiff(covered,neighbors);
end

res= keep_points(input_points,sort(keep_ind));