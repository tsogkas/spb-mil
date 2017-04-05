function string_structure = PSLz2c_find_connected_curves(current,succ,pred,imsize);
% string_structure = PSLz2c_find_connected_curves(current,succ,pred,imsize)
%
% Extract continuous contour segments using the neighborhood information of
% succ/pred.
% 
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007


%%-------------------------------------------------------------------------
%% lookup table
%% if lookup(k) = i,
%% this means that at the image index `k' the i'th skeleton point is present.
%% Otherwise (lookup(k) = 0) no point is active.
%% make also a matrix (passed) to indicate whether a curve has already
%% taken hold of the corresponding skeleton point
%%-------------------------------------------------------------------------

structure = current; expand_structure;
lookup(max(indexes))=0;

passed(indexes)  = false;
for k=1:length(indexes),
    lookup(indexes(k)) = k;
end
permt =  randperm(length(indexes));
%for k= [permt,permt,permt,permt],
%    passed(indexes(k)) = 1;
    %end
%% get those points that have both a predecessor and successor
%% the other points serve as starting points

isok =  (succ~=-1)&(pred~=-1);
starting_points = find(~isok);

for k=1:length(starting_points)    
    [passed,string_str,ener_str,scale_str] = PSLz2cI_track_curve(starting_points(k),passed,indexes,ener,scl,lookup,succ,pred,2);
    lstrings{k} = string_str;
    attribs{k}.ener  = ener_str;
    attribs{k}.scl   = scale_str;
    lst(k) = length(lstrings{k});
end

%% Check if some line points were not assigned to any curve.
%% This can happen e.g if we have a closed contour, and therefore 
%% no starting point. 
%% If this is the case, start a curve at an arbitrary point and 
%% track the curve both forward and backward.
cont =1;
while cont
    not_passed = find(passed(indexes)~=1);
    cont= ~isempty(not_passed);
    if cont
        [passed,string_str,ener_str,scale_str] = PSLz2cI_track_curve(not_passed(1),passed,indexes,ener,scl,lookup,succ,pred,2);
        lstrings{end+1} = string_str;
        lst(end+1)      = length(lstrings{end});
        attribs{end+1}.ener  = ener_str;
        attribs{end}.scl   = scale_str;        
    end
end

fields_wt = {'lstrings','lst','attribs','imsize'};
compress_structure;
string_structure = structure;