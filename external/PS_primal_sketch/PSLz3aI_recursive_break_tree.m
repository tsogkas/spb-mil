function bin_tree =PSLz3aI_recursive_break_tree(c_x,c_y,index)
% bin_tree =PSLz3aI_recursive_break_tree(c_x,c_y,index)
%
% Break a continuous curve into a set of straight line 
% tokens using the procedure described in Lowe's SCEPRO '87 paper. 
% The procedure is recursive and scale invariant.
% At each recursion level it breaks  the curve into two curve segments
% at the point where the curve difference from its straight-line
% approximation is maximal.
% Terminates when either the curve length is too short,
% or the scale-invariant merit value does not increase when further 
% breaking the curve apart. 
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007


%% bin_tree struct 
%% allow -> 1/0: split was successful/failed
%% cost  -> figure of merit at current tree
%% points: indexes of starting and ending points  


if length(c_x)<3,
    %% if the tree is too short, stop    
    bin_tree.allow = 0;
    bin_tree.merit = .1;
    bin_tree.points =  [index(1),index(end)];
else
    %% otherwise go on
    bin_tree.allow    = 1;    
    diff_x = c_x(end) - c_x(1);
    diff_y = c_y(end) - c_y(1);
    
    %%  line approximation is locus of a x + b y + c =0
    %% a,b,c are calculated from the endpoint locations
    if diff_y~=0,
        b=  - (diff_x)/(diff_y);
        nrm  = sqrt(1 + b*b);
        a = 1/nrm; b = b/nrm;
    else 
        a = -diff_y/diff_x;
        nrm = sqrt(1 + a*a);
        a = a/nrm; b = 1/nrm;
    end    
    length_here = sqrt(diff_x*diff_x + diff_y*diff_y);
    c = - a*c_x(1) - b*c_y(1);
    
    %% get the distances of all curve points from the straight line
    %% approximation, find the point with maximal distance and 
    %% get figure of merit for current straight-line approximation

    dists = abs(a*c_x + b*c_y + c);
    [md,ind] = max(dists(2:end-1));
	merit_here = (length_here)/max(md,1);
    keep_original = 1;    
    
    if md>0,

        %% go on, and do the same thing for the subsegments
        ind = ind+1;
        left_relative  = PSLz3aI_recursive_break_tree(c_x(1:ind-1),c_y(1:ind-1),index(1:ind-1));
        right_relative = PSLz3aI_recursive_break_tree(c_x(ind:end),c_y(ind:end),index(ind:end));

        merit_left  = left_relative.merit;
        merit_right = right_relative.merit;

        allow_left = left_relative.allow;
        allow_right = right_relative.allow;
        points_left = left_relative.points;
        points_right = right_relative.points;
    else
        allow_left = 0; allow_right = 0;
    end
    
    %% both leaves allow for break
    if (allow_left&allow_right),  
        %% and the merit for one of the two leaves is higher
        if max([merit_left,merit_right])>=merit_here, 
            keep_original = 0;
            bin_tree.merit    = [merit_left, merit_right];
            bin_tree.points   = [points_left(1:end), points_right(2:end)];
        end
    end
    if keep_original,
       bin_tree.merit = merit_here;
       bin_tree.points = [index(1),index(end)];
    end
end
