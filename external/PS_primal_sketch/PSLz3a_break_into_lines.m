function res = PSLz3a_break_into_lines(components,size_m,ener);
% res = PSLz3a_break_into_lines(components,size_m,size_n,ener)
% 
% Takes continuous curves and breaks them into straight line segments.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

structure = components;
expand_structure;
lines = []; merit_geom = [];
for k=1:length(lstrings),
    coords = lstrings{k};
    if length(coords)>5

        %% turn curve 1-d indexes into [c_m,c_n] coordinates
        coord_n = ceil(coords/size_m);
        coord_m = coords - (coord_n - 1)*size_m;

        lines_comp    = {};
        bin_tree = PSLz3aI_recursive_break_tree(coord_m,coord_n,1:length(coord_n));

        %% bin_tree struct
        %% cost  -> figure of merit at current tree
        %% points: indexes of starting and ending points of leaves

        returned = bin_tree.points;
        for ln=1:length(returned)-1,
            wt = returned(ln):returned(ln+1);
            lines_comp{ln} = [coord_m(wt);coord_n(wt);coords(wt);attribs{k}.ener(wt);attribs{k}.scl(wt)];
        end

        lines = [lines,lines_comp(:)'];
        merit_geom = [merit_geom,bin_tree.merit];
    end
end

fields_wt = {'lines','merit_geom'};
compress_structure;
res = structure;