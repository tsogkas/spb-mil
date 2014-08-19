function inp_points = PS2z1_choose_points(inp_points,rel);
% inp_points = PS2z1_choose_points(inp_points,rel)
% 
% Rel is a set of acceptance criteria applied to the fields of 
% inp_points. All of them are enforced to keep a subset of inp_points.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

for k=1:length(rel),
    eval(sprintf('inp_points = keep_points(inp_points,find(inp_points.%s));',rel{k}));
end
