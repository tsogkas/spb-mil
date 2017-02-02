function h = draw_arrow(x, r, orientation,outline_color)
% draw filled circles at centers x with radii r.
% x is a matrix of columns.  r is a row vector.

n = 3;					% resolution
coords_x = [0:r/5:r];
coord_x = x(1) + cos(orientation)*(coords_x);
coord_y = x(2) + sin(orientation)*(coords_x);

h = [];
% hold is needed for fill()
held = ishold;
hold on


h = [h line(coord_x,coord_y, 'Color', outline_color)];
