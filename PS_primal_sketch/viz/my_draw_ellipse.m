
function h = draw_ellipse(x, rotation,eigenvals, outline_color, fill_color)
% DRAW_ELLIPSE(x, c, outline_color, fill_color)
%   Draws ellipses at centers x with covariance matrix c.
%   x is a matrix of columns.  c is a positive definite matrix.
%   outline_color and fill_color are optional.

n = 20;					% resolution
radians = [0:(2*pi)/(n-1):2*pi];
unitC = [sin(radians); cos(radians)];

if ~exist('outline_color'),
    outline_color = 'g';
end
coord_x = unitC(1,:)*eigenvals(1);
coord_y = unitC(2,:)*eigenvals(2);
h = [];

y(1,:) = rotation(1,1)*coord_x + rotation(1,2)*coord_y + x(1);
y(2,:) = rotation(2,1)*coord_x + rotation(2,2)*coord_y + x(2);
if 1==1,
    h  =plot(y(1,:),y(2,:),outline_color);
    set(h,'linewidth',2)
else
    %y = r*unitC + repmat(x(:, i), 1, n);
    if nargin < 4
        h = [h line(y(1,:), y(2,:), 'Color', outline_color)];
    else
        h = [h fill(y(1,:), y(2,:), fill_color, 'EdgeColor', outline_color)];
    end
end