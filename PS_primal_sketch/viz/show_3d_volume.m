function mov_volume =  show_3d_volume(val_3d_array,step);
% mov_volume =  show_3d_volume(val_3d_array,step);
% 
% Receives as input a 3d image  (3rd dimension is time/scale % etc.)
% and visualizes this as a 3d cube.
% `step' is for subsampling the x/y dims (default is 2);
% 
% If mov_volume is requested, it makes a movie of this evolution.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

if nargin ==1, step =2; end
load mri;
%isocaps(flipdim(1-val_3d_array(1:step:end,1:step:end,:),1),-.1);

views = [   0.9063    0.4226         0   -0.6645;...
        -0.3696    0.7927    0.4848   -0.4539;...
        -0.2049    0.4394   -0.8746    8.9803;...
        0         0         0    1.0000];

ha = axes;
nvol = size(val_3d_array,3);
for k=2:nvol,
    isocaps(flipdim(1-val_3d_array(1:step:end,1:step:end,1:k),1),-.1);
    axis('tight'); %axis(h0);
    colormap(map)
    view(views);
    h0 = axis; axis([h0(1:5),nvol]);
    set(ha,'xtick',[],'ytick',[],'ztick',[]);
    if nargout~=0,
        mov_edge(k-1) = getframe;
    end
    pause(.01);
end