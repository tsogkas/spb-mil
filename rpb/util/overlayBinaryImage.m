function out = overlayBinaryImage(im,bin,color)
% function out = overlayBinaryImage(im,bin,color)
% 
% Overlays binary image on image, using the selected color.
% 
% INPUTS
% im:       input image (RGB or grayscale).
% bin:      binary image.
% color:    color used for overlay.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% February 2011

if nargin<3, color = 'y'; end   % default color is yellow
bin = logical(bin); % turn to binary if not already
idx = find(bin);
switch color
    case 'y'
        rgb = [1 1 0];
    case 'r'
        rgb = [1 0 0];
    case 'g'
        rgb = [0 1 0];
    case 'b'
        rgb = [0 0 1];
    case 'm'
        rgb = [1 0 1];
    case 'w'
        rgb = [1 1 1];
    otherwise
        error('Invalid color')
end

if ndims(im)==1, im = cat(3,im,zeros(size(im)),zeros(size(im))); end
ch1 = im(:,:,1); ch1(idx) = rgb(1);
ch2 = im(:,:,2); ch2(idx) = rgb(2);
ch3 = im(:,:,3); ch3(idx) = rgb(3);
out = cat(3,ch1,ch2,ch3);



