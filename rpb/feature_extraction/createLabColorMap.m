function cmap = createLabColorMap(im,nbins)
% function cmap = createLabColorMap(im,nbins)
% 
% Create LAB-space colormap for input image (RGB or grayscale).
% 
% INPUTS
% im:       original image
% nbins:    number of bin labels

if nargin<2, nbins = ones(3,1)*32; end
if numel(nbins)==1, nbins = ones(3,1)*nbins; end

gamma = 2.5;
% min and max values for a,b channels of LAB
% used to scale values into the unit interval
abmin = -73;
abmax = 95;
if ndims(im)==2, % grayscale image
    cmap = max(1,ceil(im*nbins(1)));
else % RGB image
    cmap = zeros(size(im,1),size(im,2),3);
    % convert gamma-corrected image to LAB and scale values into [0,1]
    lab = RGB2Lab(im.^gamma);
    lab(:,:,1) = lab(:,:,1) ./ 100;
    lab(:,:,2) = (lab(:,:,2) - abmin) ./ (abmax-abmin);
    lab(:,:,3) = (lab(:,:,3) - abmin) ./ (abmax-abmin);
    lab(:,:,2) = max(0,min(1,lab(:,:,2)));
    lab(:,:,3) = max(0,min(1,lab(:,:,3)));
    for i=1:3, cmap(:,:,i) = max(1,ceil(lab(:,:,i)*nbins(i))); end
end

