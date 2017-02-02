function map = plotRectangleGrid(im,scale,ratio,color)

[height,width,~] = size(im);
horizontalStep = 2*scale*ratio + 1;
verticalStep   = 2*scale + 1;
map = false(height,width);
map(:,1:horizontalStep:width) = true;
map(1:verticalStep:height,:) = true;
imshow(overlayBinaryImage(im,map,color),[]);



