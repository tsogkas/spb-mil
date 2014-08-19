%% This demo demonstrates the use of rpb symmetry detector 
% Each one of the cells should be executed separately (e.g. by using ctrl+r). 
% They are meant to give some quick examples of the various arguments that 
% can be passed into rpb.m. De-comment the cell you want to try out.


%% Add directories to path and read image
startup
im = double(imread('101087.jpg'))/255;  % read input image
% (if the berkeley segbench code is included in your path, you can use
% im = imgRead(image_iid) instead.

%% Straightforward use of rpb function

[sup,pfat] = rpb(im);           % run detector
% [sup,pfat] = rpb(im,'color'); % same as previous command
% [sup,pfat] = rpb(im,'gray');  % run detector using gray features
                                % (brightness and texture)

%% Visualizing results

figure;
subplot(2,2,1); imshow(im); title('Input image')
subplot(2,2,2); imshow(pfat); title('Responses before non-maximum suppression')
subplot(2,2,3); imshow(sup); title('rpb (after non-maximum suppression)')
subplot(2,2,4); imshow(sup>0.5); title('rpb (thresholded)')

%% Use rpb to test with a new beta vector instead of the learned ones
% You can use rpb with a new beta vector as input (the dimensions of the
% bew weight vector must agree with the dimensionality of the feature
% vector used).

% newBeta = rand(14,1);
% [sup,pfat] = rpb(im,'color',newBeta,0);    

%% Multiple tests on one image
% If you want to make multiple tests on the same image, avoid extracting 
% features all the time, by calculating them once and passing them as an
% argument to the rpb function.

% [dlc,drc,dlr] = histGradFeatures(im);
% [sup,pfat] = rpb(im,'color',[],0,dlc,drc,dlr);          

%% Use rpb to return feature array
% You can use rpb to return only the multidimensional array of concatenated
% features, avoiding the cost of evaluating the rpb response.

% [~,~,~,f] = rpb(im,'color',[],1,dlc,drc,dlr);

%% Using spectral feature
% If you have pre-calculated the spectral feature for symmetry, you can use
% it as an additional parameter.

% [sup,pfat] = rpb(im,'spectral',[],0,dlc,drc,dlr,spectralFeatureArray);          


