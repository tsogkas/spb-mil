%% This demo demonstrates the use of rpb ridge detector 

%% Add directories to path and read image
rootFolder = 'symmetry';
addpath(genpath(rootFolder));   % add necessary directories to path
im = double(imread('101087.jpg'))/255;  % read input image
% (im the berkeley segbench code is included in your path, you can use
% im = imgRead(image_iid) instead

%% Straightforward use of rpb function

[sup,pfat] = rpb(im);           % run detector
% [sup,pfat] = rpb(im,'color'); % same as previous command
% [sup,pfat] = rpb(im,'gray');  % run detector using only brightness

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

%% Use rpb to test with a new beta vector instead of the learned ones
% You can use rpb with a new beta vector as input

% [sup,pfat] = rpb(im,'color',newBeta,0);          

%% Using spectral feature
% If you have pre-calculated the spectral feature for symmetry, you can use
% it as an additional parameter.

% [sup,pfat] = rpb(im,'spectral',[],0,dlc,drc,dlr,spectralFeatureArray);          

%% Visualizing results

figure(1), imshow(im);      title('Input image')
figure(2), imshow(pfat);    title('Responses before non-maximum suppression')
figure(3), imshow(sup);     title('rpb (after non-maximum suppression)')
figure(4), imshow(sup>0.5); title('rpb (thresholded)')

