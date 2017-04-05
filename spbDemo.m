%% Demonstrates the use of the spbMIL symmetry detector 

img = im2double(imread('101087.jpg'));  % read input image

%% Straightforward use of spbMIL

spb = spbMIL(img); % run detector using "color" set of features
% spb = spbMIL(img,'featureSet','color'); % same as previous command
% spb = spbMIL(img,'featureSet','gray');  % run detector using only brightness features

%% Multiple tests on one image
% If you want to make multiple tests on the same image, avoid extracting 
% features all the time, by calculating them once and passing them as an
% argument to the rpb function.

histFeatures = computeHistogramFeatures(img);
spb = spbMIL(img,'histFeatures', histFeatures);          

%% spb fields
% The output of spbMIL contains various useful fields:
% spb.fat:        symmetry probability map (before nonmax suppresion)
% spb.thin:       thinned symmetry probability map
% spb.orientMap:  dense orientation map of symmetry responses

%% Use spbMIL to test with a new weight vector 
% model = trainMIL(...)
% spb = spbMIL(img,'w',newWeights);          

%% Using spectral feature
% If you have pre-calculated the spectral feature for symmetry, you can do:

% tmp = load(spfeatPath); spfeat = tmp.spectral;
% spb = spbMIL(img,'spectralFeature',spectral);          

%% Visualizing results

figure(1), imshow(img);             title('Input image')
figure(2), imshow(spb.pfat);        title('Symmetry probability map')
figure(3), imshow(spb.thin);        title('Spb after non-maximum suppression')
figure(4), imshow(spb.thin>0.5);    title('spb thresholded')

