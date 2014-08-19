%% Train script
% Train on the Berkeley Segmentation Dataset using Multiple Instance Learning.
% 
% iids:             training images vector
% groundPath:       path where the groundtruth ridge maps are stored
% featurePath:      path where the features used for training are stored
% trainMethod:      'mil' for Multiple Instance Learning
% numberOfSamples:  total number of samples used in training
% samplingMethod:   one of {'random','balanced','boundary'}
% featConfig:       one of {'color', 'gray', 'no_texture','no_boundary',
%                   'gray_no_boundary','fullcolor','color_spectral'}. 
%                   Determines the features used.
% boundPerc:        percentage of boundary pixels used for training (per image)
% buffer:           determines the area around groundtruth ridge pixels, from where
%                   we do not choose training samples
% 
% OUTPUTS
% beta:     The weights of the detector.
% xstd:     Standard deviation of the vector weights.
% p_bags:   The posterior probabilities calculated as an output by the
%           symmetry detector.
% lli:      The final value of the log-likelihood after the optimization.
% x:        The N x K x nInst array, carrying the features for all bags corresponding
%           to each one of the training samples. N is the number of samples, K is the
%           length of the feature vector and nInst is the number of instances
%           included in each bag.
% y:        The 0/1 labels for each sample.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% Setup

iids = imgList('train');    % list with training images
groundPath = fullfile('/home','tsogkas','symmetry_detector','Groundtruth','train');
featPath = fullfile('/home','tsogkas','symmetry_detector','Features');
savePath = fullfile('/media','Data','symmetry_detector','beta');
trainMethod = 'mil';
numberOfSamples = 500000;
samplingMethod = 'balanced';
featConfig = 'color';
buffer = 3;
boundPerc = 0.5;


%% Multiple Instance Learning

[beta, xstd] = trainSymMIL(iids,groundPath,featPath,numberOfSamples,...
        samplingMethod,featConfig,buffer,boundPerc);
save(fullfile(savePath,'learned_weight_vector.mat'))        


    
