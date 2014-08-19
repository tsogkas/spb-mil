%% Symmetry benchmark

% Script is used to test the performance of the symmetry detector on
% the full testing dataset. 
% See also function symmetryBench()

%% Setup

% Test image dataset
test_iids = imgList('test');
% Paths for groundtruth data, stored features and train data (beta
% weights). These are to be set according to your directory structure,
groundPath = fullfile('/media','Data','symmetry_detector','Groundtruth','test2');
betaPath = fullfile('/media','Data','symmetry_detector','beta');
scorePath = fullfile('/media','Data','symmetry_detector','scores','release_test');
nthresh = 50;
symmetryBench(test_iids,groundPath,betaPath,scorePath,nthresh);
symmetryBenchGraphs(scorePath)

