function paths = setPaths()
% SETPATHS Setup paths based on default directory structure.
% 
%   paths = setPaths()
% 
%   Stavros Tsogkas, <tsogkas@cs.toronto.edu>
%   Last update: February 2017

paths.spbmil.root      = fileparts(mfilename('fullpath'));

% data
paths.data             = fullfile(paths.root, 'data');
paths.bsds500          = fullfile(paths.data,'BSDS500');
paths.bsds500gt        = fullfile(paths.bsds500, 'groundtruth');
paths.bsds500im        = fullfile(paths.bsds500, 'images');
paths.symmax500        = fullfile(paths.data,'SYMMAX500');

% output/results
paths.spbmil.output    = fullfile(paths.spbmil.root, 'output');
psths.spbmil.plots     = fullfile(paths.spbmil.output, 'plots');
paths.spbmil.models    = fullfile(paths.spbmil.output,'models');
paths.spbmil.spectral  = fullfile(paths.spbmil.output, 'spectral');
