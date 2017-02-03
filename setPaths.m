function paths = setPaths()
% SETPATHS Setup paths based on default directory structure.
% 
%   paths = setPaths()
% 
%   Stavros Tsogkas, <tsogkas@cs.toronto.edu>
%   Last update: January 2017

paths.root             = fileparts(mfilename('fullpath'));

% data
paths.data             = fullfile(paths.root, 'data');
paths.output           = fullfile(paths.root, 'output');
paths.bsds500          = fullfile(paths.data,'BSDS500');
paths.bsds500gtTrain   = fullfile(paths.bsds500,'groundtruth','train');
paths.bsds500gtTest    = fullfile(paths.bsds500,'groundtruth','test');
paths.bsds500gtVal     = fullfile(paths.bsds500,'groundtruth','val');
paths.bsds500imTrain   = fullfile(paths.bsds500,'images','train');
paths.bsds500imTest    = fullfile(paths.bsds500,'images','test');
paths.bsds500imVal     = fullfile(paths.bsds500,'images','val');
paths.symmax500        = fullfile(paths.data,'SYMMAX500');
paths.symmax500gtTrain = fullfile(paths.symmax500,'train');
paths.symmax500gtTest  = fullfile(paths.symmax500,'test');
paths.symmax500gtVal   = fullfile(paths.symmax500,'val');

% models
paths.models           = fullfile(paths.output,'mil','models');


