% Parameter management script, for use with any function that
% accepts a variable list of arguments.  Note: defaults should be
% assigned to variables before running this script

if ~exist('varargin','var')
    error('The variable ''varargin'' is not defined. Please include it is the arguments of the function.')
end

for i=1:2:length(varargin)
  if ~isstr(varargin{i}) | ~exist(varargin{i},'var')
    if ~isstr(varargin{i})
      error('invalid parameter');
    else
      error(['invalid parameter ' varargin{i}]);
    end
  end
  eval([varargin{i} '= varargin{i+1};']);
end
