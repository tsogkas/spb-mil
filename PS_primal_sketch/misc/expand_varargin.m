for k_varargin=1:length(varargin)/2,
    eval(sprintf('%s = varargin{2*k_varargin};',varargin{2*k_varargin-1}));
end