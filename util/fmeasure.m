% Compute f-measure from recall and precision
% 
%   f = fmeasure(r,p)

function [f] = fmeasure(r,p)
f = 2*p.*r./(p+r+((p+r)==0));
