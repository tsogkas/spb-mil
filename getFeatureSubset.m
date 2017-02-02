function ind = getFeatureSubset(subset)
% Return indexes of feature subset
% 
%   ind = getFeatureSubset(subset)
% 
%   subset: one of {'spectral','color','gray','no_texture'}
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last update: August 2013

switch subset
    case 'spectral'
        ind = 1:14;
    case 'color'
        ind = 1:13;
    case 'gray'
        ind = [1,2,5,6,9,10,13];
    case 'no_texture'
        ind = [1:4,6:8,10:12];
    otherwise
        error('Feature subset not supported')
end