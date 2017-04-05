function spb = spbMIL(img,varargin)
% Construct features array, probability map at all scales and orientations,
% orientation map and non-maximum suppressed output, depending on the input
% parameters.
% 
%   spb = SPBMIL(img,varargin), where img is a grayscale or RGB input image
% 
%   OPTIONAL INPUT ARGUMENTS
%   featureSet:   feature configuration used to extract features and evaluate
%                 probability responses. One of {'spectral','color','gray','no_texture'}.
%   fineScale:    use features at fine scale {true}.
%   w:            use input weight vector instead of the default.
%   histFeatures: histogram feature struct if already pre-computed, to avoid
%                 unnecessary computations.
%   spectralFeature: pre-computed spectral feature.
% 
%   OUTPUTS
%   spb:            struct containing the following fields:
%   spb.fat:        symmetry probability map.
%   spb.thin:       thinned symmetry probability map
%   spb.orientMap:  dense orientation map of symmetry responses
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: February 2017


% Default parameters 
prm = { 'featureSet',        'color', ...
        'fineScale',         true,...
        'w',                 [],...
        'histFeatures',      [],...
        'spectralFeature',   []
      };
prm = parseVarargin(prm, varargin);     

% If input image is grayscale, use appropriate model
if ismatrix(img) && ~strcmp(prm('featureSet'),'gray')
    warning('The input image is grayscale but you are not using ''gray'' feature subset')
    warning('Changing feature subset used to ''gray''')
    featureSet = 'gray';
else
    featureSet = prm('featureSet');
end

if strcmp(featureSet, 'spectral') && isempty(spectral)
    error('You have not provided spectral feature');
end

% Optionally use user-specified weight vector
if isempty(prm('w'))
    w = defaultWeights(featureSet);
else
    w = prm('w');
end

% If the feature set is 'gray' convert image to grayscale
img = im2double(img);
if strcmp(featureSet, 'gray'), img = rgb2gray(img); end


% Compute histogram features
if isempty(prm('histFeatures'))  
    histFeatures = computeHistogramFeatures(img,'fine',prm('fineScale'));
else
    histFeatures = prm('histFeatures');
end
[height,width,nOrient,nScales,~] = size(histFeatures.dlc);
thetas   = histFeatures.thetas;
scales   = histFeatures.scales;
step     = histFeatures.step;
features = cat(5, ones(height,width,nOrient,nScales,1,'single'),...
    histFeatures.dlc, histFeatures.drc, histFeatures.dlr);

% Combine with spectral feature if provided
if ~isempty(prm('spectralFeature'))
    features = cat(5,features,repmat(spectral,[1,1,1,nScales]));
end
features = features(:,:,:,:,getFeatureSubset(featureSet));
assert(length(w)==size(features,5),...
    'Size of width vector does not agree with feature configuration');

% Compute dense symmetry response
features = reshape(features,height*width*nOrient*nScales,[]);
pb = reshape(1./(1+exp(-features*w)),[height,width,nOrient,nScales]);
[~,orientMap] = max(max(pb,[],4),[],3);
[~,scalesMap] = max(max(pb,[],3),[],4);
for s = 1:nScales
    scale = floor(scales(s)/step(s));
    for o = 1:nOrient
        tmp = imresize(pb(:,:,o,s),1/step(s),'bilinear');
        a = fitparab2(tmp,scale,scale,thetas(o));
        pb(:,:,o,s) = imresize(a,[height,width],'bilinear');
    end
end

% Noisy-OR response and non-maximum suppression 
p_inst      = reshape(pb,[height*width,nOrient*nScales]);
p_inst      = max(0,min(p_inst,1));
logp_bags   = sum(log(1-p_inst+eps),2);
p_bags      = reshape(1-exp(logp_bags),[height,width]);
sup         = nonmax(p_bags,thetas(orientMap));
spb.fat     = p_bags;
spb.thin    = sup;
spb.scales  = histFeatures.scales;
spb.thetas  = histFeatures.thetas;
spb.orientMap = orientMap;
spb.scalesMap = scalesMap;



% -------------------------------------------------------------------------
function w = defaultWeights(featureSet)
% -------------------------------------------------------------------------
switch featureSet
    case 'spectral'
        w = [-7.277058539808386;   0.075769718374560; 
              0.878473484015558;   1.583249092499895;
              3.899242166812703;  -0.121294943600757;
              1.452889638705204;   1.319490136714483;
              3.927183882128103;  -1.261347594455979;
             -1.003812310494814;  -1.043351555773494;
             -1.770948704033194;   0.009880627418217;];
    case 'color'
%         weights = [-6.54001569881004;-0.317591380311371;1.87794646755063;1.28269948188804;3.99614300468390;0.289616770646934;0.957724691090984;1.09607455266995;4.28043755023855;-1.01656506721445;-1.19665557212154;-0.860407782628580;-2.59753117106842;]; % old
        w = [-7.318570926307534;    0.369244954373050;
              0.937115112670649;    1.654541496325269;
              4.072846883190388;    0.017486829124990;
              1.448453574841385;    1.513514483735794;
              4.242503892502128;    -1.427126231729479;
              -1.095028987718087;   -1.132117076482029;
              -2.103062445150381;];
    case 'gray'
%         weights = [-6.44514046390470;1.20879203131223;4.37178250414567;1.48029354631580;4.33415914637308;-2.24317741581830;-2.86022220830624;]; %old
        w = [-7.188621019359869;    1.822216363721022; 
             4.201182314246688;     1.536088795176045; 
             4.519480613821219;     -2.797315940486013; 
             -2.356838424154210;];
    case 'no_texture'
%         weights = [-5.73872428102262;1.88976506325115;1.87838349589835;2.07200591874450;2.44810248820298;1.19908797114823;1.75638250148779;-2.34004274862761;-1.18088564539031;-1.44830602726111;]; % old
        w = [-6.188826728436059;    2.423153581623645; 
             1.282953171347384;     2.164423693889976; 
             2.437648610064436;     1.704396774046794; 
             1.859980407702409;     -2.524623681633751; 
             -1.107941426953260;    -1.603304887584006;];
    otherwise
        error('Feature configuration not supported')
end

