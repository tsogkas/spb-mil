function [sup,p_bags,orientMap,f] = rpb(im,featConfig,newBeta,featOnly,dlc,drc,dlr,sPb)
% function [sup,p_bags,orientMap,f] = rpb(im,featConfig,newBeta,featOnly,dlc,drc,dlr,sPb)
% 
% Construct features array, probability map at all scales and orientations,
% orientation map and non-maximum suppressed output, depending on the input
% parameters.
% 
% INPUTS
% im:           original image (grayscale or RGB).
% featConfig:   feature configuration used to extract features and evaluate
%               probability responses. One of {'spectral','color','gray','no_texture'}.
% newBeta:      use input weight vector instead of the learned one.
% featOnly:     true if you want to return only array of features and false if you
%               want to evaluate the detector response (default is false).
% dlc/drc/dlr:  histogram differences if already pre-computed, to avoid
%               unnecessary computations.
% sPb:          pre-computed spectral feature using gPb code.
% 
% OUTPUTS
% f:            array of concatenated features.
% orientMap:    orientation map of symmetry (also used for nmax suppression).
% p_bags:       array of symmetry probabilities at all scales and orientations.
% sup:          non-maximum suppressed response of symmetry.
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% February 2011

%% Parameter setup

if nargin<2, featConfig = 'color'; end
if nargin<3, newBeta = []; end
if nargin<4, featOnly = false; end
if ndims(im)==2, featConfig = 'gray'; end
if strcmp(featConfig,'gray')&&(ndims(im)==3), im = rgb2gray(im); end
if nargin<7, [dlc, drc, dlr] = histGradFeatures(im); end
[h,w,norient,nscales,nchannels] = size(dlc);
if strcmp(featConfig,'gray') && (nargin >= 7) && (nchannels>2)
    dlc = dlc(:,:,:,:,[1,4]);
    drc = drc(:,:,:,:,[1,4]);
    dlr = dlr(:,:,:,:,[1,4]);
end
if nargin < 8
    sPb = [];
else
    sPb = repmat(sPb,[1,1,1,nscales]);  % spectral feature
end

%% Construct feature array

switch featConfig
    case 'spectral'
        if isempty(sPb), error('you have not provided spectral feature');end
        f = cat(5, ones(h,w,norient,nscales,'single'), dlc, drc , dlr,sPb);
        b = [-6.51802316701700;-0.493295831296959;1.83540354693975;1.21578765355725;3.88698064420630;0.117687234039617;0.955056692421465;0.996236738670149;4.06113533145988;-0.859428201196241;-1.15913172258885;-0.766825416349586;-2.45020647533604;0.00741091440699249;];
    case {'color'}
        f = cat(5, ones(h,w,norient,nscales,'single'), dlc, drc, dlr);
        b = [-6.54001569881004;-0.317591380311371;1.87794646755063;1.28269948188804;3.99614300468390;0.289616770646934;0.957724691090984;1.09607455266995;4.28043755023855;-1.01656506721445;-1.19665557212154;-0.860407782628580;-2.59753117106842;];
    case {'gray'}
        f = cat(5, ones(h,w,norient,nscales,'single'), dlc, drc,dlr);
        b = [-6.44514046390470;1.20879203131223;4.37178250414567;1.48029354631580;4.33415914637308;-2.24317741581830;-2.86022220830624;];
    case 'no_texture'
        f = cat(5,ones(h,w,norient,nscales,'single'),dlc(:,:,:,:,1:3),...
            drc(:,:,:,:,1:3),dlr(:,:,:,:,1:3));
        b = [-5.73872428102262;1.88976506325115;1.87838349589835;2.07200591874450;2.44810248820298;1.19908797114823;1.75638250148779;-2.34004274862761;-1.18088564539031;-1.44830602726111;];
    otherwise
        error('Illegal value for string variable "featConfig"\n')
end

clear dlc dlr drc bfeat

if ~isempty(newBeta)
    if length(newBeta)==size(f,5)
        b = newBeta;
    else
        error('Size of weight vector does not agree with feature congiguration');
    end
end
                
%% Evaluate dense symmetry response

if featOnly
    sup = []; p_bags = []; orientMap = [];
else
    % Calculate probabilities
    K = size(f,5);
    reshapedf = reshape(f,[h*w*norient*nscales,K]);
    pb = 1./(1+exp(-reshapedf*b));
    pb = reshape(pb, [h,w,norient,nscales]);
    [~,orientMap] = max(max(pb,[],4),[],3); % integer index orientation map
    thetas= (0:7)*pi/8;
    scales = [4:2:14, 16:4:28, 32:8:48];
    orientMap = thetas(orientMap);
    % Savitsky-Golay filtering
    fittedPb = zeros(size(pb),'single');
    for i = 1:norient
      for j = 1:nscales
          a = fitparab2(pb(:,:,i,j),scales(j),scales(j),thetas(i));
          fittedPb(:,:,i,j) = a;
      end
    end
    % Non-maximum suppression and noisy-or response
    p_inst = reshape(fittedPb,[h*w,norient*nscales]);
    p_inst = max(0,min(p_inst,1));
    logp_bags = log(1-p_inst+eps);
    logp_bags = sum(logp_bags,2);
    p_bags = reshape(1-exp(logp_bags),[h,w]);
    sup = nonmax(p_bags,orientMap);
    % mask out 1-pixel border where nonmax suppression fails
    sup(1,:) = 0; sup(end,:) = 0; sup(:,1) = 0; sup(:,end) = 0;
end
