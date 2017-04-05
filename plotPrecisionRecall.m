function fig = plotPrecisionRecall(models,dataset,set)
% PLOTPRECISIONRECALL Plot precision-recall curves and print figure. 
% 
%   PLOTPRECISIONRECALL(models) where models can be one of the following:
% 
%   i) a cell array with the names or paths of the mat-files that contain 
%      the models to be compaired and their respective statistics.
%   ii)a cell array whose elements are structs with the trained models and
%      their statistics.
%       
%   See also: testSPB
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: March 2017 (added support for 'human' performance)

if nargin < 2, dataset = 'BMAX500'; end
if nargin < 3, set = 'val'; end

paths = setPaths();
models = loadModels(models, paths, dataset, set);
fig = setupFigure;
hpr = zeros(numel(models),1); % pr plot handles
f   = zeros(numel(models),1); % best f-measure 
txt = cell(numel(models),1);  % legend text 
for m=1:numel(models)
    [hpr(m), f(m), txt{m}] = plotpr(models{m}); 
end

% Sort by f-score, create legend and print figure
[~,inds] = sort(f,'descend');
legend(hpr(inds), 'Location','SouthWest', txt(inds)); 
if ~isdir(paths.spbmil.plots), mkdir(paths.spbmil.plots); end
export_fig('-transparent',fullfile(paths.spbmil.plots, 'pr.pdf'),fig)


% -------------------------------------------------------------------------
function [h, bestF, txt] = plotpr(model,lineWidth,markerSize)
% -------------------------------------------------------------------------
if nargin < 2, lineWidth  = 4; end
if nargin < 3, markerSize = 16; end

% Compute P,R,F
P = sum(model.stats.cntP,1) ./ max(eps, sum(model.stats.sumP,1));
R = sum(model.stats.cntR,1) ./ max(eps, sum(model.stats.sumR,1));
bestP = model.stats.odsP;
bestR = model.stats.odsR;
bestF = model.stats.odsF;

% Plot P-R curve
color = model2color(model);
if ~ismember(model.name, {'levinstein','amat','human'})
    plot(P,R,[color '-'],'LineWidth',lineWidth);
end
% Add marker for best F-measure value
h = plot(bestP,bestR,[color 'o'],'MarkerFaceColor',color,'MarkerSize',markerSize); 

% Create legend text 
txt = model2legend(model,bestF);


% -------------------------------------------------------------------------
function h = setupFigure()
% -------------------------------------------------------------------------
h = figure; clf;
hold on; box on; grid on;
set(gca,'Fontsize',24);
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',0:0.2:0.8);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Recall'); ylabel('Precision');
axis([0 1 0 0.8]); 

% Plot iso-contours
[r_gt,p_gt] = meshgrid(0:0.01:1,0:0.01:1);
f_gt        = fmeasure(r_gt,p_gt);
[C,cl]      = contour(0:0.01:1,0:0.01:1,f_gt,.1:.1:.8);
clabel(C,cl)

% -------------------------------------------------------------------------
function F = fmeasure(P,R), F = 2 .* P .* R ./ max(eps, P+R);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function c = model2color(model)
% -------------------------------------------------------------------------
if isfield(model,'opts') && isfield(model.opts, 'featureSet')
    switch model.opts.featureSet
        case 'color'
            c = 'c';
        case 'gray'
            c = 'g';
        case {'no-texture','no_texture'}
            c = 'k';
        case 'spectral' 
            c = 'm';
    end
elseif strcmpi(model.name, 'human')
    c = 'r';
elseif strcmpi(model.name, 'amat')
    c = 'b';
elseif strcmpi(model.name, 'deepskel')
    c = 'k';
else error('Unknown model type')
end

% -------------------------------------------------------------------------
function t = model2legend(model,f)
% -------------------------------------------------------------------------
if isfield(model,'opts') && isfield(model.opts, 'featureSet')
    t = sprintf('MIL-%s: F=%.2f', model.opts.featureSet, f);
elseif strcmpi(model.name, 'human')
    t = sprintf('Human: F=%.2f', f);
elseif strcmpi(model.name, 'amat')
    t = sprintf('AMAT: F=%.2f', f);
elseif strcmpi(model.name, 'deepskel')
    t = sprintf('DeepSkel: F=%.2f', f);
else error('Unknown model type')    
end

% -------------------------------------------------------------------------
function models = loadModels(models,paths,dataset,set)
% -------------------------------------------------------------------------
for m=1:numel(models)
    if ischar(models{m})
        if exist(models{m}, 'file') % absolute path
            tmp = load(models{m});
        elseif exist([models{m}, '.mat'],'file')
            tmp = load([models{m}, '.mat']);
        % partial path
        elseif exist(fullfile(paths.spbmil.models, models{m}), 'file') 
            tmp = load(fullfile(paths.models, models{m}));
        elseif exist(fullfile(paths.spbmil.models, [models{m} '.mat']), 'file')
            tmp = load(fullfile(paths.models, [models{m} '.mat']));
        else
            error([models{m} ' was not found'])
        end
        models{m} = tmp.model; % this depends on the name of the results struct
    end
    % Create convenient stats field for requested dataset and subset
    if isfield(models{m}, dataset) && isfield(models{m}.(dataset), set)
        models{m}.stats = models{m}.(dataset).(set).stats;
    else
        error(['Model has not been evaluated on ' dataset ' ' set])
    end
end

