% Read precision-recall scores from text files and plot precision-recall
% curves. 
% 
%   plotPrecisionRecall(plotGrouped)
% 
%   If plotGrouped = false, then this function plots the pr curves
%   comparing our algorithm trained with different feature sets and other
%   methods. If plotGrouped = false, it plots pr curves comparing results
%   before and after grouping using [1].
% 
%   [1]: I. Kokkinos, Highly Accurate Boundary Detection and Grouping,
%        Proc. IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), 2010
% 
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% Last update: August 2013

function plotPrecisionRecall(plotGrouped)

if nargin < 1, plotGrouped = false; end

% Load score files
scores = loadScores;

if plotGrouped
    % Grouping 50 contours
    h1 = setupFigure;
    [phColor              ,fmColor]              = plotpr(scores,'color'      ,'b',1);
    [phSpectral           ,fmSpectral]           = plotpr(scores,'spectral'   ,'k',1);
    [phColorGrouped50     ,fmColorGrouped50]     = plotpr(scores,'cgrouped50' ,'r',1);
    [phSpectralGrouped50  ,fmSpectralGrouped50]  = plotpr(scores,'sgrouped50' ,'g',1);
    legend([phColor, phColorGrouped50, phSpectral, phSpectralGrouped50],...
        'Location','SouthWest',...
        sprintf('CBT           : F=%4.3f'  ,fmColor),...
        sprintf('CBT+G_{50}   : F=%4.3f'  ,fmColorGrouped50),...    
        sprintf('CBTS          : F=%4.3f'  ,fmSpectral),...
        sprintf('CBTS+G_{50}  : F=%4.3f'  ,fmSpectralGrouped50));
    
    % Grouping 100 contours
    h2 = setupFigure;
    [phColor              ,fmColor]              = plotpr(scores,'color'      ,'b',1);
    [phSpectral           ,fmSpectral]           = plotpr(scores,'spectral'   ,'k',1);
    [phColorGrouped100    ,fmColorGrouped100]    = plotpr(scores,'cgrouped100','r',1);
    [phSpectralGrouped100 ,fmSpectralGrouped100] = plotpr(scores,'sgrouped100','g',1);
    legend([phColor, phColorGrouped100, phSpectral, phSpectralGrouped100],...
        'Location','SouthWest',...
        sprintf('CBT            : F=%4.3f'  ,fmColor),...
        sprintf('CBT+G_{100}   : F=%4.3f'  ,fmColorGrouped100),...    
        sprintf('CBTS           : F=%4.3f'  ,fmSpectral),...
        sprintf('CBTS+G_{100}  : F=%4.3f'  ,fmSpectralGrouped100));
    print(h1,'-depsc2',fullfile(scores.self,'pr_grouped50'))
    print(h2,'-depsc2',fullfile(scores.self,'pr_grouped100'))
    close(h1); close(h2);
else
    h = setupFigure;
    [phColor,     fmColor]     = plotpr(scores,'color'     ,'b');
    [phGray,      fmGray]      = plotpr(scores,'gray'      ,'m');
    [phNoTexture, fmNoTexture] = plotpr(scores,'noTexture' ,'g');
    [phSpectral,  fmSpectral]  = plotpr(scores,'spectral'  ,'r');
    [phLind,      fmLind]      = plotpr(scores,'lind'      ,'k');
    [phLevin,     fmLevin]     = plotpr(scores,'levin'     ,'c');
    legend([phSpectral,phColor,phGray,phNoTexture,phLind,phLevin],...
        'Location','SouthWest',...
        sprintf('CBTS     : F=%4.3f'  ,fmSpectral),...    % spectral
        sprintf('CBT      : F=%4.3f'  ,fmColor),...       % color
        sprintf('BT       : F=%4.3f'  ,fmGray),...        % gray
        sprintf('CB       : F=%4.3f'  ,fmNoTexture),...   % no texture
        sprintf('Lindeberg: F=%4.3f'  ,fmLind),...        % Lindeberg
        sprintf('Levinshtein: F=%4.3f',fmLevin));         % Levinshtein
    print(h,'-depsc2',fullfile(scores.self,'pr'))
    close(h)
end


% --- Plot precision-recall curve for a feature combination ---------------
function [h, F] = plotpr(scores,features,color,lineWidth,markerSize)
% h: plot handle
% F: maximum F-measure

assert(ischar(features));
assert(ischar(color));
if nargin < 4, lineWidth  = 2; end
if nargin < 5, markerSize = 8; end

pr    = dlmread(scores.(features).pr);
score = dlmread(scores.(features).final);
F     = score(4);
if ~strcmp(features,'levin')
    plot(pr(:,2),pr(:,3),[color '-'],'LineWidth',lineWidth);
end
h = plot(score(2),score(3),[color 'o'],'MarkerFaceColor',color,...
    'MarkerSize',markerSize); % add marker for best F-measure value

% --- Setup figure for plotting -------------------------------------------
function h = setupFigure()

h = figure; clf;
hold on; box on; grid on;
set(gca,'Fontsize',14);
set(gca,'XTick',[0 .25 .5 .75]);
set(gca,'YTick',[0 .25 .5 .75]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Recall'); ylabel('Precision');
axis([0 .75 0 .75]); 

% Plot iso-contours
[r_gt,p_gt] = meshgrid(0:0.01:1,0:0.01:1);
f_gt        = fmeasure(r_gt,p_gt);
[C,cl]      = contour(0:0.01:1,0:0.01:1,f_gt,.1:.1:.8);
clabel(C,cl)


% % Read precision-recall scores from text files and plot precision-recall
% % curves. 
% % 
% %   plotPrecisionRecall(plotGrouped)
% % 
% %   If plotGrouped = false, then this function plots the pr curves
% %   comparing our algorithm trained with different feature sets and other
% %   methods. If plotGrouped = false, it plots pr curves comparing results
% %   before and after grouping using [1].
% % 
% %   [1]: I. Kokkinos, Highly Accurate Boundary Detection and Grouping,
% %        Proc. IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), 2010
% % 
% % Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% % Last update: August 2013
% function plotPrecisionRecall(plotGrouped)
% 
% if nargin < 1, plotGrouped = false; end
% 
% % Load score files
% scores = loadScores;
% 
% if plotGrouped
%     % Grouping 50 contours
%     h1 = setupFigure;
%     [phColor              ,fmColor]              = plotpr(scores,'color'      ,'b',1);
%     [phSpectral           ,fmSpectral]           = plotpr(scores,'spectral'   ,'k',1);
%     [phColorGrouped50     ,fmColorGrouped50]     = plotpr(scores,'cgrouped50' ,'r',1);
%     [phSpectralGrouped50  ,fmSpectralGrouped50]  = plotpr(scores,'sgrouped50' ,'g',1);
%     legend([phColor, phColorGrouped50, phSpectral, phSpectralGrouped50],...
%         'Location','SouthWest',...
%         sprintf('CBT           : F=%4.3f'  ,fmColor),...
%         sprintf('CBT+G_{50}   : F=%4.3f'  ,fmColorGrouped50),...    
%         sprintf('CBTS          : F=%4.3f'  ,fmSpectral),...
%         sprintf('CBTS+G_{50}  : F=%4.3f'  ,fmSpectralGrouped50));
%     
%     % Grouping 100 contours
%     h2 = setupFigure;
%     [phColor              ,fmColor]              = plotpr(scores,'color'      ,'b',1);
%     [phSpectral           ,fmSpectral]           = plotpr(scores,'spectral'   ,'k',1);
%     [phColorGrouped100    ,fmColorGrouped100]    = plotpr(scores,'cgrouped100','r',1);
%     [phSpectralGrouped100 ,fmSpectralGrouped100] = plotpr(scores,'sgrouped100','g',1);
%     legend([phColor, phColorGrouped100, phSpectral, phSpectralGrouped100],...
%         'Location','SouthWest',...
%         sprintf('CBT            : F=%4.3f'  ,fmColor),...
%         sprintf('CBT+G_{100}   : F=%4.3f'  ,fmColorGrouped100),...    
%         sprintf('CBTS           : F=%4.3f'  ,fmSpectral),...
%         sprintf('CBTS+G_{100}  : F=%4.3f'  ,fmSpectralGrouped100));
%     print(h1,'-depsc2',fullfile(scores.self,'pr_grouped50'))
%     print(h2,'-depsc2',fullfile(scores.self,'pr_grouped100'))
%     close(h1); close(h2);
% else
%     h = setupFigure;
%     [phColor,     fmColor]     = plotpr(scores,'color'     ,'b');
%     [phGray,      fmGray]      = plotpr(scores,'gray'      ,'m');
%     [phNoTexture, fmNoTexture] = plotpr(scores,'noTexture' ,'g');
%     [phSpectral,  fmSpectral]  = plotpr(scores,'spectral'  ,'r');
%     [phLind,      fmLind]      = plotpr(scores,'lind'      ,'k');
%     [phLevin,     fmLevin]     = plotpr(scores,'levin'     ,'c');
%     legend([phSpectral,phColor,phGray,phNoTexture,phLind,phLevin],...
%         'Location','SouthWest',...
%         sprintf('CBTS     : F=%4.3f'  ,fmSpectral),...    % spectral
%         sprintf('CBT      : F=%4.3f'  ,fmColor),...       % color
%         sprintf('BT       : F=%4.3f'  ,fmGray),...        % gray
%         sprintf('CB       : F=%4.3f'  ,fmNoTexture),...   % no texture
%         sprintf('Lindeberg: F=%4.3f'  ,fmLind),...        % Lindeberg
%         sprintf('Levinshtein: F=%4.3f',fmLevin));         % Levinshtein
%     print(h,'-depsc2',fullfile(scores.self,'pr'))
%     close(h)
% end
% 
% 
% % --- Plot precision-recall curve for a feature combination ---------------
% function [h, F] = plotpr(scores,features,color,lineWidth,markerSize)
% % h: plot handle
% % F: maximum F-measure
% 
% assert(ischar(features));
% assert(ischar(color));
% if nargin < 4, lineWidth  = 2; end
% if nargin < 5, markerSize = 8; end
% 
% pr    = dlmread(scores.(features).pr);
% score = dlmread(scores.(features).final);
% F     = score(4);
% if ~strcmp(features,'levin')
%     plot(pr(:,2),pr(:,3),[color '-'],'LineWidth',lineWidth);
% end
% h = plot(score(2),score(3),[color 'o'],'MarkerFaceColor',color,...
%     'MarkerSize',markerSize); % add marker for best F-measure value
% 
% % --- Setup figure for plotting -------------------------------------------
% function h = setupFigure()
% 
% h = figure; clf;
% hold on; box on; grid on;
% set(gca,'Fontsize',14);
% set(gca,'XTick',[0 .25 .5 .75]);
% set(gca,'YTick',[0 .25 .5 .75]);
% set(gca,'XGrid','on');
% set(gca,'YGrid','on');
% xlabel('Recall'); ylabel('Precision');
% axis([0 .75 0 .75]); 
% 
% % Plot iso-contours
% [r_gt,p_gt] = meshgrid(0:0.01:1,0:0.01:1);
% f_gt        = fmeasure(r_gt,p_gt);
% [C,cl]      = contour(0:0.01:1,0:0.01:1,f_gt,.1:.1:.8);
% clabel(C,cl)
% 
% 
