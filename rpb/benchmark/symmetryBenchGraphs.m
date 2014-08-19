function symmetryBenchGraphs(scorePath)
% function symmetryBenchGraphs(scorePath)
%
% Create graphs, after symmetryBench has been run.
%
% See also symmetryBench.
%
% Stavros Tsogkas <stavros.tsogkas@ecp.fr>
% November 2011

%% Get file names
% MIL color
fname = fullfile(scorePath,'mil_color_score.txt');
score_mil_color = dlmread(fname); % thresh,r,p,f
% MIL gray
fname = fullfile(scorePath,'mil_gray_score.txt');
score_mil_gray = dlmread(fname); % thresh,r,p,f
% MIL no texture
fname = fullfile(scorePath,'mil_notext_score.txt');
score_mil_notext = dlmread(fname); % thresh,r,p,f
% MIL spectral
fname = fullfile(scorePath,'mil_spectral_score.txt');
score_mil_spectral = dlmread(fname); % thresh,r,p,f
% Lindeberg
fname = fullfile(scorePath,'lind_score.txt');
score_lind = dlmread(fname); % thresh,r,p,f
% Levinshtein
fname = fullfile(scorePath,'levin_score.txt');
score_levin = dlmread(fname); % thresh,r,p,f
% % Ncut score
% fname = fullfile(scorePath,'ncut_score.txt');
% score_ncut = dlmread(fname); % thresh,r,p,f
% 
%% Create the overall PR grap
h = figure; clf;
figure(h); hold on;
box on; grid on;
set(gca,'Fontsize',12);
set(gca,'XTick',[0 .25 .5 .75 1]);
set(gca,'YTick',[0 .25 .5 .75 1]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Recall'); ylabel('Precision');
% title('Precision-Recall curves for MIL and Lindeberg');
axis square; 
% axis([0 0.75 0 0.75]);
% MIL, color configuration figure
fname = fullfile(scorePath,'mil_color_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
ph_mil_color = plot(pr(:,2),pr(:,3),'b-','LineWidth',3);
% plot(score_mil_color(2),score_mil_color(3),'ro','MarkerFaceColor','r','MarkerSize',10); % add marker for best F-measure value
% MIL, gray configuration figure
fname = fullfile(scorePath,'mil_gray_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
hold on;
ph_mil_gray = plot(pr(:,2),pr(:,3),'m-','LineWidth',3);
% plot(score_mil_gray(2),score_mil_gray(3),'bo','MarkerFaceColor','b','MarkerSize',10);
% MIL, no texture configuration figure
fname = fullfile(scorePath,'mil_notext_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
hold on;
ph_mil_notext = plot(pr(:,2),pr(:,3),'g-','LineWidth',3);
% plot(score_mil_notext(2),score_mil_notext(3),'mo','MarkerFaceColor','m','MarkerSize',10);
% MIL, spectral configuration figure
fname = fullfile(scorePath,'mil_spectral_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
hold on;
ph_mil_spectral = plot(pr(:,2),pr(:,3),'r-','LineWidth',3);
% % plot(score_mil_spectral(2),score_mil_spectral(3),'go','MarkerFaceColor','g','MarkerSize',10);
% Lindeberg detector
fname = fullfile(scorePath,'lind_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
hold on;
ph_lind = plot(pr(:,2),pr(:,3),'k-','LineWidth',3);
% % plot(score_lind(2),score_lind(3),'ko','MarkerFaceColor','k','MarkerSize',10);
% Levinshtein detector
fname = fullfile(scorePath,'levin_pr.txt');
pr = dlmread(fname); % thresh,r,p,f
hold on;
% ph_levin = plot(pr(:,2),pr(:,3),'b-','LineWidth',3);
ph_levin = plot(score_levin(2),score_levin(3),'ro','MarkerFaceColor','r','MarkerSize',10);
% Ground truth score 
[r_gt,p_gt] = meshgrid(0:0.01:1,0:0.01:1);
f_gt = fmeasure(r_gt,p_gt);
[C,cl] = contour(0:0.01:1,0:0.01:1,f_gt,[0.3,0.4,0.5,0.6,0.73,0.8]);
clabel(C,cl)
hold on;

%% Add legend 
legend([ph_mil_spectral,ph_mil_color,ph_mil_gray,ph_mil_notext,ph_lind,ph_levin],...
    'Location','SouthWest',...
    ...% MIL, spectral configuration figure
    sprintf('CG+BG+TG+SP: F=%4.3f',score_mil_spectral(4)),...
    ...% MIL, color configuration figure
    sprintf('CG+BG+TG: F=%4.3f',score_mil_color(4)),...
    ...% MIL, gray configuration figure
    sprintf('BG+TG: F=%4.3f',score_mil_gray(4)),...
    ...% MIL, no texture configuration figure
    sprintf('CG+BG: F=%4.3f',score_mil_notext(4)),...
    ...% Lindeberg detector
    sprintf('Lindeberg: F=%4.3f',score_lind(4)),...
    ...% Lindeberg detector
    sprintf('Levinshtein: F=%4.3f',score_levin(4)));

% Print figures
print(h,'-depsc2',fullfile(scorePath,'pr.eps'));
print(h,'-djpeg95','-r36',fullfile(scorePath,'pr_half.jpg'));
print(h,'-djpeg95','-r0',fullfile(scorePath,'pr_full.jpg'));
close(h)

%% Not necessarily used
function prplot(h,r,p,ti)
figure(h); 
plot(r,p,'ko-');
box on;
grid on;
set(gca,'Fontsize',12);
set(gca,'XTick',[0 .25 .5 .75 1]);
set(gca,'YTick',[0 .25 .5 .75 1]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Recall');
ylabel('Precision');
title(ti);
axis square;
axis([0 1 0 1]);

