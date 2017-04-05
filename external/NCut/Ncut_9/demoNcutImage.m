% function demoNcutImage;
% demoNcutImage
% 
% demo for NcutImage
% also initialize matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
% close all
disp('Ncut Image Segmentation demo');

%% read image, change color image to brightness image, resize to 160x160
% I = imread_ncut('jpg_images/3.jpg',160,160);
% [h w] = size(I);
I = imread_ncut('~/symmetry_detector/images/BSDS300/images/test/3096.jpg');   
% I = imread_ncut('C:\Users\MAS User\Desktop\Thesis\images\BSDS300\images\train\22090.jpg',321,481);

%% display the image
% figure(1);clf; imagesc(I);colormap(gray);axis off;
% disp('This is the input image to segment, press Enter to continue...');
% pause;

%% compute the edges imageEdges, the similarity matrix W based on
%% Intervening Contours, the Ncut eigenvectors and discrete segmentation
nbSegments = 7;
disp('computing Ncut eigenvectors ...');
tic;
[SegLabel,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(I,nbSegments);
disp(['The computation took ' num2str(toc) ' seconds on the ' num2str(size(I,1)) 'x' num2str(size(I,2)) ' image']);

%% display the edges
% figure(2);clf; imagesc(imageEdges); axis off
% disp('This is the edges computed, press Enter to continue...');
% pause;

%% display the segmentation
% figure(3);clf
% bw = edge(SegLabel,0.01);
% J1=showmask(I,imdilate(bw,ones(2,2))); imagesc(J1);axis off
% disp('This is the segmentation, press Enter to continue...');
% pause;

%% display Ncut eigenvectors
figure;clf;set(gcf,'Position',[100,500,200*(nbSegments+1),200]);
[nr,nc,nb] = size(I);
for i=1:nbSegments
    subplot(1,nbSegments,i);
    imagesc(reshape(NcutEigenvectors(:,i) , nr,nc));axis('image');axis off;
end
disp('This is the Ncut eigenvectors...');
disp('The demo is finished.');

%% display W affinities for a chosen pixel and a specific radius

% binaryEdges = imageEdges>80; % note that edge strength is in [0,100]
% distMap = bwdist(binaryEdges);
% idx = find(distMap<5);
% heatMap = reshape(W(idx(1),:),h,w);
% figure, imagesc(binaryEdges) 
% hold on
% imagesc(heatMap)
% 
% 
