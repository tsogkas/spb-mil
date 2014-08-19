function tmap = createTextonMap(im,k)
% function tmap = createTextonMap(im,k)
% 
% Create texton map of an image (Berkeley Pb code)
% 
% INPUTS
% im:   original image
% k:    number of texton labels

if nargin<2, k = 64; end

no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
fname = sprintf( ...
'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
textonData = load(fname); % defines fb,tex,tsim
if size(im,3)==3, tmapim = rgb2gray(im); else tmapim = im; end
tmap = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);
