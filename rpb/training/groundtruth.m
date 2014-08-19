function gt = groundtruth(imagelist,mode,groundpath,start)

% Creates groundtruth used for symmetry detector.
% 
% imagelist:    list with image iids (for Berkeley Segmentation Dataset).
% groundpath:   path where the ground truths are going to be saved.
% mode:         'auto' uses skeltonization of image segments based on
%               BSDS300. 'manual' lets user draw a ground truth map for
%               each image.
% start:        used to continue groundtruth creation from a specific image

if nargin<2, mode = 'auto'; end
if nargin<3, groundpath = []; end
if nargin<4, start = 1; end


gt = cell(numel(imagelist),1);
if strcmp(mode,'manual')
    for i=start:numel(imagelist)
        i
        iid = imagelist(i);
        im = imgRead(iid);
        [h w ~] = size(im);
        gt{i} = false(h,w);
        gt_manual = false(h,w);
        done = false;
        while ~done
            idx = [];
            fprintf('Processing image with iid %d\n',iid)
            fig1 = figure(1); imshow(im); title(sprintf('iid = %d', iid));
            disp('"f": adds freehand contour')
            disp('"l": adds linear segment')
            disp('"r" or "d": removes last contour')
            disp('return: ends annotation for current image')
            reply = input('Choose type for next contour: ','s');
            while ~isempty(reply)
                switch reply
                    case 'f'
                        p = imfreehand('Closed',false);
                        pos = getPosition(p);
                        pos = round(pos);
                        linpos = sub2ind([h,w],pos(:,2),pos(:,1));
                        tmp = zeros(h,w); tmp(linpos) = 1;
                        tmp = bwmorph(bwmorph(bwmorph(tmp,'dilate'),'bridge'),'thin',inf);
                        linpos = find(tmp);
                        idx = [idx; linpos];
                        reply = input('Choose type for next contour: ','s');
                    case 'l'
                        p = imline;
                        pos = getPosition(p);
                        linpos = linseg([h,w],pos);
                        idx = [idx; linpos];
                        reply = input('Choose type for next contour: ','s');
                    case {'d','r'}
                        delete(p)
                        idx = idx(1:end-numel(linpos));                    
                        reply = input('Choose type for next contour: ','s');
                    otherwise
                        reply = input(['Wrong key! Press "f" for freehand contour,"l" for straight line,\n',...
                                    '"d" to delete last contour, or return to end annotation: '],'s');
                end
            end
            gt_manual(:) = false;
            gt_manual(idx) = true;
            fig2 = figure(2); imshow(overlayBinaryImage(im,gt_manual)); title(sprintf('iid = %d', iid));
            reply = input('Accept ground truth? [y/n]: ','s');  
            while ~strcmp(reply,'y') && ~strcmp(reply,'n')
                reply = input('INVALID ANSWER!: Accept ground truth? [y/n]: ','s');
            end
            done = strcmp(reply,'y');   % if not satisfied with the ground truth, start over
            close(fig2)
        end
        gt{i} = gt_manual;
        clear tmp 
        if ~isempty(groundpath)
            save(fullfile(groundpath, sprintf('groundtruth_manual_%d', iid)),'gt_manual','iid','i','idx')
        end
        close all
    end
elseif strcmp(mode,'auto')
    format = 'color';
    skelThresh = 35;
    regionPixelThreshold = 50;  % if a region has less pixels than this threshold,
                                % then we do not compute the skeleton transform
    for i=1:numel(imagelist)
        iid = imagelist(i);
        im = imgRead(iid, format);
        figure, imshow(im); title(sprintf('iid = %d', iid));
        [h w ~] = size(im);
        [segs uids] = readSegs(format, iid);
        nUsers = length(segs);  % number of different human segmentations
        nRegions = zeros(length(segs),1); % vector that contains the number of regions
                                          % in which each user segmented the image
        ridge_gt = false(h,w,nUsers);   % ridge_gt for different human segmentations 
        ridgeUnion = false(h,w);

        for iu=1:nUsers % for each subject segmentation do... 
            nRegions(iu) = max(max(segs{iu}));   % number of regions in which image is segmented
            % for each region of segmentation do...
            for ireg = 1:nRegions(iu)  % ireg==1 usually corresponds to background so we ignore it
                seg = (segs{iu}==ireg); % examine each region seperately
                regidx = find(seg);
                if (length(regidx) < regionPixelThreshold), continue; end;
                if strcmp(format,'gray')
                    region = repmat(im, [1 1 3]); % turn the gray image to RGB
                    tmp = im; tmp(regidx) = 1;
                    region(:,:,1) = tmp;   %  mark the current region pixels as red
                    tmp(regidx) = 0;
                    region(:,:,2) = tmp;
                    region(:,:,3) = tmp;
                else
                    region = im;    
                    tmp = region(:,:,1); tmp(regidx) = 1;
                    region(:,:,1) = tmp;   
                    tmp = region(:,:,2); tmp(regidx) = 0;
                    region(:,:,2) = tmp;
                    tmp = region(:,:,3); tmp(regidx) = 0;
                    region(:,:,3) = tmp;
                end

                fig = figure; imshow(region);
                title(sprintf('Region %d/%d (user %d/%d)', ireg, nRegions(iu), iu, nUsers));

                fprintf(2, 'Do you want to add this region to the skeleton transform?\n');
                reply = input('enter: skips region without adding, "y": adds region, "n" skips current image): ', 's');
                close(fig); 
                while ~strcmp(reply,'') && ~strcmp(reply,'n') && ~strcmp(reply,'y')
                    reply = input('INVALID INPUT!: enter: skips region without adding, "y": adds region, "n" skips current segmentation): ', 's');
                end
                switch reply
                    case '', continue;
                    case 'n', break;
                    case 'y',
                        % skeletonization using AAFM in Matlab (Alex Telea-Nicholas Howe)
                        skel = logical(bwmorph(skeleton(seg) > skelThresh,'skel',Inf));
                        fig = figure('WindowStyle','docked'); imshow(skel)
                        % union of skeletons for different regions
                        ridge_gt(:,:,iu) = ridge_gt(:,:,iu) | skel;
                end
            end
            ridgeUnion = ridgeUnion | ridge_gt(:,:,iu);
        end
        fig = figure; imshow(ridgeUnion);
        gt{i} = ridgeUnion;
        if ~isempty(groundpath)
            save(fullfile(groundpath, sprintf('groundtruth_seg_%d', iid)),'nUsers','nRegions','iid','ridge_gt','ridgeUnion','i')
        end
        close all
    end
else
    error(['Choose "auto" mode for annotation using the Berkeley Segmentation',...
        ' Dataset BSDS300 or "manual" mode for manual annotations']);
end
                
                
        
