function [valsum valdiff]= boundValFeat(pb, btheta, scales)

% function [valsum valdiff]= boundValFeat(pb, btheta, scales)
% 
% Function that is used to collect the feature of the boundary responce for
% all scales and orientations.
% 
% INPUT
% pb:       boundary image used for feature collection
% btheta:   the orientation vector for the boundary map in radians.
% scales:   Vector of radii for the circular disk. Scale determines at 
%           what distance from the pixel we look for boundary.
%
% OUTPUT
% valsum:    feature of sum of bilateral boundary responses
% valdiff:   feature of difference of bilateral boundary responses


% distance from boundary confirmation (we check the existence of
% boundary in an area that is part of a circular ring at orientation
% vertical to the orientation of the symmetry in pixel (x,y).
[h w] = size(pb);
if nargin<3, scales = [4:2:14, 16:4:28, 32:8:48]; end
scales = floor(scales);
thetas = (0:7)*pi/8; % orientation vector
nscales = length(scales);
nthetas = length(thetas);
valsum = zeros(h,w,nthetas,nscales,'single');
valdiff = zeros(h,w,nthetas,nscales,'single');

for isc=1:nscales
    wl = scales(isc); 
    [u,v] = meshgrid(-wl:wl,-wl:wl);
    gamma = atan2(v,u);
    cRing = ((u.^2 + v.^2 <= ceil((1.1*wl)^2)) & (u.^2 + v.^2 >= floor((0.9*wl)^2))); % circular ring mask
    cRing = imclose(cRing,strel('disk',1));

    for it=1:nthetas
        % we need accordance of the orientation of the symmetry axis and the
        % orientation of the corresponding boundary. We allow a margin of
        % pi/norient rad angle difference (default is pi/8 rad)
        temp = single(pb);
        switch it
            case 1        
                temp(btheta>thetas(2) & btheta<thetas(nthetas)) = 0;
            case nthetas
                temp(btheta>thetas(1) & btheta<thetas((nthetas-1))) = 0;
            otherwise
                temp(btheta>thetas(it+1) | btheta<thetas(it-1)) = 0;
        end
        side = 1 + (mod(gamma-thetas(it),2*pi) < pi);    % split ring into 2 sides
        side = side.*cRing;
        lmask = (side==1);
        rmask = (side==2);
        rectside = sqrt(sum(lmask(:))/3);  % straight line approximation of the circular area length
        if (thetas(it) < pi/2)
            angleMask = ((mod(gamma,pi) <= (thetas(it) + pi/2 + 0.2))...
                    & (mod(gamma,pi) >= (thetas(it) + pi/2 - 0.2)));
            mask = cRing & angleMask;
        else
            angleMask = ((mod(gamma,pi) <= (thetas(it) - pi/2 + 0.2))...
                    & (mod(gamma,pi) >= (thetas(it) - pi/2 - 0.2)));
            mask = cRing & angleMask;
        end
        lmask = single(lmask & mask);  
        rmask = single(rmask & mask);
        lbound = conv2(temp, lmask,'same');   % response only for left boundary
        rbound = conv2(temp, rmask,'same');   % response only for right boundary
        valsum(:,:,it,isc) = (lbound + rbound)/(2*rectside);
        valdiff(:,:,it,isc) = abs(lbound - rbound)/(2*rectside);
    end
end




