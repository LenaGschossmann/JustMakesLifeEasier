function [regioMask, xyz_range, threshold, pxLim] = regiongrowing_3d_Lena(I,seedPoint, threshold, pxLim, maxDist, fixedMean)
% adapted from:
% https://de.mathworks.com/matlabcentral/fileexchange/19084-region-growing

% This function performs "region growing" in an image from a specified
% seedpoint (x,y)
%
% J = regiongrowing(I,x,y,t) 
% 
% I : input image 3D
% J : logical output image of region
% x,y : the position of the seedpoint (if not given uses function getpts)
% t : maximum intensity distance (defaults to 0.2)
%
% The region is iteratively grown by comparing all unallocated neighbouring pixels to the region. 
% The difference between a pixel's intensity value and the region's mean, 
% is used as a measure of similarity. The pixel with the smallest difference 
% measured this way is allocated to the respective region. 
% This process stops when the intensity difference between region mean and
% new pixel become larger than a certain treshold (t)
%
% Example:
%
% I = im2double(imread('medtest.png'));
% x=198; y=359;
% J = regiongrowing(I,x,y,0.2); 
% figure, imshow(I+J);
%
% Author: D. Kroon, University of Twente
c = seedPoint(1); r = seedPoint(2); z = seedPoint(3);
J = false(size(I)); % Output 
I = double(I);
Isizes = [size(I,2) size(I,1) size(I,3)]; % Dimensions of input image: x, y, z
reg_mean = I(r,c,z); % The mean of the segmented region

if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('pxLim','var')==0), pxLim=numel(I); end
if(exist('maxDist','var')==0) || isempty(maxDist), maxDist=max(Isizes); end
if(exist('fixedMean','var')==0), fixedMean = false; end
if(exist('threshold','var')==0), threshold = 2; end

% Free memory to store neighbours of the (segmented) region
reg_free = 100000;
reg_list = zeros(reg_free,4);
reg_pos = 1;

% Neighbor locations (footprint)
% neigb= ...
%     [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; -1 1 0; 1 1 0; -1 -1 0; 1 1 0;
%     0 0 -1; -1 0 -1; 1 0 -1; 0 -1 -1; 0 1 -1; -1 1 -1; 1 1 -1; -1 -1 -1; 1 1 -1;
%     0 0 1; -1 0 1; 1 0 1; 0 -1 1; 0 1 1; -1 1 1; 1 1 1; -1 -1 1; 1 1 1];
neigb = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];

% Start regiogrowing until distance between region and posible new pixels become
% higher than a certain treshold
old_neigb = [r,c,z,I(r,c,z)];
num_old = 1;
cont = 1;

while cont==1 && reg_pos-1<pxLim
    % Find new neighboring pixel
    new_neigb = zeros(num_old*size(neigb, 1),4);
    neigb_pos = 1;
    for k = 1:num_old
        for j = 1:size(neigb,1)
            rn = old_neigb(k,1) + neigb(j,1);
            cn = old_neigb(k,2) + neigb(j,2);
            zn = old_neigb(k,3) + neigb(j,3);
            
            % Check if neighbour is inside or outside the image and within
            % specified max distance
            ins = (cn>=1)&&(rn>=1)&&(zn>=1)&&...
                (cn<=Isizes(1))&&(rn<=Isizes(2))&&(zn<=Isizes(3)) && ...
                sqrt(((cn-seedPoint(1))^2 + (rn-seedPoint(2))^2)+(zn-seedPoint(3))^2) <= maxDist;
            
            if (ins&&(J(rn,cn,zn)==0)) 
                new_neigb(neigb_pos,:) = [rn,cn,zn,I(rn,cn,zn)];
                J(rn,cn,zn)=true;
                neigb_pos = neigb_pos+1;
            end
        end
    end
    new_neigb(new_neigb(:,1)==0,:) = [];
    
    if ~isempty(new_neigb)
        % Add all neighbours within a specified range from mean to region
        dist = abs(new_neigb(:,4)-reg_mean);
        [indexVal, ~] = find(dist < threshold);
        reg_list(reg_pos:reg_pos+numel(indexVal)-1,:) = new_neigb(indexVal,:);
        % calculate region mean
        if ~fixedMean
            reg_mean= (reg_mean*(reg_pos-1) + sum(new_neigb(indexVal,4)))/(reg_pos-1+numel(indexVal));
        end
        reg_pos = reg_pos+numel(indexVal);
%         J(new_neigb(indexVal,1:3)) = 2;
%         [m,i] = find(J==2); numel(m)
        old_neigb = new_neigb(indexVal,:);
        num_old = size(old_neigb,1);       
    else
        cont = 0;
    end
end

% Return the segmented area as logical matrix
regioMask = J;
for cnt = 1:reg_pos-1
    regioMask(reg_list(cnt,1),reg_list(cnt,2),reg_list(cnt,3)) = true;
end

% Return col, row, depth range
xyz_range = zeros(3,2);
xyz_range(1,:) = [find(sum(sum(J,2),3), 1, 'first') find(sum(sum(J,2),3), 1, 'last')];
xyz_range(2,:) = [find(sum(sum(J,1),3), 1, 'first') find(sum(sum(J,1),3), 1, 'last')];
xyz_range(3,:) = [find(sum(sum(J,1),2), 1, 'first') find(sum(sum(J,1),2), 1, 'last')];
end