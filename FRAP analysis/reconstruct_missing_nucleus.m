function [L_nuc,L_cyt,nuc_data,cyt_data,Nucindex] = reconstruct_missing_nucleus(mask_data, label, nuc_data, I_data, I_mask, t, ...
    idx_bleach, fr, datum, filenameshort)


% === Function: reconstruct_missing_nucleus ===
% Reconstructs and tracks a bleached nucleus across frames and extracts labeled
% nuclear and cytoplasmic regions with their properties
%
% Inputs:
%   mask_data     - Cell array containing binary masks and regionprops structures for each frame
%   label         - Cell array with labeled mask images and corresponding regionprops structures
%   I_data        - Image data from the protein channel
%   I_mask        - Image data from the nuclear channel
%   t             - Vector of time points corresponding to each frame
%   idx_bleach    - Index of the bleached nucleus to be reconstructed
%   fr            - Total number of frames to process
%   filenameshort - String with short filename identifier used for saving output images
%
% Outputs:
%   L_nuc     - Store final nuclear label image
%   L_cyt     - Store final labeled cytoplasm image
%   nuc_data  - Cell array with extracted properties for each nucleus across frames
%   cyt_data  - Cell array with extracted properties for each cytoplasmic region across frames
%   Nucindex  - Matrix with index of nucleus labels per frame 

% The function includes the following major steps:
%   1. Label consistency assignment
%   2. Estimating nucleus motion via neighbors
%   3. Reconstructing the bleached nucleus forward in time
%   4. Saving visualization of updated masks
%   5. Relabeling and tracking nuclei
%   6. Generating and assigning cytoplasm masks
%   7. Extracting and labeling nuclear and cytoplasmic regions
%   8. Saving label images with annotation
 
label2 = label;
mask3_data = mask_data;

% === Assign Consistent Labels to Regions ===
for i = 1:size(label2,1)
    
    for j = 1:size(label2{i,2},1)
        
        nmask = label2{i,1};
        pxls = label2{i, 2}(j).PixelIdxList;
        val_label = mean(nmask(pxls));
        label2{i,2}(j).Label = val_label;
    end
    
    
end


% === Collect All Unique Labels ===
uniq_lbs = [];

for i = 1:size(label2,1)
    
    temp_unq_lbs = unique(label{i,1});
    
    uniq_lbs = unique([uniq_lbs;temp_unq_lbs]);
    
    
end

% === Extract Centroid Locations per Label ===
for i = 1:fr
    
    for j = 1:size(label2{i, 2},1)
        
        label_id = label2{i, 2}(j).Label;
        cent_temp = label2{i, 2}(j).Centroid;
        centroid_lbl_locs{i,label_id} = cent_temp; 
        
    end
    
end

% === Estimate Nucleus Motion Across Frames ===
for i = 1:fr-1
   
    for j = 1:size(centroid_lbl_locs,2)
        
        if (isempty(centroid_lbl_locs{i+1,j}) ~= 1 && isempty(centroid_lbl_locs{i,j}) ~= 1)
            
            movement{i,j} = centroid_lbl_locs{i+1,j} - centroid_lbl_locs{i,j};
        end
        
    end
    
end

% === Compute Average Motion of Closest Neighbors ===
indexes = persistent_closest_neighbors(nuc_data,idx_bleach,6);
indexes = indexes';
temp_sum = cell2mat(movement(:,indexes(1)));

for i = 2:size(indexes,2)
    
    sum = temp_sum + cell2mat(movement(:,indexes(i)));
    temp_sum = sum;
    
end

mean_mov = sum./size(indexes,2);

% === Reconstruct Missing Nucleus in Each Frame ===
temp_img = zeros(size(mask_data{1, 1}));
lbl_number = idx_bleach;% label number of missing nucleus in first frame
pxl = mask_data{1, 2}(lbl_number).PixelIdxList; 
temp_img(pxl) = 1;
temp_img = imerode(temp_img,strel('disk',3));

mask2_data{1,1} = mask_data{1,1};
mask2_data{1,2} = mask_data{1,2};


for i = 2:fr
    temp_label_matrix = vertcat(label2{i, 2}.Label);
    indx_label = find(temp_label_matrix == lbl_number); 
    
    if ~isempty(indx_label)
        temp_img2 = zeros(size(mask_data{1, 1}));
        pxl2 = mask_data{i, 2}(indx_label).PixelIdxList;
        temp_img2(pxl2) = 1;
        shiftx = (mean_mov(i-1,1));
        shifty = (mean_mov(i-1,2));
        shiftedImage = fraccircshift(temp_img, [shifty, shiftx]);
        mask2_data{i,1} = imbinarize(shiftedImage+mask_data{i,1}-temp_img2);
        mask2_data{i,2} = regionprops( mask2_data{i,1},'Centroid','Area','PixelIdxList');
        temp_img = shiftedImage;
        mask_data{i,1} = mask2_data{i,1};
        mask_data{i,2} = mask2_data{i,2};
        
    else
        shiftx = (mean_mov(i-1,1));
        shifty = (mean_mov(i-1,2));
        shiftedImage = fraccircshift(temp_img, [shifty, shiftx]);
        mask2_data{i,1} = imbinarize(shiftedImage+mask_data{i,1});
        mask2_data{i,2} = regionprops( mask2_data{i,1},'Centroid','Area','PixelIdxList');
        temp_img = shiftedImage;
        mask_data{i,1} = mask2_data{i,1};
        mask_data{i,2} = mask2_data{i,2};
    end
    
end

% === Update Centroids and Generate Point Clouds ===
for fr = 1:fr
    
    x = [];
    y = [];
    
    for c = 1:size(mask_data{fr,2},1)
        
        tempx = mask_data{fr,2}(c).Centroid(1);
        tempy = mask_data{fr,2}(c).Centroid(2);
        x = [x;tempx];
        y = [y;tempy];
        loc = [x,y];
        mask_data{fr,3} = loc;
        
    end

end

for fr = 1:fr
    
    for c = 1:size(mask_data{fr,2},1)
        
        loc = [mask_data{fr,3} zeros(size(mask_data{fr,2},1),1)];
        mask_data{fr,4} = pointCloud(loc);
        
    end

end

% === Extract Perimeter for Visualization and save Updated Masks with Overlaid Boundaries ===
for fr = 1:size(mask_data,1)
   
    stats = regionprops(mask_data{fr,1},I_mask(:,:,fr),'MeanIntensity','Centroid','PixelList','PixelIdxList');

    H_u = size(I_mask(:,:,fr),1);
    W_u = size(I_mask(:,:,fr),2);
    [B_u,~,~,Idx_u] = traceobject(stats,H_u,W_u);
    for j = 1:length(Idx_u)
        stats(j).PerimIdxList = Idx_u{j};
        stats(j).PerimXYList = fliplr(B_u{j});
    end
    
    
		figure('visib','off')
		
        if ~(exist('Image results','dir') == 7)
            mkdir 'Image results'
        end
        
        if ~(exist('Image results\New mask','dir') == 7)
			mkdir 'Image results\New mask'
        end
        
        baseDir = fullfile('Image results', 'New mask', datum);
        subDir = fullfile(baseDir, filenameshort);

        % Create base and subdirectory if they don't exist
        if ~exist(subDir, 'dir')
            mkdir(subDir);
        end
        
        
        
        set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
%       imshowcat(gcf,imadjust(mat2gray(I_mask(:,:,fr))),stats,(WS));
%       imshowcat(WS,imlocalbrighten(mat2gray(Igf)),imregionalmax(temp6))
% 
%       imshowcat2(gcf,imadjust(mat2gray(I_data(:,:,fr))),stats); %defualt
%       imshowcat(gcf,imadjust(mat2gray(I_mask(:,:,fr))),stats);
%       imshowcat(gcf,imadjust(mat2gray(I_mask_org(:,:,fr))),stats);
        imshowcat(gcf,imadjust(mat2gray(I_data(:,:,fr))),stats);
        title(['Time = ' num2str(t(fr)) ' seconds'],'FontSize',30)
        saveas(gcf, fullfile('Image results', 'New mask', datum, filenameshort, ...
    ['new_mask_', filenameshort(1:end-4), '_', num2strDU(fr,3), '.jpg']));
	
    
end

% === Relabel First Frame and Track Nuclei ===
label{1,1} = bwlabel(mask_data{1,1});
label{1,2} = mask_data{1,2};
mask_data{1,5} = mask_data{1,4};
uniq_lbl = unique(label{1,1});
max_lookback = 3; %default value set at 3 
% [label, uniq_lbl] = track_nuclei_labels(mask_data, label, max_lookback);
[label, uniq_lbl] = track_nuclei_labels_exact(mask_data, label, max_lookback);

% === Visualize and Save Label Overlays ===
for fr = 1:fr
    
    RGB = label2rgb(label{fr,1},'jet','w','noshuffle');
    
    figure('visib','off'),
    %figure ,   
    imageHandle = imshow(RGB, 'InitialMagnification', 'fit');

    % Get the axes handle containing the image.  Use this handle in the
    % remaining code instead of relying on gca.
    axesHandle = ancestor(imageHandle, 'axes');

    % Get the extrema points for each labeled object.
    s = regionprops(label{fr,1}, 'Extrema','Centroid');

    % Superimpose the text label at the left-most top extremum location
    % for each object.  Turn clipping on so that the text doesn't
    % display past the edge of the image when zooming.
    hold(axesHandle, 'on');

    % commented lines 637 - 654
    for nuc = 1:size(label{fr,2},1)
       id = label{fr,1}(round(label{fr,2}(nuc).Centroid(2)), round(label{fr,2}(nuc).Centroid(1)));
       x = (label{fr,2}(nuc).Centroid(1));
       y = (label{fr,2}(nuc).Centroid(2));
       text(x, y, sprintf('%d', id), ...
          'Parent', axesHandle, ...
          'Clipping', 'on', ...
          'Color', 'k', ...
          'FontSize', 15, ...
          'FontWeight', 'bold');
    end
    hold(axesHandle, 'off');
    
    title(['Time = ' num2str(t(fr)) ' seconds'])
    
    if ~(exist('Image results\Tracking','dir') == 7)
        mkdir 'Image results\Tracking'
    end

    baseDir = fullfile('Image results', 'Tracking', datum);
    subDir = fullfile(baseDir, filenameshort);

    % Create base and subdirectory if they don't exist
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end     

    saveas(gcf, fullfile('Image results', 'Tracking', datum, filenameshort, ...
['label_data_', filenameshort(1:end-4), '_', num2strDU(fr,3), '.jpg']));
	   
end

% === Generate Cytoplasm Mask using Watershed ===
for i = 1:size(mask_data,1)

    % Invert nuclear mask and apply watershed to find boundaries between regions
    temp = watershed(~mask_data{i,1});
    cyt_oppo = temp == 0;  % Watershed lines (0 values)
    
    % cyt_gbl{i,1} is a binary mask of cytoplasm (i.e., everything but watershed boundaries)
    cyt_gbl{i,1} = ~cyt_oppo;

    % Extract cytoplasmic region stats
    cyt_gbl{i,2} = regionprops(~cyt_oppo, 'Centroid', 'Area', 'PixelIdxList');

    % Erode nucleus labels slightly to make masks conservative
    L_nucout(:,:,i) = imerode(label{i,1}, strel('disk',3));
end

L_nuc = L_nucout;  % Store final nuclear label image


% === Assign each cytoplasmic region to the nearest nucleus ===
for i = 1:size(mask_data,1)

    temp = double(cyt_gbl{i,1});  % Start with binary cytoplasm mask
    lab = label{i,1};             % Nucleus labels in the same frame

    for j = 1:size(cyt_gbl{i,2},1)
        coord = cyt_gbl{i,2}(j).Centroid;
        r1 = round(coord(1));     % x
        r2 = round(coord(2));     % y

        % Assign the label of the nucleus closest to the cytoplasm centroid
        temp(cyt_gbl{i,2}(j).PixelIdxList) = lab(r2,r1);
    end

    cyt_gbl{i,1} = temp;  % Store label-assigned cytoplasm
end


% === Combine cytoplasm and mask to get labeled cytoplasm image ===
for i = 1:size(mask_data,1)

    % Invert nuclear mask to get cytoplasm
    bub = ~mask_data{i,1};

    % Erode to reduce overlap and separate touching cells
    bub = imerode(bub, strel('disk',1));

    % Multiply label image with cytoplasm binary mask to keep only cytoplasm
    L_cytout(:,:,i) = immultiply(cyt_gbl{i,1}, bub);

    % Store output
    L_cyt = L_cytout;
end


% === Create index of nucleus labels per frame ===
Nucindex = zeros(size(uniq_lbl,1)-1, size(mask_data,1));

for i = 1:size(mask_data,1)
    templ = unique(label{i,1});     % Unique labels (includes background)
    lim = length(templ) - 1;        % Number of labeled nuclei
    Nucindex(1:lim, i) = templ(2:end); % Skip background (label 0)
end


% === Extract properties for nuclei ===
for i = 1:size(I_mask,3)
    % Get basic shape info from binary mask
    nuc_data{i,1} = regionprops(imbinarize(L_nuc(:,:,i)), ...
        'Centroid','Area','PixelIdxList','PixelList');
end

for i = 1:size(I_mask,3)
    si = L_nuc(:,:,i);
    for j = 1:length(nuc_data{i,1})
        comp_image = zeros(size(si));
        comp_image(nuc_data{i,1}(j).PixelIdxList) = si(nuc_data{i,1}(j).PixelIdxList);

        % Boundary of the nucleus region for visualization or analysis
        nuc_data{i,1}(j).Boundary_pixels = bwboundaries(comp_image);

        % Assign label by checking the value at the centroid
        x_coord = round(nuc_data{i,1}(j).Centroid(2));
        y_coord = round(nuc_data{i,1}(j).Centroid(1));
        nuc_data{i,1}(j).label = si(x_coord,y_coord);
    end
end


% === Extract properties for cytoplasm ===
for i = 1:size(I_mask,3)
    cyt_data{i,1} = regionprops(imbinarize(L_cyt(:,:,i)), ...
        'Centroid','Area','PixelIdxList','PixelList');
end

for i = 1:size(I_mask,3)
    si = L_cyt(:,:,i);
    for j = 1:length(cyt_data{i,1})
        comp_image = zeros(size(si));
        comp_image(cyt_data{i,1}(j).PixelIdxList) = si(cyt_data{i,1}(j).PixelIdxList);

        % Get the boundary pixels for visualization or later analysis
        cyt_data{i,1}(j).Boundary_pixels = bwboundaries(comp_image);

        % Assign the label from the average pixel value in the region
        x_coord = round(cyt_data{i,1}(j).Centroid(2));
        y_coord = round(cyt_data{i,1}(j).Centroid(1));
        cyt_data{i,1}(j).label = mean(comp_image(cyt_data{i,1}(j).PixelIdxList));
    end
end
end