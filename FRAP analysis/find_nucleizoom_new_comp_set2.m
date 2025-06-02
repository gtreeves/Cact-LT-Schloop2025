function [L_nucout,L_cytout,nuc_data,cyt_data,Nucindex,side,mask_data,label] = ...
	find_nucleizoom_new_comp_set2(I_data,I_mask,t,scalings,datum,filenameshort,filenum,yesplot,h)


% This function segments and tracks the nuclei and the surrounding
% cytoplasm. This includes the bleached nucleus and the neighbouring nuclei
% as well.
%
% Inputs:
%   I_data        - 3D matrix (HxWxT) of raw protein channel images
%   I_mask        - 3D matrix (HxWxT) of nuclear channel images 
%   t             - Vector of timepoints corresponding to each frame
%   scalings      - Vector [scalex, scaley] for converting pixels to microns
%   datum         - String used to name/save result folders
%   filenameshort - Short filename string used for saving outputs
%   filenum       - Integer used to apply dataset-specific preprocessing
%   yesplot       - Logical flag (true/false) for plotting visual outputs
%   h             - Threshold value
%
% Outputs:
%   L_nucout      - 3D labeled nuclear mask (HxWxT)
%   L_cytout      - 3D labeled cytoplasmic mask (HxWxT)
%   nuc_data      - Cell array of regionprops for nuclei (one cell per frame)
%   cyt_data      - Cell array of regionprops for cytoplasm (one cell per frame)
%   Nucindex      - Matrix of nucleus label indices across frames
%   side          - Orientation of embryo 
%   mask_data     - Cell array of nuclear masks and point cloud data
%   label         - Cell array of labeled nuclear mask and corresponding regionprops


I_mask2 = I_mask;

% --- Determine scale and filtering based on 'filenum' ---
if ismember(filenum, [1:11])
    scale = 0.6;
    apply_median_filter = true;
elseif ismember(filenum, [12:14])
    scale = 2;
    apply_median_filter = true;
elseif ismember(filenum, [15:18, 25:28])
    scale = 1;
    apply_median_filter = true;
elseif ismember(filenum, [21])
    scale = .25;
    apply_median_filter = true;
elseif ismember(filenum, [22:24])
    scale = .35;
    apply_median_filter = true;
else
    scale = [];  % No processing applied
end

% --- Apply preprocessing if filenum is in the specified list ---
if ismember(filenum, [1:18, 21:28])
    for i = 1:size(I_data, 3)
        temp = immultiply(I_data(:,:,i), scale);             % Intensity scaling
        temp1 = imcomplement(uint8(temp));                   % Invert image

        if apply_median_filter
            temp1 = medfilt2(temp1);                         % Apply median filter
        end

        I_mask(:,:,i) = uint8(temp1);                        % Store in I_mask
    end
end


%
% Below are placeholders to make sure the function call isn't changed.
% 
if ~exist('yesplot','var')
	yesplot = false;
end
L0 = yesplot;

%
% This value is for the watershed taken on the averaged frame.
%
if ~exist('h','var') || h > 1 || h < 0
	h = 0.40;
	if h > 1 || h < 0
		warning(['h must be between zero and one. Default h = ',num2str(h),' used.'])
	end
end
h_bleach = 0.25;
[H,W,num_frames] = size(I_data);
k_bleach = zeros(num_frames,1);

% checking to see the image of I_data
%figure, imshow(uint8(I_data(:,:,1)))

% checking to see the image of I_data
%figure, imshow(uint8(I_mask(:,:,1)))

%
% Determine orientation of embryo
%
if numel(I_mask) > 1
	
	%
	% Max intensity projections
	% 
	I_mx = max(I_data,[],3);
	I_mx2 = max(I_mask,[],3);
    
    %figure, imshow(uint8(I_mx))
    %figure, imshow(uint8(I_mx2))
    
	
	%
	% Correlation coefficient: does dl intensity correlate with nuclei?
	%
	C = corrcoef(double(I_mx(:)),double(I_mx2(:)));
	if C(1,2) > 0.3
		side = 'ventral';
% 		I0 = (mean(double(I_mask(:,:,[1,10:end])),3) + ...
% 			mean(double(I_data(:,:,[1,10:end])),3))/2;
		I0 = (mean(double(I_mask),3) + mean(double(I_data),3))/2;
	elseif C(1,2) < -0.3
		side = 'dorsal';
% 		I0 = (mean(double(I_mask(:,:,[1,10:end])),3) + ...
% 			mean(double(imcomplement(I_data(:,:,[1,10:end]))),3))/2;
		I0 = (mean(double(I_mask),3) + mean(double(imcomplement(I_data)),3))/2;
	else
		side = 'lateral';
% 		I0 = mean(double(I_mask(:,:,[1,10:end])),3);
		I0 = mean(double(I_mask),3);
	end
else
	side = 'ventral';
% 	I0 = mean(double(I_data(:,:,[1,10:end])),3);
	I0 = mean(double(I_data),3);
end

%figure, imshow(mat2gray(I0))


% =========================================================================
% First, attempt to identify and protect the bleached nucleus
% =========================================================================

%
% Setting up erosion and dilation structuring elements. Want to erode more
% than we dilate.
%
d_erode = 2; % erode by this many microns
d_dilate = 0.5; % dilate by this many microns
see = strel('disk',round(2));
sed = strel('disk',round(1));

%
% Subtract the average of frames 2:10 from frame 1. This should give us an
% estimate of which nucleus is bleached.
%
I210 = uint16(mean(double(I_mask(:,:,2:10)),3));
% figure, imshow(uint8(I210)), title('I210 line 98')
Isub = imsubtract(I_mask(:,:,1),I210);
% figure, imshow(uint8(Isub)), title('Isub line 100')

%
% Blur the subtraction, erode, dilate, then blur again.
%
Isub1 = gaussFiltDU(Isub,round(0.5/scalings(1)));
I1 = imerode(Isub1,see);
Ibn = imdilate(I1,sed); % bn = "bleached nucleus"
Ibn = gaussFiltDU(Ibn,round(0.5/scalings(1)));
Ibn1 = mDU(double(Ibn));

%
% Find the brightest "nucleus" by applying a global threshold, then
% regionprops
%
theta = graythresh(uint8(255*Ibn1));
bw = Ibn1 > theta;
stats_bleachednuc = regionprops(bw,Ibn1,'MeanIntensity','PixelIdxList','PixelList');
m = [stats_bleachednuc.MeanIntensity];
[mmax,imax] = max(m);
bleachednucmask = false(H,W);
bleachednucmask(stats_bleachednuc(imax).PixelIdxList) = true;
bleachednucmask = imdilate(bleachednucmask,sed); % bn = "bleached nucleus"
bleachednuc = zeros(H,W);
bleachednuc(bleachednucmask) = Ibn1(bleachednucmask);




% =========================================================================
% Now we loop through each frame
% =========================================================================

d_nucerode = 1; % erode nuclei by this many microns
s_nuc = strel('disk',round(2));
d_cyterode = 0.5; % erode cytoplasm by this many microns, was 0.1
s_cyt = strel('disk',round(2));
L_nucout = zeros(size(I_data));
L_cytout = zeros(size(I_data));
Nucstats = cell(num_frames,1);
for i = 1:num_frames % can test changes with fewer frames like 1:30
	
	if numel(I_mask) > 1
		I0 = double(I_mask(:,:,i));
	else
		I0 = double(I_data(:,:,i));
    end
    
    
%     if numel(I_mask) > 1
% 		I0 = (mat2gray(I_mask(:,:,i)));
%         I0 = imadjust(I0,[.02,.4]);
% 	else
% 		I0 = double(I_data(:,:,i));
% 	end

	%
	% Erode, dilate, and blur.
	%
%     figure, imshowcat(imlocalbrighten(mat2gray(I0)))
	I1 = imerode(I0,see);
	I1 = imdilate(I1,sed);
	Igf = gaussFiltDU(I1,round(2));
%     figure, imshow(uint8(Igf))

	% Saturating some pixels, then do watershed. The new watershed will
	% serve as the basis for local thresholding.
	%
	Y = prctile(Igf(:),[10 90]); 
    %removing top 20 for 7,8,9,10
    % saturate 5% at each end
	Igf = (Igf - Y(1))/(Y(2) - Y(1));
	Igf(Igf > 1) = 1; Igf(Igf < 0) = 0;
    %figure, imshow((Igf))
    
    
    % for all general cases
    temp = histeq(Igf);
%     temp = (Igf);
    %figure, imshow(temp)

    temp2 = imerode(temp,strel('disk',5));
    %figure, imshow(temp2)
    temp3 = imcomplement(temp2);
    %figure, imshow(temp3)
    temp4 = imbinarize(temp3);
    %figure, imshow(temp4)
    temp5 = imfill(imcomplement(temp4),'holes');
    %figure, imshow(temp5)
    temp6 = imerode(temp5,strel('octagon',6));

    Iws = watershed(~(temp6));
    WS = Iws == 0;
    Lws = bwlabel(~WS);
    WS_plot{i} = WS;

    %figure, imshow(WS)
	
	%
	% Determine the label of the bleached nucleus in this frame
	%
	L_bn = Lws(bleachednucmask);
	k_bleach(i) = mode(L_bn);
	
	%
	% First we do a nucstats to get information about the nucleus
	% intensities in the current frame. When divided into watershed
	% sections this way, we can then normalize the whole image according to
	% its local intensity. Then, we will do a hard threshold on the whole
	% locally-normalized image to get the nucleus in each watershed section
	%
	Igf0 = gaussFiltDU(I0,round(10));
    % Igf0 = gaussFiltDU(I0,round(12)); % for for i = [1,2,3,4,6,7];
    
	stats1 = regionprops(~WS,Igf0,'Area','Maxintensity','MinIntensity','pixelvalues');
	n_nuc = length(stats1);
	I1 = double(immultiply(Igf0,~WS));
	for j = 1:n_nuc
		A = stats1(j).PixelValues;
		if i > 1 && j == k_bleach(i)
			Y = prctile(A,100*h_bleach);
		else
			Y = prctile(A,100*h);
		end
		I1(Lws == j) = I1(Lws == j)/double(Y);
		
	end
	mask_nuc = I1 >= 1;
	mask_nuc = imfill(mask_nuc,'holes'); 
	mask_nuc = imerode(mask_nuc,s_nuc); % erode to get conservative estimates
	
	%
	% Now, within each watershed area, eliminate any object that is not the
	% largest object. This is a fairly convoluted process.
	%
	nucstats = regionprops(mask_nuc,Lws,'Area','Meanintensity');
	mu = [nucstats.MeanIntensity]';
	a = [nucstats.Area]';
	v_del = false(length(mu),1);
	for j = 1:n_nuc
		k = find(mu == j);
		if ~isempty(k)
			aj = a(k);
			[~,imax] = max(aj);
			k_del = setdiff(1:length(k),imax);
			v_del(k(k_del)) = true; % objects that are not the one with largest area
		end
	end
	k_del  = find(v_del);
	L = bwlabel(mask_nuc);
	for j = 1:length(k_del)
		L(L == k_del(j)) = 0;
	end
	mask_nuc = L > 0;

	%
	% Remove nuclei that are less than 30% of the max MeanIntensity and
	% less than a certain number of microns^2 in area. To do this, we first
	% run regionprops.
	%
	L = Lws.*mask_nuc;
	nucstats = regionprops(L,Igf0,'Area','Meanintensity');
	MeanI = [nucstats.MeanIntensity];
	Area = [nucstats.Area];
	Areamin = 2.15/scalings(1)^2;%pi*((1.9-d_nucerode)/scalings(1))^2;
	Intensitymin = 0.3*max([nucstats.MeanIntensity]);
	k = find(Area < Areamin | MeanI < Intensitymin);
	
	for j = 1:length(k)
		if k(j) ~= k_bleach(i)
			L(L == k(j)) = 0;
		end
	end
	mask_nuc = L > 0;
	
	%
	% Masks and label matrices
	%
	mask_nuc0 = imdilate(mask_nuc,s_nuc); % dilate to set base mask
	mask_cyt = imerode(~mask_nuc0,s_cyt); % erode complement of base mask
	L_nuc = immultiply(Lws,mask_nuc); % label matrix
	L_cyt = immultiply(Lws,mask_cyt); % label matrix

    stats_ws_com = regionprops(imcomplement(WS),'Centroid','PixelList','PixelIdxList');


% Now working on each region and isolate individual nuclei, 
% by using binarization

num_comp = size(stats_ws_com,1);
gray_img = imadjust(mat2gray(gaussFiltDU(I0,round(.5/scalings(1)))));
% gray_img = imadjust(mat2gray(U2));
temp_final_img = false(size(I0));

for  id = 1:num_comp
        
    temp_pl = stats_ws_com(id).PixelList; %pixel list for one component
    
    temp_img = zeros(size(I0)); %temp image for one component
    
    
    %going through each pixel of one component and transferreing pixel intensity 
    % from grayscale image to a temp image
    
    for j = 1:size(temp_pl,1)
        
        temp_img(temp_pl(j,2),temp_pl(j,1)) = gray_img(temp_pl(j,2),temp_pl(j,1));
    
    end
    
    %figure, imshow(temp_img);
    % normal imbinarize
    temp_img = imbinarize(temp_img);
  
    %Cleaning the binarized image
    temp_img = imerode(temp_img,strel('disk',12)); %for 2gfp cases 
    temp_img = imdilate(temp_img,strel('disk',1));
    temp_img = bwareaopen(temp_img, 1500);
 
    
    %going through each pixel of temp image and transferreing pixel intensity 
    % to a temp final image basically stitching together different
    % binarized components
    for j = 1:size(temp_pl,1)
        
        temp_final_img(temp_pl(j,2),temp_pl(j,1)) = temp_img(temp_pl(j,2),temp_pl(j,1));
    
    end
    
    temp_final_img = bwmorph(temp_final_img,'bridge');
    temp_final_img = imfill(temp_final_img,'holes');
    
end

mask_data{i,1} = temp_final_img;
mask_data{i,2} = regionprops(temp_final_img,'Centroid','Area','PixelIdxList');


end

% === Update Centroids and Generate Point Clouds ===
for fr = 1:size(mask_data,1)
    
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

for fr = 1:size(mask_data,1)
    
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
        imshowcat(gcf,imadjust(mat2gray(I_mask2(:,:,fr))),stats);
%     imshowcat(gcf,imadjust(mat2gray(I_mask2(:,:,fr))),stats,WS_plot{fr})
        title(['Time = ' num2str(t(fr)) ' seconds'],'FontSize',15)
        saveas(gcf, fullfile('Image results', 'New mask', datum, filenameshort, ...
    ['new_mask_', filenameshort(1:end-4), '_', num2strDU(fr,3), '.jpg']));
	
end

% === Label First Frame and Track Nuclei ===
label{1,1} = bwlabel(mask_data{1,1});
label{1,2} = mask_data{1,2};
mask_data{1,5} = mask_data{1,4};

% k = [];
% dist = [];

k = {};
dist = {};

uniq_lbl = unique(label{1,1});

max_lookback = 3; %default value set at 3 

% [label, uniq_lbl] = track_nuclei_labels(mask_data, label, max_lookback);
[label, uniq_lbl] = track_nuclei_labels_exact(mask_data, label, max_lookback);

% === Visualize and Save Label Overlays ===
for fr = 1:size(label,1)
    
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
    
    title(['Time = ' num2str(t(fr)) ' seconds'],'FontSize',15)
    
    if ~(exist('Image results\Tracking','dir') == 7)
        mkdir 'Image results\Tracking'
    end

    baseDir = fullfile('Image results', 'Tracking', datum);
    subDir = fullfile(baseDir, filenameshort);

    % Create base and subdirectory if they don't exist
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end     
    
    %     saveas(gcf, ['Image results\Tracking',filesep,'label_data_',filenameshort(1:end-4),'_',num2strDU(fr,2),'.jpg']);
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
for i = 1:num_frames
    % Get basic shape info from binary mask
    nuc_data{i,1} = regionprops(imbinarize(L_nuc(:,:,i)), ...
        'Centroid','Area','PixelIdxList','PixelList');
end

for i = 1:num_frames
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
for i = 1:num_frames
    cyt_data{i,1} = regionprops(imbinarize(L_cyt(:,:,i)), ...
        'Centroid','Area','PixelIdxList','PixelList');
end

for i = 1:num_frames
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















