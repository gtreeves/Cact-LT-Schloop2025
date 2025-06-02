 function [nucstats,xnuc,ynuc,snuc,w,mask,U,CC_size] = find_nuclei_cyto(I,xp,yp,scalings,varargin)

%first set of changes for find_nuclei3
%Finds nuclei in embryo cross sections.
%
%function [nucstats,xnuc,ynuc,snuc,w,mask] = find_nuclei(I,xp,yp,scalings,varargin)
%
% This function takes an image of the nuclei, and unrolls it to find the
% location of the nuclei. These locations are mapped back onto the original
% location at the periphery of the non-unrolled embryo image.  Then we
% calculate the locations of the nuclei, "xnuc" and "ynuc".  The variable
% "snuc" corresponds to the pseudo-arc length location of the nuclei around
% the periphery of the embryo.
%
% Inputs:
%
% "I": image of the nuclei
% "xp","yp": embryo periphery as detected by "borderFinder". 
% "scalings": 1x3 vector, microns per pixel of [x y z] directions, resp.
%
% Optional argument varargin can consist of these things, in this order:
%	* "Yhatmax": number of pixels into the embryo we take our data from.
%		When we unroll the embryo, we are essentially throwing out all data
%		except a thin annulus.  Yhatmax is the width of that annulus.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yesplot": whether or not to plot the image and the boundary. The
%		behavior of the plot is contingent on passing other variables at
%		the end of varargin "filename","idx", and/or "outfilename".
%		If none of those are passed, the figure is plotted as normal. If
%		"filename" only is passed, then the outfilename will become the
%		"filenamehort" (ie, the actual lsm filename with the immediate
%		folder it was sitting in appended in front of it, separated by an
%		underscore). If both "filename" and "idx" are passed, then
%		filenameshort will be appended with a 4-digit index specified by
%		the "idx" numeric variable. If "outfilename" is passed, this
%		overrides "filename", as the outfilename is already given. If "idx"
%		is passed by itself, then it is ignored.
%		Default (for yesplot), false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "stage": what nuclear cycle is the embryo? it can be 10-14. Default,
%		14.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "filename": usually the full-path filename of the lsm file that is
%		being analyzed. filename will be converted into "filenameshort" as
%		part of the outfilename for figure export. See "yesplot" for more
%		info. Ignored if "yesplot" is false, or if "outfilename" is
%		specified. Default, non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "outidx": numerical index appended to "filenameshort" to make
%		outfilename for figure export. See "yesplot" for more info. Ignored
%		if either "yesplot" is false or if "filename" is non-existent, or
%		if "outfilename" is specified. Default, non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "outfilename": the output filename for figure export. Ignored if
%		"yesplot" is false. See "yesplot" for more info. Default,
%		non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "final_erosion": If you would like the nuclei at the end, before
%		being sent to the final "regionprops" call, to be eroded for
%		conservativeness's sake, then specify the number of pixels here to
%		make the disk-shaped structuring element.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs:
% 
% nucstats: structure containing the statistics of the nuclei, as detected
%	by "regionprops".
% "xnuc","ynuc": centroid locations of each nucleus
% "snuc": the pseudo-arclength locations of each nucleus
% "w": the length of the unrolled image.  Roughly equivalent to the
%	pseudo-perimeter of the embryo in pixels.
% "mask": binary-ish image of the nuclei.  I say "ish" because each
%	individual nucleus is actually assigned a different number.

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1; 
if nArg >= iArg && ~isempty(varargin{iArg})
	Yhatmax = varargin{iArg}; else
	Yhatmax = round(18.36/scalings(1));
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	stage = varargin{iArg}; else
	stage = 14;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	filename = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	outidx = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	outfilename = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	final_erosion = varargin{iArg}; else
	final_erosion = 0;
end%, iArg = iArg + 1;

I0 = I;

I0 = uint8(I0-(I0.*.85));

I = I0;

[H,W] = size(I);
nt = length(xp) - 1;
xc = W/2; yc = H/2;
x = (1:W); y = (1:H)';

%
% Background subtract and smooth out gaussian noise
%
I = imtophat(I,strel('disk',20));
I = gaussFiltDU(I);

%
% Unrolling
%
U = unroll2(I,xp,yp,Yhatmax); % unrolled strip of nuc layer
[h,w] = size(U);


% Expected radii of nuclei in microns
%
stages = 10:14;
% radii = [3.7,3.5,3,2.5,1.5]; % in microns, for stages 10-14
% radii2 = 1.5;
radii = [3.7,3.5,3,2.5,1.1]; % in microns, for stages 10-14
radii2 = 1.1;
nucrad = round(radii(stages == stage)/scalings(1));
% if isempty(nucrad), nucrad = round(1.5/scalings(1)); end
if isempty(nucrad), nucrad = round(1.1/scalings(1)); end
N = round(nucrad);
N2 = round(radii2/scalings(1));
% N3 = round(1.2/scalings(1));
N3 = round(1.0/scalings(1));
se1 = strel('line',round(2*N+1),90); % structuring element for 1D image, to
% determine rough estimate of where the nuclei are
se2 = strel('disk',N2); % structuing element for unrolled image, to take 
% away spurs and feathers
se3 = strel('disk',N3); % make minimalist nuclei


%
% Padding, morphological opening, then imregionalmax.
%
N = 10;
U2 = [U(:,end-(2*N-1):end),U,U(:,1:2*N)]; % padding image w periodic bc's


% Find the nuclear layer, first roughly by averaging the whole image, then
% more precisely by reregister
%
Y = mDU(sum(U,2));
hh = 0.25;
i1 = find(Y(1:end-1) < hh & Y(2:end) >= hh);
i2 = find(Y(1:end-1) > hh & Y(2:end) <= hh);
if isempty(i2), i2 = h; end
mask1 = false(size(U2));
mask1(i1:i2,:) = true;

mask2 = false(size(U2)); 
[y1,y2] = reregister(U2,0.5); % taking only the nuclei in PD axis
y1 = round(smooth(y1,N));
y2 = round(smooth(y2,N));
for i = 1:length(y1)
	mask2(y1(i):y2(i),i) = true;
end

mask = mask1 & mask2;

% Getting the ridge lines separating each nuclei

% [bgm fgm] = ridgelines(U); % ridgelines isn't as good
[bgm fgm] = ridgelines3(U); 

% % Getting the stats about each region as separated by ridge lines
stats_bgm_com = regionprops(imcomplement(bgm),'Centroid','PixelList','PixelIdxList','Area');


% Now working on each region and isolate individual nuclei, 
% by using binarization

num_comp = size(stats_bgm_com,1);
gray_img = mask2.*imadjust(mat2gray(U2));
% gray_img = imadjust(mat2gray(U2));
temp_final_img = false(size(U2));

%going through all components extracted from watershed ridge lines 
for  i = 1:num_comp
        
    temp_pl = stats_bgm_com(i).PixelList; %pixel list for one component
    
    temp_img = zeros(size(U2)); %temp image for one component
    trial = zeros(size(U2));
    
    %going through each pixel of one component and transferreing pixel intensity 
    % from grayscale image to a temp image
    
    for j = 1:size(temp_pl,1)
        
        temp_img(temp_pl(j,2),temp_pl(j,1)) = 1;
        trial(temp_pl(j,2),temp_pl(j,1)) = gray_img(temp_pl(j,2),temp_pl(j,1)); 
    
    end
    
    %figure, imshow(temp_img);
    temp_img = (temp_img-imbinarize(trial)); % converting to binary image
    %figure, imshow(temp_img);
    
    %Cleaning the binarized image
    temp_img = imerode(temp_img,strel('disk',3));
    
    %going through each pixel of temp image and transferreing pixel intensity 
    % to a temp final image basically stitching together different
    % binarized components
    for j = 1:size(temp_pl,1)
        
        temp_final_img(temp_pl(j,2),temp_pl(j,1)) = temp_img(temp_pl(j,2),temp_pl(j,1));
    
    end
    
end

stats = regionprops(temp_final_img,'Centroid','PixelList','PixelIdxList');

num_comp_final = size(stats,1);
for i = 1:num_comp_final
    CC_size(i) = size(stats(i).PixelIdxList,1);
end


%
% Getting "x" and "y" values for the (1) centroids and (2) pixel list of
% all the nuclei in the unrolled strip.  Later, these will be used to
% recontruct the segmented image in the original, non-unrolled coordinate
% system.
%
centroids = cat(1,stats(:).Centroid);
xcen = centroids(:,1); ycen = centroids(:,2);
xcen = xcen - 2*N; % shifting to account for the padding

v = xcen > w | xcen < 1; % these lines account for nuclei that are split
stats(v) = [];
xcen(v) = []; snuc = xcen;
ycen(v) = [];

n_nuc = length(stats);


%
% Now we transform the centroid coordinates from theta,Yhat coordinates to
% original (non-unrolled) image x,y coordinates.
%
xpround = roundx(xp,x); % "roundx" is a helper function.
ypround = roundx(yp,y);
ds = round(sqrt(diff(xpround).^2 + diff(ypround).^2));
s = [0;cumsum(ds)]; % pseudo arclength.

[theta,r] = cart2pol(xp-xc,yp-yc); theta(end) = theta(end) + 2*pi;
th_center = interp1(s,theta,xcen); % centers of nuclei in theta

r0 = interp1(theta,r,th_center);
r_center = -ycen*Yhatmax/h + r0; % centers of nuclei in "r", the distance
% to the center of the image

xnuc = r_center.*cos(th_center) + xc; % centroids in non-unrolled
ynuc = r_center.*sin(th_center) + yc; % (i.e., original) image coordinates.

%
% Now we do the same transformation, but on each individual pixel from the
% segmented nuclei.
%
L = zeros(size(I));
for i = 1:n_nuc
	px = stats(i).PixelList;
	xlist = px(:,1); ylist = px(:,2); % PixelList gives x first, then y.
	xlist = xlist - 2*N; % fixing the offset
	xlist = mod(xlist,w); % fixing for split nuclei
	thlist = interp1(s,theta,xlist);
	r0 = interp1(theta,r,thlist);
	rlist = -ylist*Yhatmax/h + r0;
	
	x2D = round(rlist.*cos(thlist) + xc); % pixel locations in non-unrolled
	y2D = round(rlist.*sin(thlist) + yc); % (original) image coordinates.
	idx = (x2D-1)*H + y2D; % pixel index in long-vector format
	idx(idx < 1) = []; idx(idx > H*W) = [];
	
	L(idx) = i; % make that pixel be equal to "i", for the i-th nucleus
end
% Lholes = L;

% L = imfill(L,8,'holes'); % filling in holes
% mask = imopen(L,se2);
mask = L;

%mask = imopen(L,strel('disk',2)% cleaning up weird edges of nuclei, again
% if final_erosion > 0
% 	mask = imerode(mask,strel('disk',ceil(final_erosion)));
% end


nucstats = regionprops(mask,I0,'PixelIdxList','PixelList','MeanIntensity');
% Note that here, "mask" is not a bw image, but is a "double" array.  This
% is because each nucleus is represented by a distinct integer.  This is
% important because, after transforming back to the original coordinates,
% and filling in the "holes" that may have appeared as the result of the
% "rerolling", some of the distinctly-segmented nuclei may then touch. They
% can still be distinguished because they have different integer labels.

%
% Removing "nuclei" that eroded away from the morph opening.  The "if"
% statment is if the ones that eroded away were last, indexically.  In
% those cases, nucstats won't even know they were there, so it will have
% fewer than n_nuc (the orignial number of nuclei) elements.  If the ones
% that eroded away were in the middle of the indexing, then nucstats will
% still have those nuclei, but they will be empty (so the "for" loop takes
% care of those).
%
if length(nucstats) < n_nuc
	k = length(nucstats)+1:n_nuc;
	xnuc(k) = []; ynuc(k) = []; snuc(k) = [];
end
v = false(length(nucstats),1);
for i = 1:length(nucstats)
	v(i) = length(nucstats(i).PixelIdxList) <= 1;
% 	v(i) = isempty(nucstats(i).PixelIdxList);
end
nucstats(v) = [];
xnuc(v) = []; ynuc(v) = []; snuc(v) = [];

%
% Nuclei whose intensity is 2 IQR's from the median will be considered
% spurious and be cast off.
%
nuc = [nucstats.MeanIntensity]';
Y = prctile(nuc,[25 50 75]);
v = nuc < Y(2)-2*(Y(3)-Y(1));
nucstats(v) = [];
xnuc(v) = []; ynuc(v) = []; snuc(v) = [];

%
% Plotting, if asked for. Two plots are made.
%
if yesplot% && ispc
% 	disp('job done')
	%
	% Trace perimeter of each nucleus
	%
	[B,~,~,Idx] = traceobject(nucstats,H,W);
	for j = 1:length(Idx)
		nucstats(j).PerimIdxList = Idx{j};
		nucstats(j).PerimXYList = fliplr(B{j});
    end
	
    % Trace perimeter of each nucleus in unrolled image
	%
    H_u = size(U2,1);
    W_u = size(U2,2);
    [B_u,~,~,Idx_u] = traceobject(stats,H_u,W_u);
	for j = 1:length(Idx_u)
		stats(j).PerimIdxList = Idx_u{j};
		stats(j).PerimXYList = fliplr(B_u{j});
    end
		
	% {
	%
	% Plot, either invisibly if a filename is specified, or visibly. This
	% second plot is the full cross section.
	%
	if exist('filename','var') || exist('outfilename','var')
		figure('visib','off')
		
		if ~(exist('NucImages','dir') == 7)
			mkdir NucImages
		end
	else
		figure
	end
	
% 	imshowcat(gcf,double(~~mask),mDU(double(I)),double(~~mask))
	imshowcat(gcf,mDU(double(I)),nucstats)
	set(gca,'Position',[0 0 1 1])
	
	
	% Export image/plot and close it
	
	if exist('filename','var') && ~exist('outfilename','var')
		% Extract only the file name with extension from the full path
        [~, name, ext] = fileparts(filename);
        filenameshort = [name, ext];
        
        % Set default outidx if it doesn't exist
        if ~exist('outidx','var')
            outidx = 0;
        end
			
		% Create the output directory if it doesn't exist
        outputDir = fullfile('NucImages', 'Cytoplasmic', filenameshort);
        
        % Create the folder if it doesn't exist
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        
        % Create the full filename for the output image
        outfile = fullfile(outputDir, [filenameshort(1:end-4), '_', num2strDU(outidx, 4), '.jpg']);
        
        % Save the figure
        print(gcf, outfile, '-djpeg', '-r150');
        close(gcf);
		
% 		figure('visib','off')
% 		%imshowcat(gcf,mDU(double(U2)),stats_bgm_com,stats)
%         %imshowcat(gcf,imadjust(mat2gray(U2)),bgm,fgm,stats)
%         imshowcat(gcf,mDU(double(U2)),bgm,fgm,stats)
% 		set(gca,'Position',[0 0 1 1])
% 		print(gcf,['NucImages',filesep,...
% 			'U2mn_',filenameshort(1:end-4),'_',num2strDU(outidx,4),'.jpg'],'-djpeg','-r150')
% 		close(gcf)
		
		
% 		imwrite(U2mn,['NucImages',filesep,...
% 			'U2mn_',filenameshort(1:end-4),'_',num2strDU(outidx,4),'.jpg'])
		
	elseif exist('outfilename','var')
		print(gcf,['NucImages',filesep,outfilename,'.jpg'],'-djpeg','-r150')
		close(gcf)
	end
	%}
	
% 	figure
% 	imshow(u,[])
% 	hold on
% 	plot(xnuc,ynuc,'*')
end



