function [IM,lsminf1,lsminf2,bg,bgsig] = ftn_imread(filename,subtr_bg,stabilize)
%Read in lsm file, perform background subtr if asked for
%
%function [IM,lsminf1,lsminf2] = ftn_imread(filename,subtr_bg)
%
% This function takes an lsm file (the full path to the file, and the
% filename should be passed in the character variable "filename"), reads it
% in and reads in the metadata. Then the image data is background
% subtracted with the mode of the first frame/z-slice. Alternatively, one
% can feed to background subtraction value.
%
% Inputs:
% "filename": character variable with full path to the file, including the
%	filename and extension 
% "subtr_bg": can be a logical or numeric variable. If it is false, then no
%	background subtraction will take place. If it is logical true, then
%	background subtraction will take place with the mode of each channel's
%	first frame/z-slice. If passed as a numeric, non-zero, non-unity
%	variable, then this value will be used to subtract the background. If
%	the numeric value is non-scalar, then the number of elements should
%	match the number of channels in the image data, and each value will be
%	used to subtract the background of the corresponding channel.



if strcmp(filename(end-2:end),'lsm')
	lsm =  lsmRead2(filename);
	lsminf1 = lsm.lsm;
	IM = lsm.data;
	clear lsm % to save memory.
	lsminf2 = lsminfo(filename);
elseif strcmp(filename(end-2:end),'czi')
	[IM,lsminf1,lsminf2] = openczi(filename);
else
	error('Only lsm or czi files supported')
end

%
% Check for notes that suggest only some of the frames should be used. In
% particular, the user can save a text file under the filename (including
% the full path) with a ".txt" extension instead of ".lsm" or ".czi". In
% this text file, the user can specify "i_start=<M>;" and/or "i_end=<N>;"
% where "M" is the first and "N" is the last frame the user wants included
% in the analysis.
%
W = lsminf2.DimensionX;
H = lsminf2.DimensionY;
D = lsminf2.DimensionZ;
T = lsminf2.DimensionT;
numch = lsminf2.NUMBER_OF_CHANNELS;
IM = reshape(IM,[H,W,numch,D,T]);

i_start = 1; i_end = T;
if exist([filename(1:end-4),'.txt'],'file')
	s = fileread([filename(1:end-4),'.txt']);
	eval(s)
	i_start = max(i_start,1);
	i_end = min(i_end,T);
	IM = IM(:,:,:,:,i_start:i_end);
end




%
% Extracting our region of interest.  There will be a rectangle
% that contains the embryo, outside of which all pixels will be
% zero.
%
if lsminf2.ScanInfo.USE_ROIS&&~lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS
	I = sum(sum(sum(IM,3),4),5);
	maxbw = logical(max(I)); j = find(maxbw); j1 = j(1); j2 = j(end);
	bwj = logical(I(:,j1)); i = find(bwj); i1 = i(1); i2 = i(end);
	
	IM = IM(i1:i2,j1:j2,:,:,:);
end

%
% Subtract background, if asked for.
%
num_ch = size(IM,3);
if exist('subtr_bg','var') && ~isempty(subtr_bg)
	if isnumeric(subtr_bg)
		
		if length(subtr_bg) < num_ch
			bg = subtr_bg(1)*ones(1,num_ch);
			IM = imsubtract(IM,bg(1));
			bgsig = NaN;
		else
			bg = subtr_bg(1:num_ch);
			for i = 1:num_ch
				IM(:,:,i,:,:) = imsubtract(IM(:,:,i,:,:),bg(i));
			end
			bgsig = NaN;
		end
		
	elseif islogical(subtr_bg) && subtr_bg

		%
		% We assume the background intensity level (i.e., the true black
		% level) for each channel is equal to the mode of intensities seen
		% in a slice. 
		%
		[bg,bgsig] = ftn_calcbg(IM(:,:,:,1,1),false);
		for i = 1:num_ch
			IM(:,:,i,:) = imsubtract(IM(:,:,i,:,:),bg(i));
		end
	else
		bg = NaN;
		bgsig = NaN;
	end
else
	bg = NaN;
	bgsig = NaN;
end



%
% Stabilize movement in the image time series
%
if exist('stabilize','var') && ((isnumeric(stabilize) && stabilize > 0) || ...
		(islogical(stabilize) && stabilize))
	[dx,dy] = image_movement(squeeze(sum(IM,3))); % this is probably fine because we won't saturate
	if isnumeric(stabilize)
		minframesize = stabilize;
	else
		minframesize = [];
	end
	IM = image_stabilize(IM,dx,dy,minframesize);
end



