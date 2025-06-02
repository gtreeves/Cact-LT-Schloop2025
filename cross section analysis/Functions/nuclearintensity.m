function [Y,stdY] = nuclearintensity(I,nucstats)
%Computes nuclear intensity of given channels.  Works with analyze_xs.
%
%function [Y,stdY] = nuclearintensity(I,nucstats)
%
% This function calculates the mean (and standard error of the mean)
% intensity of the nuclei in I, specified by the pixel index lists in
% "nucstats".  The structure "nucstats" comes from a call to "regionprops"
% from "find_nuclei".  
%
% The image "I" is a z-slice of an embryo cross section, and only contains
% the channels that are either nuclei or contain nuclear protein.  Other
% channels in the original embryo image are not passed to this file.   The
% output mean intensities has the "n" columns, where n = (1 + number of
% nuclear protein channels).  You add the "1" because one channel is by
% definition the nuclei, and we are interested in that channel as well.
% The order is specified by the original image; so channels are removed
% from I but the remaining ones are not shuffled around.
%
% Inputs:
%
% "I": image of z-slice of embryo cross section, including only the
%	channels that have a nuclear protein.
% "nucstats": structure from "find_nuclei" (calculated by "regionprops")
%	that contains the PixelList and PixelIdxList of each nucleus.
%
% Outputs:
% 
% "Y": n_nuc-by-numch array of mean nuclear intensities.
% "stdY": misnomer, should be semY since it's the standard error of the
%	mean nuclear intensity.
%

numch = size(I,3);
n_nuc = length(nucstats);
Y = zeros(n_nuc,numch); stdY = Y;

%
% Run a loop to go through each nucleus
%
for i = 1:length(nucstats)
	v = nucstats(i).PixelIdxList; % the pixels on nucleus "i"
	if ~isempty(v)
		for j = 1:numch
			I1 = I(:,:,j);
			Iv = double(I1(v));
			
			%
			% Avoid errors of long tails: find the mean of the middle 80%
			% of pixels
			%
			P = prctile(Iv,[10 90]);
			Iv(Iv < P(1)) = [];
			Iv(Iv > P(2)) = [];
			
			Y(i,j) = mean(Iv);
% 			Y(i,j) = median(Iv);
			stdY(i,j) = std(Iv)/sqrt(length(Iv));
			
% 			if j == 2
% 				if i == 1
% 					figure
% 				end
% 				cdfplot(Iv)
% 				hold on				
% 			end
		end
	end
end
	
	