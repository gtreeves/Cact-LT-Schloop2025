function [nucdl,cytdl,std_nucdl,std_cytdl] = ...
	dl_vs_t(I_data,L_nuc,L_cyt,Nucindex)
%
%
%
%
% This function will take a time series of data images and a time series of
% nuclear masks and extract the nuclear and cytoplasmic concentrations as a
% function of time for each nucleus and cytoplasmic island.

num_frames = size(I_data,3);
% n_nuc = max(L_nuc(:));
n_nuc = size(Nucindex,1);
nucdl = zeros(num_frames,n_nuc);
cytdl = nucdl;
std_nucdl = nucdl;
std_cytdl = nucdl;

for i = 1:num_frames
	L1 = L_nuc(:,:,i);
	L2 = L_cyt(:,:,i);
	I1 = I_data(:,:,i);
	nucindex = Nucindex(:,i);
	
	%
	% Loop through each nucleus
	%
	for j = 1:n_nuc
		label = nucindex(j);
		
		if label > 0
			Iv = double(I1(L1 == label));
			Iv = Iv(:);
			if ~isempty(Iv)
% 				Y = prctile(Iv,[10 90]); % remove top and bottom 10outliers
                Y = prctile(Iv,[1 99]); % remove top and bottom 1 outliers
				v = Iv > Y(1) & Iv < Y(2);
				nucdl(i,label) = mean(Iv(v)); % was median
				std_nucdl(i,label) = std(Iv(v))/sqrt(length(Iv(v)));
			else
				nucdl(i,j) = NaN;
				std_nucdl(i,j) = NaN;
			end
			
			Iv = double(I1(L2 == label));
			Iv = Iv(:);
			if ~isempty(Iv)
				Y = prctile(Iv,[1 99]);
				v = Iv > Y(1) & Iv < Y(2);
				cytdl(i,label) = mean(Iv(v)); % was median
				std_cytdl(i,label) = std(Iv(v))/sqrt(length(Iv(v)));
			else
				cytdl(i,j) = NaN;
				std_cytdl(i,j) = NaN;
			end
		else
			cytdl(i,j) = NaN;
			std_cytdl(i,j) = NaN;
		end
		
	end
end




