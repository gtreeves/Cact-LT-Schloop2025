function [y1,y2,umean] = reregister(u,h)
% Realigns unrolled nuclei to match heights
%
%function [y1,y2,umean] = reregister(u,h)
%
% "u": the unrolled z-slice of a cross-sectioned embryo nuclear stain
% "h": the scaled cutoff of where the edge is
%
% "y1,y2": the top and bottom indices, resp, of the nuclei for each col of
%			"u".
% "umean": the mean intensity of "u" between the indicies "y1" and "y2".
%	Only calculated if asked for.
%

%
% Initializing
%
[m,n] = size(u);
y1 = zeros(n,1);
y2 = y1;
if nargout > 2
	umean = zeros(n,1);
end

if ~exist('h','var')
	h = 0.25;
end

for j = 1:n
	t = double(u(:,j));
	
	[tmax,imax] = max(t);
	tmin = min(t);
	
	%
	% Normalizing:
	%
	t = (t - tmin)/(tmax - tmin);
	
	%
	% Finding top index (outer edge of embryo):
	%
% 	try
	i1 = find(t(2:imax) >= h & t(1:imax-1) < h);
% 	catch
% 		1
% 	end
	
	if ~isempty(i1)
		y1(j) = i1(1);
	else
		y1(j) = 1;
	end
	
	%
	% Finding bottom index (inner edge of the nuclear layer):
	%
% 	try
	i2 = find(t(imax+1:end) <= h & t(imax:end-1) > h);
% 	catch
% 		1
% 	end
	
	if ~isempty(i2)
		y2(j) = i2(end) + imax;
	else
		y2(j) = m;
	end
	
	if nargout > 2
		umean(j) = mean(double(u(y1(j):y2(j),j)));		
	end
end


