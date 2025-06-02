function [B,X,Y,Idx] = traceobject(mask,m,n)
%Finds boundary of bw object.
%
%function [B,X,Y,Idx] = traceobject(mask)
%
% This function takes a bw image "mask", which has exactly one object in
% it, and finds the pixels that perfectly trace the boundary.  This
% function is based on the image processing toolbox function
% "bwtraceboundary", but it makes the implementation easier for my
% purposes.
%
% The input "mask" can be one of two things:
% (1) A bw image with exactly one object. In this case, the interpretation
%	of this function is straightforward: the output is a set of pixels that
%	trace that object.  The inputs "m" and "n" are not used.  In this case,
%	the outputs B,X,Y are double arrays, with X being the x-coordinates of
%	the boundary, and Y the y-coords, and B = [Y X].  (Thus, "B" is more in
%	the form of array format, with rows coming first, then columns).  The
%	fourth output, "Idx" is the indices in B but in long-column format.
% (2) A structure which is the output of "regionprops" performed on a bw or
%	label image with "N" objects.  This structure must have the field
%	"PixelIdxList".  Also, you must pass the size of the original image in
%	"m" and "n" (with "m" the number of rows, "n" the number of columns).
%	Alternatively, you could put [nrows ncols] into the input "m", in which
%	case the third input is not needed.  In this case, the outputs B,X,Y
%	are cell arrays with N rows and one column, and the i-th element of
%	these corresponds to the "B","X","Y" of the i-th object.

if isstruct(mask)
	N = length(mask);
	if ~isfield(mask,'PixelIdxList')
		error('You need to have a PixelIdxList in your structure')
	end
	if ~exist('m','var') || (length(m) < 2 && ~exist('n','var'))
		error('If you pass a regionprops struct as 1st input, you also need array size.')
	end
	if length(m) == 2
		n = m(2);
		m = m(1);
	end
else
	N = 1;
end

B = cell(N,1); X = B; Y = B; Idx = B;
for i = 1:N
	if isstruct(mask)
		P = mask(i).PixelIdxList;
	else
		P = find(mask);
		[m,n] = size(mask);
	end
		
	bw = false(m,n);
	bw(P) = true;
	
	
	maxbw = max(bw);
	j0 = find(maxbw); j0 = j0(1);
	bwj0 = bw(:,j0); i0 = find(bwj0); i0 = i0(1); % now (i0,j0) is true.
	
	B{i} = bwtraceboundary(bw,[i0,j0],'N');
	X{i} = B{i}(:,2); Y{i} = B{i}(:,1);
	Idx{i} = Y{i} + m*(X{i} - 1);
end

if N == 1 && ~isstruct(mask)
	B = B{1}; X = X{1}; Y = Y{1}; Idx = Idx{1};
end