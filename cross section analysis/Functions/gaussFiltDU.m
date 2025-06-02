function J = gaussFiltDU(I,sig)
%Runs a Gaussian filter on an image "I", w/std "sig".
%
%function J = gaussFiltDU(I,sig)
%
% This function performs a Gaussian filter (or blur) on an image "I" with
% standard deviation "sig", outputting the image "J".  The image "I" is
% assumed to have "reflective" boundaries to perform the, and are padded
% with pixels accordingly, so that the Gaussian blur can end up with an
% image the same size as the orig.
%
% "I": Image of dimension 1,2, or 3.  Can be any class.
% "sig": standard deviation.  (Optional; default, 2).  This can be scalar
%	(in which case it applies to all dimensions), or a vector of 2 or 3
%	values, with each value corresp to a dimension of the array we're
%	blurring.  If the array has 3 dimensions and sig is length 2, then the
%	blur will use sig(1) on the first two dimensions and sig(2) on the
%	third dimension.
%
% "J": the blurred image.  Same size and class as "I".

if ~exist('sig','var')
	sig = 2;
end

%
% Checking out our dimensions of "I" and the number of elts of sig.
%
[m,n,o,p] = size(I);
if p > 1
	error('Error: you screwed up.')
end
if o > 1 && (m == 1 || n == 1)
	error('You made a stupidly-dimensioned matrix, meester.')
end
p = length(sig);
if o > 1 && p == 2
	sig = [sig(1) sig(1) sig(2)];
elseif o > 1 && p ~= 2
	sig = sig*ones(1,3);
elseif m > 1 && n > 1 && p == 1
	sig = sig*ones(1,2);
end
p = length(sig);

%
% Building the kernel.  The kernel is truncated to ignore elements that
% have magnitude less than tolerance.
%
tolerance = 0.001;
mu = 0; % The mean is always zero here.
ker = cell(p,1);
for i = 1:p
	if sig(i) <= 0
		ker{i} = NaN;
	else
		X = norminv(1-tolerance,mu,sig(i));
		floorx = floor(X);
		nn = -floorx:floorx;
		ker{i} = normpdf(nn,mu,sig(i));
	end
end

s = 'symmetric';

if m == 1
	J = imfilter(I,ker{1},s,'same','conv');
elseif n == 1
	J = imfilter(I,ker{1}',s,'same','conv');
elseif o == 1
	I = imfilter(I,ker{2},s,'same','conv');
	J = imfilter(I,ker{1}',s,'same','conv');
else
	I = imfilter(I,ker{2},s,'same','conv');
	I = imfilter(I,ker{1}',s,'same','conv');
	if isnan(ker{3})
		J = I;
	else
		ker{3} = reshape(ker{3},1,1,2*floorx+1);
		J = imfilter(I,ker{3},s,'same','conv');
	end
end




