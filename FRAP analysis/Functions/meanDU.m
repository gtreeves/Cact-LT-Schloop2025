function M = meanDU(x,dim)
%Calculates the mean of an array, just like "mean", but ignores NaN's.
%
%function M = meanDU(x,dim)
%
% This function calculates the mean of your array "x", but ignores NaN's
% where they occur.  If an entire row (or col, if you've set dim = 2) of
% your array is NaN, then the result will still be NaN, but otherwise any
% NaN's won't factor into the calculation.
%
% "x": the array you're trying to take the mean of.  Expect vector or 2D
%	matrix.
% "dim": Optional argument; the dimension along the array that you're
%	averaging over. Default, 1.
%
% "M": your mean.

[m,n] = size(x);
% if ~exist('dim','var') || n == 1
% 	dim = 1;
% elseif m == 1
% 	dim = 2;
% end
if ~exist('dim','var')
 	dim = 1;
end
if m == 1 && dim == 1
	dim = 2;
end

xnan = isnan(x);
sumnotxnan = sum(~xnan,dim);
sumnotxnan(sumnotxnan == 0) = NaN; % prevent divideByZero warnings
x(xnan) = 0;
M = sum(x,dim)./sumnotxnan;