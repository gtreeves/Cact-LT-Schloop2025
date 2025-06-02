function V = conseccheck(x)
%Tells you which groups of elements are consecutive
%
%function V = conseccheck(x)
%
% This function accepts a vector "x" that contains a sequence of increasing
% integers.  If not integers, the numbers will be rounded.  If not
% increasing, an error will be returned.
%
% The cell variable "V" will output which groups of indices of "x" contain
% numbers that increase by one each time.  The number of elements of "V"
% will correspond to the number of groups.
%
% Ex: say x = [5 6 7 10 11 15 16 17 18].  Then "V" will have 3 elements
% (because there are three groups of consecutive integers) and will be:
%	V = {[1 2 3];[4 5];[6 7 8 9]}
%
% A group must contain at least two integers (ie, no group sizes of one!).
%
%

x = round(x);
if sort(x) ~= x
	error('values must be monotonically increasing!')
end

v = []; V = {}; k = 1;
for i = 1:length(x)-1
	if x(i+1) - x(i) == 1
		if isempty(v)
			v = i;
		end
		v = [v i+1];
		V{k,1} = v;
	elseif ~isempty(v)
		k = k + 1;
		v = [];
	end
end

