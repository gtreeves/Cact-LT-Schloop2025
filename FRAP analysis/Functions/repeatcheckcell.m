function [V,y] = repeatcheckcell(x,nosort)
%Tells you which groups of elements in cell str vbl "x" are repeated
%
%function V = repeatcheck(x)
%
% This function accepts a cell array of strings "x". The second, optional
% output "nosort" forces the function to ignore the sorting requirement.
%
% The cell variable "V" will output which groups of indices of "x" contain
% elements that are equal.  The number of elements of "V" will correspond
% to the number of groups. 
%
% Ex: say x = {'5' '7' '11' '11' '11' '16' '17.14' '17.14'}.  Then "V" will
% have 2 elements (because there are two groups of repeated elements) and
% will be: 
%	V = {[3 4 5];[7 8]}
%
% A group must contain at least two elements (ie, no group sizes of one!).
%
% The output "y" will be the values of the repeated elements.

if ~isequal(sort(x),x) && ((exist('nosort','var') && ~nosort) || ~exist('nosort','var'))
	error('Cell array must be sorted!')
end

if ~iscellstr(x)
	error('Input must be cell array of strings!')
end
v = []; V = {}; k = 1; y = {};
for i = 1:length(x)-1
	if strcmp(x{i+1},x{i})
		if isempty(v)
			v = i;
			y = [y;x(i)];
		end
		v = [v i+1];
		V{k,1} = v;
	elseif ~isempty(v)
		k = k + 1;
		v = [];
	end
end

