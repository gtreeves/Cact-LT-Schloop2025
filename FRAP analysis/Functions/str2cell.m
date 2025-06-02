function c = str2cell(s,sep)
%Parses a string variable into a list contained within a cell array.
%
%function c = str2cell(s,sep)
%
% This function takes a row charater array, "s", and converts it into a
% cell array of strings, "c".  The substrings of "s" that get parsed into
% the different elements of "c" are separated by a special character,
% "sep".  For example, if "sep" is ',' (a comma), and "s" is a comma
% separated list, then the parts of the list (with the comma as well as the
% spaces before and after the comma removed) get parsed into the cell array
% of strings.
%
% Ex: s = 'dog, cat, mouse'
%	c = str2cell(s,',') returns c = {'dog','cat','mouse'}.
%
% In addition, any resulting empty cells get removed. 
% Ex: s = 'dog,,, cat, mouse'
%	c = str2cell(s,',') returns c = {'dog','cat','mouse'}
%		instead of c = {'dog','','','cat','mouse'}
%
% If "sep" is not passed, a comma is used as default.



%
% Input checks
%
if ~ischar(s) && ~iscellstr(s)
	error('Input "s" must be a character array or cell array of strings')
end
[m,n] = size(s);
if ischar(s) && (isempty(s) || m > 1)
	error('Character array input "s" must be non-empty with one row')
elseif iscellstr(s) && n > 1
	error('Cell string array input "s" must be non-empty with one column')
end
if ~exist('sep','var')
	sep = ',';
end
nsep = length(sep);

%
% Converting
%
C = cell(m,1);
for j = 1:m
	if iscellstr(s)
		s0 = s{j};
	else
		s0 = s;
	end
	n = length(s0);
	v = strfind(s0,sep);
	N = length(v) + 1;
	v = [0 v n+1];
	c = cell(1,N);
	for i = 1:N
		if i == 1
			k = 1;
		else
			k = nsep;
		end
		s1 = strtrim(s0(v(i)+k:v(i+1)-1));
% 		s1 = s0(v(i)+1:v(i+1)-1);
% 		s1 = deblank(s1); % nifty function to remove trailing whitespace
% 		s1 = fliplr(s1);
% 		s1 = deblank(s1);
% 		s1 = fliplr(s1);
		c{i} = s1;
	end
	
	v = cellfun(@isempty,c);
	c(v) = [];
	
	C{j} = c;
end

if iscellstr(s)
	c = C;
end

