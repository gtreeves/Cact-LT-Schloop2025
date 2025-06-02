function tf = isrowvec(x)
%tells whether "x" is a row vector
%
%function tf = isrowvec(x)

N = size(x);
if length(N) == 2 && N(1) == 1 && N(2) > 1
	tf = true;
else
	tf = false;
end