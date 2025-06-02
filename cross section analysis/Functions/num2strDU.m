function s = num2strDU(n,digits)
%Makes a string out of an integer, possibly with leading zeros.
%
%function s = num2strDU(n,digits)
%
% This function is an extension of "num2str" which changes a number into
% the corresponding string.  In this function, we are especially concerned
% with taking integer numbers and converting them into a string with a
% fixed number of elements.  For example, if the target number of elements
% in the string is 3, and the integer is "n = 7", then the string will be
% "s = '007'".
%
% "n": number to be converted to string.  Must be a nonnegative integer.
% "digits": the number of digits the output string should have.  This will
%	determine the number of zeros to lead with (if any).  Must be a number
%	greater than or equal to the number of digits in "n", except in the
%	special case that "digits" is zero, in which case the  program returns
%	an empty string for "s". 
%
% "s": the output string, perhaps with leading zeros.

DU = n == round(n);
if ~DU || length(n) > 1 || n < 0
	error('You need to pay better attention to how to input "n".')
end

s = num2str(n);
d2 = length(s);
difc = digits - d2;

if digits == 0
	s = '';
elseif difc < 0
	error('You need to pay better attention to how to input "digits".')
elseif difc > 0
	Z = repmat('0',1,difc);
	s = [Z,s];
end