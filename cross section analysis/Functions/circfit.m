function   [xc,yc,R,a] = circfit(x,y)
%Given points in vectors x,y, finds best-fit circle.
%
%function [xc,yc,R,a] = circfit(x,y)
%
% This function fits a circle  in x,y plane in a more accurate (less prone
% to ill condition ) procedure than circfit2 but using more memory.
%
% "x,y": column vecs where (x(i),y(i)) is a measured point
%
% "xc,yc": estimate of the center of the circle.
% "R": estimate of the radius of the circle.
% "a": optional output is the vector of coeficient "a" describing the
%	circle's equation:
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 

x = x(:); y = y(:);
a = -[x y ones(size(x))]\(x.^2 + y.^2);
xc = -a(1)/2;
yc = -a(2)/2;
R  =  sqrt((a(1)^2 + a(2)^2)/4 - a(3));