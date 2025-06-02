function [t,raw,s] = domainMeas(I,xp,yp,varargin)
%Finds domains of gene expression in a cross-sectioned embryo.
%
%function [t,raw,s] = domainMeas(I,xp,yp,varargin)
%
% This function takes a truecolor cross-sectional image of an embryo and
% finds the fluorescence intensity (gene expression) as a function of
% fractions of total circumference.
%
% Inputs:
%
% "I": can be either a uint8/uint16 grayscale image, or string of filename
% "xp,yp": the points along the perimeter of the embryo (outputs of
%	"borderFinder.m".)
% 
% Optional argument varargin can consist of the following things:
%	* "Yhatmax":  the distance, in pixels) into the embryo we go when
%		creating the trapezoidal region that we measure the intensity in.
%		Default, 30 pixels.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ns": choice for number of bins in "s" (pseudoarclength) when
%		measuring the intensity around the periphery.  Must be an even
%		number.  Default, 300.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs:
%
% "t": npts-by-nChannels array of smoothened fluorescent intensities.
% "raw": same as "t", but the raw (non-smoothened) data.
% "s": arclength, going from -1 to +1, with "npts" points.
% 


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Yhatmax = varargin{iArg}; else
	Yhatmax = 30;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	ns = varargin{iArg}; else
	ns = 300;
end%, iArg = iArg + 1;


%
% Reading in "I"
%
if ischar(I)
	I = imread(I);
end

[m,n,o] = size(I);

%
% Perimeter of the embryo
%
dxp = diff(xp); dyp = diff(yp); ds = sqrt(dxp.^2 + dyp.^2); 
perim = sum(ds); 

%
% Increasing the density of points around the embryo
%
if isodd(ns), ns = ns + 1; end
S = perim/ns*ones(ns,1); S = [0;cumsum(S)]; % global pseudo arclength

x = zeros(ns+1,1); y = zeros(ns+1,1);
x(1) = xp(1); y(1) = yp(1);
Sr = 0; % running total of pseudoarclength.
i = 1;
for j = 1:ns-1
	D = S(j+1) - Sr;
	D2 = D + sqrt((x(j)-xp(i))^2 + (y(j)-yp(i))^2);
	if D2 >= ds(i) && i < length(ds)
		D2 = D - sqrt((xp(i+1)-x(j))^2 + (yp(i+1)-y(j))^2);
		i = i + 1;
	end
	d = D2/ds(i);
	x(j+1) = d*dxp(i) + xp(i);
	y(j+1) = d*dyp(i) + yp(i);
	
	Sr = Sr + sqrt((x(j+1) - x(j))^2 + (y(j+1) - y(j))^2);	
end
x(end) = x(1); y(end) = y(1);
dx = diff(x); dy = diff(y); ds = sqrt(dx.^2 + dy.^2); 
perim2 = sum(ds); ss = 2*[0;cumsum(ds)]/perim2 - 1;

%
% Defining image coordinates
%
x_im = (1:n); y_im = (1:m)';
X = repmat(x_im,m,1); Y = repmat(y_im,1,n);

%
% Building the raw data
%
raw = zeros(ns,o);

for i = 1:ns
	
	%
	% Defining window that we care about (don't want to do all these
	% calculations on the whole image...what a waste of time that would be!
	%
	[jx,jx] = roundx(x(i),x_im);
	[iy,iy] = roundx(y(i),y_im);
	
	j1 = max(jx-(Yhatmax+10),1);
	j2 = min(jx+(Yhatmax+10),n);
	i1 = max(iy-(Yhatmax+10),1);
	i2 = min(iy+(Yhatmax+10),m);
	
	%
	% Transforming coordinates on our little local rectangle.
	%
	theta = atan2(y(i+1)-y(i),x(i+1)-x(i));
	Xhat = (X(i1:i2,j1:j2) - x(i))*cos(theta) + ...
		(Y(i1:i2,j1:j2) - y(i))*sin(theta);
	Yhat = -(X(i1:i2,j1:j2) - x(i))*sin(theta) + ...
		(Y(i1:i2,j1:j2) - y(i))*cos(theta);
	
	xhat = (x(i+1) - x(i))*cos(theta) + ...
		(y(i+1) - y(i))*sin(theta); % xhat should be ~ds.
% 	yhat = -(x(i+1) - x(i))*sin(theta) + ...
% 		(y(i+1) - y(i))*cos(theta); % yhat should be zero.
	
	%
	% Finding the trapezoidal region that we will include in our
	% measurement of the raw data at this point s(i).  To do this, we need
	% to approximate the normals to our points (0,0) and (xhat,yhat).
	%
	if i == 1
		x0 = x(ns); y0 = y(ns);
	else
		x0 = x(i-1); y0 = y(i-1);
	end
	if i == ns
		x2 = x(2); y2 = y(2);
	else
		x2 = x(i+2); y2 = y(i+2);
	end
	xhat0 = (x0 - x(i))*cos(theta) + (y0 - y(i))*sin(theta);
	yhat0 = -(x0 - x(i))*sin(theta) + (y0 - y(i))*cos(theta);
	xhat2 = (x2 - x(i))*cos(theta) + (y2 - y(i))*sin(theta);
	yhat2 = -(x2 - x(i))*sin(theta) + (y2 - y(i))*cos(theta);
	
	m1 = -yhat0/(xhat - xhat0);
	m2 = yhat2/xhat2;	
	
	u = Yhat < Yhatmax & Yhat >= 0;
	v = Xhat > -m1*Yhat & Xhat < -m2*Yhat + xhat;
	bw = u & v;
	
	for j = 1:o
		I1 = I(i1:i2,j1:j2,j);
		I1 = I1(bw);
		raw(i,j) = mean(I1(:));
	end

end

%
% Smoothing. 
%
t = zeros(ns+1,o);
for j = 1:o
	I2 = [raw(1:end-1,j);raw(:,j);raw(2:end,j)];
	
	smth1 = smooth(I2,10);	
	t(:,j) = smth1(ns:2*ns);
end

raw = raw([1:end,1],:);

s = linspace(-1,1,ns+1)';
t = interp1(ss,t,s);
raw = interp1(ss,raw,s);












