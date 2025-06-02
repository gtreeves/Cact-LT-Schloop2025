function [xp,yp] = borderFinder(I,varargin)
%Finds the border that surrounds the cross-sectioned embryo.
%
%function [xp,yp] = borderFinder(I,varargin)
%
% This function takes an image of am embryo cross-section and finds points
% on the outer edge of the embryo.  First, the image is gaussian
% filtered.  Then, taking 6-degree slices of the image (like a pizza),
% average intensity of the pixels within the slice is found as a function
% of distance from the center of the image.  The radius where the intensity
% drops to 25% maximal is where we say the edge of the embryo is.  We
% attempt some filtering of outliers at the end of the function. Single
% outliers can be detected by drastic changes in curvature.  Others can be
% corrected by fitting a circle to the embryo and removing those that are
% beyond two time the inner quartile range from the mean.
%
% "I": can be either a uint8 grayscale image, or string of filename.
%
% Optional argument varargin can consist of these things, in this order:
%	* "h": height of intensity cutoff.  Default, 0.25.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yesplot": whether or not to plot the image and the boundary. The
%		behavior of the plot is contingent on passing other variables at
%		the end of varargin "filename","idx", and/or "outfilename".
%		If none of those are passed, the figure is plotted as normal. If
%		"filename" only is passed, then the outfilename will become the
%		"filenamehort" (ie, the actual lsm filename with the immediate
%		folder it was sitting in appended in front of it, separated by an
%		underscore). If both "filename" and "idx" are passed, then
%		filenameshort will be appended with a 4-digit index specified by
%		the "idx" numeric variable. If "outfilename" is passed, this
%		overrides "filename", as the outfilename is already given. If "idx"
%		is passed by itself, then it is ignored.
%		Default (for yesplot), false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "nt": choice for number of bins in theta.  Will determine length of
%		outputs "xp" and "yp" (which will be nt+1).  Default, 60.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "xc": x-location of the center of the object in the image.  Default,
%		not specified, in which case it will be calculated as the image
%		width over 2.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yc": y-location of the center of the object in the image.  Default,
%		not specified, in which case it will be calculated as the image
%		height over 2.
%	* "filename": usually the full-path filename of the lsm file that is
%		being analyzed. filename will be converted into "filenameshort" as
%		part of the outfilename for figure export. See "yesplot" for more
%		info. Ignored if "yesplot" is false, or if "outfilename" is
%		specified. Default, non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "idx": numerical index appended to "filenameshort" to make
%		outfilename for figure export. See "yesplot" for more info. Ignored
%		if either "yesplot" is false or if "filename" is non-existent, or
%		if "outfilename" is specified. Default, non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "outfilename": the output filename for figure export. Ignored if
%		"yesplot" is false. See "yesplot" for more info. Default,
%		non-existent. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%		
%
%
% "xp,yp": the points along the perimeter of the embryo.


warning off MATLAB:divideByZero

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	h = varargin{iArg}; else
	h = 0.25;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nt = varargin{iArg}; else
	nt = 60;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	xc = varargin{iArg}; else
	xc = NaN;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yc = varargin{iArg}; else
	yc = NaN;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	filename = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	idx = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	outfilename = varargin{iArg}; else
	%
end%, iArg = iArg + 1;

%
% Reading in "I"
%
I0 = I;
if ischar(I)
	I = imread(I);
end
o = size(I,3);
if o > 1
	I = sum(double(I),3);
end

[m,n] = size(I);

%
% Gaussian filter, w/std 3.  Then, background subtract with a disk of
% radius 25 pixels.
%
sig = 3;
I = gaussFiltDU(I,sig);
I = imtophat(I,strel('disk',25));

%
%  Now we divide the image into small regions in polar coordinates.
%
if isnan(xc)
	xc = round(n/2);
end
if isnan(yc)
	yc = round(m/2);
end
x = (1:n) - xc;  y = (1:m)' - yc;
X = repmat(x,m,1); Y = repmat(y,1,n);

[Theta,R] = cart2pol(X,Y);
clear X Y

r_out = zeros(nt,1);
theta = linspace(-pi,pi,nt+1)';

%
% Here we loop through each "pizza-slice" of theta.  In each slice, we try
% to predict the location of the edge of the embryo.  In some cases, this
% is a non-trivial problem.  In particular, if there is a lot of noise/haze
% outside the embryo, there is not good contrast between the embryo and
% outside the embryo, isolated bright spots inside the embryo (interior
% nuclei eg), isolated bright spots outside the embryo (dust particle, eg),
% or non-continuous contrast at the perimeter of the embryo (sparse nuclei,
% for example).
%
for i = 1:nt
	
	%
	% I1 is a column vector of image intensity values that lie within our
	% current slice.  The col vec r1 coresponds to the radii of these
	% intensity values.  Then, we sort these vectors together such that the
	% values of r1 are ascending, and smooth.
	%
	I1 = I(Theta > theta(i) & Theta < theta(i+1));
	r0 = R(Theta > theta(i) & Theta < theta(i+1));
	
	[r0,ir0] = sort(r0);
	r0 = round(r0);
	r = (r0(1):r0(end))'; nr = length(r);
	I1 = I1(ir0);
	npt = length(I1); % number of points for raw data.
	n0 = floor(npt/nr); % properties of npt/nt.  This helps us determine the
	% shape and locations of "ones" in the matrix "M", which is the linear
	% transform from a string of individual pixel intensity values to those
	% averaged into bins.
	
	rowind = r0 - r0(1) + 1;
	colind = (1:npt)';
	M = sparse(rowind,colind,1,nr,npt); sumM = sum(M,2);
	Y0 = (M*double(I1))./sumM;
	Y0(isnan(Y0)) = 0;
	Y1 = imtophat(Y0,strel('line',100,90));
	
	%
	% Transforming the data to get between zero and one.  Then, we throw
	% out the background (that is, intensities that are below 0.05).
	%
	Y1min = min(Y1);
	Y1max = max(Y1);
	y1 = (Y1 - Y1min)/(Y1max - Y1min); % normalized intensity
	
	isig = find(y1 > 0.05); % location of signal.
	ibkgrnd = isig(end); r_end = r(ibkgrnd);
	
	%
	% Now we have an upper limit to the "r" of the embryo, r_end.  We
	% assume the edge of the embryo is in the interval (r_end/3,r_end).
	%
	[iY1,iY1] = roundx(r_end/3,r); % location of smallest r.
	[maxY,iY2] = max(Y1(iY1:ibkgrnd)); iY2 = iY2 + (iY1 - 1);
	Y2 = Y1(iY2:ibkgrnd);
	
	%
	% Normalizing the intensities again, this time only taking the chunk
	% that we care about.
	%
	Y2min = Y1min;
	Y2max = max(Y2);
	y2 = (Y2 - Y2min)/(Y2max - Y2min); % normalized intensity #2
	i_out = find(y2(1:end-1) >= h & y2(2:end) < h);
	if ~isempty(i_out)
		i_out = i_out(1);
	else
		i_out = length(y2);
	end
	
	r_out(i) = r(iY2-1+i_out); % the radius of the edge of the embryo.
end

%
% Trying to discard outliers.  Who knows if this'll work or not.  It should
% work for single outliers.  I don't think it would work for mulitple
% outliers, or outlier regions.
%
dtheta = theta(2) - theta(1);
rtt = zeros(size(r_out));
for k = 2:length(r_out)-1
	rtt(k) = (r_out(k-1) - 2*r_out(k) + r_out(k+1))/dtheta^2;
end
rtt(1) = (r_out(end) - 2*r_out(1) + r_out(2))/dtheta^2;
rtt(end) = (r_out(end-1) - 2*r_out(end) + r_out(1))/dtheta^2;

Y = prctile(rtt,[25 50 75]);
k = find(abs(rtt) > abs(3*(Y(3)-Y(1)) - Y(2)));


if ~isempty(k)
	K = conseccheck(k);
	for i = 1:length(K)
		theta2 = theta(1:end-1);
		r_out2 = r_out;
		k1 = k(K{i});
		if length(k1) > 1
			theta2(k1(2:end-1)) = [];
			r_out2(k1(2:end-1)) = [];
			r_out(k1(2:end-1)) = spline(theta2,r_out2,theta(k1(2:end-1)));
		end
	end
end

%
% The points we predict are on the outside of the embryo.
%
dr = 0; % an arbitrary offset.
xp = (r_out([1:end,1])-dr).*cos(theta) + xc; % I used to round these.
yp = (r_out([1:end,1])-dr).*sin(theta) + yc;

%
% Fitting points to a circle.
%
% xp = data.Xp{7}; yp = data.Yp{7};
[xc1,yc1,R] = circfit(xp(1:end-1),yp(1:end-1));
[th,r] = cart2pol(xp(1:end-1)-xc1,yp(1:end-1)-yc1); % new xc and yc now.

%
% The outliers are in the index variable "k"
%
Y = prctile(r-R,[25 50 75]);
IQR = (Y(3)-Y(1))/2; % inner quartile range.  Points more than twice this
% away from the mean will be modified.
k = abs((r-R)-Y(2)) > 2*IQR;

%
% Redefining th,r according to the original xc,yc
%
[th,r] = cart2pol(xp(1:end-1)-xc,yp(1:end-1)-yc);

%
% Circularizing our points
%
k2 = [k;k;k];
th2 = [th-2*pi;th;th+2*pi];
r_out2 = [r;r;r];

th2(k2) = [];
r_out2(k2) = [];
r(k) = spline(th2,r_out2,th(k));
xp = r.*cos(th) + xc;
yp = r.*sin(th) + yc;
xp = xp([1:end,1]); yp = yp([1:end,1]);


%
% Plotting, if you want.
%
if yesplot && ispc
	
	%
	% Plot, either invisibly if a filename is specified, or visibly.
	%
	if exist('filename','var') || exist('outfilename','var')
		figure('visib','off')
		
		if ~(exist('BorderImages','dir') == 7)
			mkdir BorderImages
		end
	else
		figure
	end
	
	set(gcf,'paperpositionmode','auto')
	imshowcat(gcf,I);
	hold on
	plot(xp,yp,'r o -')
	plot(r_out([1:end,1]).*cos(theta)+xc,r_out([1:end,1]).*sin(theta)+yc,'y')
	set(gca,'Position',[0 0 1 1])
	
	%
	% Export image/plot and close it
	%
	if exist('filename','var') && ~exist('outfilename','var')
		kslash = union(strfind(filename,'\'),strfind(filename,'/'));
		filename(kslash) = '_';
		filenameshort = filename(kslash(end-1)+1:end);
				
		if ~exist('idx','var')
			idx = 0;
		end
		
		print(gcf,['BorderImages',filesep,filenameshort(1:end-4),'_',num2strDU(idx,4),'.jpg'],'-djpeg','-r150')
		close(gcf)
		
	elseif exist('outfilename','var')
		print(gcf,['BorderImages',filesep,outfilename,'.jpg'],'-djpeg','-r150')
		close(gcf)
	end
	
	% 	c1628(I)
	% 	hold on
	% 	plot(xp,yp,'r o -')
	% 	plot(r_out([1:end,1]).*cos(theta)+xc,r_out([1:end,1]).*sin(theta)+yc,'y')
	
end


