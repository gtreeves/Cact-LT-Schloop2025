function U = unroll2(I,xp,yp,varargin)
%Finds domains of gene expression in a cross-sectioned embryo.
%
%function U = unroll2(I,xp,yp,varargin)
%
% This function takes a truecolor cross-sectional image of an embryo and
% unrolls it into a straight line.  This function is "2" because, unlike
% "unroll.m", this function requires both an inner and outer boundary for
% the embryo's nuclei.
%
% "I": can be either a uint8/uint16 grayscale image, or string of filename
% "xp,yp": points around the outer perimeter of the embryo.
% 
% Optional argument varargin can consist of the following things:
%	(1) "Yhatmax": the depth into the embryo that we keep. Default, 60 pxl.
%
% "U": the unrolled image.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Yhatmax = varargin{iArg}; else
	Yhatmax = 60;
end%, iArg = iArg + 1;

%
% Reading in "I"
%
if ischar(I)
	I = imread(I);
end

o = size(I,3);
if o > 1
	I1 = sum(double(I),3);
else
	I1 = double(I);
end

% %
% % First we must extract our region of interest.  There will be a rectangle
% % that contains the embryo, outside of which all pixels will be zero.
% %
% maxbw = logical(max(I1)); j = find(maxbw); j1 = j(1); j2 = j(end);
% bwj = logical(I1(:,j1)); i = find(bwj); i1 = i(1); i2 = i(end);
% I = I(i1:i2,j1:j2,:); clear I1
[m,n,o] = size(I);

ns = length(xp)-1;
x_im = (1:n); y_im = (1:m)';

%
% Finding "x_in" and "y_in", the inner boundary of the embryo, by moving
% the periphery in by "Yhatmax" pixels.
%
xc = round(n/2); yc = round(m/2);
[theta,r] = cart2pol(xp-xc,yp-yc);
r_in = r - Yhatmax;
[x_in,y_in] = pol2cart(theta,r_in);
x_in = x_in + xc; y_in = y_in + yc;

x_out = roundx(xp,x_im); y_out = roundx(yp,y_im);
x_in = roundx(x_in,x_im); y_in = roundx(y_in,y_im);
dx = diff(x_out); dy = diff(y_out); ds = round(sqrt(dx.^2 + dy.^2));

%
% Each trapezoidaloid that we perform the keystone transform on has four
% vertices, and the x and y coordinates for these are in the rows of X_t
% and Y_t (in image coordinates).
%
X_t = [x_out(1:end-1),x_out(2:end),x_in(1:end-1),x_in(2:end)];
Y_t = [y_out(1:end-1),y_out(2:end),y_in(1:end-1),y_in(2:end)];
% Above: the four vertices in image coordinates.

%
% Each trapezoidaloid has a bounding box.  We worry about this because,
% every time we perform a transformation on a trapezoidaloid, we don't want
% to do the transformation on the whole image, since we only interested in
% the data in the trapezoidaloid.  So we create a new image, which is the
% old image cropped to the bounding box of the trapezoidaloid.  But when we
% do this, we have to redefine our coordinates.  Below are the bounding box
% and vertices in cropped coordinates.
%
xL = floor(min(X_t,[],2)); yL = floor(min(Y_t,[],2));
xU = ceil(max(X_t,[],2)); yU = ceil(max(Y_t,[],2));

X_t1 = X_t - repmat(xL,1,4) + 1;
Y_t1 = Y_t - repmat(yL,1,4) + 1; 

%
% Cropping, then performing the keystone/rotation transformation.
%
% s = [0;cumsum(ds)];
w = sum(ds);
U = zeros(Yhatmax,w,o);
ustart = 1;
for i = 1:ns
	
	%
	% Cropping.
	%
	I1 = I(yL(i):yU(i),xL(i):xU(i),:,:);
% 	[m1,n1,o1] = size(I1);
	
	%
	% The four points that define the transformation
	%		
	inpts = [X_t1(i,:)',Y_t1(i,:)'];
	outpts = [1,1
			  ds(i),1
			  1,Yhatmax
			  ds(i),Yhatmax];
	
	%
	% The actual transformation
	%
	try
	
	tform2 = maketform('projective',inpts,outpts);
	catch
		1
	end
	It = imtransform(I1,tform2,'XYScale',1,...
		'XData',outpts([1,2],1)','YData',outpts([1,3],2)');	
	
	ufinish = ustart+ds(i)-1; % this ensures the correct width of the
	U(:,ustart:ufinish,:) = It; % unroled strip, and is better for MEMORY 
	ustart = ufinish + 1;% purposes than the commented-out command below.
% 	U = [U,It];	
end










