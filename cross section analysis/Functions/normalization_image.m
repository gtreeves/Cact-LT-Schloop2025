function [normimg,t,raw] = normalization_image(I,xp,yp,xc,yc,R,Yhatmax)
%
%
%
%
% This function takes an image with a cross section of an embryo that has
% been nuclear stained (or stained with some other stain that is supposed
% to be uniform, so that any variations in intensity of the image are the
% result in variations in illumination/capture of photons) and creates a
% normalization image. This image will smooth the variation around the
% embryo's perimeter and fill the rest of the image with it (radially).
% This image will have a max of 1. That way, dividing a data image by the
% normalization image will (supposedly) create a uniformly bright data
% image. It's worked ok in the past (see Liberman et al., 2009).

%
% Parsing inputs
%
if ~exist('Yhatmax','var')
	Yhatmax = 30;
end
if ~exist('xc','var') || ~exist('yc','var') || ~exist('r','var') 
	[xc,yc,R] = circfit(xp,yp);
end
ns = length(xp) - 1;




% =========================================================================
% First step is to do something very similar to domainMeas.
% =========================================================================

%
% Defining image coordinates
%
[m,n] = size(I);
x_im = (1:n); y_im = (1:m)';
X = repmat(x_im,m,1); Y = repmat(y_im,1,n);


%
% Put together the points for the calculations of the inner vertices (two
% inner vertices for each trapezoidaloid). First we shift our perimeter
% points so that we can easily calc distances between next-nearest
% neighbors.
%
xp0 = [xp(ns); xp(1:ns-1)];
yp0 = [yp(ns); yp(1:ns-1)];
xp2 = [xp(2:ns); xp(1)];
yp2 = [yp(2:ns); yp(1)];
dx2 = xp2 - xp0;
dy2 = yp2 - yp0;

%
% Unit vectors pointing in the direction of the local norm to our perimeter
% points. The reason why these vectors are [-dy dx] is because that is the
% minus one over the slope vector. Then we add Yhatmax times these unit
% vectors to get a set of inner points, A.
%
r1 = [-dy2 dx2]; r1 = r1./repmat(sqrt(r1(:,1).^2 + r1(:,2).^2),1,2);
r2 = -r1;

A1 = [xp(1:ns) yp(1:ns)] + Yhatmax*r1; % don't know which is the right direction.
A2 = [xp(1:ns) yp(1:ns)] + Yhatmax*r2;

%
% Since we didn't know which was the right direction, we'll test that to
% see which "inner" points are closer to the center of the image.
%
D1 = A1 - repmat([xc yc],ns,1); d1 = sqrt(D1(:,1).^2 + D1(:,2).^2);
D2 = A2 - repmat([xc yc],ns,1); d2 = sqrt(D2(:,1).^2 + D2(:,2).^2);

% take whichever one had smaller distance
A = A1.*repmat((d1 < d2),1,2) + A2.*repmat((d2 < d1),1,2); 
x_in = A(:,1); y_in = A(:,2);
x_in = [x_in(1:end); x_in(1)];
y_in = [y_in(1:end); y_in(1)];


%
% For loop to go around the perimeter of the embryo and build the raw data.
% Each point in the raw data corresponds to the average intensity within
% the corresponding trapezoidaloid.
%
raw = zeros(ns,1);
dx = diff(xp); dy = diff(yp);
dx_in = diff(x_in); dy_in = diff(y_in);
M3_out = false(m,n,ns); M4_out = false(m,n,ns);
for i = 1:ns
	
	%
	% Need to define the conditions of each side of our trapezoidaloid.
	% First condition is the perimeter of the embryo. This should have the
	% same sign as when the line is compared to the center.
	%
	sgn = sign((yc - yp(i))*dx(i) - (xc - xp(i))*dy(i));
	M1 = sgn*((Y - yp(i))*dx(i) - (X - xp(i))*dy(i)) > 0;
	
	%
	% Second condition is the inner border. This should have the opposite
	% sign from when the line is compared to the center.
	%
	sgn = sign((yc - y_in(i))*dx_in(i) - (xc - x_in(i))*dy_in(i));
	M2 = sgn*((Y - y_in(i))*dx_in(i) - (X - x_in(i))*dy_in(i)) < 0;
	
	%
	% Third condition: i-side of the trapezoidaloid between the perimeter
	% and inner. This should have the same sign as when the line is
	% compared to the perimeter point xp(i+1).
	%
	sgn = sign((yp(i+1) - yp(i))*(xp(i)-x_in(i)) - (xp(i+1) - xp(i))*(yp(i)-y_in(i)));
	M3 = sgn*((Y - yp(i))*(xp(i)-x_in(i)) - (X - xp(i))*(yp(i)-y_in(i))) > 0;
	M3_out(:,:,i) = M3; % store for later
	
	%
	% Fourth condition: (i+1)-side of the trapezoidaloid between the
	% perimeter and inner. This should have the same sign as when the line
	% is compared to the perimeter point xp(i).
	%
	sgn = sign((yp(i) - yp(i+1))*(xp(i+1)-x_in(i+1)) - (xp(i) - xp(i+1))*(yp(i+1)-y_in(i+1)));
	M4 = sgn*((Y - yp(i+1))*(xp(i+1)-x_in(i+1)) - (X - xp(i+1))*(yp(i+1)-y_in(i+1))) > 0;
	M4_out(:,:,i) = M4; % store for later
	
	%
	% All four constraints together to make the raw data.
	%
	bw = M1 & M2 & M3 & M4;	
	I1 = I(bw);
	raw(i) = mean(double(I1(:)));
end



%
% Smoothing the raw data and normalizing to the mean.
%
I2 = [raw(1:end-1);raw;raw(2:end)];
smth1 = smooth(I2,10);	
t = smth1(ns:2*ns);
t = t/mean(t);



% =========================================================================
% Next, we transform this into variation along the entire image.
% =========================================================================

normimg = ones(m,n);

%
% Each trapezoidaloid will be extended out from the embryo and also into
% the embryo towards the center. To do this, we "remember" the direction
% vectors that we used to make the inner boundary. (We also have to
% remember which of the two antiparallel directions was correct.)
%
r = r1.*repmat((d1 < d2),1,2) + r2.*repmat((d2 < d1),1,2); 
A = [xp(1:ns) yp(1:ns)] + 0.8*R*r;
x_in2 = A(:,1); y_in2 = A(:,2);
x_in2 = [x_in2(1:end); x_in2(1)];
y_in2 = [y_in2(1:end); y_in2(1)];
dx_in2 = diff(x_in2); dy_in2 = diff(y_in2);

%
% Now, another loop around the perimeter of the embryo.
%
% figure
% imshow(I,[]);
% hold on
for i = 1:ns
	
	%
	% Here, we reuse the third and fourth constraints from the first for
	% loop. Those were stored in the arrays M3_out and M4_out. Then, we
	% need a different version of the second constraint, but this time with
	% an x_in and y_in that are much closer to the center of the embryo
	% (and may overlap with other x_in,y_in's, but that should be ok).
	%
	% This second constraint should have the same sign as when the line is
	% compared to the perimeter.
	%
	sgn = sign((yp(i) - y_in2(i))*dx_in2(i) - (xp(i) - x_in2(i))*dy_in2(i));
	M2 = sgn*((Y - y_in2(i))*dx_in2(i) - (X - x_in2(i))*dy_in2(i)) > 0;
	
	%
	% Now, combine the three constraints:
	%
	bw = M2 & M3_out(:,:,i) & M4_out(:,:,i);
	normimg(bw) = t(i);
	
% 	plot([xp(i:i+1);x_in(i+1:-1:i);xp(i)],[yp(i:i+1);y_in(i+1:-1:i);yp(i)],'y o -')
% 	plot([x_in(i:i+1);x_in2(i+1:-1:i);x_in(i)],[y_in(i:i+1);y_in2(i+1:-1:i);y_in(i)],'c o -')
	
end















