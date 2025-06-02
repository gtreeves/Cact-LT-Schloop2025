function [A,B,mu,sig] = fit_gaussian2(x,y,sig0)
%Fits x,y data to a gaussian.
%
%function [A,B,mu,sig] = fit_gaussian2(x,y,sig0)
%
% This function is "2" because the original was written to work with
% analyze_xs and specifically works with a circular function on the domain
% (-1,1]. So, for example, this works for the Dl gradient or pMad gradient.
% (We didn't need the von mises distribution because those Gaussians don't
% extend so far as to wrap around.) However, this function is just for
% normal Gaussians on an infinite domain not restricted to (-1,1].
%
% This function accepts as input x,y, and output the parameters of the
% gaussian:
% Y = A*exp(-(x-mu).^2/2/sig^2) + B;
%
% The third input is optional and represents an estimate of sigma.



%
% Getting an estimate of what and where the max is
%
ysmooth = smooth(x,y); % smoothing
[ymax,imax] = max(ysmooth); ymin = min(ysmooth);
xmid = x(imax);

%
% The initial guesses and upper and lower bounds of each parameter
%
A = ymax - ymin; AL = 0.1*A; AU = 10*A;
B = ymin; BL = 0; BU = A + ymin;
mu = xmid; muL = x(1); muU = x(end);
if ~exist('sig0','var')
	sig0 = (x(end)-x(1))/3;
	sigL = 0;
	sigU = Inf;
else
	sigL = 0.5*sig0;
	sigU = 2*sig0;
end

%
% Defining options and the equation to fit to
%
% f = fittype('A*exp(-(x-mu)^2/2/sigma^2)+B');
% opts = fitoptions('Method','NonlinearLeastSquares',...
% 	'Startpoint',[A B M mu sig0],...
% 	'Lower',[AL BL ML muL sigL],...
% 	'Upper',[AU BU MU muU sigU]);


fhandle = @(p,x)p(1)*exp(-(x-p(3)).^2/2/p(4)^2)+p(2);
StartPoint = [A B mu sig0];
Lower = [AL BL muL sigL];
Upper = [AU BU muU sigU];

tols = 1e-16;
options = optimset('Display','off',...
	'MaxFunEvals',5e2,...
	'MaxIter',5e2,...
	'TolX',tols,...
	'TolFun',tols,...
	'TolCon',tols ,...
	'UseParallel',false);

%
% The actual fit
%
p = lsqcurvefit(fhandle,StartPoint,x,y,Lower,Upper,options);
A = p(1);
B = p(2); 
mu = p(3);
sig = p(4);







