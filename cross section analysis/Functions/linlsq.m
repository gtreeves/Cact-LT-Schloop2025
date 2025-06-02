function [m,b,stats,h] = linlsq(x,y,varargin)
%Performs simple linear least squares (y = mx + b)
%
% function [m,b,stats,h] = linlsq(x,y,varargin)
%
% This function performs simple linear least squares on a data set x,y.
% The model is y = m*x + b, where "m" is the slope and "b" is the
% intercept.  The inputs "x" and "y" must be col vecs of the same length.
%
% Inputs:
%	"x","y": col vecs of your data, must be same length.  "x" is the
%		explanatory, or independent vbl, and "y" is the dependent vbl.
%
% Optional argument varargin can consist of these things, in this order:
%	* "sig": the sqrt of the variance of the data in y (ie, the
%		"errorbars").  This can either be a single value or a col vec same
%		length as x and y.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "plotstats": do you want the statistics printed directly on the plot?
%		Default, false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "COLOR": Color parameter passed to the plot.  Can be in the form of
%		an abbreviation (such as 'b' for blue) or an RGB colorspace row
%		vector (such as [0 0 1] for blue).  Default, 'b'.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "Bergmann": whether or not you want to plot a "Bergmann-type" plot
%		(see Fig 3 from laChapelle and Bergmann, 2010).  Default, "false".
%	* "XLIM": if you have pre-specified x limits to the output plot.
%
% Outputs:
%	"m": best-fit slope
%	"b": best-fit intercept
%	"stats": structure with the following fields:
%		"m","b": obvious
%		"b0": right now, ybar/2.
%		"rsquare": the r-squared value, or coeff of var explained.
%		"pval_m": p-value for test that "m" isn't different from zero.
%		"pval_b": p-value for test that "b" isn't different from zero.
%		"pval_b0": p-value for test that "b" isn't in [-b0,b0].
%		"p_b0minus": probability that b < -b0
%		"p_b0plus": probability that b > b0
%		"ci95_m": 1x2 vec showing 95% confidence interval for "m"
%		"ci95_b": 1x2 vec showing 95% confidence interval for "b"
%		"df": number of degrees of freedom (that's n-2).
%		"stdparam": 1x2 vector, the standard deviation of [m,b], and is
%			also the denominator of the t-statistic for [m,b].  This is
%			useful to either construct an arbitrary-percent C.I., or to
%			recover the t-distribution that describes your params.  For
%			example, for m, if the tstatdenom is D, then the 70% confidence
%			interval for m is "m +/- t*D", where "t" is the value where the
%			tcdf is equal to 1-(1-70%)/2 or 85%, with (n-2) degrees of
%			freedom.
%		"S": Bergmann's "scaling coefficient", equal to m*xbar/ybar.
%		"pval_S0": p-value for test that "S" isn't different from zero.
%		"pval_S1": p-value for test that "S" isn't different from one.
%		"ci95_S": 1x2 vec showing 95% confidence interval for "S"
%		"std_S": the denomenator of the t-statistic for S.
%	"h": plot handles.
%


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	sig = varargin{iArg}; else
	sig = 1;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	plotstats = varargin{iArg}; else
	plotstats = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	COLOR = varargin{iArg}; else
	COLOR = 'b';
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Bergmann = varargin{iArg}; else
	Bergmann = false; % should be default false.
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	XLIM = varargin{iArg}; else
	%
end%, iArg = iArg + 1;

%
% least squares solution for m,b:
%
xbar = mean(x); ybar = mean(y);
Sxx = sum((x - xbar).^2);
Syy = sum((y - ybar).^2);
Sxy = sum((x - xbar).*(y - ybar));
m = Sxy/Sxx;
b = ybar - m*xbar;
yhat = m*x + b;
SSE = sum((y - yhat).^2);
rsquare = 1 - SSE/Syy;

S = m*xbar/ybar;

%
% t-tests for both "m" and "b" to be different from zero, as well as for S
% being different from zero and from one.  Also, now added t-test for b not
% being in the interval [-b0,b0].
%
n = length(y);
df = n - 2;

stdparam(1) = sqrt(1/Sxx*SSE/df);
tstat = m/stdparam(1);
pval_m = 2*tcdf(-abs(tstat),df);

sxx = sum(x.^2);
stdparam(2) = sqrt(sxx/n/Sxx*SSE/df);
tstat = b/stdparam(2);
pval_b = 2*tcdf(-abs(tstat),df);

std_S = stdparam(1)*xbar/ybar;
tstat = S/std_S;
pval_S0 = 2*tcdf(-abs(tstat),df);
tstat = (S-1)/std_S;
pval_S1 = 2*tcdf(-abs(tstat),df);

b0 = ybar/2;
tminus = (b+b0)/stdparam(2);
tplus = (b-b0)/stdparam(2);
p_b0minus = 1 - tcdf(tminus,n-2);
p_b0plus = tcdf(tplus,n-2);
pval_b0 = p_b0plus + p_b0minus;

%
% 95% confidence interval (also, this shows you how to construct a
% confidence interval).
%
alpha = 0.05;
t = -tinv(alpha/2,df);
dm = t*stdparam(1);
db = t*stdparam(2);
ci95_m = m + dm*[-1 1];
ci95_b = b + db*[-1 1];

dS = t*std_S;
ci95_S = S + dS*[-1 1];

%
% 68% conf int
%
alpha68 = 2*(1 - normcdf(1,0,1));
t68 = -tinv(alpha68/2,df);
dm68 = t68*stdparam(1);
db68 = t68*stdparam(2);
ci68_m = m + dm68*[-1 1];
ci68_b = b + db68*[-1 1];

dS68 = t68*std_S;
ci68_S = S + dS68*[-1 1];

%
% Closing remarks:
%
stats.m = m;
stats.b = b;
stats.b0 = b0;

stats.rsquare = rsquare;
stats.pval_m = pval_m;
stats.pval_b = pval_b;
stats.pval_b0 = pval_b0;
stats.p_b0minus = p_b0minus;
stats.p_b0plus = p_b0plus;
stats.ci95_m = ci95_m;
stats.ci95_b = ci95_b;
stats.ci68_m = ci68_m;
stats.ci68_b = ci68_b;
stats.df = df;
stats.stdparam = stdparam;

stats.S = S;
stats.pval_S0 = pval_S0;
stats.pval_S1 = pval_S1;
stats.ci95_S = ci95_S;
stats.ci68_S = ci68_S;
stats.std_S = std_S;

%
% Plotting, if asked for
%
if exist('yesplot','var') && ((islogical(yesplot) && yesplot) || ishandle(yesplot))
	if Bergmann
		x = x/xbar-1;
		y = y/ybar-1;
		yhat = yhat/ybar-1;
		if ~isequal(sig,1)
			sig = sig/ybar;
		end
		xbar = 0;
		ybar = 0;
		m = S;
		b = 0; % b = 1 - S;
		dm68 = dS68;
		dm = dS;
	end
	
	if isnumeric(yesplot)
		figure(yesplot)
	elseif ishandle(yesplot)
		figure(yesplot.Number)
	else
		figure
	end
	if ~isequal(sig,1) 
		h1 = errorbar(x,y,sig,'.');
	else
		h1 = plot(x,y,'.');
	end
	set(h1,'Color',COLOR,'Markersize',10)
	hold on
	
	%
	% Setting limits
	%
	if Bergmann
		if ~exist('XLIM','var')
			XLIM = get(gca,'XLim');
			XLIM(1) = min([-0.2,XLIM(1)]);
			XLIM(2) = max([XLIM(2),0.2]);
		end
		xlim(XLIM)
		YLIM = get(gca,'YLim');
		YLIM(1) = min([-0.2,YLIM(1)]);
		YLIM(2) = max([YLIM(2),0.2]);
		ylim(YLIM)
	else
		if ~exist('XLIM','var')
			XLIM = get(gca,'XLim');
			XLIM(1) = min([0,XLIM(1)]);
		end
		xlim(XLIM)
		YLIM = get(gca,'YLim');
		YLIM(1) = min([0,YLIM(1)]);
		ylim(YLIM)
	end
	
	%
	% Plotting best-fit line on top
	%
% 	[x,i] = sort(x);
% 	X = [XLIM(1);x;XLIM(2)];
% 	Y = [b+m*XLIM(1);yhat(i);b+m*XLIM(2)];
	Y = b+m*XLIM;
	h2 = plot(XLIM,Y);
	set(h2,'Color',COLOR,'Linewidth',1)
	
	%
	% Plotting the 68% error-bar on the slope.
	%
	m1 = m + dm68; m2 = m - dm68;
	Y1 = m1*(XLIM-xbar)+ybar;
	Y2 = m2*(XLIM-xbar)+ybar;
% 	h3 = plot(XLIM,[Y1;Y2]',':');
% 	set(h3,'Color',COLOR)
	h3 = fill([XLIM fliplr(XLIM)],[Y1 fliplr(Y2)],[0 173 238]/255);
	
	%
	% Filling the 95% confidence region on the slope
	%
	m1 = m + dm; m2 = m - dm;
	Y1 = m1*(XLIM-xbar)+ybar;
	Y2 = m2*(XLIM-xbar)+ybar;
	h4 = fill([XLIM fliplr(XLIM)],[Y1 fliplr(Y2)],[166 206 238]/255);%,'FaceColor',[0.5 1 0.5]);
	
	ax_ch = get(gca,'Children');
	ax_ch1 = [h1;h2;h3;h4];
	ax_ch(1:length(ax_ch1)) = ax_ch1;
	set(gca,'Children',ax_ch)
	
	%
	% If we're doing Bergmann plots, we now plot the diagonal line on top.
	%
	if Bergmann
		h5 = plot(XLIM,XLIM,'k :');
	else
% 		h5 = plot(XLIM,ybar/xbar*XLIM,'k :');
	end
	
	if plotstats
	s0 = sprintf('r^{2} coeff: %g\n',rsquare);
	s1 = sprintf('p-val for m: %g\n',pval_m);
	s2 = sprintf('p-val for b: %g\n',pval_b);
	s3 = sprintf('95%% conf int for m: [%6.2f %6.2f]\n',ci95_m);
	s4 = sprintf('95%% conf int for b: [%6.2f %6.2f]\n',ci95_b);
	
	htb = annotation('textbox',[0.4 0.17 0.500 0.2000]);
	set(htb,'String',[s0,s1,s2,s3,s4],'Fontsize',12,...
		'HorizontalAlignment','right','LineStyle','none')
	h = [h1;h2;htb];
	else
		h = [h1;h2];
	end
else
	h = [];
end
