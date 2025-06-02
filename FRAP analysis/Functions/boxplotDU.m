function varargout = boxplotDU(X,varargin)
%Easy boxplot
%
%function varargout = boxplotDU(X,varargin)
%
% "X" is am n-by-1 cell variable, and each element of the cell is its own
% set of data. 
%
% Optional argument varargin can consist of the following things:
%	* "G":  The variable "G" is a cell string array listing of group
%		names, and it must have the same number of elements as "X".  The
%		groups will be listed in the order they appear in "X" and "G".
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "plotdots": if logical "true", then in addition to the regular
%		box-and-whisker plot, we will also plot on top of that each
%		individual data point.  If numerical, it will be treated as "true"
%		and the numerical value will be taken as the "jitter".  Default,
%		false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "Colors": if you'd like the boxes to take on certain colors, pass
%		this argument.  This argument can be: (1) a 1-by-n character
%		array of one-letter abbrev. of colors (i.e., 'rgbcmyk'), (2)an
%		n-by-3 matrix of RGB values, or (3) a cell variable with n
%		elements, where each element is either a one-letter color abbrev.
%		or a 1-by-3 row vec of RGB values. In all cases, n = number of
%		elements in X. Default, the first seven rows of the "lines"
%		colormap.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "Markersize": if you are plotting dots and you'd like them to have a
%		certain markersize, then pass this argument.  This argument should
%		be numerical. Default, 6.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "plotviolin": if logical "true", then in addition to the regular
%		box-and-whisker plot, we will also plot on behind that the violin
%		plot-type histogram.  If numerical, it will be treated as "true"
%		and the numerical value will be taken as the "histogram type". See
%		distributionPlot for more info. Default, false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "newfigure": if logical "true", then the plot will be made in a new
%		figure. Default, true. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	G = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	plotdots = varargin{iArg}; else
	plotdots = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Colors = varargin{iArg}; else
	Colors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840]; % lines colormap
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Markersize = varargin{iArg}; else
	Markersize = 6;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	plotviolin = varargin{iArg}; else
	plotviolin = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	newfigure = varargin{iArg}; else
	newfigure = true;
end%, iArg = iArg + 1;

%
% Check input
%
if ~iscell(X)
	X = {X}; % in the case of making a boxplot of one set of data
end

[m,n] = size(X);
if m > 1 && n > 1
	error('Your cell array "X" must be non-singleton in exactly one dimension.')
end

%
% Check to see if the "grouping" cell variable exists. This variable
% contains the titles of the different categories.
%
if ~exist('G','var')
	G = cell(length(X),1);
	for i = 1:length(X)
		G{i} = num2str(i);
	end
end

%
% If our X variable has any empty elements, we'll remove them.
%
v = cellfun(@isempty,X);
X(v) = [];
G(v) = [];
if exist('Colors','var')
	if ischar(Colors) || iscell(Colors)
		Colors(v) = [];
	elseif isnumeric(Colors)
		Colors(v,:) = [];
	end
end

%
% If the grouping variable has repeated elements, we will re-group them.
%
[~,y] = repeatcheckcell(sort(G));
K = [];
if ~isempty(y)
	for i = 1:length(y)
		k = find(strfindDU(G,y{i}));
		k = k(:);
		
		if isrowvec(X)
			X(k(1)) = [X(k(1)) X(k(2:end))];
		else
			X(k(1)) = [X(k(1)); X(k(2:end))];
		end
		
		K = [K;k(2:end)];
	end
	X(K) = [];
	G(K) = [];
	
	if exist('Colors','var')
		if ischar(Colors) || iscell(Colors)
			Colors(K) = [];
		elseif isnumeric(Colors)
			Colors(K,:) = [];
		end
	end
end

if newfigure
	figure
end
if plotviolin
	if isnumeric(plotviolin)
		distributionPlot(X,'color',[0.8 0.8 0.8],'histopt',plotviolin,'showmm',2)
	else
		distributionPlot(X,'color',[0.8 0.8 0.8],'showmm',2)
	end
	hold on
end


if plotdots
	for i = 1:length(X)
		if isnumeric(plotdots)
			jitter = min([plotdots,0.5]);
			epsln = rand(length(X{i}),1) - 0.5;
		else
			jitter = 0;
			epsln = zeros(length(X{i}),1);
		end
		if exist('Colors','var')
			if ischar(Colors)
				plot(i*ones(length(X{i}),1)+epsln*jitter,X{i},'.','Color',Colors(i),'Markersize',Markersize)
			elseif isnumeric(Colors)
				plot(i*ones(length(X{i}),1)+epsln*jitter,X{i},'.','Color',Colors(i,:),'Markersize',Markersize)
			elseif iscell(Colors)
				plot(i*ones(length(X{i}),1)+epsln*jitter,X{i},'.','Color',Colors{i},'Markersize',Markersize)			
			end
		else
			plot(i*ones(length(X{i}),1)+epsln*jitter,X{i},'.','Markersize',Markersize)
		end
		hold on
	end
end

if isrowvec(X), X = X'; end
L = cell(size(X));
for i = 1:length(X)
	L{i} = i*ones(size(X{i}));
end
x = cell2mat(X);
L = cell2mat(L);

if exist('Colors','var')
	h = boxplot(x,G(L),'colors',Colors,'Symbol','+');
else
	h = boxplot(x,G(L));
end
if nargout > 0
	varargout = {h};
end