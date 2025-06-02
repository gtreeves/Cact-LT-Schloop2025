function epsDU(fh,filename,fontsize)
%Saves matlab figure to eps with better standards.
%
%function epsDU(fh,filename)
%
% This function is an alternative to Matlab's "figure" function.  I am
% writing this because, since r2104b, some of the defaults for figures have
% become a bit annoying.  
%
% Inputs: fh = figure handle, usually "gcf", and "filename" is the name of
% 	the file without extension. NOTE: if you "accidentally" don't pass the
% 	figure handle, and instead just the filename as the only argument, the
% 	function will know what to do. The figure that will be saved is the
% 	current figure.
%
% This function fixes:
%	* the default color for x,y,z axes (and also xlabel, ylabel, etc) is
%		now 75% K.
%	* I like to have 'paperpositionmode' set to 'auto'
%	* I like to have 'fontsize' set to 24

if ~ishandle(fh)
	if ~exist('fontsize','var') && exist('filename','var')
		fontsize = filename;
	end
	filename = fh;
	fh = gcf;
end
if length(filename) > 4 && (strcmp(filename(end-3:end),'.eps') || strcmp(filename(end-3:end),'.jpg'))
	filename = filename(1:end-4);
end
if ~exist('fontsize','var')
	fontsize = 24;
end

%
% Focus on the given figure handle, if it is not already current.
%
if any(findall(0,'Type','Figure') == fh) % check if figure is already open
	s = get(fh,'visib');
	
	%
	% If the figure is open and the visibility is off, then check if
	% desired figure is current. If so, then don't do anything. If not,
	% then we have to make it current, which will turn its visibility on.
	%
	if strcmp(s,'off') && ~isequal(fh,gcf)
		figure(fh) % this flashes the figure
		set(fh,'visib','off')
	end
else % open new version of the figure if it's not already open
	error('no such figure exists')
end

% set(gcf,'paperpositionmode','auto')
% set(gca,'xcolor','k','ycolor','k','zcolor','k')
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0])
set(gca,'FontSize',fontsize)
set(1,'renderer','painters')
% saveas(fh,[filename,'.fig'])

% filename = ['/home/babel/gtreeves/Office/SourceEPS/',filename,'.eps'];
print(gcf,[filename,'.eps'],'-depsc')%,'-tiff')
print(gcf,[filename,'.jpg'],'-djpeg','-r300')





