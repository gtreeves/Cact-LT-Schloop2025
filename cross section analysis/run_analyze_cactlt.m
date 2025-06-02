function soln = run_analyze_cactlt(filename,filenum,data_ch,varargin)
%Runs "analyze_bleach" on each lsm file in dir given in "pth".
%
%function soln = run_analyze_bleach(filename,genotype,data_ch,varargin)
%
% Optional argument varargin can consist of these things, in this order:
%	* "mask_ch": numeric input that tells the function which channel has
%		the image data from which to derive the mask.  This can be used to
%		calculate the autocorrelation function for any arbitrarily-shaped
%		region, and perhaps even disjoint regions.  If passed as a scalar
%		numerical zero, or logical false, then the mask will not be used
%		and the autocorrelation function for the entire frame will be
%		calculated. Default, 2.
%		If this is not specified, but you still want to
%		specify other arguments, put empty brackets -- [] -- in place of
%		this argument.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "outfilename": If you want the output save files to have a name other
%		than just "run_analyze_xsCurrent.mat", it appends this 
%		character variable onto the end of the name.  Must be a character
%		variable, otherwise defaults.  Default, empty string, ''.  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else 
	mask_ch = 2;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	outfilename = varargin{iArg}; else
 	outfilename = '';
end%, iArg = iArg + 1;

if ischar(filename) && ~isdir(filename)
	filenames = {filename};
elseif iscellstr(filename)
	filenames = filename;
elseif isdir(filename)
	pth = filename;
	filenames = readdir2(pth,'lsm');
else
	%
end
nCZI = length(filenames);



%
% Looping through and extract data from each file.
%
j = 1; LE = {}; FE = {}; e = 1;
for i = 1:nCZI
	filename = filenames{i};
	
	try
		
        % first run the analysis for the cytoplasmic regions
        cyt_data = analyze_cactlt_cyto(filename,data_ch,mask_ch,yesplot);
 
        % next run the analysis for the nuclear regions
		data = analyze_cactlt_nuc(filename,filenum,cyt_data,data_ch,mask_ch,yesplot);
		
		%
		% Immediate info about fault in analyze_cactlt. Error will come
		% when trying to fit to a bleach recovery model. That error will be
		% caught by the try/catch block.
		%
		if ~isstruct(data)
			kslash = strfind(filename,filesep);
			pth = filename(1:kslash(end));
			filenameshort = filename(kslash(end)+1:end);
			
			disp([data,' in:'])
			disp(pth)
			disp(filenameshort)
			disp(' ')
		end
		
	
		% Save data file near the original czi file
		
		save([filename(1:end-4),'_data.mat'],'data');
		
		% Save the current step in the analysis
        
		soln(j,1) = data;
		save(['run_analyze_cactlt_Current',outfilename],'soln')
		disp(['j = ',num2str(j)])
		j = j + 1;
		
	catch errorinfo
		
		% Information passed when one of the files throws an error
		
		lastE = errorinfo;
		LE{end+1} = lastE;
		FE{end+1} = filename;
		fprintf('%s in %s\n',LE{end}.message,filename)
		save(['run_analyze_bleach_Error',outfilename],'LE','FE')
		disp(['e = ',num2str(e)])
		e = e + 1;		
	end
	
end


if ~exist('soln','var')
	soln = 'no successful runs completed';
end




