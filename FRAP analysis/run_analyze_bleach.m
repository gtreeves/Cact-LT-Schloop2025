function soln = run_analyze_bleach(filename,genotype,filenum,data_ch,varargin)
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
nLSM = length(filenames);



%
% Looping through and extract data from each file.
%
j = 1; LE = {}; FE = {}; e = 1;
for i = 1:nLSM
	filename = filenames{i};
	
	try
		
		data = analyze_bleach_mat_cact(filename,genotype,filenum,data_ch,mask_ch,yesplot);
		
        
        %
		% Immediate info about fault in analyze_bleach. Error will come
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
		
		        
        % Get the folder containing the file
        fileFolder = fileparts(filename);
        
        % Get all folder parts in the path
        folders = strsplit(fileFolder, filesep);
        
        % Ensure there are enough folders to extract from
        if numel(folders) >= 1
            % The folder right before the .mat file (i.e., second last folder)
            datum_candidate = folders{end-1};
            
            % Validate it contains a date pattern
            if ~isempty(regexp(datum_candidate, '\d{4}-\d{2}-\d{2}', 'once'))
                datum = datum_candidate;
            else
                datum = '';
                warning('Expected date-containing folder not found.');
            end
        else
            datum = '';
            warning('Folder structure too shallow to extract datum.');
        end
        data1 = fit_bleach_mat_cact(filename,data,datum);
		data = data1; % so that "data" is not erased in debugging		
		
        soln(j,1) = data;
        
		
		% Save data file near the original lsm file
        [~, filenameshort, ~] = fileparts(filename);
        
        % Determine category based on keywords in the path or filename
        
        % Strip leading/trailing quotes and convert to lowercase
        clean_filename = lower(strip(filename, ''''));

        % Check for substring and assign group folder accordingly
        if contains(clean_filename, '2x egfp')
            group_folder = '2 GFP';
        elseif contains(clean_filename, 'free')
            group_folder = 'Free GFP';
        else
            group_folder = 'Normal gfp';
        end
        
        % Extract the folder containing the file (e.g., '2023-07-06 free GFP FRAP')
        fileFolder = fileparts(filename);
        folderParts = strsplit(fileFolder, filesep);
        
        if numel(folderParts) >= 1
            parentFolder = folderParts{end-1};  % This gives '2023-07-06 free GFP FRAP'
        else
            error('Folder structure is not deep enough.');
        end
        
        % Construct the save directory path
        save_dir = fullfile('czi mat results', group_folder, parentFolder, filenameshort);
        
        % Create the directory if it doesn't exist
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        % Save the file in the constructed path
        save(fullfile(save_dir, [filenameshort '_results.mat']), 'soln');
		
		% Save the current step in the analysis
		save(['run_analyze_bleach_Current',outfilename],'soln')
		disp(['j = ',num2str(j)])
		j = j + 1;
		
	catch errorinfo
		
		%
		% Information passed when one of the files throws an error
		%
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




