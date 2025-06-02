function [filenames,filenames_short] = readdir2(pth,ext)
%Outputs list of filenames found in "pth" with extension "ext".
%
%function [filenames,filenames_short] = readdir2(pth,ext)
%
% This function is "2" because I replaced using "strtok" with "deblank" and
% "strfind" on the "ispc" to handle the case where the filename has spaces
% in it.
%
% UPDATE on 2018-04-15: changed both windows and unix to use "dir"
% function, which is cleaner.
%
% Inputs:
% 
% "pth": path to the directory you want (character variable)
% "ext": optional; character variable of extension of files you want.  If
%	you don't give it this variable, the default extension is "lsm".  
%	
%	NOTE: If you give '-' as the value of "ext", then the function will
%	look for files (or folders) without an extension (or, more precisely,
%	for files or folders that do not have any dots in their names).
%
%	NOTE2: If you give '\' or '/' as the value of "ext", then the function
%	will look just for folders.
%
% Outputs:
% 
% Both "filenames" and "filenames_short" are cell arrays that contain the
% files in your folder "pth" with extension "ext".  "filenames" includes
% the full path, and "filenames_short" only the filename itself.

if ~exist('ext','var')
	ext = 'lsm'; % default extension is "lsm"
end
if strcmp(ext,char(92)) || strcmp(ext,char(47))
	wantdir = true;
else
	wantdir = false;
end
if ~strcmp(ext(1),'.') && ~strcmp(ext,'-')
	ext = ['.',ext];
end

filenames = {};
filenames_short = {};
if ischar(pth)
	if ~strcmp(pth(end),filesep) % if "pth" does not end in a filesep
		pth(end+1) = filesep; % then add the filesep.
	end
	w = dir(pth);    

    %
    % "w" is a structure vbl with info/metadata about each entity in the
    % folder "pth".
    %
    k = 1;
    n_ext = length(ext);
    for i = 3:length(w)
        filename = w(i).name;
		
		%
		% Get a list of directories
		%
		if wantdir
			if w(i).isdir
				filenames{k,1} = [pth,filename,filesep];
				filenames_short{k,1} = [filename,filesep];
				k = k + 1;
			end
			
		%
		% Get files (but not folders) with no extension (presumably). As
		% long as it has no dot in it.
		%
		elseif strcmp(ext,'-')
			if ~w(i).isdir && isempty(strfind(filename,'.'))
				filenames{k,1} = [pth,filename];
				filenames_short{k,1} = filename;
				k = k + 1;
			end
		
		%
		% Get files with the specified extension
		%
		else
			if ~w(i).isdir && length(filename)>= n_ext && strcmp(filename(end-(n_ext-1):end),ext)
				filenames{k,1} = [pth,filename];
				filenames_short{k,1} = filename;
				k = k + 1;
			end
		end
    end
	
	
	
end

filenames = sort(filenames);
filenames_short = sort(filenames_short);








