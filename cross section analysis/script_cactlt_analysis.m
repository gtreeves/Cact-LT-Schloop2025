% script_cactlt_analysis
%
% The old function was "1" because we skipped some czi files and this one
% is a lot better organized. The files processed will be the same as can be
% found in an excel table "Bleaching_stats.xlsx".
%
% This function is also possible because of the re-structuring of the
% folders in which the files are stored. Bad images are sent to the "Bad"
% folder, and zstacks are found in the "zstack" folder.
%
% This processes the directories in pths, below, which should be
% photobleaching time series. 

clear
close all
load Mat/cactlt_filenames % after running "script_gathercactltfiles"


% Loop through each directory

soln = [];
for ii = 1:length(filenames)
	
	filename = filenames{ii};
    
    filenum = ii;
	
	data = run_analyze_cactlt(filename,filenum,1,2);
	if ~isstruct(data)
		disp([data,' in:'])
		disp(filename)
		disp(' ')
	else
		soln = [soln;data];
	end
	disp(['ii = ',num2str(ii),' out of ',num2str(length(filenames))])
end




