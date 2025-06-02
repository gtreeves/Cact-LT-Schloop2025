% script_mats_analysis
%
% The old function was "1" because we skipped some lsm files and this one
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
yesplot = false;
load bleaching_filenames % after running "script_gatherbleachfiles"


%
% Loop through each directory
%
Soln_1cpt = [];
for ii = 1:length(filenames)
	
	filename = filenames{ii};
	genotype = 'GFP';
    filenum = ii;
	
	data = run_analyze_bleach(filename,genotype,filenum,1,2);
	if ~isstruct(data)
		disp([data,' in:'])
		disp(filename)
		disp(' ')
	else
		Soln_1cpt = [Soln_1cpt;data];
	end
	disp(['ii = ',num2str(ii),' out of ',num2str(length(filenames))])
end


save Mat/1componentfits Soln_1cpt



