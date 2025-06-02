% script_fit_twocompFRAP_cact
% 
% This script loads the previously-calculated one-component fits to the
% data and fits the two-component model to our Cact-LT frap data. Note that
% it is not important that we are loading one-component fits, we just need
% the data (not the fits), but it was expedient to load the fits too. Note
% also that we are fitting the 2-cpt model only to the 1x1x data, so we
% will skip the free GFP and the 2xGFP embryos.

clear
close all
yesplot = true;

load Mat/1componentfits Soln_1cpt

%
% Now that we have all the mat files loaded, do the fit on each one
%
Soln_2cpt = [];
for ii = 1:length(Soln_1cpt)

	%
	% The two component model is only appropriate for 1x1x embryos, so we
	% skip the others
	%
	if contains(Soln_1cpt(ii).matfilename,'Free GFP') || contains(Soln_1cpt(ii).matfilename,'2 GFP')
		disp(['ii = ',num2str(ii),' of ',num2str(length(Soln_1cpt))])
		continue
	end
	data = fit_FRAP_2cpt(Soln_1cpt(ii),yesplot);
	Soln_2cpt = [Soln_2cpt;data];

	disp(['ii = ',num2str(ii),' of ',num2str(length(Soln_1cpt))])

end

save Mat/2componentfits Soln_2cpt

%}
