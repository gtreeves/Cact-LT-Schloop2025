function data = fit_FRAP_2cpt(data)
%Fits our secondary data (curves vs t) to a 2-cpt model for Cact-LT
%
%function data =  fit_FRAP_2cpt(data)
%
% This function fits a model of two exponentially-decaying functions to the
% bleached nucleus data, with the cytoplasmic concentration as a forcing
% function. This function uses features of the data, such as the initial
% intensity, the final intensities of the nucleus and the cytoplasm, and
% the known equilibrium constant (KnucG) for free GFP control embryos to
% constrain the parameters of these two models. This leaves us with three
% remaining parameters: phi, which is the fraction of nuclear GFP that is
% bound to Cact-LT, kout, and koutG.
%
% It should also be noted that we correct for bleaching by dividing each
% time point of the bleached nucleus by the average of all other nuclei,
% then multiplying by the average's first time point. We do the same for
% the bleached nucleus's cytoplasm (but with the average of all the other
% cytoplasmic compartments).

smoothwindow = 5;

%
% Equilibrium from free GFP controls
%
KnucG = 1.2;

%
% Extracting data
%
[t,ynuc_smooth,ycyt_smooth,ynuc,ycyt] = smooth_ynuc(data,smoothwindow);

%
% Correct for autofluorescence
%
load Mat/autofluorescence_FRAP nuc_autofluorescence_FRAP cyt_autofluorescence_FRAP
ynuc = ynuc - nuc_autofluorescence_FRAP; ynuc(ynuc < 0) = 0;
ycyt = ycyt - cyt_autofluorescence_FRAP; ycyt(ycyt < 0) = 0;
ynuc_smooth = ynuc_smooth - nuc_autofluorescence_FRAP; ynuc_smooth(ynuc_smooth < 0) = 0;
ycyt_smooth = ycyt_smooth - cyt_autofluorescence_FRAP; ycyt_smooth(ycyt_smooth < 0) = 0;

%
% Estimate the steady states
%
v = t > 360; % s
Cnucbar = mean(ynuc(v));
Ccytbar = mean(ycyt(v));

%
% Initial time point
%
c0 = ynuc(1);
	

% =========================================================================
% Fit our nuclear data to the solved differential equation
% =========================================================================

%
% Make options
%

if contains(data.matfilename,'Normal')
	phi = [0.5 0.1 0.7];
else
	phi = [0.2 0.05 0.5];
end

Lower = [0.0001 0.05 phi(2)];
Upper = [0.1 0.08 phi(3)];
StartPoint = [0.003 0.07 phi(1)];

tols = 1e-16;
options = optimset('Display','off',...
	'MaxFunEvals',5e2,...
	'MaxIter',5e2,...
	'TolX',tols,...
	'TolFun',tols,...
	'TolCon',tols ,...
	'UseParallel','never');

%
% Do the fit.  If there are outliers, remove them and re-perform
% the fit.  This only works if there are only a handful of
% outliers.  And we only do this once.
%
fhandle = @ftn_fit_FRAP_2cpt;
for jj = 1:2

	%
	% Do the fit
	%
	[P,resnorm,residual,~,~,~,J] = ...
		lsqcurvefit(fhandle,StartPoint,t,ynuc,...
		Lower,Upper,options,ycyt_smooth,c0,Cnucbar,Ccytbar,KnucG);
	
	%
	% Outlier analysis
	%
	m = mean(residual);
	s = std(residual);
	v = residual > m + 3*s | residual < m - 3*s;
	if jj == 1 && sum(v) > 0
		t = t(~v);
		ynuc_smooth = ynuc_smooth(~v);
		ycyt_smooth = ycyt_smooth(~v);
		ynuc = ynuc(~v);
		ycyt = ycyt(~v);
	else
		continue
	end
end

%
% Final parameter relationships
%
kout = P(1);
koutG = P(2);
phi = P(3);
phicyt = 1 - (1 - phi)*Cnucbar/Ccytbar/KnucG;
Knuc = phi*Cnucbar/phicyt/Ccytbar;
kin = Knuc*kout;
kinG = KnucG*koutG;

%
% Reconstruct final results
%
yfit = fhandle(P,t,ycyt_smooth,c0,Cnucbar,Ccytbar,KnucG);

%
% Calc goodness of fit
%
R2 = 1 - resnorm/norm(ynuc_smooth-mean(ynuc_smooth))^2;

%
% Make errorbars
%
nu = length(J)-length(P); % DOF
rmse = norm(residual(:))/sqrt(nu);
C = (J'*J) \ eye(length(P)); % C = inv(X'*X);
se = sqrt(diag(C))*rmse;
covar_matrix = C*rmse^2;

% 68% errorbar (one sigma)
alpha = 1 - 0.68;
delta1 = se*tinv(1-alpha/2,nu);
delta1(1:2) = delta1(1:2)*60; % convert to min^-1

% 95% errorbar
alpha = 0.05;
delta2 = se*tinv(1-alpha/2,nu);
delta2(1:2) = delta2(1:2)*60; % convert to min^-1


%
% Output
%
data.yfit = yfit;
data.Cnucbar = Cnucbar;
data.Ccytbar = Ccytbar;
data.c0 = c0;
data.phi = phi;
data.phicyt = phicyt;
data.k_in = kin*60; % min^-1
data.k_out = kout*60; % min^-1
data.K_nuc = Knuc;
data.kinG = kinG*60; % min^-1
data.koutG = koutG*60; % min^-1
data.KnucG = KnucG;
data.r2 = R2;
data.delta1 = delta1;
data.delta2 = delta2;
data.covar = covar_matrix;




