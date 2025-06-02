function data = fit_bleach_mat_cact(filename,data,datum,experimenttype,yesplot)
%Fits our secondary data (curves vs t) to a bleach recovery model
%
%function data = fit_bleach(data,experimenttype)
%
% This function fits the bleach data to an exponentially-decaying model,
% but we have the cytoplasmic concentration as a forcing function. To make
% sure the fitted nuc data does not have artificial jiggling from the noise
% in cytoplasmic data acquisition, we smooth cytoplasmic dl data for each
% bleached nucleus.

if ~exist('experimenttype','var')
	experimenttype = 'bleach';
end
if ~exist('yesplot','var')
	yesplot = false;
end

vsep = strfind(filename,filesep);
vsep1 = vsep(end);
% vsep2 = vsep(end-1);
filenameshort = filename(vsep1+1:end);

%
% Extracting
%

t = data.t(2:end);
ynuc = data.y_nuc(2:end);
ycyt = data.y_cyt(2:end);
nuc_data = data.nuc_data;
cyt_data = data.cyt_data;
locs = data.locs;


mean_oth_nuc = mean(fillmissing(data.y_nucall(:,locs),'linear',1),2);
mean_oth_cyt = mean(fillmissing(data.y_cytall(:,locs),'linear',1),2);


% changing the multiplication to initial value of mean
ynuc = (data.y_nuc(2:end)./mean_oth_nuc(2:end)).*mean_oth_nuc(1);
ycyt = (data.y_cyt(2:end)./mean_oth_cyt(2:end)).*mean_oth_cyt(1);



%
% Smooth the cyt data
%
smoothwindow = 10;
ycyt_fit = smooth(t,ycyt,smoothwindow);
	

% =========================================================================
% Fit our nuclear data to the solved differential equation
% =========================================================================

%
% Make options
% p = [c0,k_in,k_out]
Lower = [min(ynuc(:))/2 0 0];
Upper = [2*max(ynuc(:)) Inf Inf];
StartPoint = [ynuc(1) 0.01 0.01];

tols = 1e-16;
% options = optimset('Display','off',...
% 	'MaxFunEvals',5e2,...
% 	'MaxIter',5e2,...
% 	'TolX',tols,...
% 	'TolFun',tols,...
% 	'TolCon',tols ,...
% 	'UseParallel','always');

options = optimset('Display','off',...
	'MaxFunEvals',5e2,...
	'MaxIter',5e2,...
	'TolX',tols,...
	'TolFun',tols,...
	'TolCon',tols);

%
% Do the fit.  If there are outliers, remove them and re-perform
% the fit.  This only works if there are only a handful of
% outliers.  And we only do this once.
%
fhandle = @ftn_FRAP; % y = A*(1 - exp(-k*t)) + B;
for jj = 1:2
	ynuc_smooth = smooth(t,ynuc,smoothwindow);
	[P,resnorm,residual,~,~,~,J] = ...
		lsqcurvefit(fhandle,StartPoint,t,ynuc_smooth,...
		Lower,Upper,options,ycyt_fit,experimenttype);
	
	m = mean(residual);
	s = std(residual);
	v = residual > m + 3*s | residual < m - 3*s;
	if jj == 1 && sum(v) > 0
		t = t(~v);
		ynuc = ynuc(~v);
		ycyt_fit = ycyt_fit(~v);
        ycyt=ycyt(~v);
	else
		continue
	end
end
ynuc_smooth = smooth(t,ynuc,smoothwindow);
yfit = ynuc_smooth - residual;
% y = ftn_FRAP(P,t,t,ycyt_fit,'bleach');
c0 = P(1);
k_in = P(2);
k_out = P(3);

%
% Calc goodness of fit
%
R2 = 1 - resnorm/norm(ynuc-mean(ynuc))^2;

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

% 95% errorbar
alpha = 0.05;
delta2 = se*tinv(1-alpha/2,nu);


%
% Output
%


figure('visib','off')
% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
plot(t,yfit,'linewidth',2)
hold on
plot(t,ycyt_fit,'linewidth',2)
hold on
plot(t,ynuc,'m--o')
hold on
plot(t(1:end),ycyt(1:end),'r--o')
xlabel(['k_{out} = ',num2str(k_out),...
		',   k_{in} = ',num2str(k_in),',   r^2 = ',num2str(R2)])
legend('nuclear int.','cyt. int.','actual nuc values','actual cyt values')
%saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'recovery_',filenameshort(1:end-4),'.fig']);
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'recovery1_norm_',filenameshort(1:end-4),'.jpg']);
saveas(gcf, fullfile('Image results', 'Neighbor nuc stats', datum, filenameshort, ...
['recovery1_norm_', filenameshort(1:end-4),'.jpg']));    

figure('visib','off'),
subplot(2,1,1)
plot(t,yfit,'linewidth',2)
hold on
plot(t,ynuc,'m--o')
xlabel(['k_{out} = ',num2str(k_out),...
		',   k_{in} = ',num2str(k_in),',   r^2 = ',num2str(R2)])
legend('nuclear int.','actual nuc values')	

subplot(2,1,2)
plot(t,ycyt_fit,'linewidth',2)
hold on
plot(t(1:end),ycyt(1:size(t,1)),'r--o')
xlabel(['k_{out} = ',num2str(k_out),...
		',   k_{in} = ',num2str(k_in),',   r^2 = ',num2str(R2)])
legend('cyt. int.','actual cyt values')	
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'recovery2_norm_',filenameshort(1:end-4),'.jpg']);
saveas(gcf, fullfile('Image results', 'Neighbor nuc stats', datum, filenameshort, ...
['recovery2_norm_', filenameshort(1:end-4),'.jpg']));    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.yfit = yfit;
data.c0 = c0;
data.k_in = k_in;
data.k_out = k_out;
data.r2 = R2;
data.delta1 = delta1;
data.delta2 = delta2;
data.covar = covar_matrix;

%
% plotting invisibly, if asked for
%
if yesplot
	filename = data.filename;
	if exist([filename(1:end-4),'.fig'],'file')
		open([filename(1:end-4),'.fig'])
	else
		figure
	end
	
	set(gcf,'visib','off','paperpositionmode','auto','pos',[20 355 1478 360]);
	
	subplot(1,3,3)
	hold on
	plot(t,[yfit,ycyt_fit],'linewidth',2)
	xlabel(['k_{out} = ',num2str(k_out),...
		',   k_{in} = ',num2str(k_in),',   r^2 = ',num2str(R2)])

	print(gcf,[filename(1:end-4),'.jpg'],'-djpeg','-r300')
	saveas(gcf,[filename(1:end-4),'.fig'])
	close
end


