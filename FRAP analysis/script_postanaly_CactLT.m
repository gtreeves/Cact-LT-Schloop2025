% script_postanaly_CactLT
%
% What to do after we are done with the analysis..make boxplot, etc
%
% Running this script should reproduce the FRAP plots from the paper

clear
close all
plotviolin = false;
boxplotcolors = [[0 0.4471 0.7412];[0 0.4471 0.7412];
	[0.8510 0.3255 0.0980];[0.8510 0.3255 0.0980];
	[0.9294 0.6941 0.1255];[0.9294 0.6941 0.1255];
	[125,46,143]/255;[125,46,143]/255;
	[0.4660 0.6740 0.1880]; [0.4660    0.6740    0.1880];
	];


load Mat/Cnucbar_out Cnucbar_out Ccytbar_out
load Mat/1componentfits Soln_1cpt
load Mat/2componentfits Soln_2cpt


%
% Extract variables
%
% for 1 component
matfiles1 = {Soln_1cpt.matfilename}';
kin1 = [Soln_1cpt.k_in]';
kout1 = [Soln_1cpt.k_out]';
Knuc1 = [Soln_1cpt.K_nuc]';
r21 = [Soln_1cpt.r2]';
delta11 = [Soln_1cpt.delta1]';
dlogkout1 = delta11(:,2)./kout1;

% for 2 component
matfiles2 = {Soln_2cpt.matfilename}';
kin2 = [Soln_2cpt.k_in]';
kout2 = [Soln_2cpt.k_out]';
Knuc2 = [Soln_2cpt.K_nuc]';
phi2 = [Soln_2cpt.phi]';
phicyt2 = [Soln_2cpt.phicyt]';
kinG2 = [Soln_2cpt.kinG]';
koutG2 = [Soln_2cpt.koutG]';
NCR = Cnucbar_out./Ccytbar_out;
NCR2 = [Soln_2cpt.Cnucbar]'./[Soln_2cpt.Ccytbar]';
r22 = [Soln_2cpt.r2]';
delta12 = [Soln_2cpt.delta1]';
dlogkout2 = delta12(:,1)./kout2;

%
% Load-in Sharva's analysis of Dl-GFP from Etika, and extract those vbls
%
load('Mat/2024-06-14 import_export_normalized rates')
kin_Dl = k2_in_d'*60;
kout_Dl = k2_out_d'*60;

%
% Filtering
%
v1x = contains(matfiles1,'Normal GFP'); % 1x1x embryos
v2G = contains(matfiles1,'2 GFP');
vfree = contains(matfiles1,'free GFP');
cvmax1 = 0.3; 
cvmax2 = 1.5; 
r2min = 0.4; % poor-fit filter
vr1 = r21 > r2min;
vr2 = r22 > r2min;
vcv1 = dlogkout1 < cvmax1;
vcv2 = dlogkout2 < cvmax2;


%% ========================================================================
% Boxplot of kin,kout by genotype
% =========================================================================
% {


vplot1 = vr1 & vcv1;
vplot2 = vr2 & vcv2;
G_list = {'kin1x','kout1x','kinG','koutG','kin2x','kout2x','kin2c','kout2c','kinDl','koutDl'};
C = {kin1(v1x & vplot1) kout1(v1x & vplot1) ...
	kin1(vfree & vplot1) kout1(vfree & vplot1) ...
	kin1(v2G & vplot1) kout1(v2G & vplot1) ...
	kin2(vplot2) kout2(vplot2) ...
	kin_Dl kout_Dl ...
	};
figure('pos',[390         282        1035         423])
boxplotDU(C,G_list,0.2,boxplotcolors,12,plotviolin,false);
set(gca,'yscale','log')
ylabel('k_{in}  or  k_{out} [min^{-1}]')

ylim([0.01 100])
xlim([0.5 length(C)+0.5])
set(gca,'ytick',10.^([-2:2]))

% save Mat/kin_kout_Cact_w_Dl G_list C

for i = 1:length(C)
	fprintf([G_list{i},' = %f +/- %f, n = %i\n'],mean(C{i}),std(C{i}),length(C{i}))
end
fprintf('\n')

%}

%% ========================================================================
% Boxplot of phinuc,phicyt (2-cpt model for 1x1x embryos)
% =========================================================================
% {


vplot2 = vr2 & vcv2;
G_list = {'psinuc','psicyt'};
C2 = {phi2(vplot2) phicyt2(vplot2) ...
	};
boxplotDU(C2,G_list,0.2,zeros(2,3),12,plotviolin);
ylabel('\psi_{nuc}  or  \psi_{cyt}') % supposed to be psi
ylim([0 1])

set(gca,'ytick',0:0.2:1)
% title('\phi')

fprintf('phi(1x1x) = %f +/- %f, n = %i\n',mean(C2{1}),std(C2{1}),length(C2{1}))
fprintf('phicyt(1x1x) = %f +/- %f, n = %i\n',mean(C2{2}),std(C2{2}),length(C2{2}))
fprintf('\n')

%}



%% ========================================================================
% Boxplot of Knuc by genotype
% =========================================================================
% {


vplot1 = vr1 & vcv1;
vplot2 = vr2 & vcv2;
G_list = {'Knuc1x','KnucG','Knuc2x','Knuc2c'};
C3 = {Knuc1(v1x & vplot1) ...
	Knuc1(vfree & vplot1) ...
	Knuc1(v2G & vplot1) ...
	Knuc2(vplot2) ...
	};
boxplotDU(C3,G_list,0.2,boxplotcolors(1:2:end,:),12,plotviolin);
ylabel('K_{nuc}')
ylim([0 1.4])

for i = 1:length(C3)
	fprintf([G_list{i},' = %f +/- %f, n = %i\n'],mean(C3{i}),std(C3{i}),length(C3{i}))
end
fprintf('\n')

set(gca,'ytick',0:0.2:1.4)

%
% Checking on propagation of error
%
for i = 1:length(C3)
	kin = mean(C{2*i-1});
	dkin = std(C{2*i-1});
	kout = mean(C{2*i});
	dkout = std(C{2*i});
	cstd = sqrt(cov(C{2*i-1},C{2*i}));
	cstd = cstd(1,2);
	dKnuc = kin/kout.*sqrt((dkin/kin).^2 + (dkout/kout).^2 - 2*cstd^2/(kin*kout));

	fprintf(['For ',G_list{i},', std = %f; prop err = %f\n'],std(C3{i}),dKnuc)
end
fprintf('\n')

%}

%% ========================================================================
% Boxplot of NCR by genotype
% =========================================================================
% {


vplot1 = vr1 & vcv1;
vplot2 = vr2 & vcv2;
G_list = {'NCR1x','NCRG','NCR2x','NCR2c'};
C4 = {NCR(v1x & vplot1) ...
	NCR(vfree & vplot1) ...
	NCR(v2G & vplot1) ...
	NCR(vplot2) ...
	};
boxplotDU(C4,G_list,0.2,boxplotcolors(1:2:end,:),12,plotviolin);
ylabel('NCR')
ylim([0 1.4])

for i = 1:length(C4)
	fprintf([G_list{i},' = %f +/- %f, n = %i\n'],mean(C3{i}),std(C3{i}),length(C3{i}))
end
fprintf('\n')

set(gca,'ytick',0:0.2:1.4)

%}



%% ========================================================================
% Plot of recovery curve for one specific embryo
% =========================================================================
% {

%
% Pick the right file
%
matfilenames1 = {Soln_1cpt.matfilename}';
matfilenames2 = {Soln_2cpt.matfilename}';
for i = 3%1:length(Soln2)
	
	matfilename2 = matfilenames2{i};

	idx1 = contains(matfilenames1,matfilename2);
	idx2 = i;

	%
	% Extract data from the 1-cpt model
	%
	smoothwindow = 10;
	[t1,ynuc_smooth1,ycyt_smooth1,ynuc1,ycyt1] = smooth_ynuc(Soln_1cpt(idx1),smoothwindow);
	P1 = [Soln_1cpt(idx1).c0 Soln_1cpt(idx1).k_in Soln_1cpt(idx1).k_out];
	yfit1 = ftn_FRAP(P1,t1,ycyt_smooth1);

	P2 = [Soln_2cpt(idx2).k_out Soln_2cpt(idx2).koutG Soln_2cpt(idx2).phi];
	yfit2 = ftn_fit_FRAP_2cpt(P2,t1,ycyt_smooth1,ynuc1(1),...
		Soln_2cpt(idx2).Cnucbar,Soln_2cpt(idx2).Ccytbar,Soln_2cpt(idx2).KnucG);

	%
	% Plotting
	%
	figure('visib','off')
	t = [Soln_1cpt(idx1).t(1); t1]; % pre-bleach data point
	y = [Soln_1cpt(idx1).y_nuc(1); ynuc1];
	ys = [Soln_1cpt(idx1).y_nuc(1); ynuc_smooth1];

	plot(t,y,'.','Markersize',6)
	hold on
	plot(t1,yfit1)
	plot(t1,yfit2)


end

%}


