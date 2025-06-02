% script_postanaly_cact
%
% This script loads the tertiary data from RICS (the fitted parameters) and
% generates plots

clear
close all
% load Mat/2022-07-27_Solnnew
load Mat/2023-02-13_Solnall_Solnonly
Soln1 = Soln;
load Mat/2023-05-31_10-04-30_2xCact-LT_GFP_Solnonly % 2x LT
Soln2 = Soln;
load Mat/2023-05-31_11-27-47_Cact-LT_2xGFP_Solnonly % 2x GFP
Soln3 = Soln;
load Mat/2023-06-07_23-27-44_free_GFP_all_Solnonly % free GFP only, with hist
Soln = [Soln1;Soln2;Soln3;Soln];

% load Mat/2023-05-31_2c_diffusion_fit_all phi_out R2
load Mat/2024-01-22_23-41-00_2c_diffusion_fit_all phi_out R2
nsoln = length(Soln);

%
% Correct for phi_out and R2, which were loaded in separately
%
if length(phi_out) < nsoln
	phi_out(end+1:nsoln) = 0;
	R2(end+1:nsoln) = 0;
end

%
% Autofluorescence
%
load Mat\autofluorescence_RICS nuc_autofluorescence_RICS cyt_autofluorescence_RICS

plotviolin = false;
boxplotcolors = [[0 0.4471 0.7412];[0.8510 0.3255 0.0980];[0.9294 0.6941 0.1255];[125,46,143]/255;
	[0 0.4471 0.7412];[0.8510 0.3255 0.0980];[0.9294 0.6941 0.1255]];

%%

% {

w00 = 0.25;
wz = 3*w00;

fnames = fieldnames(Soln);
v = strfindDU(fnames,'fit');
fnames = fnames(v);
append = cell(1,5);
for i = 1:5
	append{i} = fnames{i}(4:end);
end
idx = 1:length(append);

%
% Remove the empty data sets
%
filenames = {Soln.filename}';
filenamesshort1 = {Soln.filenameshort}';
v = cellfun(@isempty,filenames);
Soln(v) = [];
filenames(v) = [];
phi_out(v) = [];
R2(v) = [];


P = {Soln.pth}';
G = {Soln.genotype}';
Itot = [Soln.totsignal]';
Inuc = [Soln.nucsignal]' - nuc_autofluorescence_RICS; % autofluorescence subtract
Icyt = [Soln.cytsignal]' - cyt_autofluorescence_RICS; % autofluorescence subtract
sigtot = sqrt([Soln.totvar]');
signuc = sqrt([Soln.nucvar]');
sigcyt = sqrt([Soln.cytvar]');
NCRI = Inuc./Icyt; % nuclear to cytoplasmic ratio of intensity

%
% Create filenamesshort, including end node folder, and pthshort, which is
% just the end node folder
%
pthshort = cell(nsoln,1);
filenamesshort = cell(nsoln,1);
for i = 1:nsoln
	filename = filenames{i};
	kslash = strfind(filename,filesep);
	filename(kslash) = '_';
	filenamesshort{i} = filename(kslash(end-1)+1:end);
	pthshort{i} = filename(kslash(end-1)+1:kslash(end)-1);
end


%
% Extract tertiary data
%
w01 = zeros(nsoln,length(append));
w02 = zeros(nsoln,length(append));
D = zeros(nsoln,length(append));
A = zeros(nsoln,length(append));
dA = zeros(nsoln,length(append));
B = zeros(nsoln,length(append));
dB = zeros(nsoln,length(append));
dD = zeros(nsoln,length(append));
dphi = zeros(nsoln,length(append));
dphicyt = zeros(nsoln,length(append));
R20 = zeros(nsoln,length(append)); % R2 values for the PSF step (fast dir)
R21 = zeros(nsoln,length(append)); % R2 values for the 1-cpt fit
R22 = zeros(nsoln,length(append)); % R2 values for the 2-cpt fit
for j = 1:length(append)
	A(:,j) = [Soln.(['A',append{j}])]';
	w01(:,j) = [Soln.(['w01',append{j}])]';
	R20(:,j) = [Soln.(['R2_PSF',append{j}])]';

	
	a = [Soln.(['fit',append{j}])]';
	E_PSF = [a.errorbar68_PSF]';
	dA(:,j) = E_PSF(:,1);

	if strcmp(append{j},'cc')
	else
		D(:,j) = [Soln.(['D',append{j}])]';
		w02(:,j) = [Soln.(['w02',append{j}])]';
		R21(:,j) = [Soln.(['R2_1',append{j}])]';
		R22(:,j) = [Soln.(['R2_2',append{j}])]';

		E_1 = [a.errorbar68_1]';
		E_2 = [a.errorbar68_2]';
		dB(:,j) = E_PSF(:,2);
		dD(:,j) = E_1(:,3);
		dphi(:,j) = E_2(:,4);

	end
end
phi = [Soln.phinuc]';
phicyt = [Soln.phicyt]';
dxcc = [Soln.dxcc]';
dycc = [Soln.dycc]';

%
% Concentrations in particles per micron cubed, then convert to nM
%
gamma = sqrt(2)/4; % Brown2008 has this as sqrt(2)/4
V_psf = pi^1.5*w00^2*wz;
V_eff = V_psf/gamma;
C = 1./(V_eff*A)/ 0.602;
dC = C.*dA./A; % Propagation of error: dC = C*(dA/A)
NCRC = C(:,2)./C(:,3); % nuclear to cytoplasmic ratio

%
% Gathering data on the microscope parameters. This is important because we
% didn't hold LP or taup/tauL/dr consistent.
%
metadata = [Soln.metadata]';
lsminf2 = [metadata.lsminf2]';
data_ch = [metadata.data_ch]';
mask_ch = [metadata.mask_ch]';
cziinfo = [lsminf2.cziinfo]';
LP = [cziinfo.Laserpower]';
LP = reshape(LP,2,nsoln);
LP = LP(1,:)'.*(data_ch == 1) + LP(2,:)'.*(data_ch == 2);
DG = [cziinfo.Detectorgain]';
DG = reshape(DG,2,nsoln);
DG = DG(1,:)'.*(data_ch == 1) + DG(2,:)'.*(data_ch == 2);
taup = [metadata.taup]';
tauL = [metadata.tauL]';
dr = [metadata.dr]';
H = [metadata.H]';
W = [metadata.W]';

sig0 = [metadata.std_bg]';
sig0 = reshape(sig0,2,nsoln);
sig0 = sig0(1,:)'.*(data_ch == 1) + sig0(2,:)'.*(data_ch == 2);

%
% Calc variances, epsilon, S-factor 
%
signtot = sqrt(A(:,1).*Itot.^2);
signnuc = sqrt(A(:,2).*Inuc.^2); 
signcyt = sqrt(A(:,3).*Icyt.^2);
sigdtot = sqrt(sigtot.^2 - signtot.^2);
sigdnuc = sqrt(signuc.^2 - signnuc.^2);
sigdcyt = sqrt(sigcyt.^2 - signcyt.^2);
sigstot = sqrt(sigdtot.^2 - sig0.^2);
sigsnuc = sqrt(sigdnuc.^2 - sig0.^2);
sigscyt = sqrt(sigdcyt.^2 - sig0.^2);
epcyt = signcyt.^2./sigscyt.^2;%./LP*0.01;
epnuc = signnuc.^2./sigsnuc.^2;%./LP*0.01;
Stot = sigstot.^2./Itot;
Snuc = sigsnuc.^2./Inuc;
Scyt = sigscyt.^2./Icyt;


%
% Calc apparent number, N, from the RICS data
%
Nnuc_RICS = Inuc.^2./(signuc.^2 - sig0.^2);
Ncyt_RICS = Icyt.^2./(sigcyt.^2 - sig0.^2);


%
% Filters
%
gofcutoff = 0.85;
vA = dA./A < 0.05; % errorbars in amplitude must be 5% or less
ve = dD <= 4;%1.2;%3.5;%1.15; % goodness-of-fit cutoff as assayed by errorbar
vR2 = R20 > 0.00 & R21 > gofcutoff & R22 > 0.00; % goodness-of-fit cutoff
vR2_phi = R2 > gofcutoff;


genotype{1} = 'Cact-LT/CyO; mat-GFP/H2A-RFP';
genotype{2} = 'cact-LT;H2A-RFP/mat-mGFP';
genotype{3} = 'cact-LT/H2A-RFP;mat-mGFP';
genotype{4} = {'H2A-RFP/CyO;mat-mGFP/mat-GFP';'H2A-RFP/+;mat-GFP/+'};
wt = strfindDU(G,genotype{1});
LT_2x = strfindDU(G,genotype{2});
GFP_2x = strfindDU(G,genotype{3});
GFP_only = strfindDU(G,genotype{4}{1}) | strfindDU(G,genotype{4}{2});
Gtype = {wt LT_2x GFP_2x GFP_only};
G_list = {'wt' '2x LT' '2x GFP' 'free GFP'};

G2 = G;
for i = 1:length(G_list)
	G2(Gtype{i}) = G_list(i);
end

G_list2 = G_list; G_list2{1} = '1x1x';

%}

%% ========================================================================
% Read-in the embryos to exclude from the excel-annotated data. These
% embryos had bad masks, old nuclei, or, in one case, looked like it was nc
% 15 somehow.
% =========================================================================
%{

exclude = cell(1,length(Gtype));
v_ex = false(nsoln,1);
for i = 1:length(Gtype)
	[numdata,txtdata] = xlsread('RICS_stats_Cact.xlsx',G_list{i});
	exclude{i} = numdata(:,23);
	v_ex(Gtype{i}) = ~isnan(numdata(:,23));
end
v_ex = v_ex | tauL > 4; % also exclude embryos with too-large tauL

save Mat/cact_filters v_ex vA vR2 vR2_phi wt

%}

load Mat/cact_filters v_ex vA vR2 vR2_phi wt



%% ========================================================================
% Boxplots of Dnuc,Dcyt, for 1x,1x embryos, then by genotype
% =========================================================================
% {


vplot = vR2(:,2) & vR2(:,3) & ~v_ex;
C0 = {D(wt & vplot,2) D(wt & vplot,3)};
boxplotDU(C0,{'nuc','cyt'},0.2,boxplotcolors,12,plotviolin);
title('D')


vplot = vR2(:,2) & ~v_ex;
C1 = {D(wt & vplot,2) D(LT_2x & vplot,2) D(GFP_2x & vplot,2) D(GFP_only & vplot,2)};
boxplotDU(C1,G_list,0.2,boxplotcolors,12,plotviolin);
title('D_{nuc}')
disp('Nuclear: average D, then s.d.')
disp(cellfun(@mean,C1))
disp(cellfun(@std,C1))


vplot = vR2(:,3) & ~v_ex;
C2 = {D(wt & vplot,3) D(LT_2x & vplot,3) D(GFP_2x & vplot,3) D(GFP_only & vplot,3)};
boxplotDU(C2,G_list,0.2,boxplotcolors,12,plotviolin);
title('D_{cyt}')
disp('Cytoplasmic: average D, then s.d.')
disp(cellfun(@mean,C2))
disp(cellfun(@std,C2))


for i = 1:length(C1)
	fprintf(['n (',G_list2{i},', nuclear) = %d\n'],length(C1{i}))
	fprintf(['n (',G_list2{i},', cytoplasmic) = %d\n'],length(C2{i}))
end


%}


%% ========================================================================
% Plots of average ACF by genotype
% =========================================================================
% {

load Mat/ACFslowdir Gnuc Gcyt
etalimit = size(Gnuc,1);
Gnuc0 = A(:,2)';
Gnuchat = Gnuc./Gnuc0;
Gnuchat(1,:) = 1;
Gcyt0 = A(:,3)';
Gcythat = Gcyt./Gcyt0;
Gcythat(1,:) = 1;
vbad = false;

vplot = vR2(:,2) & vR2(:,3) & ~v_ex;


%
% all embryos
%
for i = 1:length(Gtype)
	
	figure
	h1 = plot(0:(etalimit-1),Gnuchat(:,Gtype{i} & vplot & ~vbad),'k');
	hold on
	h2 = plot(0:(etalimit-1),Gcythat(:,Gtype{i} & vplot & ~vbad),'b');
	legend([h1(1) h2(1)],'nuclear','cytoplasmic')
	ylabel('G/G0')
	xlabel('\Deltay')
	title(G_list2{i})

end


%
% avg +/- sd
%
for i = 1:length(Gtype)

	%
	% Make average curves
	%
	Gnucavg = mean(Gnuchat(:,Gtype{i} & vplot & ~vbad),2);
	Gcytavg = mean(Gcythat(:,Gtype{i} & vplot & ~vbad),2);

	% ----------------------------------------------
	% Below I used to fit. Now just plot the avg D
	% ----------------------------------------------
	%
	% Set up the startpoint and boundaries
	%
	A1 = 1; % A is the amplitude.
	B1 = 0; % "B" is the background
	D1 = 1; % micron^2/s
	DU = 30;
	DL = 1e-4;
	w0 = median(w00);

	Lower = [A1-1e-6 -1e-6 DL w0*(1-1e-6)];
	Upper = [A1+1e-6 1e-6 DU w0*(1+1e-6)];
	StartPoint = [A1 B1 D1 w0];

	%
	% Set up fittype and options.
	%
	tols = 1e-16;
	options = optimset('Display','off',...
		'MaxFunEvals',5e2,...
		'MaxIter',5e2,...
		'TolX',tols,...
		'TolFun',tols,...
		'TolCon',tols ,...
		'UseParallel',false);
	fhandle = @ftn_fit_1comp_diffusion;
	XI = {0 (0:(etalimit-1))'};

	%
	% perform the fittings
	%
	p = lsqcurvefit(fhandle,StartPoint,XI,Gnucavg(2:end),Lower,Upper,options,...
		median(dr),median(wz),median(taup)*1e-6,median(tauL)*1e-3,true);
	Fnuc = ftn_fit_1comp_diffusion(p,XI,median(dr),median(wz),...
		median(taup)*1e-6,median(tauL)*1e-3,false);
	Dnuc = p(3);


	p = lsqcurvefit(fhandle,StartPoint,XI,Gcytavg(2:end),Lower,Upper,options,...
		median(dr),median(wz),median(taup)*1e-6,median(tauL)*1e-3,true);
	Fcyt = ftn_fit_1comp_diffusion(p,XI,median(dr),median(wz),...
		median(taup)*1e-6,median(tauL)*1e-3,false);
	Dcyt = p(3);

	%
	% Fit 1-cpt model to average curves
	%
	figure
	h1 = errorbar(0:(etalimit-1),Gnucavg,stdDU(Gnuchat(:,Gtype{i} & vplot & ~vbad),2),'.k');
	hold on
	plot(0:(etalimit-1),Fnuc,'k')

	h2 = errorbar(0:(etalimit-1),Gcytavg,stdDU(Gcythat(:,Gtype{i} & vplot & ~vbad),2),'db');
	plot(0:(etalimit-1),Fcyt,'b')

	legend([h1 h2],'nuclear','cytoplasmic')
	ylabel('G/G0')
	xlabel('\Deltay')
	title(G_list2{i})
	YLIM = ylim;
	ylim([0 YLIM(2)])

	disp(['nuclear diffusivity for ',G_list2{i},': ',num2str(Dnuc)])
	disp(['cytoplasmic diffusivity for ',G_list2{i},': ',num2str(Dcyt)])

end

%}




%% ========================================================================
% Boxplot of phi by genotype
% =========================================================================
% {

vplot = vR2(:,2) & vR2_phi & ~v_ex;
C3 = {phi_out(wt & vplot) phi_out(LT_2x & vplot) phi_out(GFP_2x & vplot)};
CC = {NCRI(wt & vplot) NCRI(LT_2x & vplot) NCRI(GFP_2x & vplot)};
boxplotDU(C3,G_list,0.2,boxplotcolors,12,plotviolin);
ylim([0 1])
set(gca,'ytick',0:0.2:1)
title('\phi')

for i = 1:length(C3)

	fprintf(['phi(',G_list{i},') = %f +/- %f (n = %d)\n'],...
		mean(C3{i}),std(C3{i}),length(C3{i}))
end

%}









