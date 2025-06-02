function [data,F] = ftn_gradfit7(data,filenum,cyt_data,varargin)
%Fits Dorsal/Nuc data (from lish_border) to a time/z-series of Gaussians
%
%function [data,F] = lish_gradfit(data,varargin)
%
% This function reads in the data output from "lish_border.m" and fits
% them to determine a,b,sig, and mu for each timepoint and z-slice (see
% Liberman et al, 2009).  The best fit mu for the whole time/z-series
% becomes "s_mid", which is considered the ventral midline for each
% timepoint. Then these four fields (A,B,Sig,s_mid) are appended to the end
% of "data", which is saved back to the mat file called:
%	[data.pth,'data.mat']
%
% "data": structure containing the data and metadata.
%
% Optional argument varargin can consist of these things, in this order:
%	(1) "yesplot": whether you want to plot the outcome.  Default, "false".
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% "data": the output structure is the same as the input, plus the fit
%	parameters
% "F": the movie of the plots, only made if asked for by "yesplot"

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end%, iArg = iArg + 1;

% tstart = data.tstart;
% zstart = data.zstart;
% tend = data.tend;
% zend = data.zend;
% T = tend - tstart + 1;
% if data.thstart ~= data.thend, error('not ready for theta stack'), end
% D = zend - zstart + 1;
T = data.D;
D = 1;
dz = 0;
bcg = data.bcg;

S = data.S;
RatioS = data.R2;
S_cyt = cyt_data.S_cyt;
RatioS_cyt = cyt_data.R2;

dt = data.dt; % timestep in seconds.
z0 = data.z0; % how deep the first slice is
A = zeros(T,D); B = A; Mu = A; Sig = A; Rsquare = A;
dA = A; dB = A; dSig = A;
st = A; mo = A;
Rsmooth = cell(T,D);
x = linspace(-1,1);
t = linspace(0,(T-1)*dt,T)'/60; % time vector in minutes
z = linspace(z0,z0+(D-1)*dz,D)'; % depth in microns

if data.filenameshort == "Time course nc10 or 11 to gastrulation.czi"
        s_mid = .3727; % Time course nc10 or 11 to gastrulation.czi
        
    elseif data.filenameshort == "nc 12 to gastrulation.czi"
        s_mid = .8178; % nc 12 to gastrulation
        
    elseif data.filenameshort == "Embryo 1 nc12-gastrulation.czi"
        s_mid = .3594; % Embryo 1 nc12-gastrulation
        
    else 
        s_mid = .0087; % Embryo 2 nc11-gastrulation
       
end

r_mean_cyt = [];
r_mean_nuc = [];

%
% Plotting and making movie, if asked for.
%
	for i = 1:T
		for j = 1:D
			
			%
			% Extract and plot points
			%
			s_cyt = mod(S_cyt{i,j}-s_mid+1,2) - 1;
            [s_cyt,isort] = sort(s_cyt);
            r_cyt = RatioS_cyt{i,j};
            r_cyt = r_cyt(isort);
            
            % mirroring the pts from lhs(i.e s<0 onto the right side to check for slope...
            ... this is meant to see if there is bilateral symmetry)
            indxs = find(s_cyt<0); % indexes of all the pts below s<0
            midpt = max(indxs); % the last index of where s is < 0

            s_cyt_temp = zeros(2*midpt,1);
            s_cyt_temp(1:midpt) = s_cyt(1:midpt);
            s_cyt_temp(midpt+1:end) = sort(-s_cyt(1:midpt));

%             s_cyt = s_cyt_temp;
            s_cyt = s_cyt_temp(midpt+1:end);

            r_cyt_temp = zeros(2*midpt,1);
            r_cyt_temp(1:midpt) = r_cyt(1:midpt);
            r_cyt_temp(midpt+1:end) = flipud(r_cyt(1:midpt));

%             r_cyt = r_cyt_temp;
            r_cyt = r_cyt_temp(midpt+1:end);
           
            
            Y_cyt = prctile(r_cyt,[20 80]); % remove top and bottom 5 outliers
%             Y_cyt = [225-bcg 275-bcg];
            v = r_cyt > Y_cyt(1) & r_cyt < Y_cyt(2);
            r_new_cyt = r_cyt(v);
            s_new_cyt = s_cyt(v);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            s = mod(S{i,j}-s_mid+1,2) - 1;
            [s,isort] = sort(s);
            r = RatioS{i,j};
            r = r(isort);
            
            
            % mirroring the pts from lhs(i.e s<0 onto the right side to check for slope...
            ... this is meant to see if there is bilateral symmetry)
            indxs = find(s<0); % indexes of all the pts below s<0
            midpt = max(indxs); % the last index of where s is < 0

            s_temp = zeros(2*midpt,1);
            s_temp(1:midpt) = s(1:midpt);
            s_temp(midpt+1:end) = sort(-s(1:midpt));

%             s = s_temp;
            s = s_temp(midpt+1:end);

            r_temp = zeros(2*midpt,1);
            r_temp(1:midpt) = r(1:midpt);
            r_temp(midpt+1:end) = flipud(r(1:midpt));

%             r = r_temp;
            r = r_temp(midpt+1:end);
           
            

            Y = prctile(r,[20 80]); % remove top and bottom 5 outliers
%             Y = [125-bcg 200-bcg];
            v = r > Y(1) & r < Y(2);
            r_new = r(v);
            s_new = s(v);
            
            [slope1,intercept1,stats_cyt,h_p1] = linlsq(s_new_cyt,r_new_cyt);
            [slope2,intercept2,stats_nuc,h_p2] = linlsq(s_new,r_new);

        end
%         
        r_mean_cyt(i) = mean([intercept1,slope1+intercept1]);
        r_mean_nuc(i) = mean([intercept2,slope2+intercept2]);
        sig_nuc(i) = std(r_new);
        sig_cyt(i) = std(r_new_cyt);
        
    end
    
colors1 = cbrewer('seq', 'PuRd', 9); 
colors2 = cbrewer('seq', 'PuBuGn', 9);

if yesplot
    
    % --- Create Output Folder Structure ---
    figBaseFolder = 'Figures';
    figSubFolder = fullfile(figBaseFolder, data.filenameshort);
    
    if ~exist(figSubFolder, 'dir')
        mkdir(figSubFolder);
    end
    
    
    figure
	set(gcf,'Position',[100 100 400 300],...
		'Paperpositionmode','auto','Color',[1 1 1])
	plot(t,r_mean_cyt,'Color',colors1(7, :),'LineWidth',2)
    hold on 
    plot(t,r_mean_nuc,'Color',colors2(7, :),'LineWidth',2)
    xlabel('Time(in min)')
    ylabel('mean intensity')
    legend('mean cyt int.','mean nuc int.','Location','best')
    
    % --- Save the Figure ---
    figName = fullfile(figSubFolder, [data.filenameshort(1:end-4), '_mean_intensities']);
    saveas(gcf, [figName, '.jpg']);
    saveas(gcf, [figName, '.fig']);
	
    figure
	set(gcf,'Position',[100 100 400 300],...
		'Paperpositionmode','auto','Color',[1 1 1])
	plot(t,r_mean_cyt,'Color',colors1(7, :),'LineWidth',2)
    hold on 
    plot(t,r_mean_nuc,'Color',colors2(7, :),'LineWidth',2)
    hold on
    patch([t' fliplr(t')], [(r_mean_cyt-sig_cyt) fliplr((r_mean_cyt+sig_cyt))],colors1(5, :), 'FaceAlpha', .2)
    hold on
    patch([t' fliplr(t')], [(r_mean_nuc-sig_nuc) fliplr((r_mean_nuc+sig_nuc))],colors2(5, :), 'FaceAlpha', .2)
    xlabel('Time(in min)')
    ylabel('mean intensity')
    legend('mean cyt int.','mean nuc int.','\mu_{cyt}\pm\sigma','\mu_{nuc}\pm\sigma','Location','best')
    
    % --- Save the Figure ---
    figName = fullfile(figSubFolder, [data.filenameshort(1:end-4), '_mean_intensities_std']);
    saveas(gcf, [figName, '.jpg']);
    saveas(gcf, [figName, '.fig']);
    
    F(i,j) = getframe(gcf); 
    
    % Define the output folder path
    baseFolder = fullfile('Figures', 'Averaged embryo figures and data');
    
    % Create the folder if it doesn't exist
    if ~exist(baseFolder, 'dir')
        mkdir(baseFolder);
    end
    
    eval(['data' num2str(filenum) '.t' ' = t;']);
    eval(['data' num2str(filenum) '.mean_nuc' ' = r_mean_nuc;']);
    eval(['data' num2str(filenum) '.mean_cyt' ' = r_mean_cyt;']);
    eval(['data' num2str(filenum) '.std_nuc' ' = sig_nuc;']);
    eval(['data' num2str(filenum) '.std_cyt' ' = sig_cyt;']);
    
    % Create folder
    baseFolder = fullfile('Figures', 'Averaged embryo figures and data');
    if ~exist(baseFolder, 'dir')
        mkdir(baseFolder);
    end
    
    % Save it
    eval(['save(fullfile(baseFolder, ''data' num2str(filenum) '.mat''), ''data' num2str(filenum) ''');']);
    
else
	F = 'no plot asked for';	
end

end

