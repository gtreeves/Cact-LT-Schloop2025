function [data,F] = ftn_gradfit8_plot_slope(data,cyt_data,varargin)
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

% Plotting and making movie, if asked for.

pval_nuc = [];
pval_cyt = [];

if yesplot
	figure
% 	set(gcf,'Position',[100 100 400 300],...
% 		'Paperpositionmode','auto','Color',[1 1 1])
	YL = 1.2*max(A(:)+B(:));
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
            
            x_s_int = 0:.05:1;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [slope1,intercept1,stats_cyt,h_p1] = linlsq(s_new_cyt,r_new_cyt);
            [slope2,intercept2,stats_nuc,h_p2] = linlsq(s_new,r_new);
            
            plot(s_cyt,r_cyt,'b.')
            hold on
            plot(s_new_cyt,r_new_cyt,'r*')
            hold on
            c_cyt = polyfit(s_cyt,r_cyt,1);
            % Evaluate fit equation using polyval
            y_est_cyt = polyval(c_cyt,s_cyt);
            % Add trend line to plot
%             plot(s_cyt,y_est_cyt,'LineWidth',2)
            plot(x_s_int,(slope1.*x_s_int)+intercept1,'LineWidth',2)
            hold on
            plot(s,r,'m.')
            hold on
            plot(s_new,r_new,'g*')
            hold on
            c = polyfit(s,r,1);
            % Evaluate fit equation using polyval
            y_est = polyval(c,s);
            % Add trend line to plot
%             plot(s,y_est,'LineWidth',2)
            plot(x_s_int,(slope2.*x_s_int)+intercept2,'LineWidth',2)
            ylim([0 500])
%             [slope1,intercept1,stats_cyt,h_p1] = linlsq(s_cyt,r_cyt);
%             [slope2,intercept2,stats_nuc,h_p2] = linlsq(s,r);
%             xlabel(['slope_{cyt} = ',num2str(stats_cyt.m),...
% 		',   pval_{slope-cyt} = ',num2str(stats_cyt.pval_m),...
%         ',   slope_{nuc} = ',num2str(stats_nuc.m),...
% 		',   pval_{slope-nuc} = ',num2str(stats_nuc.pval_m)])
            legend('out-cyt','in-cyt','fitted line-cyt',...
                'out-nuc','in-nuc','fitted line-nuc',...
                'Location','best','FontSize',8)
            legend('boxoff')

			%
			% Add time annotation
			%
% 			timestr = sprintf('%6.2f min, %6.2f microns',t(i),z(j));
			timestr = sprintf('%6.2f min',t(i));
% 			htb = annotation('textbox',[0.5450 0.8105 0.3500 0.1000]);
			htb = annotation('textbox',[0.3050 0.8105 0.6000 0.1000]);
			set(htb,'String',timestr,'Fontsize',12,...
				'HorizontalAlignment','right','LineStyle','none')
			
			%
			% Other annotations
			%
			set(gca,'Fontsize',12)
% 			xlabel('s [fractional half circumference]')
            xlabel([',   slope_{nuc} = ',num2str(stats_nuc.m),...
                    ',   pval_{slope-nuc} = ',num2str(stats_nuc.pval_m)])
%             title(['slope_{cyt} = ',num2str(stats_cyt.m),...
%                     ',   pval_{slope-cyt} = ',num2str(stats_cyt.pval_m),...
%                     ',   slope_{nuc} = ',num2str(stats_nuc.m),...
%                     ',   pval_{slope-nuc} = ',num2str(stats_nuc.pval_m)],'FontSize',5)
            title(['slope_{cyt} = ',num2str(stats_cyt.m),...
                    ',   pval_{slope-cyt} = ',num2str(stats_cyt.pval_m)])
%             title(['frame = ',num2str(i)])
			ylabel('r [relative intensity]')
			
			hold off			
% 			xlim([-1,1])	
            xlim([0,1])
			
            all_y_vals = [r_cyt(:); r(:); r_new_cyt(:); r_new(:); ...
              (slope1.*x_s_int(:) + intercept1); ...
              (slope2.*x_s_int(:) + intercept2)];

            % Remove NaNs or Infs
            all_y_vals = all_y_vals(isfinite(all_y_vals));

            % If we have valid values, use them to set limits
            if ~isempty(all_y_vals)
                ymin = 0;
                ymax = max(all_y_vals) * 1.1;  % add 10% headroom
                if ymax <= ymin
                    ymax = ymin + 1;  % fallback to avoid ylim error
                end
                ylim([ymin ymax])
            else
                ylim([0 500]);  % fallback if something went wrong
            end
			
			% Cleaning up
			
			F(i,j) = getframe(gcf); %#ok<AGROW>
			delete(htb)
        end
        
        r_mean_cyt(i) = mean([intercept1,slope1+intercept1]);
        pval_nuc(i) = stats_nuc.pval_m;
        pval_cyt(i) = stats_cyt.pval_m;
        slope_nuc(i) = slope2;
        slope_cyt(i) = slope1;
        
    end
    
% pval_nuc = pval_nuc';
% pval_cyt = pval_cyt';
% slope_nuc = slope_nuc';
% slope_cyt = slope_cyt';

r_mean_cyt_all_time = mean(r_mean_cyt);
    
colors1 = cbrewer('seq', 'PuRd', 9); 
colors2 = cbrewer('seq', 'PuBuGn', 9);



% Identify outliers based on p-values
cyt_outliers = pval_cyt < 0.001; % Outliers for cyt where pval < 0.02
nuc_outliers = pval_nuc < 0.001; % Outliers for nuc where pval < 0.02

% --- Create Output Folder Structure ---
figBaseFolder = 'Figures';
figSubFolder = fullfile(figBaseFolder, data.filenameshort);

if ~exist(figSubFolder, 'dir')
    mkdir(figSubFolder);
end

figure
set(gcf, 'Position', [100 100 400 300], ...
    'Paperpositionmode', 'auto', 'Color', [1 1 1])
% Plot normal points for cyt
plot(t(~cyt_outliers), slope_cyt(~cyt_outliers)./r_mean_cyt_all_time, 'o', ...
    'Color', colors1(5, :), 'LineWidth', 2)
hold on
% Highlight outliers for cyt
plot(t(cyt_outliers), slope_cyt(cyt_outliers)./r_mean_cyt_all_time, 'rx', ...
    'LineWidth', 2)
% Plot normal points for nuc
plot(t(~nuc_outliers), slope_nuc(~nuc_outliers)./r_mean_cyt_all_time, '*', ...
    'Color', colors2(5, :), 'LineWidth', 2)
% Highlight outliers for nuc
plot(t(nuc_outliers), slope_nuc(nuc_outliers)./r_mean_cyt_all_time, 'b+', ...
    'LineWidth', 2)

if data.filenameshort == "Time course nc10 or 11 to gastrulation.czi"
        xline(4,'r','LineWidth', 1.4);
        xline(16,'r','LineWidth', 1.4);
        xline(34,'r','LineWidth', 1.4);
        xline(8,'g','LineWidth', 1.4);
        xline(20,'g','LineWidth', 1.4);
        xline(36,'g','LineWidth', 1.4);
        
    elseif data.filenameshort == "nc 12 to gastrulation.czi"
        xline(6,'r','LineWidth', 1.4);
        xline(30,'r','LineWidth', 1.4);
        xline(8,'g','LineWidth', 1.4);
        xline(32,'g','LineWidth', 1.4);
    
    elseif data.filenameshort == "Embryo 1 nc12-gastrulation.czi"
        xline(6,'r','LineWidth', 1.4);
        xline(24,'r','LineWidth', 1.4);
        xline(8,'g','LineWidth', 1.4);
        xline(26,'g','LineWidth', 1.4);
        
    else 
        xline(10,'r','LineWidth', 1.4);
        xline(28,'r','LineWidth', 1.4);
        xline(12,'g','LineWidth', 1.4);
        xline(30,'g','LineWidth', 1.4);
 end

xlabel('Time(in min)')
ylabel('Slope')
% legend('cyt (n.s.)', 'cyt (pval<.02)', 'nuc (n.s.)', 'nuc (pval<.02)', ...
%     'Location', 'best')
legend('cyt (n.s.)', 'cyt (pval<.001)', 'nuc (n.s.)', 'nuc (pval<.001)', ...
    'Location', 'best')

% --- Save the Figure ---
figName = fullfile(figSubFolder, [data.filenameshort(1:end-4), '_slope_plot']);
saveas(gcf, [figName, '.jpg']);
saveas(gcf, [figName, '.fig']);


  
else
	F = 'no plot asked for';	
end


%
% Concluding remarks
%
data.t = t;
data.z = z;
data.s_mid = s_mid;

try
	save([data.pth,data.prefix,'data.mat'],'data');
end

