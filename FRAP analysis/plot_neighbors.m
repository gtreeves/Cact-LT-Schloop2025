function [mean_oth_nuc,mean_oth_cyt] = plot_neighbors(nucdl,cytdl,idx,locs,t,datum,filenameshort,filenum)

% saves the image results for the 'plot of nucall and cytall intensities',
% 'plot of mean of other nucleus,cytoplasm and bleached nucleus ,cytoplasm'
% and 'diff. between mean of other nuclei and bleached nucleus' and 'diff.
% between mean of other cytoplasm and bleached cytoplasm'

% function [mean_oth_nuc,mean_oth_cyt] = plot_neighbors(nucdl,cytdl,idx,locs,t,filenameshort)

% nucdl contains the intensities of all the nuclei for all the frames. It is
% a (frame x num of nuc) matrix

% cytdl contains the intensities of all the cytoplasm for all the frames. It is
% a (frame x num of cyt) matrix

% idx is the identity/label of the bleached nucleus

% locs has the label numbers of the closest specified
% number of neighbors

% t is vector conatining the time when each frame was captured

nucdl(nucdl == 0) = NaN;
cytdl(cytdl == 0) = NaN;

nucdl= fillmissing(nucdl,'nearest',1);
cytdl=fillmissing(cytdl,'nearest',1);


soln_temp.y_nucall = nucdl;
soln_temp.y_cytall = cytdl;
soln_temp.y_nuc = nucdl(:,idx);
soln_temp.y_cyt = cytdl(:,idx);
soln_temp.t = t;


% soln.y_nuc(soln.y_nuc == 0) = NaN
% soln.y_cyt(soln.y_cyt == 0) = NaN

%
soln_temp.y_nucall(soln_temp.y_nucall == 0) = NaN;
soln_temp.y_cytall(soln_temp.y_cytall == 0) = NaN;


% removing the bleached nucleus column for nucall and cytall
% soln.y_nucall1(:,idx) = [];
% soln.y_cytall1(:,idx) = [];


figure('visib','off')
% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
p1 = plot(soln_temp.t,fillmissing(soln_temp.y_nucall(:,locs),'nearest',1),'r', 'LineWidth', 1);
hold on
p2 = plot(soln_temp.t,fillmissing(soln_temp.y_cytall(:,locs),'nearest',1),'b', 'LineWidth', 1);
hold on
p3 = plot(soln_temp.t,fillmissing(soln_temp.y_nuc,'nearest',1),'g', 'LineWidth', 2);
hold on
p4 = plot(soln_temp.t,fillmissing(soln_temp.y_cyt,'nearest',1),'m', 'LineWidth', 2);
hold on
title('plot of nucall and cytall intensities')
legend([p1(1) p2(1) p3 p4],{'nucall int','cytall int','blc nuc','blc cyt'})
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'nucall_cytall_',filenameshort(1:end-4),'.fig']);
if ~(exist('Image results\Neighbor nuc stats','dir') == 7)
        mkdir 'Image results\Neighbor nuc stats'
    end

    baseDir = fullfile('Image results', 'Neighbor nuc stats', datum);
    subDir = fullfile(baseDir, filenameshort);

    % Create base and subdirectory if they don't exist
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end     
saveas(gcf, fullfile('Image results', 'Neighbor nuc stats', datum, filenameshort, ...
['nucall_cytall_', filenameshort(1:end-4),'.jpg']));
    

mean_oth_nuc = mean(fillmissing(soln_temp.y_nucall(:,locs),'nearest',1),2);
mean_oth_cyt = mean(fillmissing(soln_temp.y_cytall(:,locs),'nearest',1),2);

std_oth_nuc = std(fillmissing(soln_temp.y_nucall(:,locs),'nearest',1),0,2);
std_oth_cyt = std(fillmissing(soln_temp.y_cytall(:,locs),'nearest',1),0,2);

%ser_nuc = std_oth_nuc/sqrt(length(soln.y_nucall));
%ser_cyt = std_oth_cyt/sqrt(length(soln.y_nucall));

ts = norminv([0.025  0.975]); 

ci_nuc = mean_oth_nuc+ts.*std_oth_nuc;
ci_cyt = mean_oth_cyt+ts.*std_oth_cyt;

figure('visib','off')
% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
plot(soln_temp.t,mean_oth_nuc,'r', 'LineWidth', 2)
hold on
plot(soln_temp.t,fillmissing(soln_temp.y_nuc,'nearest',1),'m', 'LineWidth', 2)
hold on
patch([soln_temp.t' fliplr(soln_temp.t')], [ci_nuc(:,1)' fliplr(ci_nuc(:,2)')], 'r', 'FaceAlpha', .2)
hold on
plot(soln_temp.t,fillmissing(soln_temp.y_cyt,'nearest',1),'b', 'LineWidth', 2)
hold on
plot(soln_temp.t,mean_oth_cyt,'g', 'LineWidth', 2)
hold on
patch([soln_temp.t' fliplr(soln_temp.t')], [ci_cyt(:,1)' fliplr(ci_cyt(:,2)')], 'g', 'FaceAlpha', .2)
hold on
title('plot of mean of other nuc,cyt and bleached nuc,cyt')
legend('mean nuc','bleach nuc','95% conf nuc','bleach cyt','mean cyt','95% conf cyt')
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'mean_other_nuc_cyt_',filenameshort(1:end-4),'.fig']);
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'mean_other_nuc_cyt_',filenameshort(1:end-4),'.jpg']);
saveas(gcf, fullfile('Image results', 'Neighbor nuc stats', datum, filenameshort, ...
['mean_other_nuc_cyt_', filenameshort(1:end-4), '_','.jpg'])); 

diff_nuc = mean_oth_nuc - fillmissing(soln_temp.y_nuc,'linear',1);
diff_cyt = mean_oth_cyt - fillmissing(soln_temp.y_cyt,'linear',1);

figure('visib','off')
% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
subplot(2,1,1)
plot(soln_temp.t',diff_nuc','r', 'LineWidth', 2)
hold on
plot(soln_temp.t',(0.*soln_temp.t'),'g', 'LineWidth', 2)
patch([soln_temp.t' fliplr(soln_temp.t')], [zeros(size(diff_nuc')) fliplr(diff_nuc')], 'c', 'FaceAlpha', .2)
title('diff. between mean of other nuc and bleached nuc')

subplot(2,1,2)
plot(soln_temp.t',diff_cyt','r', 'LineWidth', 2)
hold on
plot(soln_temp.t',(0.*soln_temp.t'),'g', 'LineWidth', 2)
patch([soln_temp.t' fliplr(soln_temp.t')], [zeros(size(diff_cyt')) fliplr(diff_cyt')], 'c', 'FaceAlpha', .2)
title('diff. between mean of other cyt and bleached cyt')
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'diff_between_int_',filenameshort(1:end-4),'.fig']);
% saveas(gcf, ['Image results\Neighbor nuc stats',filesep,'diff_between_int_',filenameshort(1:end-4),'.jpg']);
saveas(gcf, fullfile('Image results', 'Neighbor nuc stats', datum, filenameshort, ...
['diff_between_int_', filenameshort(1:end-4), '_','.jpg'])); 
 


