function [t,ynuc_smooth,ycyt_smooth,ynuc,ycyt] = smooth_ynuc(soln,smoothwindow)

soln.y_nucall(soln.y_nucall == 0) = NaN;
soln.y_cytall(soln.y_cytall == 0) = NaN;


% removing the bleached nucleus column for nucall and cytall
% soln.y_nucall1(:,idx) = [];
% soln.y_cytall1(:,idx) = [];


% %
% % GTR: count how many NaN's are in each column (ie, for each nucleus or
% % cytoplasm). These would be NaN if they didn't show up for a given
% % frame. They were originally zeros but got converted to NaN's above.
% %
% total_nan_y_nucall = sum(isnan(soln.y_nucall));
% total_nan_y_cytall = sum(isnan(soln.y_cytall));
%
% % max_fr_missing = floor(.1*size(soln.y_nuc,1));
% max_fr_missing = 0; % GTR: "max frames missing"...not allowing for any
%
% %
% % The next several lines look like they're devoted to finding the
% % variable "locs", but in reality, there is a field stored in "soln"
% % where "locs" is defined. But in a nutshell, we're figuring out how
% % many nuclei to use in our average of all other nuclei. "locs" defines
% % the indices/ID#'s of those nuclei. Since "max_fr_missing" has been
% % set to zero, this means we only take nuclei that are present the
% % whole time.
% %
% % However, as I mentioned, there is already a field in "soln" that
% % defines "locs". So most of this code is obsolete, so I will comment
% % it out.
% %
% % I assume that, since there are six indices in that field, that means
% % we're looking at nearest neighbors.
% %
% loc_nuc = find(total_nan_y_nucall<= max_fr_missing);
% loc_cyt = find(total_nan_y_cytall<= max_fr_missing);
%
% if size(loc_cyt,2)<= size(loc_nuc,2)
%
%     locs = loc_cyt;
% else
%     locs = loc_nuc;
% end

locs = soln.locs; % six nearest neighbors?

% locs = intersect(locs,locs_n);

mean_oth_nuc = mean(fillmissing(soln.y_nucall(:,locs),'movmedian',10),2);
mean_oth_cyt = mean(fillmissing(soln.y_cytall(:,locs),'movmedian',10),2);

std_oth_nuc = std(fillmissing(soln.y_nucall(:,locs),'movmedian',10),0,2);
std_oth_cyt = std(fillmissing(soln.y_cytall(:,locs),'movmedian',10),0,2);

%ser_nuc = std_oth_nuc/sqrt(length(soln.y_nucall));
%ser_cyt = std_oth_cyt/sqrt(length(soln.y_nucall));

ts = norminv([0.025  0.975]);

ci_nuc = mean_oth_nuc+ts.*std_oth_nuc;
ci_cyt = mean_oth_cyt+ts.*std_oth_cyt;


% diff_nuc = mean_oth_nuc - fillmissing(soln.y_nuc,'linear',1);
% diff_cyt = mean_oth_cyt - fillmissing(soln.y_cyt,'linear',1);


% ---------------------------------------------------------------------
% Normalize data for bleaching -- Sharva's way (nearest neighbors?)
% ---------------------------------------------------------------------

% estimating the recovaery rates for cact data after normalizing

% for new fitting
t = soln.t(2:end);
ynuc = (soln.y_nuc(2:end)./mean_oth_nuc(2:end)).*mean_oth_nuc(2);
ycyt = (soln.y_cyt(2:end)./mean_oth_cyt(2:end)).*mean_oth_cyt(2);

%
% GTR added: Removing NaNs from our bleached nucleus
%
v = isnan(ynuc) | isnan(ycyt);
t(v) = [];
ynuc(v) = [];
ycyt(v) = [];


% ---------------------------------------------------------------------
% Normalize data for bleaching -- my way (all nuclei)
% ---------------------------------------------------------------------
%{
ynuc1 = soln.y_nuc(2:end);
ycyt1 = soln.y_cyt(2:end);

ynucall = soln.y_nucall(2:end);
ycytall = soln.y_cytall(2:end);
v = false(1,size(ynucall,2));
for i = 1:size(ynucall,2)
	v(i) = isequalwithequalnans(ynucall(:,i),ynuc1);
end
ynucall(:,v) = [];
ycytall(:,v) = [];

ynucall(ynucall == 0) = NaN;
ynucavg = meanDU(ynucall,2);
ycytall(ycytall == 0) = NaN;
ycytavg = meanDU(ycytall,2);

% %
% % Removing NaNs from our bleached nucleus
% %
% v = isnan(ynuc) | isnan(ycyt);
% t(v) = [];
% ynuc(v) = [];
% ycyt(v) = [];

%
% Normalizing for bleaching
%
% ynuc = ynuc./ynucavg.*ynucavg(1);
ynuc1 = ynuc1./ycytavg.*ycytavg(1);
ycyt1 = ycyt1./ycytavg.*ycytavg(1);

%}


% %
% % Smooth the nuc and cyt data
% %
% ycyt_fit = smooth(t,ycyt,smoothwindow);
% % ycyt_fit(1) = ycyt(1);
% % ycyt_fit = smoothdata(ycyt,'sgolay',10); % for 2 gfp and free gfp cases
% ynuc_smooth = smooth(t,ynuc,smoothwindow);
% % ynuc_smooth(1) = ynuc(1);

%
% Experimental: smooth the data ignoring the first "smoothwindow" points
%
ycyt_smooth = ycyt;
ycyt_smooth = smooth(t,ycyt,smoothwindow);
ycyt_smooth(1:smoothwindow) = ycyt(1:smoothwindow);
ynuc_smooth = smooth(t,ynuc,smoothwindow);
ynuc_smooth(1:smoothwindow) = ynuc(1:smoothwindow);
















