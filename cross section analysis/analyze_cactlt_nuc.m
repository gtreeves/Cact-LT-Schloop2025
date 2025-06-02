function data = analyze_cactlt_nuc(filename,filenum,cyt_data,data_ch,varargin)
%Stabilizes image & extracts nuc & cyt time course; finds bleached nucleus
%
%function data = analyze_bleach(filename,yesplot)
%
% "filename": filename (incl path) to the czi file that contains the
%	bleaching time course.
% "genotype": string containing the genotype of the embryo in question
% "data_ch": the numerical tag for the channel where the bleaching data are
%
% Optional argument varargin can consist of these things, in this order:
%	* "mask_ch": numeric input that tells the function which channel has
%		the image data from which to derive the mask.  This can be used to
%		calculate the autocorrelation function for any arbitrarily-shaped
%		region, and perhaps even disjoint regions.  If passed as a scalar
%		numerical zero, or logical false, then the mask will not be used
%		and the autocorrelation function for the entire frame will be
%		calculated. Default, 2.
%		If this is not specified, but you still want to
%		specify other arguments, put empty brackets -- [] -- in place of
%		this argument.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else 
	mask_ch = 2;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end%, iArg = iArg + 1;


vsep = strfind(filename,filesep);
vsep1 = vsep(end);
% vsep2 = vsep(end-1);
filenameshort = filename(vsep1+1:end);
pthname = filename(1:vsep1);


% lsm = lsmRead2(filename);
[IM,lsminf1,lsminf2] = openczi(filename);
% lsminf1 = lsm.lsm;
% lsminf2 = lsminfo(filename);
calib_mean = mean(mean(double(IM(:,:,1))));

scalings = lsminf2.VOXELSIZES*1e6; % microns per pixel.

% t = lsminf2.TimeStamps.TimeStamps;
% t = t - t(2);

r=bfGetReader(filename);
nSeries = r.getSeriesCount();
m = r.getMetadataStore();
nImages = m.getImageCount();
nPlanes = m.getPlaneCount(0);

t1 = zeros(nPlanes/2,1);
t2 = zeros(nPlanes/2,1);
for i=0:nPlanes-1
    
    if rem(i,2) == 1
        
        t1((i+1)/2,1) = m.getPlaneDeltaT(0,(i)).value;
         
    else
        
        t2((i+2)/2,1) = m.getPlaneDeltaT(0,(i)).value;
                    
    end
    
end

t = t1;
dt = t(2) - t(1); % Assuming the time interval is constant;

rho = round(scalings(2)/scalings(1)); % hope we don't have to round much.
LP = [lsminf2.ScanInfo.POWER{1} lsminf2.ScanInfo.POWER{2}];

% IM = lsm.data;
[H,W,num_ch,num_frames] = size(IM);
numch = num_ch;
D = num_frames;



if num_frames < 5
	% skip this file if not enough frames
	data = 'Too few frames';
	return
end

%default
% data_ch = 2;
% mask_ch = 1;

% for Cact-LT czi data
data_ch = 1;
mask_ch = 2;

if num_ch >= 2
	I_dl = squeeze(IM(:,:,data_ch,:));
	I_nuc = squeeze(IM(:,:,mask_ch,:)); % red channel is best
else
	I_dl = squeeze(IM(:,:,data_ch,:));
	I_nuc = I_dl; % but some time series don't have it
end


if H < 50 || W < 50
	% skip this file if, after stabilizing, there is not enough image left
	data = 'Nuclei moved too much';
	return
end

%
% Call function to segment the nuclei.  The output of this function
% includes L, which is essentially a label array for every frame.
%

% figure,imshow(uint8(I_dl(:,:,4)))
% figure,imshow(uint8(I_nuc(:,:,4)))

ring_width = 18.36; % microns
nuc_ch = 2; % nuclear channel. EDIT THIS FOR EVERY DIFFERENT FILE RUN
np_ch = 1; % nuclear protein channel
nt = 60;
ns = 300;
z0 = 100; % how deep the slice is
h = 0.25;
size_cc = [];
Soln = [];


Yhatmax = round(ring_width/scalings(1)); 
% This is the width of the ring
% (radius_outer - radius_inner) in pixels.  Default width: 18.36 microns.


% First we find the periphery of the embryo and detect the intensity of
% each channel in quadrilaterals as you go around the periphery.  Next, we
% perform calculations on the nuclei/nuclear proteins/nuclear dots, if
% these things are present.  To do these things, we go through a loop that
% takes each slice individually.
%		
arc = zeros(D,1);
Xp = cell(D,1); Yp = cell(D,1);
X1 = cell(D,1); Y1 = X1; S = X1;
Nuc = X1; Std_nuc = X1;
Nuc2 = X1; Std_nuc2 = X1;
Nuc_protein = X1; Std_nuc_protein = X1;
Nuc_protein2 = X1; Std_nuc_protein2 = X1;
Intron = X1; Std_intron = X1;
U1 = X1;
T_out = zeros(ns+1,numch,D);
% 		Raw = zeros(ns,numch,D);
raw_normimg = zeros(nt,D);


background = ftn_calcbg(I_dl);
bcg = mean(background);

for i = 1:D%size(I_dl,3)
    
    I1 = I_dl(:,:,i);
    [xp,yp] = borderFinder(I1,h,yesplot,nt,[],[],filename,i);
    
    Xp{i} = xp; 
    Yp{i} = yp;
    arc(i) = sum(sqrt(diff(xp).^2 + diff(yp).^2));
    
     
    stage = '';
    yesplot = true;
    
    % for nuclei
    if filenameshort == "Time course nc10 or 11 to gastrulation.czi"
        if i<20
            lims = [50 100];
        else
            lims = [10 100];
        end
    elseif filenameshort == "nc 12 to gastrulation.czi"
        if i<8
            lims = [50 100];
        else
            lims = [10 100];
        end        
    elseif filenameshort == "Embryo 1 nc12-gastrulation.czi"
        if i<10
            lims = [60 100];
        elseif (10<=i)&&(i<=15)
            lims = [40 100];
        else
            lims = [10 100];
        end
    else 
        if i<12
            lims = [60 100];
        elseif (12<=i)&&(i<=16)
            lims = [50 100];
        else
            lims = [10 100];
        end
        
    end
    
    [nucstats,xnuc,ynuc,snuc,w,mask,u1,CC_size] = ...
        find_nuclei_nuc(I_nuc(:,:,i),I_dl(:,:,i),lims,...
        xp,yp,scalings,Yhatmax,yesplot,stage,filename,i,[],4);
    
%     % for cytoplasm
%     [nucstats,xnuc,ynuc,snuc,w,mask,u1,CC_size] = ...
%         find_nuclei3_c2(I_nuc(:,:,i),...
%         xp,yp,scalings,Yhatmax,yesplot,stage,filename,i,[],4);
    
    U1{i} = mask;
    X1{i} = xnuc;
    Y1{i} = ynuc;
    S{i} = 2*snuc/w - 1;
    %size = [];
    size_cc = horzcat(size_cc,CC_size);
    num_nuc{i} = CC_size;
    
    
    % Create normalization image
    %
    [xc,yc,R] = circfit(xp,yp);
    [normimg,t,raw] = normalization_image(I1,xp,yp,xc,yc,R,Yhatmax);
    raw_normimg(:,i) = raw;

    
    % Trace perimeter of each nucleus
    %
    [B,~,~,Idx] = traceobject(nucstats,H,W);
    for j = 1:length(Idx)
        nucstats(j).PerimIdxList = Idx{j};
        nucstats(j).PerimXYList = fliplr(B{j});
    end

    
    % Invisibly plot, then export each nucleus
    %
    figure('visib','off')
    set(gcf,'paperpositionmode','auto')
%     imshowcat(gcf,mDU(double(I1)),nucstats)
%     imshowcat(gcf,mDU(double(imadjust(I1))))
    imshowcat(gcf,mDU(double(imadjust(I1))),nucstats)
%     imshowcat(gcf,mDU(double(I1)),double(~~mask))
    set(gca,'Position',[0 0 1 1])

    %noi = noise_reduction(I(:,:,np_ch),xp,yp,Yhatmax,filenameshort1,i);
    
    % Obtaining the intensity of the nuclear protein channels, as well
    % as the nuclear channel itself, for each nucleus.
    %
    [np,stdnp] = nuclearintensity(I1,nucstats);
    [np2,stdnp2] = nuclearintensity(double(I1)./...
        repmat(normimg,1,1,numch),nucstats);
    Nuc{i} = np(:,1)-bcg;
    Std_nuc{i} = stdnp(:,1);
    Nuc2{i} = np2(:,1)-bcg;
    Std_nuc2{i} = stdnp2(:,1);
    
    if ~isempty(np_ch)
        Nuc_protein{i} = np(:,2:end);
        Std_nuc_protein{i} = stdnp(:,2:end);
        Nuc_protein2{i} = np2(:,2:end)-bcg;
        Std_nuc_protein2{i} = stdnp2(:,2:end);
    elseif i == D
        Nuc_protein = NaN; Std_nuc_protein = NaN;
        Nuc_protein2 = NaN; Std_nuc_protein2 = NaN;
    end

	disp(['borderfinder and nucfinder i = ',num2str(i),' out of ',num2str(D)])
         

end



% Transforming data on the nuclear protein into the scaled data. To
% do this, we have to perform several normalizations. Some are
% within the for loop, some without.
%		
nuc = cell2mat(Nuc);
Timecourse_norm = median(nuc)/LP(1)/100; % Dividing by this
% scalar corrects for an overall brighter timecourse in the
% nuclei. This could be an issue if, say, one embryo was taken
% at z = 160 and another at z = 140, by accident. The overall
% nuclear intensity should be the same across every embryo
% imaged.
Calib_norm = calib_mean/1.75e4;
Calib_norm = 1;
% Dividing by this scalar
% corrects for day-to-day laser power changes. The calib_mean
% is the mean intensity of the calibration image taken on a
% given day.

%
% Normalization for the case in which we pre-normalized by our
% normalization image.
%
Nuc_norm2 = 1;

std_nuc = cell2mat(Std_nuc);
R = cell(D,1); Std_R = cell(D,1);
R2 = cell(D,1);

for i = 1:D

%
% Below is the typical, nucleus-by-nucleus calculation.
%
    Nuc_norm = Nuc{i}./median(Nuc{i}); % Dividing by this corrects
% for changes in intensity around the circumference of the
% embryo, at frame i. It is the nuclear intensity, which is
% supposed to be uniform, normalized by the median nuclear
% intensity (so that dividing by this vector does not bring our
% nuc protein intensity down to order 1.

    R{i} = Nuc_protein{i}/LP(2) ./ Nuc_norm / Timecourse_norm / Calib_norm;

    Std_R{i} = R{i}.*sqrt((Std_nuc_protein{i}./Nuc_protein{i}).^2 + ...
    repmat((Std_nuc{i}./Nuc{i}).^2,1,size(Nuc_protein{i},2)));

%
% We also created a version of our nuclear protein measurement
% that came from our original image pre-normalized by our
% normalization image. This nucprotein measurement doesn't have
% to be later normalized by the individual nuclei to correct
% for variations in intensity along the DV axis. So this
% calculation will take that into account. We assume the
% variability here comes almost 100% from the Nuc_protein part,
% so we don't have to calculate a "Std_R2".
%
%     R2{i} = Nuc_protein2{i}/LP(2) ./ Nuc_norm2 / Timecourse_norm / Calib_norm;
    R2{i} = Nuc_protein2{i};
end

%
% Normalization for the domainMeas case
%
T_out1 = squeeze(T_out(:,nuc_ch,:));
T_out2 = squeeze(T_out(:,np_ch,:));
Timecourse_norm3 = mean(T_out1(:))/LP(1)/100;
Nuc_norm3 = T_out1./repmat(mean(T_out1),ns+1,1);
R3 = T_out2/LP(2) ./ Nuc_norm3 / Timecourse_norm3 / Calib_norm;


% metadata

data.filenameshort = filenameshort;
data.filename = filename;
data.metadata.lsminf1 = lsminf1;
data.metadata.lsminf2 = lsminf2;
data.LP = lsminf2.ScanInfo.POWER;
data.H = H;
data.W = W;
data.num_frames = num_frames;
data.num_ch = num_ch;

data.D = D;
data.dt = dt;
data.z0 = z0;

data.Lbar = mean(arc)/2*scalings(1);
data.s_mid = NaN;
% data.s = s;
data.t_out = T_out;
data.R3 = R3;
data.S = S;
data.R = R;
data.Std_R = Std_R;
data.R2 = R2;
data.Std_R2 = Std_nuc_protein2;
data.metadata.LP = LP;
data.metadata.calib_mean = calib_mean;
data.metadata.calib_LP = lsminf2.ScanInfo.POWER;
data.metadata.w = mean(arc);
data.metadata.arc = arc;
data.metadata.L = arc/2*scalings(1);
data.metadata.Xp = Xp;
data.metadata.Yp = Yp;
data.metadata.X = X1;
data.metadata.Y = Y1;
data.metadata.S = S;
data.metadata.Nuc = Nuc;
data.metadata.Std_nuc = Std_nuc;
data.metadata.Nuc_protein = Nuc_protein;
data.metadata.Std_nuc_protein = Std_nuc_protein;
data.metadata.Nuc_protein2 = Nuc_protein2;
data.metadata.Std_nuc_protein2 = Std_nuc_protein2;
data.bcg = bcg;

data.nucprotein_names = {};
data.A = NaN;
data.B = NaN;
data.M = NaN;
data.Sig = NaN;
data.dA = NaN;
data.dB = NaN;
data.dM = NaN;
data.dSig = NaN;
data.gof = NaN;

save([data.filename(1:end-4),'_data.mat'],'data');

% Now we try to plot it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These will generate video and figures for slope and p-vals across the dorsa-ventral
% axis

yesplot = true;
[data,F] = ftn_gradfit8_plot_slope(data,cyt_data,yesplot);

% --- Extract file path and short name (remove .czi) ---
filename_full = char(data.filename);                % Ensure it's a character vector
[fpath, ~, ~] = fileparts(filename_full);           % Extract the folder path
filenameshort_noext = data.filenameshort(1:end-4);  % Remove extension (e.g., '.czi')

% --- Save .mat files to the same directory as the input file ---
datafile = fullfile(fpath, [filenameshort_noext, '_data.mat']);
save(datafile, 'data');

Soln = [Soln; data];

Ffile = fullfile(fpath, [filenameshort_noext, '_F.mat']);
save(Ffile, 'F');

% --- Generate Output Video Frames ---
F1 = ftn_output2(data, F, IM);

% --- Save video under Figures/filenameshort/ ---
figBaseFolder = 'Figures';
figSubFolder = fullfile(figBaseFolder, data.filenameshort);
if ~exist(figSubFolder, 'dir')
    mkdir(figSubFolder);
end

vidfile = fullfile(figSubFolder, [filenameshort_noext, '_slope_vid.avi']);
v_out = VideoWriter(vidfile);
v_out.FrameRate = 2;
open(v_out)
writeVideo(v_out, F1)
close(v_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These will generate video slope intercept values at locations s = o and 1

yesplot = true;
[~,F_intercepts] = ftn_gradfit6(data,cyt_data,yesplot);

% --- Generate Output Video Frames ---
F_intercepts1 = ftn_output2(data, F_intercepts, IM);

% --- Save video under Figures/filenameshort/ ---
figBaseFolder = 'Figures';
figSubFolder = fullfile(figBaseFolder, data.filenameshort);
if ~exist(figSubFolder, 'dir')
    mkdir(figSubFolder);
end

vidfile = fullfile(figSubFolder, [filenameshort_noext, '_slope_intercepts_vid.avi']);
v_out = VideoWriter(vidfile);
v_out.FrameRate = 2;
open(v_out)
writeVideo(v_out, F_intercepts1)
close(v_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These will generate figures of the mean nuclear and cytoplasmic intensities 
%with standard deviation across frames. Also save data rewuired for
%generating the average embryo later on

yesplot = true;
[~,~] = ftn_gradfit7(data,filenum,cyt_data,yesplot);

end





