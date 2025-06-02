% script for obtaining the Averaged plot of intensity for 4 embryo cross sections

clc 
clear all

% Assuming data1, data2, ..., dataN contain 'x' for time and 'y' for intensity

% Load your data here. For example:
load('data1.mat') % embryo 'Time course nc10 or 11 to gastrulation'
load('data2.mat') % embryo 'nc 12 to gastrulation'
load('data3.mat') % embryo 'Embryo 1 nc12-gastrulation'
load('data4.mat') % embryo 'Embryo 2 nc11-gastrulation'

N = 4;

data1.t = data1.t(6:end);
data1.mean_nuc = data1.mean_nuc(6:end);
data1.mean_cyt = data1.mean_cyt(6:end);
data1.std_nuc = data1.std_nuc(6:end);
data1.std_cyt = data1.std_cyt(6:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtaining the amplitude sclaing coeff for curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find the indices of the last highest peaks for each dataset for nuclei
[~, peakIdx1] = max(data1.mean_nuc(1:end)); % Assuming the peak will be in the last 37 points
[~, peakIdx2] = max(data2.mean_nuc(1:end));
[~, peakIdx3] = max(data3.mean_nuc(1:end)); % Assuming the peak will be in the last 37 points
[~, peakIdx4] = max(data4.mean_nuc(1:end));

% using the peak intensity of data 2 as reference to scale up the other
% data
referencePeakIntensity = data2.mean_nuc(peakIdx2);


% scaling factor for others 
scale1 = referencePeakIntensity/data1.mean_nuc(peakIdx1);
scale2 = referencePeakIntensity/data2.mean_nuc(peakIdx2);
scale3 = referencePeakIntensity/data3.mean_nuc(peakIdx3);
scale4 = referencePeakIntensity/data4.mean_nuc(peakIdx4);


% Find the indices of the lowest peaks for each dataset for cytoplasm
[~, peakcIdx1] = min(data1.mean_cyt(1:end)); % Assuming the peak will be in the last 37 points
[~, peakcIdx2] = min(data2.mean_cyt(1:end));
[~, peakcIdx3] = min(data3.mean_cyt(1:end)); % Assuming the peak will be in the last 37 points
[~, peakcIdx4] = min(data4.mean_cyt(1:end));

% using the peak intensity of data 2 as reference to scale up the other
% data for cytoplasm
referencePeakcIntensity = data2.mean_cyt(peakcIdx2);


% scaling factor for other cytoplasm
scale1c = referencePeakcIntensity/data1.mean_cyt(peakcIdx1);
scale2c = referencePeakcIntensity/data2.mean_cyt(peakcIdx2);
scale3c = referencePeakcIntensity/data3.mean_cyt(peakcIdx3);
scale4c = referencePeakcIntensity/data4.mean_cyt(peakcIdx4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the scaled and offset curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Find the two highest peaks from the end towards the start for the reference dataset
% (assuming data2 is the reference)
[peaksRef, locsRef] = findpeaks(flip(data1.mean_nuc));
[~, highestPeaksIdx] = maxk(peaksRef, 2); % Find the two highest peaks in reverse order
locsRef = length(data1.mean_nuc) - locsRef(highestPeaksIdx) + 1; % Convert indices back to original order
highPeakTimesRef = data1.t(locsRef);

% Sort the times of the two highest peaks for consistency
highPeakTimesRef = sort(highPeakTimesRef, 'ascend');

figure('Visible', 'off'),
% Loop through each dataset
for i = [2:N]  % Assuming you have N datasets, and data2 is the reference
    % Find peaks from the end towards the start
    
    
    [peaks, locs] = findpeaks(flip(eval(sprintf('data%d.mean_nuc', i))));
    [~, highestPeaksIdx] = maxk(peaks, 2); % Find the two highest peaks in reverse order
    locs = length(eval(sprintf('data%d.mean_nuc', i))) - locs(highestPeaksIdx) + 1; % Convert indices back
    highPeakTimes = eval(sprintf('data%d.t(locs)', i));
    
    % Sort the times of the two highest peaks for consistency
    highPeakTimes = sort(highPeakTimes, 'ascend');
    
    % Calculate the scale and offset to align the peaks
    scale = (highPeakTimesRef(2) - highPeakTimesRef(1)) / (highPeakTimes(2) - highPeakTimes(1));
    offset = highPeakTimesRef(1) - (highPeakTimes(1) * scale);
    
    % Apply the transformation to the time data
    eval(sprintf('data%d.t = data%d.t * scale + offset;', i, i));
    
    % Plot the transformed data to verify alignment
    plot(eval(sprintf('data%d.t', i)), eval(sprintf('scale%d.*data%d.mean_nuc',i, i))); hold on;
end

% Plot the reference data (data2 in this case)
plot(data1.t, scale1.*data1.mean_nuc, 'LineWidth', 2);

% Add labels and legend
xlabel('Time (min)');
ylabel('Mean Intensity');
legendArray = arrayfun(@(x) sprintf('Data Set %d', x), [2:N], 'UniformOutput', false);
legend([legendArray, 'Reference Data 1'], 'Location', 'best');

% Release the plot hold
hold off;


% for cytoplasm
figure('Visible', 'off'),

for i = [2:N]  % Assuming you have N datasets, and data2 is the reference
    % Find peaks from the end towards the start
    
    
    [peaks, locs] = findpeaks(flip(eval(sprintf('data%d.mean_nuc', i))));
    [~, highestPeaksIdx] = maxk(peaks, 2); % Find the two highest peaks in reverse order
    locs = length(eval(sprintf('data%d.mean_nuc', i))) - locs(highestPeaksIdx) + 1; % Convert indices back
    highPeakTimes = eval(sprintf('data%d.t(locs)', i));
    
    % Sort the times of the two highest peaks for consistency
    highPeakTimes = sort(highPeakTimes, 'ascend');
    
    % Calculate the scale and offset to align the peaks
    scale = (highPeakTimesRef(2) - highPeakTimesRef(1)) / (highPeakTimes(2) - highPeakTimes(1));
    offset = highPeakTimesRef(1) - (highPeakTimes(1) * scale);
    
    % Apply the transformation to the time data
    eval(sprintf('data%d.t = data%d.t * scale + offset;', i, i));
    
    % Plot the transformed data to verify alignment
    plot(eval(sprintf('data%d.t', i)), eval(sprintf('scale%d.*data%d.mean_cyt',i,i))); 
    hold on;
end

% Plot the reference data (data2 in this case)
plot(data1.t,  scale1c.*data1.mean_cyt, 'LineWidth', 2);

% Add labels and legend
xlabel('Time (min)');
ylabel('Mean Intensity');
legendArray = arrayfun(@(x) sprintf('Data Set %d', x), [2:N], 'UniformOutput', false);
legend([legendArray, 'Reference Data 1'], 'Location', 'best');
% Release the plot hold
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the scaled and offset curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nuc

figure('Visible', 'off'),
plot(data1.t, scale1.*data1.mean_nuc,'r', 'LineWidth', 2)
hold on
plot(data2.t, scale2.*data2.mean_nuc,'b', 'LineWidth', 1)
hold on
plot(data3.t, scale3.*data3.mean_nuc,'g', 'LineWidth', 1)
hold on
plot(data4.t, scale4.*data4.mean_nuc,'m', 'LineWidth', 1)
hold off
xlabel('Time (min)');
ylabel('Mean Intensity');
legend('Data 1 (Ref.)', 'Data 2', 'Data 3', 'Data 4')
title('scaled and offset curves of nuc data')


% for cyt

figure('Visible', 'off'),
plot(data1.t, scale1.*data1.mean_cyt,'r', 'LineWidth', 2)
hold on
plot(data2.t, scale2.*data2.mean_cyt,'b', 'LineWidth', 1)
hold on
plot(data3.t, scale3.*data3.mean_cyt,'g', 'LineWidth', 1)
hold on
plot(data4.t, scale4.*data4.mean_cyt,'m', 'LineWidth', 1)
hold off
xlabel('Time (min)');
ylabel('Mean Intensity');
legend('Data 1 (Ref.)', 'Data 2', 'Data 3', 'Data 4')
title('scaled and offset curves of cyt data')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolating the dataset for the same time points and plotting them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nuc
interpData1 = scale1.*data1.mean_nuc;
interpData2 = interp1(data2.t, scale2.*data2.mean_nuc, data1.t, 'nearest', 'extrap');
interpData3 = interp1(data3.t, scale3.*data3.mean_nuc, data1.t, 'linear', 'extrap');
interpData4 = interp1(data4.t, scale4.*data4.mean_nuc, data1.t, 'linear', 'extrap');

figure('Visible', 'off'),
plot(data1.t, interpData1,'r', 'LineWidth', 2)
hold on
plot(data1.t, interpData2,'b', 'LineWidth', 1)
hold on
plot(data1.t, interpData3,'g', 'LineWidth', 1)
hold on
plot(data1.t, interpData4,'m', 'LineWidth', 1)
hold off
legend('Data 1', 'Data2', 'Data3', 'Data4')
title('interpolated mean curves of nuc data')

% for cyt

interpcData1 = scale1.*data1.mean_cyt;
interpcData2 = interp1(data2.t, scale2.*data2.mean_cyt, data1.t, 'nearest', 'extrap');
interpcData3 = interp1(data3.t, scale3.*data3.mean_cyt, data1.t, 'nearest', 'extrap');
interpcData4 = interp1(data4.t, scale4.*data4.mean_cyt, data1.t, 'nearest', 'extrap');

figure('Visible', 'off'),
plot(data1.t, interpcData1,'r', 'LineWidth', 2)
hold on
plot(data1.t, interpcData2,'b')
hold on
plot(data1.t, interpcData3,'g')
hold on
plot(data1.t, interpcData4,'m')
hold off
legend('Data 1', 'Data2', 'Data3', 'Data4')
title('interpolated mean curves of cyt data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating the mean of the nuc and cyt curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nuc
nuc_data = [interpData1',interpData2,interpData3,interpData4];
mean_nuc = mean(nuc_data,2);

%for cyt
cyt_data = [interpcData1',interpcData2,interpcData3,interpcData4];
mean_cyt = mean(cyt_data,2);

figure('Visible', 'off'),

plot(data1.t,mean_nuc,'r');
hold on
plot(data1.t,mean_cyt,'b');
legend('mean nuc','mean cyt')
title('average curves for nuc and cyt of 4 embryo cross sections')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating the std of the mean nuc and mean cyt curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for nuc
interpStd1 = scale1.*data1.std_nuc;
interpStd2 = interp1(data2.t, data2.std_nuc, data1.t, 'nearest', 'extrap');
interpStd3 = interp1(data3.t, scale3.*data3.std_nuc, data1.t, 'linear', 'extrap');
interpStd4 = interp1(data4.t, scale4.*data4.std_nuc, data1.t, 'linear', 'extrap');

figure('Visible', 'off'),
plot(data1.t, interpStd1,'r', 'LineWidth', 2)
hold on
plot(data1.t, interpStd2,'b', 'LineWidth', 1)
hold on
plot(data1.t, interpStd3,'g', 'LineWidth', 1)
hold on
plot(data1.t, interpStd4,'m', 'LineWidth', 1)
hold off
legend('Data 1', 'Data2', 'Data3', 'Data4')
title('interpolated standard deviations of nuc data')

% for cyt

interpcStd1 = scale1.*data1.std_cyt;
interpcStd2 = interp1(data2.t, scale2.*data2.std_cyt, data1.t, 'nearest', 'extrap');
interpcStd3 = interp1(data3.t, scale3.*data3.std_cyt, data1.t, 'nearest', 'extrap');
interpcStd4 = interp1(data4.t, scale4.*data4.std_cyt, data1.t, 'nearest', 'extrap');

figure('Visible', 'off'),
plot(data1.t, interpcStd1,'r', 'LineWidth', 2)
hold on
plot(data1.t, interpcStd2,'b')
hold on
plot(data1.t, interpcStd3,'g')
hold on
plot(data1.t, interpcStd4,'m')
hold off
legend('Data 1', 'Data2', 'Data3', 'Data4')
title('interpolated standard deviations of cyt data')


% for nuc
sq_nuc = ((interpStd1'./interpData1').^2)+((interpStd2./interpData2).^2)+...
    ((interpStd3./interpData3).^2)+((interpStd4./interpData4).^2);

std_nuc = (mean_nuc.*sqrt(sq_nuc))/4;

% for cyt
sq_cyt = ((interpcStd1'./interpcData1').^2)+((interpcStd2./interpcData2).^2)+...
    ((interpcStd3./interpcData3).^2)+((interpcStd4./interpcData4).^2);

std_cyt = (mean_cyt.*sqrt(sq_cyt))/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the mean and std
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors1 = cbrewer('seq', 'PuRd', 9); 
colors2 = cbrewer('seq', 'PuBuGn', 9);


% figure,
% set(gcf,'Position',[100 100 400 300],...
% 		'Paperpositionmode','auto','Color',[1 1 1])
% plot(data1.t,mean_cyt,'Color',colors1(7, :),'LineWidth',2)
% hold on 
% plot(data1.t,mean_nuc,'Color',colors2(7, :),'LineWidth',2)
% hold on
% patch([data1.t' fliplr(data1.t')], [(mean_cyt-std_cyt)' fliplr((mean_cyt+std_cyt)')],colors1(5, :), 'FaceAlpha', .2)
% hold on
% patch([data1.t' fliplr(data1.t')], [(mean_nuc-std_nuc)' fliplr((mean_nuc+std_nuc)')],colors2(5, :), 'FaceAlpha', .2)
% xlabel('Time(in min)')
% ylabel('mean intensity')
% legend('mean cyt int.','mean nuc int.','\mu_{cyt}\pm\sigma','\mu_{nuc}\pm\sigma','Location','best')
% title('Averaged plot of 4 embryo cross sections')
% 
%  baseFolder = fullfile('Figures', 'Averaged embryo figures and data');
%  if ~exist(baseFolder, 'dir')
%      mkdir(baseFolder);
%  end


% Create the figure
figure('Visible', 'off'); 
set(gcf,'Position',[100 100 400 300], ...
        'Paperpositionmode','auto','Color',[1 1 1])

% Plot mean intensity curves
plot(data1.t, mean_cyt, 'Color', colors1(7, :), 'LineWidth', 2)
hold on
plot(data1.t, mean_nuc, 'Color', colors2(7, :), 'LineWidth', 2)
hold on

% Plot shaded regions for standard deviation
patch([data1.t' fliplr(data1.t')], ...
      [(mean_cyt - std_cyt)' fliplr((mean_cyt + std_cyt)')], ...
      colors1(5, :), 'FaceAlpha', .2)

patch([data1.t' fliplr(data1.t')], ...
      [(mean_nuc - std_nuc)' fliplr((mean_nuc + std_nuc)')], ...
      colors2(5, :), 'FaceAlpha', .2)

% Labels and legend
xlabel('Time (in min)')
ylabel('Mean intensity')
legend('mean cyt int.', 'mean nuc int.', ...
       '\mu_{cyt} \pm \sigma', '\mu_{nuc} \pm \sigma', 'Location', 'best')
title('Averaged plot of 4 embryo cross sections')

% Create the folder if it doesn't exist
baseFolder = fullfile('Figures', 'Averaged embryo figures and data');
if ~exist(baseFolder, 'dir')
    mkdir(baseFolder);
end

% Define output file name
figName = fullfile(baseFolder, 'Averaged plot of 4 embryo cross sections');

% Save as JPG and FIG
saveas(gcf, [figName, '.jpg']);
saveas(gcf, [figName, '.fig']);

% Close figure to free memory
close(gcf);


 








