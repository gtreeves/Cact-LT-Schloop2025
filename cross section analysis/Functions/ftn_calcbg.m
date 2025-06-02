function [bg,sig0,lsminf2] = ftn_calcbg(filename,yesplot)
%
% This function is designed to find the background calibration files and
% analyze them to get the background level and standard deviation, as well
% as the so-called "S-factor," which is the relationship between the
% variance and the intensity.
%
% Input can either be a filename or an image that you're calculating the
% background of.
%
% Optional input "yesplot": if false (default), then no plotting happens,
%	and if true, then a plot is made and saved in the working directory
%	under the filename (short) with "_bg.jpg" at the end. If no filename is
%	given (that is, if the input "filename" is really an image), then a
%	clock is the name of the output file.

if ischar(filename)
	[IM,~,lsminf2] = openczi(filename);
else
	IM = filename;
	lsminf2 = NaN;
end

	
% cziinf = lsminf2.cziinfo.cziinf.Channel;
% fnames = fieldnames(cziinf);
% v1 = strfindDU(fnames,'Wavelength1');
% v2 = strfindDU(fnames,'Wavelength2');
% if any(v1) && any(v2)
% 	lambda = [cziinf.Wavelength1 cziinf.Wavelength2];
% 	data_ch = find(lambda > 450 & lambda < 500);
% 	mask_ch = find(lambda > 500 & lambda < 600);
% else
% 	data_ch = 1;
% end

% v_ch = [data_ch mask_ch];
% bg = zeros(length(v_ch),1);
A = size(IM);
%numch = A(3);
numch = size(IM,3);
bg = zeros(1,numch);
for o = 1:numch

	I0 = double(IM(:,:,o,1)); % assume frame 1 is representative

	%
	% Extracting the slice and finding the "offset"; that is, the
	% background, or "dark current," which is characterized as the most
	% populated peak.
	%
	I0 = I0(:);
	m = max(I0);
	[n,x] = hist(I0(:),0:m);
	[A,k] = max(n(2:end-1)); % cut off ends because of possibility of saturation
	k = k + 1;
	% 	bg = x(k);
	[~,ksig] = min(abs(n(k+1:end-1)-exp(-0.5)*A)); ksig = ksig + k; % estimate 1 sigma
% 	[~,k14] = min(abs(n(k+1:end-1)-exp(-2)*A)); k14 = k14 + k; % estimate 2 sigmas
	xsig = x(ksig);
	sig = xsig-x(k);

	%
	% Now fit to a gaussian
	%
	[A,B,mu,sig] = fit_gaussian2(x(2:ksig)',n(2:ksig)',sig);
	bg(o) = mu;
	sig0(o) = sig;


% 	%
% 	% There is some interval around the background that contains some
% 	% variance. To figure out this variance, we will fit a Gaussian to that
% 	% interval around the background. So first we find the interval, then
% 	% fit to Gaussian
% 	%
% 	vminus = find(n(1:k) <= 1); % looking for gray levels achieved by only zero or one pixel
% 	vplus = find(n(k+1:end) <= 1) + k;
% 	if isempty(vminus), vminus = 1; end
% 	if isempty(vplus), vplus = m; end
% 	vminus = vminus(end); vplus = vplus(1);
% 
% 	V = I0 >= vminus & I0 <= vplus;
% 	Y = I0(V);
% 	bg(o) = mean(Y);

% 	if v_ch(o) == data_ch
% 		sig0 = std(Y);
% 	end

end

if ischar(filename)
	save([filename(1:end-4),'_bg.mat'],'bg','sig0','lsminf2');
end
if exist('yesplot','var') && yesplot
	if ~ischar(filename)

		clk1 = clock;
		clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
			num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];
		filename = [clk,'_bg.jpg'];
	else
		filename = [filename(1:end-4),'_bg.jpg'];
	end

	figure('visib','off')
	bar(x(1:201),n(1:201))
	hold on
	plot(x(1:201),A*exp(-(x(1:201)-mu).^2/2/sig^2)+B,'linewidth',2)
	print(gcf,filename,'-djpeg','-r300')
	close(gcf)

end









