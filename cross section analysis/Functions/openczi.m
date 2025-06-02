function [IM,lsminf1,lsminf2] = openczi(filename)
%
% This script is an attempt to figure out how to open czi files using the
% bio-formats package from ImageJ, specifically "bfopen".
%
% See:
% https://www.openmicroscopy.org/site/support/bio-formats5.5/users/matlab/index.html
%
% Also, note that may need to increase java memory in the future
%
% Input "filename" is a string containing the filename, including path, of
% the czi file you want to open.
%
% Outputs are designed to interface with analyze_xs:
%	"IM": the 4D image data: H-by-W-by-num_ch-by-D, where H and W are the
%		height and width of your frames, num_ch is the number of channels
%		in your data, and D is the number of z-slices.
%	"lsminf1": deprecated structure that used to contain metadata when read
%		from lsm files. Now is empty, but continues to be a placeholder.
%	"lsminf2": structure containing the metadata. It may appear to be set
%		up in an illogical way, but it is organized to match how the
%		function "lsminfo" spits out metadata for lsm files, to ensure
%		back-compatibility with lsm files.

% nseries = size(imgdata,1); % assume only one series for now

%
% Check to see if we have a czi file:
%
if ~strcmp(filename(end-2:end),'czi')
	error('You must only pass czi files to "openczi"')
end

imgdata = bfopen(filename);
imgdata(3:4) = [];


% =========================================================================
% Put together metadata. Need to convert from a java hashtable into a
% Matlab structure variable.  Also need to make it back-compatible with
% reading lsm files.
% =========================================================================
jdata = imgdata{1,2};
jentries = jdata.entrySet.toArray;
nentries = length(jentries);

metadata = cell(nentries,2);
for i = 1:nentries
	metadata{i,1} = jentries(i).getKey;
	metadata{i,2} = jentries(i).getValue;
end

%
% Now that we have all of the metadata extracted from the java hashtable,
% we need to organize it so it makes sense to the human eye.
%
cziinf = struct;
for i = 1:nentries
	a = metadata{i,1}; % fieldname
	val = metadata{i,2}; % value...but is contained in a string.
	try
		val = eval(val);
	catch
		
	end
	
	b = str2cell(a,'|');
	b{end}(end-1) = '_';
	
	c = [b;repmat({'.'},1,length(b))];
	s = strcat(c{:});
	s(end) = [];
	s(s == ' ' | s == '-' | s == '#') = [];
	eval(['cziinf.',s,' = val;']);
	1;
	
	
end


%
% OK, now that the metadata are stored in a structure, let's mine the data
% for what we need. 
% 
% NOTE: for now, we are missing number of time points, and time stamps.
% Also, note that some of these fields might not exist if we are looking at
% a z-stack or a bleach box or something else.
%

% Imagesize
W = cziinf.GlobalInformation.Image.SizeX_1;
H = cziinf.GlobalInformation.Image.SizeY_1;
if ~isfield(cziinf.GlobalInformation.Image,'SizeZ_1')
	D = 1;
else
	D = cziinf.GlobalInformation.Image.SizeZ_1;
end
num_ch = cziinf.GlobalInformation.Image.SizeC_1;
if ~isfield(cziinf.GlobalInformation.Image,'SizeT_1')
	T = 1;
else
	T = cziinf.GlobalInformation.Image.SizeT_1;
end
cziinfo.DimensionXYZ = [H W D]';
cziinfo.num_ch = num_ch;
cziinfo.DimensionT = T;

% Scalings. Number of entries can depend on whether you have a zstack or
% not. Presumably, although this isn't tested, if you have a line-scan,
% maybe you'll only have one entry?
fieldnamesscaling = fieldnames(cziinf.GlobalScaling.Distance);
k = strfind(fieldnamesscaling,'Value');
nscalings = sum(~cellfun(@isempty,k));
scalings = zeros(1,nscalings);
for i = 1:nscalings
	scalings(i) = cziinf.GlobalScaling.Distance.(['Value_',num2str(i)]);
end
cziinfo.Scalings = scalings;
% cziinfo.Scalings = [
% 	cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup.ScalingX_1;
% 	cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup.ScalingY_1;
% 	cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup.ScalingZ_1;
% 	]';

% Laserpower, Detector Gain, Amplifier Offset, Pixeltime, Linetime,
% Frametime (all depend on num_ch).
%
% Laserpower is not always stored in the same place.
if isfield(cziinf.GlobalExperiment.AcquisitionBlock,'Laser')
	Laserinfo = ...
		cziinf.GlobalExperiment.AcquisitionBlock.Laser;
	sLP = 'LaserPower_';
	f = 1;
else
	Laserinfo = ...
		cziinf.GlobalInformation.Instrument.LightSource;
	sLP = 'Power_';
	f = 100;
end
% Detector gain is also not always stored in the same place.
DGinfo = cziinf.GlobalInformation.Image.Channel;
if isfield(DGinfo,'Gain_1')
	sDG = 'Gain_';
else
	sDG = 'Voltage_';
end
speedinfo = cziinf.GlobalInformation.Image.Channel.LaserScanInfo; 
LP = zeros(1,cziinfo.num_ch);
DG = zeros(1,cziinfo.num_ch);
AO = zeros(1,cziinfo.num_ch);
PT = zeros(1,cziinfo.num_ch);
LT = zeros(1,cziinfo.num_ch);
FT = zeros(1,cziinfo.num_ch);
for i = 1:cziinfo.num_ch
	LP(i) = Laserinfo.([sLP,num2str(i)])/f;
	DG(i) = DGinfo.([sDG,num2str(i)]);
	AO(i) = DGinfo.(['Offset_',num2str(i)]);
	PT(i) = speedinfo.(['PixelTime_',num2str(i)]);
	%LT(i) = speedinfo.(['LineTime_',num2str(i)]);
	FT(i) = speedinfo.(['FrameTime_',num2str(i)]);
end
cziinfo.Laserpower = LP;
cziinfo.Detectorgain = DG;
cziinfo.Amplifieroffset = AO;
cziinfo.Pixeltime = PT*1e6;
cziinfo.Linetime = LT;
cziinfo.Frametime = FT;
% XY zoom
cziinfo.ZoomXY = [
	cziinf.GlobalInformation.Image.Channel.LaserScanInfo.ZoomX_1;
	cziinf.GlobalInformation.Image.Channel.LaserScanInfo.ZoomY_1;
	]';

% Offset on zoom
cziinfo.OffsetXY = [
	cziinf.GlobalInformation.Image.Channel.LaserScanInfo.SampleOffsetX_1;
	cziinf.GlobalInformation.Image.Channel.LaserScanInfo.SampleOffsetY_1;
	]';

% rotation
cziinfo.Rotation = cziinf.GlobalInformation.Image.Channel.LaserScanInfo.SampleRotation_1;

% position within the microscope's internal parameters?
% Some czi files don't have this.
if isfield(cziinf.GlobalInformation.Image,'S')
	cziinfo.PositionXYZ = [
		cziinf.GlobalInformation.Image.S.Scene.Position.X_1;
		cziinf.GlobalInformation.Image.S.Scene.Position.Y_1;
		cziinf.GlobalInformation.Image.S.Scene.Position.Z_1;
		]';
else
	cziinfo.PositionXYZ = [];
end

% ROI stuff
if isfield(cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup,'UseRois_1')
	cziinfo.UseRois = ...
		cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup.UseRois_1;
	cziinfo.FitFramesizeToRoi = ...
		cziinf.GlobalExperiment.AcquisitionBlock.AcquisitionModeSetup.FitFramesizeToRoi_1;
else
	cziinfo.UseRois = [];
	cziinfo.FitFramesizeToRoi = true;
end

cziinfo.cziinf = cziinf; % whole honkin' thing

%
% Now re-work the structure to make it back-compatible with the original
% formulation of analyze_xs
%
lsminf2.NUMBER_OF_CHANNELS = cziinfo.num_ch;
lsminf2.VOXELSIZES = cziinfo.Scalings;
lsminf2.ScanInfo.USE_ROIS = cziinfo.UseRois;
lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS = ~cziinfo.FitFramesizeToRoi;
lsminf2.DimensionZ = cziinfo.DimensionXYZ(3);
lsminf2.ScanInfo.PIXEL_TIME = num2cell(cziinfo.Pixeltime);
lsminf2.TimeInterval = cziinfo.Frametime(1);
lsminf2.ScanInfo.POWER = num2cell(cziinfo.Laserpower);
lsminf2.ScanInfo.DETECTOR_GAIN = num2cell(cziinfo.Detectorgain);
lsminf2.cziinfo = cziinfo;

lsminf1 = [];


% =========================================================================
% Organize image data
% =========================================================================

IM = cat(3,imgdata{1,1}{:,1});
IM = reshape(IM,[H,W,num_ch,D*T]); % no timeseries yet





