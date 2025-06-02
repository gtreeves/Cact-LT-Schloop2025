function F = ftn_output2(data,F,IM)
%Makes movie of dl-GFP channel and graph.
%
%function data = movie_output2a(filename,varargin)
%
% This function reads in the movie structure "F" from:
%	[filename(1:end-4),'_graph.mat']
% It also reads in the images from:
%	[filename(1:end-4),'_images.mat']
% That is, unless the variable "imgs" is passed in "varargin".  
%
% "filename": string containing the name of the lsm file that the images we
%	analyzing came from, including path to that file (ex:
%	'Images/NUC12A.lsm'). 
%
% Optional argument varargin can consist of these things, in this order:
%	* "F": Structure output from "movie_output1," and it is in the format
%		of "getframe".
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.

% %
% % Unpacking varargin.
% %
% nArg = size(varargin,2); iArg = 1; 
% if nArg >= iArg && ~isempty(varargin{iArg})
% 	F = varargin{iArg}; else 
% 	%
% end, iArg = iArg + 1;
% if nArg >= iArg && ~isempty(varargin{iArg})
% 	IM = varargin{iArg}; else 
% 	%
% end%, iArg = iArg + 1;

if isfield(data,'metadata') && isfield(data.metadata,'np_ch')
	data_ch = data.metadata.np_ch;
else
	data_ch = 1; % arb
end
if ~exist('F','var')
	F = dosage_moviedraw(data);
end

if ~exist('IM','var')
	IM = dosage_imread(data.filename);
end

nmovie = 384; % height of final movie
dH = round(0.25*nmovie/2); % padding above and below needed for the graph.
T = data.D;
s_mid = data.s_mid;
theta = (s_mid*180+90)*pi/180; % the rotation angle.
H1 = data.H; W1 = data.W;

% Xc = data.Xc; Yc = data.Yc; 
Xp = data.metadata.Xp; Yp = data.metadata.Yp;
% Xc = cellfun(@mean,Xp);
% Yc = cellfun(@mean,Yp);

R = data.metadata.arc/2/pi;
tmesh = data.t;


%
% Getting background subtraction and most intense image
%
% bg = zeros(T,1); 
maxint = zeros(T,1);
for i = 1:T
	M1 = IM(:,:,:,i);
	M3 = M1(:,:,data_ch);
% 	v = M3(1,:);
% 	bg(i) = mean(v);
	maxint(i) = prctile(double(M3(:)),99); % we allow oversaturation
end
% bg = max(bg);
% bg = data.metadata.bg(data_ch); % it was already background subtracted
maxint = max(maxint);
r = round(max(R)); 
% H = 2*max(R) + 50; W = H; % the "50" is for padding.
H = 2*round(r) + 50; W = H; % the "50" is for padding.
xc = round(W/2); yc = round(H/2); 


for i = 1:T
	
	I_dorsal1 = IM(:,:,data_ch,i); 
	
	%
	% Location of embryo in this time point wrt original image
	%
	[Xc,Yc,R] = circfit(Xp{i}(1:end-1),Yp{i}(1:end-1));
	
	%
	% Now we register the images with respect to all timepoints.
	%
	I_dorsal = zeros(H,W);
% 	I_dorsal1 = IM(:,:,1,i); 
% 	I_nuc = I_dorsal;
% 	I_nuc1 = IM{i,2};
% 	IM{i,1} = []; IM{i,2} = [];
	j1 = max(round(Xc - R),1);
% 	j2 = min(round(Xc(i) + R(i)),N(i));
	j2 = min(round(Xc + R),W1);
	i1 = max(round(Yc - R),1);
% 	i2 = min(round(Yc(i) + R(i)),M(i));
	i2 = min(round(Yc + R),H1);
	
	j01 = round(xc - R);
	j02 = round(xc + R);
	i01 = round(yc - R);
	i02 = round(yc + R);
	if j2-j1 ~= j02-j01
		j02 = j02 - (j02-j01 - (j2-j1));
	end
	if i2-i1 ~= i02-i01
		i02 = i02 - (i02-i01 - (i2-i1));
	end
	I_dorsal(i01:i02,j01:j02) = double(I_dorsal1(i1:i2,j1:j2));
% 	I_nuc(i01:i02,j01:j02) = double(I_nuc1(i1:i2,j1:j2));
	


	%
	% Now we should have an embryo that's in the center of the image.
	% Rotating the embryo so that the vm is down, and centering it in the
	% rotated image.  We also stretch the embryo to better fill the
	% resulting movie dimensions, and perform global intensity
	% normalization. 
	%
	Irot = imrotate(I_dorsal,theta*180/pi);
	[m,n] = size(Irot);
	xc1 = round(n/2); yc1 = round(m/2);
% 	Irot = imsubtract(Irot,bg);
% 	Irot = uint8(255*Irot/(maxint-bg));
	Irot = uint8(255*Irot/maxint); % it was already background subtracted
	
% 	M3 = uint8(zeros(H,W,3));
% 	M4 = M2(y1:y2,x1:x2,j);
% 	M4 = imsubtract(M4,bg);
% 	M4 = uint8(255*double(imresize(M4,[H W]))/(maxint-bg));
% 	M3 = M4;	
	
	
	I1 = cat(3,Irot,Irot,Irot);
	I1 = I1(yc1-r-24:yc1+r+24,xc1-r-24:xc1+r+24,:);
	I1 = imresize(I1,[nmovie nmovie]);
	I2 = F(i).cdata;
	I2 = imresize(I2,[round(0.75*nmovie) nmovie]);
	Z = uint8(zeros(dH,nmovie,3));
	I2 = [Z;I2;Z]; %#ok<AGROW>
	H2 = size(I2,1);
	if H2 > nmovie
		I2(1:(H2-nmovie),:,:) = [];
	elseif size(I2,1) < nmovie
		I2 = [Z(1:(nmovie-H2),:,:);I2]; %#ok<AGROW>
	end
	
	I = [I1,I2];
	F(i).cdata = I;
	F(i).colormap = [];
end










