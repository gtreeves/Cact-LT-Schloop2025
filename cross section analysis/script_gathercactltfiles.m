% script_gathercactltfiles
%
% Put all cross-section medea files into one excel sheet.


clear all
close all


%for czi files
% chage this to the location where the data is stored
basepth2 = 'C:\Users\sharv\OneDrive\Desktop\Allison Cact-LT czi codes for paper\Data\';

[pths2,pthsshort2] = readdir2(basepth2,'/');

pths = pths2; 


% Look into our list of folders and extract all czi files.

Filenames = cell(length(pths),1);
Filenamesshort = cell(length(pths),1);
for i = 1:length(pths)
	a = dir(pths{i});
	f = {a.name}';
	vlsm = false(length(f),1);
	for j = 1:length(f)
		if length(f{j}) > 3
			vlsm(j) = strcmp(f{j}(end-2:end),'czi');
		end
	end
	v = vlsm;
	Filenames{i} = [repmat(pths(i),sum(v),1) f(v)];
	Filenamesshort{i} = f(v);
	1;
end
filenames1 = cat(1,Filenames{:});
filenamesshort = cat(1,Filenamesshort{:});

nfiles = length(filenames1);
filenames = cell(nfiles,1);
for i = 1:nfiles
	filenames{i} = [filenames1{i,1},filesep,filenames1{i,2}];
end

% Make short version of the folder names, which include the date

pthsshort = cell(nfiles,1);
for i = 1:nfiles
	filename = filenames{i};
	kslash = strfind(filename,filesep);
	pthsshort{i} = filename(kslash(end-1)+1:kslash(end)-1);
end

save Mat/cactlt_filenames filenames pthsshort filenamesshort pths





