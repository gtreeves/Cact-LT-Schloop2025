% reads all the czi files that are broken into different parts and save the
% data as a .mat file

% run openFilesWithSamePrefix_savedata('C:\Users\sharv\OneDrive\Desktop\Cact frap czi codes for paper\part data');

function openFilesWithSamePrefix_savedata(directory)
    % This function recursively searches a directory for .czi microscopy files.
    % It groups files with the same 15-character prefix, combines them,
    % extracts image and time metadata, and saves the output as a .mat file
    % in a structured folder under "mat files for data".

    % Step 1: List contents of the input directory
    contents = dir(directory);

    % Initialize containers to store file info
    fileNames = {};
    prefixes = {};
    filePaths = {};

    % Step 2: Loop over all items in the directory
    for i = 1:numel(contents)
        % Skip '.' and '..'
        if strcmp(contents(i).name, '.') || strcmp(contents(i).name, '..')
            continue;
        end

        itemPath = fullfile(directory, contents(i).name);

        % If it's a .czi file, extract its prefix and record its info
        if ~isfolder(itemPath) && endsWith(contents(i).name, '.czi', 'IgnoreCase', true)
            prefix = contents(i).name(1:15);  % Use first 15 characters to group
            fileNames{end+1} = contents(i).name;
            filePaths{end+1} = itemPath;
            prefixes{end+1} = prefix;
        end

        % If it's a folder, recurse into it
        if isfolder(itemPath)
            openFilesWithSamePrefix_savedata(itemPath);
        end
    end

    % Step 3: Identify unique prefixes (i.e., unique groups of .czi files)
    uniquePrefixes = unique(prefixes);

    % Step 4: Loop through each group of files with matching prefix
    for i = 1:numel(uniquePrefixes)
        prefix = uniquePrefixes{i};
        matchIdx = strcmp(prefixes, prefix);
        matchingFiles = fileNames(matchIdx);
        matchingPaths = filePaths(matchIdx);

        % Sort files alphabetically to ensure consistent order
        [~, sortIdx] = sort(matchingFiles);
        matchingFiles = matchingFiles(sortIdx);
        matchingPaths = matchingPaths(sortIdx);

        % Step 5: Loop over each file in the current prefix group
        for j = 1:numel(matchingPaths)
            filePath = matchingPaths{j};
            disp(['Processing file: ' filePath]);

            % Extract clean filename (e.g., "Embryo 1 Spot 1") by removing "5sec"
            [~, filebase, ~] = fileparts(filePath);
            tokens = regexp(filebase, '(Embryo\s*\d+\s*Spot\s*\d+)', 'tokens', 'once');
            if ~isempty(tokens)
                filenameshort = tokens{1};  % Use cleaned version
            else
                filenameshort = filebase;   % Use full name as fallback
                warning('Could not clean filenameshort. Using full filename.');
            end

            % Get the folder name like "2023-04-21 LT-Cact FRAP"
            parentFolder = getParentFolderName(filePath);

            % First file in the group
            if j == 1
                [IM2, lsminf1, lsminf2] = openczi(filePath);
                scalings = lsminf2.VOXELSIZES * 1e6;  % Convert to microns

                % Extract time metadata
                r = bfGetReader(filePath);
                m = r.getMetadataStore();
                nPlanes = m.getPlaneCount(0);
                t1 = zeros(nPlanes/2, 1);
                t2 = zeros(nPlanes/2, 1);
                for k = 0:nPlanes-1
                    if rem(k,2)==1
                        t1((k+1)/2,1) = m.getPlaneDeltaT(0,k).value;
                    else
                        t2((k+2)/2,1) = m.getPlaneDeltaT(0,k).value;
                    end
                end
                t = t1;
                t_prev = t - t(2);  % Time offset
                IM = IM2;  % Initialize image stack

            else
                % For additional files in the group
                [IM3, lsminf1, lsminf2] = openczi(filePath);
                r = bfGetReader(filePath);
                m = r.getMetadataStore();
                nPlanes = m.getPlaneCount(0);
                t1 = zeros(nPlanes/2, 1);
                t2 = zeros(nPlanes/2, 1);
                for k = 0:nPlanes-1
                    if rem(k,2)==1
                        t1((k+1)/2,1) = m.getPlaneDeltaT(0,k).value;
                    else
                        t2((k+2)/2,1) = m.getPlaneDeltaT(0,k).value;
                    end
                end
                t_curr = t1 - t1(2);

                % Combine image stacks and time
                IM = cat(4, IM3, IM);  % Append to existing stack
                t = [t_curr; t_prev + t_curr(end) + 2*abs(t_prev(1))];
                
                % use this if the file name is reverse
%                 IM = cat(4,IM2,IM3);
%                 t = [t_prev ; t_curr+t_prev(end)+2*abs(t_curr(1))];

                % Save combined data only after last file in group
                if j == numel(matchingPaths)
                    data.IM = IM;
                    data.t = t;
                    data.lsminf1 = lsminf1;
                    data.lsminf2 = lsminf2;
                    data.scalings = scalings;

                    % Construct save path
                    saveDir = fullfile('mat files for data', parentFolder);
                    if ~exist(saveDir, 'dir')
                        mkdir(saveDir);
                    end

                    savePath = fullfile(saveDir, [filenameshort '.mat']);
                    save(savePath, 'data');
                end
            end
        end
    end
end

function parentFolder = getParentFolderName(filePath)
    % Returns the immediate parent folder of a file
    [parentPath, ~, ~] = fileparts(filePath);         % remove file name
    [~, parentFolder] = fileparts(parentPath);        % extract last folder
end
