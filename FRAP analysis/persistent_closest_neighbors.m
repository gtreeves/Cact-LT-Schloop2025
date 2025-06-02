function locs = persistent_closest_neighbors(nuc_data, idx_bleach, neighbors)
% persistent_closest_neighbors:
% Returns the labels of the closest `neighbors` to the bleached nucleus
% that are present in all frames of the dataset.
%
% Inputs:
%   - nuc_data   : cell array of nuclear data per frame
%   - idx_bleach : label of bleached nucleus
%   - neighbors  : number of neighbors to return (default = 6)
%
% Output:
%   - locs       : label IDs of closest stable neighbors across all frames

    if nargin < 3
        neighbors = 6;
    end

    num_frames = size(nuc_data, 1);

    % Step 1: Collect label sets per frame
    label_sets = cell(num_frames, 1);
    for fr = 1:num_frames
        label_sets{fr} = [nuc_data{fr,1}.label];
    end

    % Step 2: Find labels common to all frames
    common_labels = label_sets{1};
    for fr = 2:num_frames
        common_labels = intersect(common_labels, label_sets{fr});
    end

    % Step 3: Remove the bleached nucleus from the list
    common_labels(common_labels == idx_bleach) = [];

    % Step 4: Find centroids of common labels in frame 1
    centroids = [];
    labels = [];
    for i = 1:length(nuc_data{1,1})
        if ismember(nuc_data{1,1}(i).label, common_labels)
            centroids = [centroids; nuc_data{1,1}(i).Centroid];
            labels = [labels; nuc_data{1,1}(i).label];
        end
    end

    % Step 5: Find the centroid of the bleached nucleus in frame 1
    bleach_idx = find([nuc_data{1,1}.label] == idx_bleach, 1);
    if isempty(bleach_idx)
        error('Bleached nucleus not found in frame 1.');
    end
    bleach_centroid = nuc_data{1,1}(bleach_idx).Centroid;

    % Step 6: Compute distances and get closest N labels
    [~, dist] = dsearchn(bleach_centroid, centroids);
    [~, sorted_idx] = sort(dist);
    locs = labels(sorted_idx(1:min(neighbors, length(sorted_idx))));

end
