function [label, uniq_lbl] = track_nuclei_labels_exact(mask_data, label, max_lookback)
% Faithful recreation of the original tracking and label assignment logic.
% Tracks nuclei over time with optional lookback to past N frames (max_lookback).
%
% Inputs:
%   mask_data    - structure containing nuclei info per frame
%   label        - cell array of label matrices per frame
%   max_lookback - how many frames to look back (1 = only i-1)
%
% Outputs:
%   label        - updated label matrices per frame
%   uniq_lbl     - updated list of unique label values

uniq_lbl = unique(label{1,1});

for i = 2:size(mask_data,1)
    fixed = mask_data{i-1,4};
    moving = mask_data{i,4};

    % Register current frame to previous
    [tform, ~, ~] = pcregistericp(moving, fixed);
    movingReg = pctransform(moving, tform);
    mask_data{i,5} = movingReg;

    % Get nearest-neighbor correspondences
    [k, dist] = dsearchn(fixed.Location(:,1:2), movingReg.Location(:,1:2));
    temp = double(mask_data{i,1});
    used = [];

    % One-to-one match with acceptable distance
    if (length(k) == length(unique(k))) && max(dist) < 30
        nonrep = find(0 < dist & dist < 30);
        for lennr = 1:numel(nonrep)
            idx = k(nonrep(lennr));
            cord = round(mask_data{i-1,3}(idx, :));
            prev_img = label{i-1,1};
            temp(mask_data{i,2}(nonrep(lennr)).PixelIdxList) = prev_img(cord(2), cord(1));
        end
        uniq_lbl = unique([uniq_lbl; unique(temp)]);

    else
        % Handle both matched and unmatched
        nonrep = find(0 < dist & dist < 30);
        rep = find(dist >= 30);
        for lennr = 1:numel(nonrep)
            idx = k(nonrep(lennr));
            cord = round(mask_data{i-1,3}(idx, :));
            val = label{i-1,1}(cord(2), cord(1));
            used = [used val];
            temp(mask_data{i,2}(nonrep(lennr)).PixelIdxList) = val;
        end

        for lenr = 1:numel(rep)
            assigned = false;
            val_candidates = [];
            dist_candidates = [];

            for b = 2:min(max_lookback, i-1)
                [kb, db] = dsearchn(mask_data{i-b,3}, mask_data{i,3}(rep(lenr),:));
                if db < 30
                    cordb = round(mask_data{i-b,3}(kb,:));
                    valb = label{i-b,1}(cordb(2), cordb(1));
                    val_candidates(end+1) = valb;
                    dist_candidates(end+1) = db;
                else
                    val_candidates(end+1) = -1; % invalid
                    dist_candidates(end+1) = Inf;
                end
            end

            % Evaluate fallback logic
            valid = val_candidates(val_candidates > 0 & dist_candidates < 30);
            if isempty(valid)
                new_val = max(uniq_lbl) + 1;
                temp(mask_data{i,2}(rep(lenr)).PixelIdxList) = new_val;
                uniq_lbl = [uniq_lbl; new_val];
            else
                % Sort candidates by distance
                [sorted_dists, idx] = sort(dist_candidates);
                val1 = val_candidates(idx(1));
                val2 = length(idx) > 1 && val_candidates(idx(2));

                if ismember(val1, used) && (isempty(val2) || ismember(val2, used))
                    new_val = max(uniq_lbl) + 1;
                    temp(mask_data{i,2}(rep(lenr)).PixelIdxList) = new_val;
                    uniq_lbl = [uniq_lbl; new_val];
                elseif ~ismember(val1, used)
                    temp(mask_data{i,2}(rep(lenr)).PixelIdxList) = val1;
                    used = [used val1];
                elseif ~isempty(val2) && ~ismember(val2, used)
                    temp(mask_data{i,2}(rep(lenr)).PixelIdxList) = val2;
                    used = [used val2];
                else
                    new_val = max(uniq_lbl) + 1;
                    temp(mask_data{i,2}(rep(lenr)).PixelIdxList) = new_val;
                    uniq_lbl = [uniq_lbl; new_val];
                end
            end
        end
    end

    label{i,1} = temp;
    label{i,2} = mask_data{i,2};
end
end
