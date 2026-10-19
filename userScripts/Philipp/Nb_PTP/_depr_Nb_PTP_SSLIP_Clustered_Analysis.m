%% Nb_PTP_SSLIP_Clustered_Analysis.m
% Script to cluster closely related slip systems (same b, similar n)
% and plot their combined activity on a CI backdrop.

clear; close all; clc;

%% 1. Initialization & Directories
plotx2east; plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir  = fullfile(baseDir, 'Processed_Output');

% --- Clustering Parameters ---
dic_prefix = 'B0005';
angle_threshold = 17 * degree; % Threshold to cluster slip planes (degrees)
activity_threshold = 0.00005; % Minimum mean activity to consider a system for clustering
use_trace_angle = true; % Set to true to cluster by trace angle, false for plane normal angle
use_single_color_range = true; % Toggle to use a single color range for all plots
single_color_range = [0, 0.125]; % Default single color range
compare_prefixes = {'B0007', 'B0003'};
max_compare_clusters = 1; % Maximum cluster pairs to compare (e.g., 1;1, 2;2, 3;3)
% -----------------------------


%% 2. Load Data
% Load Undeformed EBSD for CI backdrop
undeformed_dic_mat = fullfile(outDir, 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');
if ~exist(undeformed_dic_mat, 'file')
    error('Undeformed EBSD data not found at %s', undeformed_dic_mat);
end
fprintf('Loading undeformed EBSD with DIC data...\n');
undefData = load(undeformed_dic_mat, 'newEBSD');
ebsdUndef = undefData.newEBSD;
ebsdUndef = ebsdUndef.gridify;

% Compute tuned CI
blackPoint = 0.05; whitePoint = 0.98; gammaVal = 1.5;
ci_raw = ebsdUndef.ci;
ci_norm = (ci_raw - min(ci_raw(:))) ./ (max(ci_raw(:)) - min(ci_raw(:)));
ci_tuned = (ci_norm - blackPoint) ./ (whitePoint - blackPoint);
ci_tuned(ci_tuned < 0) = 0; ci_tuned(ci_tuned > 1) = 1;
ci_tuned = ci_tuned .^ gammaVal;
cibw = [ci_tuned(:), ci_tuned(:), ci_tuned(:)];

% Load SSLIP Data
sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', dic_prefix));
if ~exist(sslip_mat_path, 'file')
    error('SSLIP data not found for prefix %s.', dic_prefix);
end
fprintf('Loading SSLIP data for %s...\n', dic_prefix);
load(sslip_mat_path, 'ebsdID', 'sSLocalRef');

% --- Apply Specimen Mask ---
mask_file_undef = fullfile(baseDir, 'insitu Specimen 1 0 Map Data 1_BW263_specimen_mask.mat');
if exist(mask_file_undef, 'file')
    fprintf('Applying specimen mask to remove vacuum artifacts...\n');
    load(mask_file_undef, 'specimenlocation');
    % Map the 1D original mask to the aligned/warped dataset using 'oldId'
    if isfield(ebsdID.prop, 'oldId')
        old_id = ebsdID.prop.oldId;
        valid_id = ~isnan(old_id);
        mask_aligned = false(size(old_id));
        mask_aligned(valid_id) = specimenlocation(old_id(valid_id));
        ebsdID(~mask_aligned) = 'notIndexed';
        ebsdUndef(~mask_aligned) = 'notIndexed';
    else
        warning('ebsdID lacks ''oldId''. Cannot apply specimen mask.');
    end
else
    warning('No specimen mask found for Increment 0. Proceeding with unmasked data.');
end

% Reconstruct grains
grainThresh = 5 * degree;
minSize = 1000;
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

% Assume the two largest grains are our bicrystal
[~, sorted_grain_idx] = sort(grainsUndef.numPixel, 'descend');
top2_idx = sorted_grain_idx(1:min(2, length(grainsUndef)));
[~, left_to_right] = sort(grainsUndef(top2_idx).centroid.x);
bicrystal_idx = top2_idx(left_to_right);

g1_id = grainsUndef(bicrystal_idx(1)).id;
g2_id = grainsUndef(bicrystal_idx(2)).id;

mask_g1 = (ebsdID.grainId == g1_id);
mask_g2 = (ebsdID.grainId == g2_id);

ori1 = grainsUndef(g1_id).meanOrientation;
ori2 = grainsUndef(g2_id).meanOrientation;

sSLocal_1 = ori1 * sSLocalRef;
sSLocal_2 = ori2 * sSLocalRef;

slip_mags_g1 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g1(:))), 2, 'omitnan');
slip_mags_g2 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g2(:))), 2, 'omitnan');

%% 3. Clustering Logic
clusters_g1 = cluster_and_plot(g1_id, slip_mags_g1, sSLocal_1, mask_g1, ebsdID, ebsdUndef, cibw, grainsUndef, dic_prefix, outDir, angle_threshold, activity_threshold, sSLocalRef, use_trace_angle, use_single_color_range, single_color_range);
clusters_g2 = cluster_and_plot(g2_id, slip_mags_g2, sSLocal_2, mask_g2, ebsdID, ebsdUndef, cibw, grainsUndef, dic_prefix, outDir, angle_threshold, activity_threshold, sSLocalRef, use_trace_angle, use_single_color_range, single_color_range);

fprintf('\nFinished processing all clusters.\n');


%% ---------------------------------------------------------
% % Manually specify which clusters to plot side-by-side in section 4:
% manual_g1_cluster_idx = 1; % E.g., Cluster 1 from Grain 1
% manual_g2_cluster_idx = 1; % E.g., Cluster 2 from Grain 2
% % ---------------------------------------------------------

%% 4. Manual Combined Cluster Plotting

% if ~isempty(clusters_g1) && ~isempty(clusters_g2) && manual_g1_cluster_idx <= length(clusters_g1) && manual_g2_cluster_idx <= length(clusters_g2)
%     fprintf('\n--- Generating Manual Combined Plot for G1-C%d and G2-C%d ---\n', manual_g1_cluster_idx, manual_g2_cluster_idx);
%     fig_manual = figure('Name', 'Manual Combined Cluster Plot', 'Visible', 'off');
% 
%     % Plot CI Backdrop
%     plot(ebsdUndef, cibw, 'micronbar', 'off');
%     hold on;
% 
%     % Plot black backgound
%     plot(ebsdUndef, zeros(size(cibw)), 'micronbar', 'off');
%     warning('CI background replaced with black background.')
%     hold on;
% 
%     % Combine slip data for G1
%     combined_slip_g1 = zeros(1, length(ebsdID));
%     for idx = clusters_g1{manual_g1_cluster_idx}
%         vals = abs(ebsdID.prop.slipIDcor(idx, mask_g1));
%         combined_slip_g1(mask_g1) = combined_slip_g1(mask_g1) + vals(:)';
%     end
% 
%     % Combine slip data for G2
%     combined_slip_g2 = zeros(1, length(ebsdID));
%     for idx = clusters_g2{manual_g2_cluster_idx}
%         vals = abs(ebsdID.prop.slipIDcor(idx, mask_g2));
%         combined_slip_g2(mask_g2) = combined_slip_g2(mask_g2) + vals(:)';
%     end
% 
%     % Plot Combined Slip Activity
%     plot_data = nan(1, length(ebsdID));
%     plot_data(mask_g1) = combined_slip_g1(mask_g1);
%     plot_data(mask_g2) = combined_slip_g2(mask_g2);
%     plot(ebsdID, plot_data, 'micronbar', 'off'); %, 'FaceAlpha', 0.65);
% 
%     mtexColorMap viridis;
% 
%     plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
% 
%     % Plot traces
%     main_sys_g1 = clusters_g1{manual_g1_cluster_idx}(1);
%     main_sys_g2 = clusters_g2{manual_g2_cluster_idx}(1);
%     hold on
%     quiver(grainsUndef(g1_id), sSLocal_1(main_sys_g1).trace, 'color', 'r', 'linewidth', 2);
%     hold on
%     quiver(grainsUndef(g2_id), sSLocal_2(main_sys_g2).trace, 'color', 'r', 'linewidth', 2);
% 
%     title(sprintf('Combined Clusters: G%d-C%d & G%d-C%d', g1_id, manual_g1_cluster_idx, g2_id, manual_g2_cluster_idx));
% 
%     % Build string for G1
%     planes_g1 = '';
%     for idx = clusters_g1{manual_g1_cluster_idx}
%         planes_g1 = [planes_g1, char(sSLocalRef(idx).n.round), ', '];
%     end
%     planes_g1(end-1:end) = []; % Remove trailing comma
% 
%     % Build string for G2
%     planes_g2 = '';
%     for idx = clusters_g2{manual_g2_cluster_idx}
%         planes_g2 = [planes_g2, char(sSLocalRef(idx).n.round), ', '];
%     end
%     planes_g2(end-1:end) = [];
% 
%     info_str = sprintf('G%d Cluster %d:\n  b=[%s]\n  involved n: %s\n\nG%d Cluster %d:\n  b=[%s]\n  involved n: %s', ...
%         g1_id, manual_g1_cluster_idx, char(sSLocalRef(main_sys_g1).b.round), planes_g1, ...
%         g2_id, manual_g2_cluster_idx, char(sSLocalRef(main_sys_g2).b.round), planes_g2);
% 
%     % Add annotation box
%     old_units = get(gca, 'Units');
%     set(gca, 'Units', 'normalized');
%     dim = get(gca, 'Position');
%     set(gca, 'Units', old_units);
%     dim(1) = dim(1) + 0.02; % Inset slightly from bottom left
%     dim(2) = dim(2) + 0.02;
%     dim(3) = 0.1; dim(4) = 0.1;
%     annotation('textbox', dim, 'String', info_str, ...
%         'Interpreter', 'none', 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', ...
%         'FitBoxToText', 'on', 'FontSize', 10, 'Margin', 5, 'VerticalAlignment', 'bottom');
% 
%     if use_single_color_range
%         setColorRange(single_color_range);
%     else
%         setColorRange([0,0.025]);
%     end
%     mtexColorbar('title', 'Combined Slip Activity (\gamma)');
% 
%     if use_trace_angle, ang_suffix = 't'; else, ang_suffix = 'n'; end
%     fig_path = fullfile(outDir, sprintf('SSLIP_Combined_G%d-C%d_G%d-C%d_%s_%s.jpg', g1_id, manual_g1_cluster_idx, g2_id, manual_g2_cluster_idx, dic_prefix, ang_suffix));
%     exportScaledFigure(fig_manual, fig_path, 'target_px', 3000, 'noDateString');
%     fprintf('Saved combined manual plot to %s\n', fig_path);
%     close(fig_manual);
% else
%     fprintf('\nSkipping Manual Plot: Specified clusters not found.\n');
% end

%% 5. Compare Same Clusters Across Different Increments

% Combine the base dic_prefix with the ones to compare
all_prefixes = [{dic_prefix}, compare_prefixes];

% Determine how many cluster pairs we can actually compare
num_compare = min([max_compare_clusters, length(clusters_g1), length(clusters_g2)]);

if num_compare > 0
    fprintf('\n--- Comparing Clusters Across Increments ---\n');
    for p_idx = 1:length(all_prefixes)
        cmp_prefix = all_prefixes{p_idx};
        cmp_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', cmp_prefix));
        
        if ~exist(cmp_mat_path, 'file')
            warning('SSLIP data not found for prefix %s. Skipping.', cmp_prefix);
            continue;
        end
        
        fprintf('Loading SSLIP data for %s...\n', cmp_prefix);
        cmp_data = load(cmp_mat_path, 'ebsdID');
        cmp_ebsdID = cmp_data.ebsdID;
        
        % Apply mask if available
        if exist('mask_aligned', 'var')
            cmp_ebsdID(~mask_aligned) = 'notIndexed';
        end
        
        for c_idx = 1:num_compare
            fprintf('  -> Plotting G1-C%d & G2-C%d for %s\n', c_idx, c_idx, cmp_prefix);
            
            fig_cmp = figure('Name', sprintf('Combined Cluster %s', cmp_prefix), 'Visible', 'off');
            
            % CI Backdrop
            plot(ebsdUndef, cibw, 'micronbar', 'off');
            hold on;
            plot(ebsdUndef, zeros(size(cibw)), 'micronbar', 'off');
            hold on;
            
            % Combine slip data for G1
            combined_slip_g1_cmp = zeros(1, length(cmp_ebsdID));
            for idx = clusters_g1{c_idx}
                vals = abs(cmp_ebsdID.prop.slipIDcor(idx, mask_g1));
                combined_slip_g1_cmp(mask_g1) = combined_slip_g1_cmp(mask_g1) + vals(:)';
            end
            
            % Combine slip data for G2
            combined_slip_g2_cmp = zeros(1, length(cmp_ebsdID));
            for idx = clusters_g2{c_idx}
                vals = abs(cmp_ebsdID.prop.slipIDcor(idx, mask_g2));
                combined_slip_g2_cmp(mask_g2) = combined_slip_g2_cmp(mask_g2) + vals(:)';
            end
            
            % Plot Combined Slip Activity
            plot_data_cmp = nan(1, length(cmp_ebsdID));
            plot_data_cmp(mask_g1) = combined_slip_g1_cmp(mask_g1);
            plot_data_cmp(mask_g2) = combined_slip_g2_cmp(mask_g2);
            
            plot(cmp_ebsdID, plot_data_cmp, 'micronbar', 'off');
            mtexColorMap viridis;
            plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
            
            % Plot traces
            main_sys_g1 = clusters_g1{c_idx}(1);
            main_sys_g2 = clusters_g2{c_idx}(1);
            hold on
            quiver(grainsUndef(g1_id), sSLocal_1(main_sys_g1).trace, 'color', 'r', 'linewidth', 2);
            hold on
            quiver(grainsUndef(g2_id), sSLocal_2(main_sys_g2).trace, 'color', 'r', 'linewidth', 2);
            
            title(sprintf('%s - Combined Clusters: G%d-C%d & G%d-C%d', cmp_prefix, g1_id, c_idx, g2_id, c_idx));
            
            % Annotation
            planes_g1 = '';
            for idx = clusters_g1{c_idx}
                planes_g1 = [planes_g1, char(sSLocalRef(idx).n.round), ', '];
            end
            if length(planes_g1) >= 2, planes_g1(end-1:end) = []; end
            
            planes_g2 = '';
            for idx = clusters_g2{c_idx}
                planes_g2 = [planes_g2, char(sSLocalRef(idx).n.round), ', '];
            end
            if length(planes_g2) >= 2, planes_g2(end-1:end) = []; end
            
            info_str = sprintf('G%d Cluster %d:\n  b=[%s]\n  involved n: %s\n\nG%d Cluster %d:\n  b=[%s]\n  involved n: %s', ...
                g1_id, c_idx, char(sSLocalRef(main_sys_g1).b.round), planes_g1, ...
                g2_id, c_idx, char(sSLocalRef(main_sys_g2).b.round), planes_g2);
            
            old_units = get(gca, 'Units');
            set(gca, 'Units', 'normalized');
            dim = get(gca, 'Position');
            set(gca, 'Units', old_units);
            dim(1) = dim(1) + 0.02; dim(2) = dim(2) + 0.02;
            dim(3) = 0.1; dim(4) = 0.1;
            annotation('textbox', dim, 'String', info_str, ...
                'Interpreter', 'none', 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', ...
                'FitBoxToText', 'on', 'FontSize', 10, 'Margin', 5, 'VerticalAlignment', 'bottom');
            
            if use_single_color_range
                setColorRange(single_color_range);
            else
                setColorRange([0,0.025]);
            end
            mtexColorbar('title', 'Combined Slip Activity (\gamma)');
            
            if use_trace_angle, ang_suffix = 't'; else, ang_suffix = 'n'; end
            fig_path = fullfile(outDir, sprintf('SSLIP_Combined_G%d-C%d_G%d-C%d_%s_%s.jpg', g1_id, c_idx, g2_id, c_idx, cmp_prefix, ang_suffix));
            exportScaledFigure(fig_cmp, fig_path, 'target_px', 3000, 'noDateString');
            close(fig_cmp);
        end
    end
    fprintf('Finished comparing clusters across increments.\n');
end

%% Helper Functions
function clusters = cluster_and_plot(grain_id, slip_mags, sSLocal, mask, ebsdID, ebsdUndef, cibw, grainsUndef, dic_prefix, outDir, angle_threshold, activity_threshold, sSLocalRef, use_trace_angle, use_single_color_range, single_color_range)
    
    fprintf('\n--- Clustering for Grain %d ---\n', grain_id);
    
    % Only consider systems with activity > threshold
    active_idx = find(slip_mags > activity_threshold);
    [~, sort_order] = sort(slip_mags(active_idx), 'descend');
    sorted_active_idx = active_idx(sort_order);
    
    assigned = false(size(sSLocal));
    clusters = {};
    
    for i = 1:length(sorted_active_idx)
        sys1_idx = sorted_active_idx(i);
        if assigned(sys1_idx)
            continue;
        end
        
        current_cluster = sys1_idx;
        assigned(sys1_idx) = true;
        
        % Look for other systems to cluster with this one
        for j = i+1:length(sorted_active_idx)
            sys2_idx = sorted_active_idx(j);
            if assigned(sys2_idx)
                continue;
            end
            
            % Check if they have the same Burgers vector (or very close)
            b1 = sSLocal(sys1_idx).b;
            b2 = sSLocal(sys2_idx).b;
            b_angle = angle(b1, b2, 'noSymmetry', 'antipodal');
            
            if b_angle < 1e-4
                if use_trace_angle
                    t1 = sSLocal(sys1_idx).trace;
                    t2 = sSLocal(sys2_idx).trace;
                    compare_angle = angle(t1, t2);
                else
                    n1 = sSLocal(sys1_idx).n;
                    n2 = sSLocal(sys2_idx).n;
                    compare_angle = angle(n1, n2, 'noSymmetry', 'antipodal');
                end
                
                if compare_angle <= angle_threshold
                    current_cluster = [current_cluster, sys2_idx];
                    assigned(sys2_idx) = true;
                end
            end
        end
        clusters{end+1} = current_cluster;
    end
    
    % Plot the clusters
    if isempty(clusters)
        fprintf('No active clusters found for Grain %d.\n', grain_id);
        return;
    end
    
    num_clusters = length(clusters);
    plot_cols = min(3, num_clusters);
    plot_rows = ceil(num_clusters / plot_cols);
    
    fig = figure('Name', sprintf('Clustered SSLIP Grain %d - %s', grain_id, dic_prefix), 'Visible', 'off');
    newMtexFigure('layout', [plot_rows, plot_cols], 'figSize', 'large');
    
    for k = 1:num_clusters
        nextAxis(ceil(k/plot_cols), mod(k-1, plot_cols)+1);
        cluster_indices = clusters{k};
        
        % Combine slip data
        combined_slip = zeros(1, length(ebsdID));
        for idx = cluster_indices
            vals = abs(ebsdID.prop.slipIDcor(idx, mask));
            combined_slip(mask) = combined_slip(mask) + vals(:)';
        end
        
        % Plot CI Backdrop
        plot(ebsdUndef, cibw, 'micronbar', 'off');
        hold on;

        plot(ebsdUndef, zeros(size(cibw)), 'micronbar', 'off');
        warning('CI background replaced with black background.')
        hold on;
        
        % Plot Combined Slip Activity
        plot_data = nan(size(combined_slip));
        plot_data(mask) = combined_slip(mask);
        plot(ebsdID, plot_data, 'micronbar', 'off'); %, 'FaceAlpha', 0.65);
        hold on
        
        mtexColorMap viridis;
        plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
        
        % Plot main trace and b
        main_sys_idx = cluster_indices(1);
        hold on
        quiver(grainsUndef(grain_id), sSLocal(main_sys_idx).trace, 'color', 'r', 'linewidth', 2);
        b_vec = sSLocal(main_sys_idx).b.normalize;
        b_vec.antipodal = true;
        hold on
        quiver(grainsUndef(grain_id), b_vec, 'color', 'w', 'linewidth', 2);
        
        % Prepare Annotation Box
        % Find the rank of the main system globally
        [~, global_sort] = sort(slip_mags, 'descend');
        
        info_str = sprintf('Cluster Total Mag: %.4f\n', sum(slip_mags(cluster_indices)));
        info_str = [info_str, sprintf('Included Systems:\n')];
        
        for idx = cluster_indices
            sys_mag = slip_mags(idx);
            global_rank = find(global_sort == idx);
            ref_b = char(sSLocalRef(idx).b.round);
            ref_n = char(sSLocalRef(idx).n.round);
            
            % Compute angle to the main system
            if idx == main_sys_idx
                ang_str = 'Main';
            else
                if use_trace_angle
                    val = angle(sSLocal(main_sys_idx).trace, sSLocal(idx).trace) / degree;
                    ang_str = sprintf('%.1f deg (trace)', val);
                else
                    val = angle(sSLocal(main_sys_idx).n, sSLocal(idx).n, 'noSymmetry', 'antipodal') / degree;
                    ang_str = sprintf('%.1f deg (normal)', val);
                end
            end
            
            info_str = [info_str, sprintf('Rank %d (Sys %d): Mag %.4f\n  b=[%s] n=(%s)\n  Angle: %s\n', ...
                global_rank, idx, sys_mag, ref_b, ref_n, ang_str)];
        end
        
        title(sprintf('Cluster %d', k));
        
        % Add annotation box inside the axes boundaries
        old_units = get(gca, 'Units');
        set(gca, 'Units', 'normalized');
        dim = get(gca, 'Position');
        set(gca, 'Units', old_units);
        
        dim(3) = 0.1; dim(4) = 0.1; % FitBoxToText will shrink it around the text automatically
        annotation('textbox', dim, 'String', info_str, ...
            'Interpreter', 'none', 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', ...
            'FitBoxToText', 'on', 'FontSize', 8, 'Margin', 5, 'VerticalAlignment', 'top');
            
        fprintf('Cluster %d: Systems [%s], Total Mag: %.4f\n', k, num2str(cluster_indices), sum(slip_mags(cluster_indices)));
    end
    
    if use_single_color_range
        setColorRange(single_color_range);
    else
        setColorRange('equal');
    end
    mtexColorbar('title', 'Combined Slip Activity (\gamma)');
    
    if use_trace_angle
        ang_suffix = 't';
    else
        ang_suffix = 'n';
    end
    fig_path = fullfile(outDir, sprintf('SSLIP_Clustered_%s_Grain%d_%s.jpg', dic_prefix, grain_id, ang_suffix));
    exportScaledFigure(fig, fig_path, 'target_px', 3000, 'noDateString');
    fprintf('Saved clustered plot to %s\n', fig_path);
    close(fig);
end
