%% Nb_PTP_TemporalEvolution.m
% Script to track the temporal evolution of slip system activity across all 
% DIC increments and compare it against the lattice rotation (misorientation) 
% from TKD/GND data as a measure of SSLIP reliability.

clear; close all; clc;

%% 1. Initialization & Directories
if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir  = fullfile(baseDir, 'Processed_Output');
undeformed_dic_mat = fullfile(outDir, 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');

% Generate list of all 14 DIC prefixes
dic_prefixes = arrayfun(@(x) sprintf('B%04d', x), 1:14, 'UniformOutput', false);

%% 2. Find Top 4 Systems based on B0005 Benchmark
ref_prefix = 'B0005';
ref_sslip_file = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', ref_prefix));

if ~exist(ref_sslip_file, 'file')
    error('Reference SSLIP data (%s) not found. Please run the analysis for B0005 first.', ref_prefix);
end

fprintf('Loading reference benchmark data from %s...\n', ref_prefix);
load(ref_sslip_file, 'ebsdID', 'sSLocalRef');

% Reconstruct grains to identify the bicrystal
grainThresh = 5 * degree;
minSize = 1000;
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

% Ensure Grain 1 is always the "left" grain and Grain 2 is the "right" grain
[~, sorted_grain_idx] = sort(grainsUndef.numPixel, 'descend');
top2_idx = sorted_grain_idx(1:2);
[~, left_to_right] = sort(grainsUndef(top2_idx).centroid.x);
bicrystal_idx = top2_idx(left_to_right);

g1_id = grainsUndef(bicrystal_idx(1)).id;
g2_id = grainsUndef(bicrystal_idx(2)).id;

mask_g1 = (ebsdID.grainId == g1_id);
mask_g2 = (ebsdID.grainId == g2_id);

ori_undef_g1 = grainsUndef(g1_id).meanOrientation;
ori_undef_g2 = grainsUndef(g2_id).meanOrientation;

%% --- PLOT GRAIN IDENTIFICATION MAP ---
fprintf('\nPlotting and saving grain identification map...\n');
fig_grains = figure;
fig_grains.Name = 'Grain Identification';
plot(ebsdID, ebsdID.orientations, 'micronbar', 'off');
hold on;
plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);

text(grainsUndef(g1_id), sprintf('Grain %d\n(Left)', g1_id), 'color', 'k', 'fontsize', 14, 'fontweight', 'bold', 'backgroundcolor', 'w', 'edgecolor', 'k', 'margin', 2);
text(grainsUndef(g2_id), sprintf('Grain %d\n(Right)', g2_id), 'color', 'k', 'fontsize', 14, 'fontweight', 'bold', 'backgroundcolor', 'w', 'edgecolor', 'k', 'margin', 2);

title(sprintf('Bicrystal Identification (Prefix: %s)', ref_prefix));
fig_grains_path = fullfile(outDir, 'Grain_Identification_Map.jpg');
exportScaledFigure(fig_grains, fig_grains_path, 'target_px', 2500, 'noDateString');
fprintf('Saved grain identification figure to %s\n\n', fig_grains_path);

% Find top 4 dominant slip systems per grain in B0005
slip_mags_g1_ref = mean(abs(ebsdID.prop.slipIDcor(:, mask_g1(:))), 2, 'omitnan');
slip_mags_g2_ref = mean(abs(ebsdID.prop.slipIDcor(:, mask_g2(:))), 2, 'omitnan');

[~, sort_g1] = sort(slip_mags_g1_ref, 'descend');
[~, sort_g2] = sort(slip_mags_g2_ref, 'descend');

top4_g1 = sort_g1(1:4);
top4_g2 = sort_g2(1:4);

fprintf('Benchmark Top 4 systems for Grain %d: %s\n', g1_id, num2str(top4_g1'));
fprintf('Benchmark Top 4 systems for Grain %d: %s\n', g2_id, num2str(top4_g2'));

%% --- SCHMID FACTOR ANALYSIS ---
% Assuming uniaxial tension along Y-direction based on fullAnalysis script
sigma = stressTensor.uniaxial(yvector);
sSLocal_g1 = ori_undef_g1 * sSLocalRef;
sSLocal_g2 = ori_undef_g2 * sSLocalRef;

SF_g1 = abs(sSLocal_g1.SchmidFactor(sigma));
SF_g1 = SF_g1(:); % Force column vector for consistent sorting
SF_g2 = abs(sSLocal_g2.SchmidFactor(sigma));
SF_g2 = SF_g2(:);

[SF_g1_sorted, SF_g1_rank] = sort(SF_g1, 'descend');
[SF_g2_sorted, SF_g2_rank] = sort(SF_g2, 'descend');

fprintf('\n--- Schmid Factor Analysis for Top 4 Systems (B0005) ---\n');
for k = 1:4
    sys_idx = top4_g1(k);
    sf_val = SF_g1(sys_idx);
    sf_rank = find(SF_g1_rank == sys_idx, 1);
    delta_sf = SF_g1_sorted(1) - sf_val;
    fprintf('Grain %d - Active Rank %d (Sys %d): SF = %.3f | SF Rank = %d (%d systems higher) | Delta SF to Max = %.3f\n', ...
        g1_id, k, sys_idx, sf_val, sf_rank, sf_rank - 1, delta_sf);
end
fprintf('\n');
for k = 1:4
    sys_idx = top4_g2(k);
    sf_val = SF_g2(sys_idx);
    sf_rank = find(SF_g2_rank == sys_idx, 1);
    delta_sf = SF_g2_sorted(1) - sf_val;
    fprintf('Grain %d - Active Rank %d (Sys %d): SF = %.3f | SF Rank = %d (%d systems higher) | Delta SF to Max = %.3f\n', ...
        g2_id, k, sys_idx, sf_val, sf_rank, sf_rank - 1, delta_sf);
end
fprintf('------------------------------------------------------\n\n');

%% 3. Iterate through all increments
num_steps = length(dic_prefixes);
pct_activity_g1 = nan(num_steps, 4);
pct_activity_g2 = nan(num_steps, 4);
abs_activity_g1 = nan(num_steps, 4);
abs_activity_g2 = nan(num_steps, 4);
misorientation_g1 = nan(num_steps, 1);
misorientation_g2 = nan(num_steps, 1);
total_slip_activity_g1 = nan(num_steps, 1);
total_slip_activity_g2 = nan(num_steps, 1);
residual_eff_g1 = nan(num_steps, 1);
residual_eff_g2 = nan(num_steps, 1);

% Lazily load undeformed EBSD only if we need to calculate missing SSLIP data
ebsdUndef_loaded = false;

for i = 1:num_steps
    prefix = dic_prefixes{i};
    fprintf('\n--- Processing Increment %s (%d/%d) ---\n', prefix, i, num_steps);
    
    sslip_file = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', prefix));
    
    % --- CALCULATE MISSING SSLIP DATA ---
    if ~exist(sslip_file, 'file') 
        error('\n--- No .mat file found for increment  %s (%d/%d). Please run the full analysis script first. ---\n', prefix, i, num_steps)
    else
        fprintf('  -> Loading existing SSLIP data...\n');
        loop_data = load(sslip_file, 'ebsdID', 'sSLocalRef');
        ebsdID_loop = loop_data.ebsdID;
        sSLocalRef_loop = loop_data.sSLocalRef;
    end
    
    % Reconstruct grains for the current loop to ensure consistent "left"/"right" assignment
    [grainsUndef_loop, ebsdID_loop.grainId] = calcGrains(ebsdID_loop, 'minPixel', minSize, 'angle', grainThresh);
    ebsdID_loop(grainsUndef_loop(grainsUndef_loop.numPixel<=minSize)) = 'notIndexed';
    [grainsUndef_loop, ebsdID_loop.grainId] = calcGrains(ebsdID_loop, 'minPixel', minSize, 'angle', grainThresh);

    [~, sorted_grain_idx_loop] = sort(grainsUndef_loop.numPixel, 'descend');
    top2_idx_loop = sorted_grain_idx_loop(1:2);
    [~, left_to_right_loop] = sort(grainsUndef_loop(top2_idx_loop).centroid.x);
    bicrystal_idx_loop = top2_idx_loop(left_to_right_loop);

    g1_id_loop = grainsUndef_loop(bicrystal_idx_loop(1)).id;
    g2_id_loop = grainsUndef_loop(bicrystal_idx_loop(2)).id;

    mask_g1_loop = (ebsdID_loop.grainId == g1_id_loop);
    mask_g2_loop = (ebsdID_loop.grainId == g2_id_loop);

    % Find corresponding slip systems in current increment
    top4_g1_loop = zeros(4,1);
    top4_g2_loop = zeros(4,1);
    
    for k = 1:4
        for j = 1:length(sSLocalRef_loop)
            if angle(sSLocalRef(top4_g1(k)).b, sSLocalRef_loop(j).b, 'noSymmetry') < 1e-5 && angle(sSLocalRef(top4_g1(k)).n, sSLocalRef_loop(j).n, 'noSymmetry') < 1e-5
                top4_g1_loop(k) = j;
            end
            if angle(sSLocalRef(top4_g2(k)).b, sSLocalRef_loop(j).b, 'noSymmetry') < 1e-5 && angle(sSLocalRef(top4_g2(k)).n, sSLocalRef_loop(j).n, 'noSymmetry') < 1e-5
                top4_g2_loop(k) = j;
            end
        end
    end
    
    if any(top4_g1_loop == 0) || any(top4_g2_loop == 0)
        warning('Could not match all top 4 slip systems for increment %s. Fallback to original indices.', prefix);
        top4_g1_loop(top4_g1_loop == 0) = top4_g1(top4_g1_loop == 0);
        top4_g2_loop(top4_g2_loop == 0) = top4_g2(top4_g2_loop == 0);
    end
    
    % --- EXTRACT ACCOMMODATION PERCENTAGE ---
    slip_mags_g1 = mean(abs(ebsdID_loop.prop.slipIDcor(:, mask_g1_loop(:))), 2, 'omitnan');
    slip_mags_g2 = mean(abs(ebsdID_loop.prop.slipIDcor(:, mask_g2_loop(:))), 2, 'omitnan');
    
    total_slip_g1 = sum(slip_mags_g1);
    total_slip_g2 = sum(slip_mags_g2);
    
    total_slip_activity_g1(i) = total_slip_g1;
    total_slip_activity_g2(i) = total_slip_g2;
    
    if isfield(ebsdID_loop.prop, 'residualEeff')
        residual_eff_g1(i) = sum(ebsdID_loop.prop.residualEeff(mask_g1_loop(:)), 'omitnan');
        residual_eff_g2(i) = sum(ebsdID_loop.prop.residualEeff(mask_g2_loop(:)), 'omitnan');
    else
        warning('residualEeff not found for increment %s.', prefix);
    end

    if total_slip_g1 > 0
        pct_activity_g1(i, :) = (slip_mags_g1(top4_g1_loop) / total_slip_g1) * 100;
        abs_activity_g1(i, :) = slip_mags_g1(top4_g1_loop);
    end
    if total_slip_g2 > 0
        pct_activity_g2(i, :) = (slip_mags_g2(top4_g2_loop) / total_slip_g2) * 100;
        abs_activity_g2(i, :) = slip_mags_g2(top4_g2_loop);
    end
    
    % --- EXTRACT LATTICE ROTATION (SSLIP ERROR) ---
    gnd_file = fullfile(outDir, sprintf('GND_Data_%s.mat', prefix));
    if exist(gnd_file, 'file')
        fprintf('  -> Found TKD data. Calculating grain rotation with respect to the undeformed orientation (used for SSLIP)...\n');
        gnd_data = load(gnd_file, 'grainsDef');
        if isfield(gnd_data, 'grainsDef')
            % Find misorientation for Grain 1
            [min_angle_g1, ~] = min(angle(gnd_data.grainsDef.meanOrientation, ori_undef_g1));
            if min_angle_g1 < 15 * degree
                misorientation_g1(i) = min_angle_g1 ./ degree;
            end
            
            % Find misorientation for Grain 2
            [min_angle_g2, ~] = min(angle(gnd_data.grainsDef.meanOrientation, ori_undef_g2));
            if min_angle_g2 < 15 * degree
                misorientation_g2(i) = min_angle_g2 ./ degree;
            end
        end
    end
end

% Save output array so you don't have to rerun the analysis
save(fullfile(outDir, 'TemporalEvolution_Data.mat'), 'dic_prefixes', 'pct_activity_g1', ...
    'pct_activity_g2', 'abs_activity_g1', 'abs_activity_g2', 'total_slip_activity_g1', ...
    'total_slip_activity_g2', 'residual_eff_g1', 'residual_eff_g2', ...
    'misorientation_g1', 'misorientation_g2', 'top4_g1', ...
    'top4_g2', 'sSLocalRef', 'g1_id', 'g2_id', 'num_steps');

%% 4. Plotting
if ~exist('pct_activity_g1', 'var')
    data_path = fullfile(outDir, 'TemporalEvolution_Data.mat');
    if exist(data_path, 'file')
        fprintf('\nLoading previously saved data from %s...\n', data_path);
        load(data_path);
    else
        error('Required data not found in workspace or in %s. Please run the previous sections first.', data_path);
    end
end

% Fallbacks in case an older .mat file was loaded
if ~exist('num_steps', 'var') && exist('dic_prefixes', 'var'), num_steps = length(dic_prefixes); end
if ~exist('g1_id', 'var'), g1_id = 1; g2_id = 2; end

fprintf('\nGenerating dashboard...\n');
fig = figure;
fig.Name = 'Temporal Evolution of SSLIP and Misorientation';
newMtexFigure('layout', [2, 4], 'figSize', 'large');

x_axis = 1:num_steps;

% Generate legend strings based on B0005 systems
leg_g1 = arrayfun(@(x) sprintf('Sys %d: b=[%s] n=(%s)', x, char(sSLocalRef(x).b.round), char(sSLocalRef(x).n.round)), top4_g1, 'UniformOutput', false);
leg_g2 = arrayfun(@(x) sprintf('Sys %d: b=[%s] n=(%s)', x, char(sSLocalRef(x).b.round), char(sSLocalRef(x).n.round)), top4_g2, 'UniformOutput', false);

nextAxis(1, 1); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Top 4 SSLIP Systems', g1_id), leg_g1, pct_activity_g1, 'Relative Slip Activity (%)'); xtickangle(45);
nextAxis(2, 1); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Top 4 SSLIP Systems', g2_id), leg_g2, pct_activity_g2, 'Relative Slip Activity (%)'); xtickangle(45);

nextAxis(1, 2); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Top 4 SSLIP Systems (Abs)', g1_id), leg_g1, abs_activity_g1, 'Absolute Slip Activity (\gamma)'); set(gca, 'YScale', 'log'); xtickangle(45);
nextAxis(2, 2); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Top 4 SSLIP Systems (Abs)', g2_id), leg_g2, abs_activity_g2, 'Absolute Slip Activity (\gamma)'); set(gca, 'YScale', 'log'); xtickangle(45);

% nextAxis(1, 3); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Total Slip Activity', g1_id), {'Total Slip \Sigma|\gamma|'}, total_slip_activity_g1, 'Total Slip \Sigma|\gamma|'); xtickangle(45);
% nextAxis(2, 3); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Total Slip Activity', g2_id), {'Total Slip \Sigma|\gamma|'}, total_slip_activity_g2, 'Total Slip \Sigma|\gamma|'); xtickangle(45);

nextAxis(1, 3); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Sum of Residual Eeff', g1_id), {'\Sigma residualEeff'}, residual_eff_g1, 'Sum of Residual Eeff'); xtickangle(45);
nextAxis(2, 3); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Sum of Residual Eeff', g2_id), {'\Sigma residualEeff'}, residual_eff_g2, 'Sum of Residual Eeff'); xtickangle(45);

nextAxis(1, 4); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Rotation From Undeformed (Proxy for SSLIP accuracy.)', g1_id), {'Misorientation from Undeformed'}, misorientation_g1, 'Rotation Angle (\circ)'); xtickangle(45);
nextAxis(2, 4); plot_panel(x_axis, dic_prefixes, sprintf('Grain %d: Rotation From Undeformed (Proxy for SSLIP accuracy.)', g2_id), {'Misorientation from Undeformed'}, misorientation_g2, 'Rotation Angle (\circ)'); xtickangle(45);

drawNow(gcm);
fig_path = fullfile(outDir, 'TemporalEvolution_SlipAndMisorientation.jpg');
exportScaledFigure(gcf, fig_path, 'target_px', 3500);
fprintf('Saved dashboard figure to %s\n', fig_path);

%% 5. Plot SSLIP Activity Maps for All Increments
fprintf('\n--- Generating SSLIP Activity Maps ---\n');

% Ensure grain reconstruction thresholds exist if data was loaded from .mat
if ~exist('minSize', 'var'), minSize = 1000; end
if ~exist('grainThresh', 'var'), grainThresh = 5 * degree; end

for i = 1:num_steps
    prefix = dic_prefixes{i};
    fprintf('  -> Processing maps for %s (%d/%d)...\n', prefix, i, num_steps);
    
    sslip_file = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', prefix));
    if ~exist(sslip_file, 'file')
        warning('SSLIP data not found for %s. Skipping.', prefix);
        continue;
    end
    
    loop_data = load(sslip_file, 'ebsdID', 'sSLocalRef');
    ebsdID_loop = loop_data.ebsdID;
    sSLocalRef_loop = loop_data.sSLocalRef;
    
    % Reconstruct grains to determine left/right assignments
    [grainsUndef_loop, ebsdID_loop.grainId] = calcGrains(ebsdID_loop, 'minPixel', minSize, 'angle', grainThresh);
    ebsdID_loop(grainsUndef_loop(grainsUndef_loop.numPixel<=minSize)) = 'notIndexed';
    [grainsUndef_loop, ebsdID_loop.grainId] = calcGrains(ebsdID_loop, 'minPixel', minSize, 'angle', grainThresh);

    [~, sorted_grain_idx_loop] = sort(grainsUndef_loop.numPixel, 'descend');
    top2_idx_loop = sorted_grain_idx_loop(1:2);
    [~, left_to_right_loop] = sort(grainsUndef_loop(top2_idx_loop).centroid.x);
    bicrystal_idx_loop = top2_idx_loop(left_to_right_loop);

    g1_id_loop = grainsUndef_loop(bicrystal_idx_loop(1)).id;
    g2_id_loop = grainsUndef_loop(bicrystal_idx_loop(2)).id;

    % Find corresponding slip systems in current increment
    top4_g1_loop = zeros(4,1);
    top4_g2_loop = zeros(4,1);
    
    for k = 1:4
        for j = 1:length(sSLocalRef_loop)
            if angle(sSLocalRef(top4_g1(k)).b, sSLocalRef_loop(j).b, 'noSymmetry') < 1e-5 && angle(sSLocalRef(top4_g1(k)).n, sSLocalRef_loop(j).n, 'noSymmetry') < 1e-5
                top4_g1_loop(k) = j;
            end
            if angle(sSLocalRef(top4_g2(k)).b, sSLocalRef_loop(j).b, 'noSymmetry') < 1e-5 && angle(sSLocalRef(top4_g2(k)).n, sSLocalRef_loop(j).n, 'noSymmetry') < 1e-5
                top4_g2_loop(k) = j;
            end
        end
    end
    
    if any(top4_g1_loop == 0) || any(top4_g2_loop == 0)
        top4_g1_loop(top4_g1_loop == 0) = top4_g1(top4_g1_loop == 0);
        top4_g2_loop(top4_g2_loop == 0) = top4_g2(top4_g2_loop == 0);
    end
    
    plot_sslip_activity_maps(ebsdID_loop, grainsUndef_loop, g1_id_loop, g2_id_loop, top4_g1_loop, top4_g2_loop, sSLocalRef_loop, prefix, outDir);
end

%% 6. Plot GND Activity Maps for All Increments
fprintf('\n--- Generating GND Activity Maps ---\n');
% Find Top 4 GND Systems based on B0005 Benchmark (if exists)
ref_gnd_file = fullfile(outDir, sprintf('GND_Data_%s.mat', ref_prefix));
if exist(ref_gnd_file, 'file')
    fprintf('Loading reference benchmark GND data from %s...\n', ref_prefix);
    gnd_ref_data = load(ref_gnd_file, 'ebsdDeformed', 'rho', 'factor', 'dS', 'grainsDef');
    
    gnd_density_sum = sum(abs(gnd_ref_data.rho), 1, 'omitnan');
    [~, sorted_gnd_idx] = sort(gnd_density_sum, 'descend');
    top4_gnd = sorted_gnd_idx(1:4);
    dS_ref = gnd_ref_data.dS;
    
    fprintf('Benchmark Top 4 GND systems: %s\n', num2str(top4_gnd));
    
    for i = 1:num_steps
        prefix = dic_prefixes{i};
        gnd_file = fullfile(outDir, sprintf('GND_Data_%s.mat', prefix));
        if ~exist(gnd_file, 'file')
            fprintf('  -> GND data not found for %s. Skipping.\n', prefix);
            continue;
        end
        
        fprintf('  -> Processing GND maps for %s (%d/%d)...\n', prefix, i, num_steps);
        loop_gnd = load(gnd_file, 'ebsdDeformed', 'rho', 'factor', 'dS', 'grainsDef');
        
        % Match GND systems across increments
        top4_loop = zeros(4,1);
        
        for k = 1:4
            for j = 1:length(loop_gnd.dS)
                if angle(dS_ref(top4_gnd(k)).b, loop_gnd.dS(j).b, 'noSymmetry') < 1e-5 && ...
                   angle(dS_ref(top4_gnd(k)).l, loop_gnd.dS(j).l, 'noSymmetry') < 1e-5 && ...
                   dS_ref(top4_gnd(k)).isEdge == loop_gnd.dS(j).isEdge
                    top4_loop(k) = j;
                end
            end
        end
        
        if any(top4_loop == 0)
            top4_loop(top4_loop == 0) = top4_gnd(top4_loop == 0);
        end
        
        plot_gnd_activity_maps(loop_gnd.ebsdDeformed, loop_gnd.grainsDef, loop_gnd.rho, loop_gnd.factor, ...
            top4_loop, loop_gnd.dS, prefix, outDir);
    end
else
    fprintf('Reference GND data (%s) not found. Cannot determine benchmark systems.\n', ref_gnd_file);
end

%% Helper Functions
function plot_panel(x_axis, dic_prefixes, title_str, leg_str, y_data_left, y_label_left)
    plot(x_axis, y_data_left, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    ylabel(y_label_left);
    
    title(title_str);
    xlabel('DIC Increment');
    xticks(x_axis);
    xticklabels(dic_prefixes);
    grid on;
    legend(leg_str, 'Location', 'best', 'Interpreter', 'none');
end

function plot_sslip_activity_maps(ebsdID, grainsUndef, g1_id, g2_id, top4_g1, top4_g2, sSLocalRef, prefix, outDir)
    fig_map = figure;
    fig_map.Name = sprintf('SSLIP Activity Maps - %s', prefix);
    newMtexFigure('layout', [2, 2], 'figSize', 'large');

    mask_g1 = (ebsdID.grainId == g1_id);
    mask_g2 = (ebsdID.grainId == g2_id);
    
    % Find global maximum activity for consistent color scaling across all 4 ranks
    global_max = 0;
    for k = 1:4
        d1 = abs(ebsdID.prop.slipIDcor(top4_g1(k), mask_g1(:)));
        d2 = abs(ebsdID.prop.slipIDcor(top4_g2(k), mask_g2(:)));
        global_max = max([global_max, max(d1(:)), max(d2(:))]);
    end
    global_max = max(1e-4, global_max); % Avoid 0 limit if perfectly empty

    for k = 1:4
        row_idx = ceil(k/2);
        col_idx = mod(k-1, 2) + 1;
        nextAxis(row_idx, col_idx);
        
        plot_data = nan(1, length(ebsdID));
        sys_g1 = top4_g1(k);
        sys_g2 = top4_g2(k);
        
        plot_data(mask_g1(:)) = abs(ebsdID.prop.slipIDcor(sys_g1, mask_g1(:)));
        plot_data(mask_g2(:)) = abs(ebsdID.prop.slipIDcor(sys_g2, mask_g2(:)));
        
        plot(ebsdID, plot_data, 'micronbar', 'off');
        mtexColorMap viridis;
        hold on;
        plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
        
        % Format slip system strings
        sS_g1 = sSLocalRef(sys_g1);
        sS_g2 = sSLocalRef(sys_g2);
        str_g1 = sprintf('G1 Sys %d: b=[%s] n=(%s)', sys_g1, char(sS_g1.b.round), char(sS_g1.n.round));
        str_g2 = sprintf('G2 Sys %d: b=[%s] n=(%s)', sys_g2, char(sS_g2.b.round), char(sS_g2.n.round));
        
        title(sprintf('Rank %d Activity\n%s\n%s', k, str_g1, str_g2), 'Interpreter', 'none', 'FontSize', 10);
    end
    
    setColorRange([0, global_max]);
    mtexColorbar('title', 'Absolute Slip Activity (\gamma)');
    drawNow(gcm);
    
    fig_path = fullfile(outDir, sprintf('SSLIP_Top4_Activity_%s.jpg', prefix));
    exportScaledFigure(fig_map, fig_path, 'target_px', 2500, 'noDateString');
    close(fig_map); % Keep workspace clean
end

function plot_gnd_activity_maps(ebsdDeformed, grainsDef, rho, factor, top4_gnd, dS, prefix, outDir)
    fig_map = figure;
    fig_map.Name = sprintf('GND Activity Maps - %s', prefix);
    newMtexFigure('layout', [2, 2], 'figSize', 'large');

    for k = 1:4
        row_idx = ceil(k/2);
        col_idx = mod(k-1, 2) + 1;
        nextAxis(row_idx, col_idx);
        
        sys_idx = top4_gnd(k);
        dS_max = dS(sys_idx);
        gnd_data = factor * abs(rho(:, sys_idx));
        
        plot(ebsdDeformed, log10(gnd_data), 'micronbar', 'off');
        hold on;
        plot(grainsDef.boundary, 'linewidth', 1.5, 'linecolor', 'w');
        mtexColorMap parula;
        
        typeStr = {'Edge', 'Screw'};
        title(sprintf('Rank %d GND System:\n%s, b=[%s] l=[%s]', k, typeStr{2 - dS_max.isEdge}, char(dS_max.b.round), char(dS_max.l.round)));
    end
    
    setColorRange([14.5, 16.5]); % Limits mapped for typical log10(m^-2) densities
    mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
    drawNow(gcm);
    
    fig_path = fullfile(outDir, sprintf('GND_Top4_Activity_%s.jpg', prefix));
    exportScaledFigure(fig_map, fig_path, 'target_px', 2500, 'noDateString');
    close(fig_map); % Keep workspace clean
end
