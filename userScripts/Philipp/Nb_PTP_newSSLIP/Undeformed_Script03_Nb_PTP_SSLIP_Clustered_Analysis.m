%% Nb_PTP_SSLIP_Clustered_Analysis.m
% New script for analyzing bicrystal slip activities using modern SSLIP bundling.
% Features a comparative mode where systems are bundled on a reference increment 
% and the same bundles are tracked across subsequent comparison increments.

clear; close all; clc;

%% 1. Initialization & Directories
plotx2east; plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

if exist('a00_Dataset_Registry', 'file')
    [registry, defaults] = a00_Dataset_Registry();
    baseDir = defaults.dir;
    outDir  = fullfile(baseDir, defaults.outDir_undeformed);
else
    % Fallback if registry not found
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP_newSSLIP');
    outDir  = fullfile(baseDir, 'Processed_Output_Undeformed');
end

plotDir = fullfile(outDir, 'Plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% --- Configuration ---
ref_prefix = 'B0005';
cmp_prefixes = {'B0007', 'B0003'}; % You can add more prefixes here

% --- Plotting Toggles ---
opt_plot = struct();
opt_plot.plot_bicrystal = false;   % Toggle combined bicrystal plots
opt_plot.plot_grain1    = false;   % Toggle single Grain 1 plots (Grain 2 grayed out)
opt_plot.plot_grain2    = true;   % Toggle single Grain 2 plots (Grain 1 grayed out)
opt_plot.withQuivers    = false;   % Generate plots with quiver traces/Burgers vectors
opt_plot.withoutQuivers = true;   % Generate plots without quiver overlays

opt_bundle = struct();
opt_bundle.activityPercentile = 85;
opt_bundle.activityThreshold  = 0; 
opt_bundle.bundleTraceAngleThresh  = 17 * degree;
opt_bundle.num_top_systems = 3;
opt_bundle.custom_clims = {[0 0.125], [0 0.25]}; % Output multiple plots with these hardcoded limits
opt_bundle.quiver_length_factor = 0.2; % Scaling factor for the trace/burgers quiver length

% Toggle clustering logic as requested (3 && 4: Burgers vector + Trace angle)
opt_bundle.bundleKinematic = false;
opt_bundle.bundleCommonSlipSystems = false;
opt_bundle.bundleBurgers = true;
opt_bundle.bundleNormalTrace = true;
opt_bundle.bundleNormal = false; % Or true, depending on exact preference

% Manual Bundling Override
opt_bundle.use_manual_override = true;
opt_bundle.manual_bundles_g1 = {
    [18, 30, 46, 3, 9, 27, 31], ...
    [37, 5, 21, 33, 40, 14, 7, 26, 48, 17], ...
    [47, 16, 11]
};
opt_bundle.manual_bundles_g2 = {
    [27, 3, 15, 41, 30, 18, 2, 46, 38, 9, 19, 31], ...
    [11], ...
    [47, 16, 28, 25]
};
% ---------------------

% Search for the newest alignment file
search_pattern = fullfile(outDir, '*_undeformedEBSD_withAllDIC.mat');
files = dir(search_pattern);
if isempty(files)
    error('Master Aligned EBSD data not found matching %s. Please run Undeformed_Script01_Alignment.m first.', search_pattern);
end
[~, idx] = sort([files.datenum], 'descend');
undeformed_dic_mat = fullfile(outDir, files(idx(1)).name);

if length(files) > 1
    warning('Multiple alignment files found. Loading the most recent one:\n%s', undeformed_dic_mat);
end

fprintf('Loading undeformed EBSD with DIC data...\n');
undefData = load(undeformed_dic_mat, 'newEBSD', 'grains');
ebsdUndef = undefData.newEBSD.gridify;
grainsUndef = undefData.grains;

% Assume the two largest grains are our bicrystal
[~, sorted_grain_idx] = sort(grainsUndef.numPixel, 'descend');
top2_idx = sorted_grain_idx(1:min(2, length(grainsUndef)));
[~, left_to_right] = sort(grainsUndef(top2_idx).centroid.x);
bicrystal_idx = top2_idx(left_to_right);

g1_id = grainsUndef(bicrystal_idx(1)).id;
g2_id = grainsUndef(bicrystal_idx(2)).id;
ori1 = grainsUndef(g1_id).meanOrientation;
ori2 = grainsUndef(g2_id).meanOrientation;

%% 3. Loop over all increments (Reference first)
all_prefixes = [{ref_prefix}, cmp_prefixes];

for p_idx = 1:length(all_prefixes)
    current_prefix = all_prefixes{p_idx};
    
    % Load SSLIP Data
    sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', current_prefix));
    if ~exist(sslip_mat_path, 'file')
        warning('SSLIP data not found for prefix %s. Skipping.', current_prefix);
        continue;
    end
    fprintf('\n======================================================\n');
    fprintf('Processing Increment: %s\n', current_prefix);
    fprintf('======================================================\n');
    load(sslip_mat_path, 'ebsdID', 'sSLocalRef');
    
    % % Legacy data fix: Ensure slipIDcor is [Pixels x Systems]
    % if size(ebsdID.prop.slipIDcor, 1) == length(sSLocalRef)
    %     ebsdID.prop.slipIDcor = ebsdID.prop.slipIDcor';
    % end
    
    % Sync grains to current ebsdID
    ebsdID.grainId = ebsdUndef.grainId;
    mask_g1 = (ebsdID.grainId == g1_id);
    mask_g2 = (ebsdID.grainId == g2_id);

    sSLocal_1 = ori1 * sSLocalRef;
    sSLocal_2 = ori2 * sSLocalRef;
    
    if p_idx == 1
        % We are at the reference increment. Compute bundles!
        if isfield(opt_bundle, 'use_manual_override') && opt_bundle.use_manual_override
            fprintf('--- USING MANUAL BUNDLE OVERRIDE ---\n');
            bundles_g1 = opt_bundle.manual_bundles_g1;
            bundles_g2 = opt_bundle.manual_bundles_g2;
            
            % Compute fallback metric values as max absolute activity across the respective grain
            gamma_abs_g1 = abs(ebsdID.prop.slipIDcor(mask_g1, :));
            gamma_abs_g2 = abs(ebsdID.prop.slipIDcor(mask_g2, :));
            metric_g1 = max(gamma_abs_g1, [], 1);
            metric_g2 = max(gamma_abs_g2, [], 1);
            
            % Identify leftovers (everything not explicitly included in a manual bundle)
            all_sys = 1:size(ebsdID.prop.slipIDcor, 2);
            leftover_g1 = setdiff(all_sys, [bundles_g1{:}]);
            leftover_g2 = setdiff(all_sys, [bundles_g2{:}]);
        else
            fprintf('--- Bundling Reference Grain 1 Systems ---\n');
            ebsdID_g1 = ebsdID; ebsdID_g1(~mask_g1) = 'notIndexed';
            [bundles_g1, metric_g1, leftover_g1] = bundleSSLIPSystems(ebsdID_g1, sSLocal_1, opt_bundle);

            fprintf('--- Bundling Reference Grain 2 Systems ---\n');
            ebsdID_g2 = ebsdID; ebsdID_g2(~mask_g2) = 'notIndexed';
            [bundles_g2, metric_g2, leftover_g2] = bundleSSLIPSystems(ebsdID_g2, sSLocal_2, opt_bundle);
        end

        % Enforce max bundles
        n_top = min([length(bundles_g1), length(bundles_g2), opt_bundle.num_top_systems]);
        % Export CSV Reports for both grains using the shared utility (BEFORE consolidating leftovers for plotting!)
        report_g1 = buildReportData(bundles_g1, leftover_g1, metric_g1, sSLocal_1, ori1, n_top);
        writeSSLIP_Report(fullfile(outDir, sprintf('Bicrystal_Grain1_Report_%s.csv', current_prefix)), report_g1);
        
        report_g2 = buildReportData(bundles_g2, leftover_g2, metric_g2, sSLocal_2, ori2, n_top);
        writeSSLIP_Report(fullfile(outDir, sprintf('Bicrystal_Grain2_Report_%s.csv', current_prefix)), report_g2);
        fprintf('Exported CSV reports for both grains.\n');

        % Consolidate any remaining bundles into leftover activity for plotting purposes
        for idx = n_top+1:length(bundles_g1)
            leftover_g1 = [leftover_g1, bundles_g1{idx}];
        end
        bundles_g1(n_top+1:end) = [];

        for idx = n_top+1:length(bundles_g2)
            leftover_g2 = [leftover_g2, bundles_g2{idx}];
        end
        bundles_g2(n_top+1:end) = [];
        
        fprintf('Bundles cached from reference increment %s.\n', current_prefix);
    else
        fprintf('Using cached bundles from reference increment for comparative plot.\n');
    end

    %% 4. Prepare Plotting Data
    num_subplots = n_top + 1;
    
    % Build combined activity fields for all subplots upfront
    combined_activities = nan(length(ebsdID), num_subplots);
    
    for t_i = 1:num_subplots
        if t_i <= n_top
            act_g1 = zeros(length(ebsdID), 1);
            for mem = 1:length(bundles_g1{t_i})
                act_g1 = act_g1 + abs(ebsdID.prop.slipIDcor(:, bundles_g1{t_i}(mem)));
            end
            
            act_g2 = zeros(length(ebsdID), 1);
            for mem = 1:length(bundles_g2{t_i})
                act_g2 = act_g2 + abs(ebsdID.prop.slipIDcor(:, bundles_g2{t_i}(mem)));
            end
        else
            act_g1 = zeros(length(ebsdID), 1);
            for mem = 1:length(leftover_g1)
                act_g1 = act_g1 + abs(ebsdID.prop.slipIDcor(:, leftover_g1(mem)));
            end
            
            act_g2 = zeros(length(ebsdID), 1);
            for mem = 1:length(leftover_g2)
                act_g2 = act_g2 + abs(ebsdID.prop.slipIDcor(:, leftover_g2(mem)));
            end
        end
        combined_activities(mask_g1(:), t_i) = act_g1(mask_g1(:));
        combined_activities(mask_g2(:), t_i) = act_g2(mask_g2(:));
    end

    % Determine Color Limits
    if isfield(opt_bundle, 'custom_clims') && ~isempty(opt_bundle.custom_clims)
        clims_to_plot = opt_bundle.custom_clims;
    else
        validData = combined_activities(~isnan(combined_activities) & ~isinf(combined_activities));
        if isempty(validData)
            clims_to_plot = {[0, 0.1]};
        else
            dyn_clim = [0, prctile(validData, 99)];
            if dyn_clim(2) <= 0, dyn_clim(2) = 0.05; end
            clims_to_plot = {dyn_clim};
        end
    end
    
    % Background drawing options
    opt_ctx = struct('plotGrainContext', true, 'grains', grainsUndef, 'grainLineColor', 'k', 'grainLineWidth', 2);
    
    % Centroids for grain quiver placements
    cx1 = grainsUndef(g1_id).centroid.x;
    cy1 = grainsUndef(g1_id).centroid.y;
    cx2 = grainsUndef(g2_id).centroid.x;
    cy2 = grainsUndef(g2_id).centroid.y;

    % Determine which quiver modes to render based on toggles
    quiver_modes = [];
    if isfield(opt_plot, 'withQuivers') && opt_plot.withQuivers, quiver_modes(end+1) = true; end
    if isfield(opt_plot, 'withoutQuivers') && opt_plot.withoutQuivers, quiver_modes(end+1) = false; end
    if isempty(quiver_modes), quiver_modes = [true]; end

    %% 4. Combined Bicrystal Plotting
    if opt_plot.plot_bicrystal
        for c_idx = 1:length(clims_to_plot)
            current_clim = clims_to_plot{c_idx};
            
            for plot_quivers = quiver_modes
                if plot_quivers
                    q_str = 'withQuivers';
                else
                    q_str = 'noQuivers';
                end
                
                fig = figure('Name', sprintf('Combined SSLIP Bicrystal Analysis - %s (clim: %.3f-%.3f) [%s]', current_prefix, current_clim(1), current_clim(2), q_str), 'Visible', 'off');
                newMtexFigure('layout', [1, num_subplots], 'figSize', 'huge');
                
                for t_i = 1:num_subplots
                    nextAxis;
                    plotSSLIP_Context(opt_ctx);
                    
                    hold on;
                    plot(ebsdID, combined_activities(:, t_i), 'micronbar', 'off');
                    mtexColorMap parula; try mtexColorMap viridis; catch, end;
                    setColorRange(current_clim);
                    
                    % Redraw grain boundaries ON TOP so they survive the activity field plotting!
                    hold on;
                    plot(opt_ctx.grains.boundary, 'linecolor', opt_ctx.grainLineColor, 'linewidth', opt_ctx.grainLineWidth);
                    
                    if t_i <= n_top
                        hold on;
                        rep_g1 = bundles_g1{t_i}(1);
                        rep_g2 = bundles_g2{t_i}(1);
                        
                        if plot_quivers
                            % Centered Quiver Traces (match plotSSLIP_SlipActivities aesthetics exactly)
                            quiver_scale = opt_bundle.quiver_length_factor * (max(ebsdID.x,[],'all') - min(ebsdID.x,[],'all'));
                            
                            % Quivers for G1
                            quiver(cx1, cy1, quiver_scale * sSLocal_1(rep_g1).trace.x, 0.5*quiver_scale * sSLocal_1(rep_g1).trace.y, 'color', 'r', 'linewidth', 5, 'ShowArrowHead', 'off');
                            quiver(cx1, cy1, quiver_scale * sSLocal_1(rep_g1).b.normalize.x, quiver_scale * sSLocal_1(rep_g1).b.normalize.y, 'color', 'r', 'linewidth', 2);
                            
                            % Quivers for G2
                            quiver(cx2, cy2, quiver_scale * sSLocal_2(rep_g2).trace.x, 0.5*quiver_scale * sSLocal_2(rep_g2).trace.y, 'color', 'r', 'linewidth', 5, 'ShowArrowHead', 'off');
                            quiver(cx2, cy2, quiver_scale * sSLocal_2(rep_g2).b.normalize.x, quiver_scale * sSLocal_2(rep_g2).b.normalize.y, 'color', 'r', 'linewidth', 2);
                        end
                        
                        % Annotation String for Title
                        stress = stressTensor.uniaxial(xvector);
                        sf_g1 = abs(sSLocal_1(rep_g1).SchmidFactor(stress));
                        sf_g2 = abs(sSLocal_2(rep_g2).SchmidFactor(stress));
                        
                        sys_rep_g1 = inv(ori1) * sSLocal_1(rep_g1);
                        sys_rep_g2 = inv(ori2) * sSLocal_2(rep_g2);
                        
                        n1 = round(sys_rep_g1.n); b1 = round(sys_rep_g1.b, 'uvw');
                        n2 = round(sys_rep_g2.n); b2 = round(sys_rep_g2.b, 'uvw');
                        
                        sys_str_g1 = sprintf('(%.0f %.0f %.0f)[%.0f %.0f %.0f]', n1.h, n1.k, n1.l, b1.u, b1.v, b1.w);
                        sys_str_g2 = sprintf('(%.0f %.0f %.0f)[%.0f %.0f %.0f]', n2.h, n2.k, n2.l, b2.u, b2.v, b2.w);
                        
                        % Put the clean text into the title directly (no crash-prone annotation boxes!)
                        title_str = sprintf('[Bundle %d]\nG1 (%d sys): %s | SF=%.2f\nG2 (%d sys): %s | SF=%.2f', ...
                            t_i, length(bundles_g1{t_i}), sys_str_g1, sf_g1, ...
                            length(bundles_g2{t_i}), sys_str_g2, sf_g2);
                        title(title_str, 'Interpreter', 'none');
                    else
                        title(sprintf('[Leftover Activity]\nG1: %d sys | G2: %d sys', length(leftover_g1), length(leftover_g2)));
                    end
                end
                
                mtexColorbar('title', 'Combined Slip Activity (\gamma)');
                
                % Export
                clim_str = sprintf('clim_%.3f_%.3f', current_clim(1), current_clim(2));
                clim_str = strrep(clim_str, '.', 'p');
                fig_path = fullfile(plotDir, sprintf('Combined_Bicrystal_SSLIP_%s_%s_%s.jpg', current_prefix, clim_str, q_str));
                exportScaledFigure(fig, fig_path, 'target_px', 3000, 'noDateString');
                close(fig);
                
                fprintf('Exported plot to: %s\n', fig_path);
            end
        end
    end

    %% 5. Single Grain Clustered Analysis (Other grain remains in gray context)
    grain_plot_configs = {};
    if opt_plot.plot_grain1
        grain_plot_configs{end+1} = struct('name', 'Grain1', 'id', g1_id, 'mask', mask_g1, 'bundles', {bundles_g1}, 'leftovers', leftover_g1, 'sSLocal', sSLocal_1, 'ori', ori1, 'cx', cx1, 'cy', cy1);
    end
    if opt_plot.plot_grain2
        grain_plot_configs{end+1} = struct('name', 'Grain2', 'id', g2_id, 'mask', mask_g2, 'bundles', {bundles_g2}, 'leftovers', leftover_g2, 'sSLocal', sSLocal_2, 'ori', ori2, 'cx', cx2, 'cy', cy2);
    end

    for g_idx = 1:length(grain_plot_configs)
        cfg = grain_plot_configs{g_idx};
        n_bundles_g = length(cfg.bundles);
        num_subplots_g = n_bundles_g + 1;

        % Build single-grain activity fields (other grain pixels remain NaN -> gray background)
        grain_activities = nan(length(ebsdID), num_subplots_g);
        for t_i = 1:num_subplots_g
            if t_i <= n_bundles_g
                act_g = zeros(length(ebsdID), 1);
                for mem = 1:length(cfg.bundles{t_i})
                    act_g = act_g + abs(ebsdID.prop.slipIDcor(:, cfg.bundles{t_i}(mem)));
                end
            else
                act_g = zeros(length(ebsdID), 1);
                for mem = 1:length(cfg.leftovers)
                    act_g = act_g + abs(ebsdID.prop.slipIDcor(:, cfg.leftovers(mem)));
                end
            end
            grain_activities(cfg.mask(:), t_i) = act_g(cfg.mask(:));
        end

        for c_idx = 1:length(clims_to_plot)
            current_clim = clims_to_plot{c_idx};

            for plot_quivers = quiver_modes
                if plot_quivers
                    q_str = 'withQuivers';
                else
                    q_str = 'noQuivers';
                end

                fig = figure('Name', sprintf('SSLIP %s Analysis - %s (clim: %.3f-%.3f) [%s]', cfg.name, current_prefix, current_clim(1), current_clim(2), q_str), 'Visible', 'off');
                newMtexFigure('layout', [1, num_subplots_g], 'figSize', 'huge');

                for t_i = 1:num_subplots_g
                    nextAxis;
                    plotSSLIP_Context(opt_ctx);

                    hold on;
                    plot(ebsdID, grain_activities(:, t_i), 'micronbar', 'off');
                    mtexColorMap parula; try mtexColorMap viridis; catch, end;
                    setColorRange(current_clim);

                    % Redraw grain boundaries ON TOP
                    hold on;
                    plot(opt_ctx.grains.boundary, 'linecolor', opt_ctx.grainLineColor, 'linewidth', opt_ctx.grainLineWidth);

                    if t_i <= n_bundles_g
                        rep_g = cfg.bundles{t_i}(1);

                        if plot_quivers
                            quiver_scale = opt_bundle.quiver_length_factor * (max(ebsdID.x,[],'all') - min(ebsdID.x,[],'all'));
                            quiver(cfg.cx, cfg.cy, quiver_scale * cfg.sSLocal(rep_g).trace.x, 0.5*quiver_scale * cfg.sSLocal(rep_g).trace.y, 'color', 'r', 'linewidth', 5, 'ShowArrowHead', 'off');
                            quiver(cfg.cx, cfg.cy, quiver_scale * cfg.sSLocal(rep_g).b.normalize.x, quiver_scale * cfg.sSLocal(rep_g).b.normalize.y, 'color', 'r', 'linewidth', 2);
                        end

                        stress = stressTensor.uniaxial(xvector);
                        sf_g = abs(cfg.sSLocal(rep_g).SchmidFactor(stress));

                        sys_rep_g = inv(cfg.ori) * cfg.sSLocal(rep_g);
                        ng = round(sys_rep_g.n); bg = round(sys_rep_g.b, 'uvw');
                        sys_str_g = sprintf('(%.0f %.0f %.0f)[%.0f %.0f %.0f]', ng.h, ng.k, ng.l, bg.u, bg.v, bg.w);

                        title_str = sprintf('[Bundle %d]\n%s (%d sys)\n%s | SF=%.2f', ...
                            t_i, cfg.name, length(cfg.bundles{t_i}), sys_str_g, sf_g);
                        title(title_str, 'Interpreter', 'none');
                    else
                        title(sprintf('[Leftover Activity]\n%s: %d sys', cfg.name, length(cfg.leftovers)));
                    end
                end

                mtexColorbar('title', 'Combined Slip Activity (\gamma)');

                % Export
                clim_str = sprintf('clim_%.3f_%.3f', current_clim(1), current_clim(2));
                clim_str = strrep(clim_str, '.', 'p');
                fig_path = fullfile(plotDir, sprintf('SingleGrain_%s_SSLIP_%s_%s_%s.jpg', cfg.name, current_prefix, clim_str, q_str));
                exportScaledFigure(fig, fig_path, 'target_px', 3000, 'noDateString');
                close(fig);

                fprintf('Exported single grain plot to: %s\n', fig_path);
            end
        end
    end
end

fprintf('\nAll specified increments processed successfully.\n');

%% Helper Function
function reportData = buildReportData(bundles, leftover_idx, metric_vals, sSLocal, ori, n_top)
    reportData = struct('BundleID', {}, 'IsLeftover', {}, 'SystemID', {}, 'Activity', {}, 'SchmidFactor', {}, 'H11', {}, 'H12', {}, 'H21', {}, 'H22', {}, 'IsRepresentative', {}, 'n_x', {}, 'n_y', {}, 'n_z', {}, 'b_x', {}, 'b_y', {}, 'b_z', {});
    
    rep_sys_to_write = [];
    rep_b_ids = string.empty;
    rep_is_leftover = logical.empty;
    rep_is_rep = logical.empty;
    
    for t_i = 1:length(bundles)
        for mem = 1:length(bundles{t_i})
            rep_sys_to_write(end+1) = bundles{t_i}(mem);
            rep_b_ids(end+1) = string(t_i);
            rep_is_leftover(end+1) = (t_i > n_top);
            rep_is_rep(end+1) = (mem == 1);
        end
    end
    for i = 1:length(leftover_idx)
        rep_sys_to_write(end+1) = leftover_idx(i);
        rep_b_ids(end+1) = "Leftover";
        rep_is_leftover(end+1) = true;
        rep_is_rep(end+1) = false;
    end
    
    % Get macroscopic stress tensor for SF calculation (assume uniaxial X for now)
    stress = stressTensor.uniaxial(xvector);
    
    for idx = 1:length(rep_sys_to_write)
        sys = rep_sys_to_write(idx);
        if sys <= length(sSLocal)
            sf = abs(sSLocal(sys).SchmidFactor(stress));
            % Recover crystal coordinates using orientation
            sys_rep = inv(ori) * sSLocal(sys);
            sys_rep_n = round(sys_rep.n);
            sys_rep_b = round(sys_rep.b, 'uvw');
            
            nx = sys_rep_n.h; ny = sys_rep_n.k; nz = sys_rep_n.l;
            bx = sys_rep_b.u; by = sys_rep_b.v; bz = sys_rep_b.w;
            
            H = sSLocal(sys).deformationTensor.matrix;
            h11 = H(1,1); h12 = H(1,2); h21 = H(2,1); h22 = H(2,2);
        else
            sf = NaN;
            nx = NaN; ny = NaN; nz = NaN; bx = NaN; by = NaN; bz = NaN;
            h11 = NaN; h12 = NaN; h21 = NaN; h22 = NaN;
        end
        
        reportData(idx).BundleID = rep_b_ids(idx);
        reportData(idx).IsLeftover = rep_is_leftover(idx);
        reportData(idx).SystemID = sys;
        reportData(idx).Activity = metric_vals(sys);
        reportData(idx).SchmidFactor = sf;
        reportData(idx).H11 = h11;
        reportData(idx).H12 = h12;
        reportData(idx).H21 = h21;
        reportData(idx).H22 = h22;
        reportData(idx).IsRepresentative = rep_is_rep(idx);
        reportData(idx).n_x = nx;
        reportData(idx).n_y = ny;
        reportData(idx).n_z = nz;
        reportData(idx).b_x = bx;
        reportData(idx).b_y = by;
        reportData(idx).b_z = bz;
    end
end
