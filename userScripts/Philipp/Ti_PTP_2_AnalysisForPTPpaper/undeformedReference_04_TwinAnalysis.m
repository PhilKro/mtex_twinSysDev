%% undeformedReference_04_TwinAnalysis.m
% Performs DIC-based analysis of {10-12} deformation twinning variants in HCP Titanium.
% Integrates crystallographic orientation identification with 2D DIC Green strain tensors
% to analyze twin activity, variant identification, and shear projection factors across loading steps.
% Outputs results to: Processed_Output_Undeformed/Twin_Analysis_Results.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
if ~isempty(thisDir), addpath(thisDir); end
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\geometry');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 0. User Configuration & Execution Toggles
cfg = struct();

% --- Increment Processing Selection ---
% Specify indices of DIC increments to process (e.g., [4] or [1, 3, 4]).
% Set to empty [] to automatically process ALL discovered increments (1:length(prefix_list)).
cfg.process_increments        = [4]; 

% --- Analysis & Plotting Mode Toggles ---
cfg.enable_gradient_maps      = false;  % Plot 1: Full-map experimental displacement gradients (du/dx, du/dy, dv/dx, dv/dy)
cfg.enable_residual_norm_diff = false; % Plot 2: Residual norm difference maps to each primary variant
cfg.enable_best_fit_primary   = true; % Plot 3: Best Fit Variant from direct norm difference (BestFitPrimary_FromDispGradient)
cfg.best_fit_primary_thresh   = 0.2;  % Maximum residual norm difference threshold
cfg.best_fit_primary_mag_cutoff = 0.4; % Minimum projection factor magnitude (x_proj) below which variant is not plotted
cfg.enable_shape_mismatch     = false;  % Plot 4: Active variants from directional shape-mismatch minimization (Shape_Mismatch_Variants)
cfg.enable_projection_factor  = false;  % Plot 5: Scalar projection factor x map (Scalar_Projection_Factor)
cfg.enable_all_variants_breakdown = false; % Plot 6 & 7: Multi-panel comparison of shape mismatch and projection factor for all individual variants (1 to N)
cfg.enable_var3_vs_var6_comparison = false; % Plot 8: Direct comparison between Variant 3 and Variant 6

% --- Selection Metric Weights ---
cfg.weight_series_enable = false;          % Generate a series of active variant maps for varying magnitude weights
cfg.weight_shape = 1.0; % Weight for shape mismatch cost (minimize shape mismatch)
cfg.weight_mag   = 0.5; % Weight for projection factor deviation cost (minimize |x - 1.0|)
cfg.weight_series_mag_weights = 0.2:0.2:1.0; % The magnitude weights to loop over (shape weight stays at 1.0)

% --- Tilt Sensitivity Series Configuration ---
cfg.tilt_series.enable        = false;    % Toggle tilt sensitivity series analysis ON/OFF
cfg.tilt_series.range_deg     = 0:1:1;  % Range and step size for specimen tilt perturbation (degrees)
cfg.tilt_series.magCutoff     = 0.05;    % Minimum projection factor x required for identifying a valid twin pixel under tilt
cfg.tilt_series.shapeCutoff   = 0.50;    % Maximum allowable shape mismatch required for identifying a valid twin pixel under tilt
cfg.tilt_series.vector        = yvector; % Rotation axis vector for specimen tilt perturbation (e.g., xvector, yvector, zvector)
cfg.tilt_series.axis_label    = 'Yaxis'; % Identifier string embedded into output folders and PNG filenames (e.g., 'Xaxis', 'Yaxis')

%% 1. Load Registry & Locate Undeformed Reference Dataset
[registry, defaults] = a00_Dataset_Registry();
inDir  = fullfile(defaults.dir, defaults.outDir_undeformed);
outDir = fullfile(defaults.dir, defaults.outDir_undeformed, 'Twin_Analysis_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Find the latest _undeformedEBSD_withAllDIC.mat file
matFiles = dir(fullfile(inDir, '*_undeformedEBSD_withAllDIC.mat'));
if isempty(matFiles)
    error('No aligned reference .mat file found in %s. Run undeformedReference_03_IMPORTER first.', inDir);
end
[~, idx] = max([matFiles.datenum]);
matname = matFiles(idx).name;
fprintf('Loading undeformed reference dataset:\n  %s\n', fullfile(inDir, matname));
dataLoaded = load(fullfile(inDir, matname));
newEBSD = dataLoaded.newEBSD;
if isfield(dataLoaded, 'options')
    options = dataLoaded.options;
else
    options = struct();
end

%% 2. Retrieve or Reconstruct Grains
grainThresh = 5*degree; 
minSize = 300; % Threshold in pixels
alpha = 1;

if isfield(dataLoaded, 'grains')
    grains = dataLoaded.grains;
    fprintf('  -> Successfully loaded precomputed, FZ-projected grains directly from reference dataset.\n');
else
    fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
    thresh_ci = defaults.thresh_ci;
    newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
    [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
end

CS = newEBSD.CS;

%% 3. Identify Parent Grain and Define {10-12} Twin Systems for Titanium
fprintf('Identifying parent grain matrix and initializing {10-12} twin systems...\n');

% Identify the largest grain as the primary parent matrix for twin analysis
[~, parentId] = max(grains.numPixel);
parentGrain = grains(parentId);
tolerance = 8 * degree;
parentGrains = grains(angle(parentGrain.meanOrientation, grains.meanOrientation) < tolerance);
fprintf('  -> Identified Parent Matrix: Grain %d (%d px)\n', parentGrain.id, parentGrain.numPixel);

% Create twinSystem object for {10-12} tensile twinning in Titanium
tS_1012 = twinSystem.hexagonal_1012(CS);

% Symmetrise to obtain all unique twin variants (6 primary variants for {10-12})
tS_sym = tS_1012.symmetrise;

figure;
plot(grains);
hold on
plot(parentGrain, 'FaceColor','r')
text(parentGrain, 'Parent Grain')
close all


%% 4. Theoretical Twin Displacement Gradient Tensors & Similarity Check
fprintf('Computing Theoretical Twin Displacement Gradient Tensors & running Similarity Check...\n');
pVars = 1:numel(tS_sym);
N = numel(pVars);
colors_p = lines(N);

H = parentGrain.meanOrientation * tS_sym.displacementGradient;
H_id = H(pVars);
H_id2d = H_id.M(1:2,1:2,:);
Hv_flat = reshape(H_id2d, 4, N);

% Calculate K1 shear plane traces for plotting
pT_k1 = parentGrain.meanOrientation .* tS_sym.k1;
pT_k1_trace = cross(pT_k1, zvector);
pT_k1_trace.antipodal = true;

% Perform similarity check on 2D projections of theoretical twin displacement gradient tensors
H_check2d = reshape(H_id.M(1:2, 1:2, :), 4, []);
H_check2d_norm = H_check2d ./ vecnorm(H_check2d, 2, 1);
simMatrix = H_check2d_norm' * H_check2d_norm;
simThresh = 0.95;
[simRow, simCol] = find(abs(triu(simMatrix, 1)) > simThresh);

if ~isempty(simRow)
    warning('Very similar theoretical twin displacement gradient tensors detected (similarity > %.2f)! This might cause variant misidentification.', simThresh);
    for k = 1:length(simRow)
        fprintf('  -> Variant %d and Variant %d: similarity = %.3f\n', pVars(simRow(k)), pVars(simCol(k)), simMatrix(simRow(k), simCol(k)));
    end
else
    fprintf('  -> All theoretical twin displacement gradient tensors are sufficiently distinct in 2D shape (similarity <= %.2f).\n', simThresh);
end

%% 5. Loop over DIC Increments for Displacement-Gradient-Based Variant Analysis
all_fields = fields(newEBSD.prop);
displacement_fields = all_fields(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true));
prefix_list = unique(cellfun(@(x) x(1:find(x=='_', 1)-1), displacement_fields, 'UniformOutput', false));
fprintf('\nFound %d DIC increments to evaluate for twin deformation mechanics.\n', length(prefix_list));

if isempty(cfg.process_increments)
    inc_indices = 1:length(prefix_list);
else
    inc_indices = intersect(cfg.process_increments, 1:length(prefix_list));
end
fprintf('Evaluating %d selected increment(s): %s\n', length(inc_indices), mat2str(inc_indices));

for idx_inc = inc_indices
    prefix = prefix_list{idx_inc};
    fprintf('\n========================================================================\n');
    fprintf('Evaluating Twinning Mechanics for Increment [%s] (%d/%d)\n', prefix, idx_inc, length(prefix_list));
    fprintf('========================================================================\n');
    
    % --- Retrieve or Calculate Experimental 2D Displacement Gradient Tensor (H_exp) ---
    f_dudx = find(startsWith(all_fields, prefix, 'IgnoreCase', true) & contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    f_dudy = find(startsWith(all_fields, prefix, 'IgnoreCase', true) & contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    f_dvdx = find(startsWith(all_fields, prefix, 'IgnoreCase', true) & contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    f_dvdy = find(startsWith(all_fields, prefix, 'IgnoreCase', true) & contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    
    if ~isempty(f_dudx) && ~isempty(f_dudy) && ~isempty(f_dvdx) && ~isempty(f_dvdy)
        dudx_data = newEBSD.prop.(all_fields{f_dudx});
        dudy_data = newEBSD.prop.(all_fields{f_dudy});
        dvdx_data = newEBSD.prop.(all_fields{f_dvdx});
        dvdy_data = newEBSD.prop.(all_fields{f_dvdy});
        fprintf('  -> Using pre-computed DIC displacement gradient fields: %s, %s, %s, %s\n', all_fields{f_dudx}, all_fields{f_dudy}, all_fields{f_dvdx}, all_fields{f_dvdy});
    else
        error('  -> Displacement gradient fields not found directly. Computing gradients via Gaussian smoothed displacements...\n');
        ux_fname = [prefix, '_Displacements_x_displacement'];
        uy_fname = [prefix, '_Displacements_y_displacement'];
        if ~isfield(newEBSD.prop, ux_fname) || ~isfield(newEBSD.prop, uy_fname)
            warning('Skipping increment %s: Missing displacement fields.', prefix);
            continue;
        end
        ebsd_grid = newEBSD.gridify;
        grid_size = size(ebsd_grid);
        
        % Replace 0 correlation / unindexed edge dropouts with NaN
        U_raw = double(ebsd_grid.prop.(ux_fname)); U_raw(U_raw == 0) = NaN;
        V_raw = double(ebsd_grid.prop.(uy_fname)); V_raw(V_raw == 0) = NaN;
        U_grid = reshape(U_raw * ebsd_grid.dPos, grid_size);
        V_grid = reshape(V_raw * ebsd_grid.dPos, grid_size);
        
        % Apply SSLIP Gaussian filter (filt_std = 4) to smooth away DIC stepsize / kron staircase edges
        filt_opt.filt_std = 4;
        try
            sm_data = filterDisplacements(U_grid, V_grid, filt_opt);
            U_grid = sm_data.U;
            V_grid = sm_data.V;
            fprintf('      Applied Gaussian displacement filter (filt_std = 4) before gradient calculation to eliminate DIC stepsize edges.\n');
        catch
            warning('filterDisplacements helper function not found in path. Falling back to raw displacements.');
        end
        
        [Hxx_grid, Hxy_grid] = gradient(U_grid, ebsd_grid.dPos, ebsd_grid.dPos);
        [Hyx_grid, Hyy_grid] = gradient(V_grid, ebsd_grid.dPos, ebsd_grid.dPos);
        
        grid_ids  = ebsd_grid.id(:);
        valid_idx = ~isnan(grid_ids) & grid_ids > 0;
        valid_grid_ids = grid_ids(valid_idx);
        [~, loc] = ismember(newEBSD.id, valid_grid_ids);
        
        Hxx_flat = Hxx_grid(:); dudx_data = Hxx_flat(valid_idx); dudx_data = dudx_data(loc);
        Hxy_flat = Hxy_grid(:); dudy_data = Hxy_flat(valid_idx); dudy_data = dudy_data(loc);
        Hyx_flat = Hyx_grid(:); dvdx_data = Hyx_flat(valid_idx); dvdx_data = dvdx_data(loc);
        Hyy_flat = Hyy_grid(:); dvdy_data = Hyy_flat(valid_idx); dvdy_data = dvdy_data(loc);
    end
    
    nPoints = numel(dudx_data);
    H_exp = zeros(2, 2, nPoints);
    H_exp(1,1,:) = dudx_data(:);
    H_exp(1,2,:) = dudy_data(:);
    H_exp(2,1,:) = dvdx_data(:);
    H_exp(2,2,:) = dvdy_data(:);
    
    % Identify valid shear domains and filter out severe DIC decorrelation / edge anomalies
    exx_approx = dudx_data(:);
    eyy_approx = dvdy_data(:);
    exy_approx = 0.5 .* (dudy_data(:) + dvdx_data(:));
    geff = sqrt(((exx_approx - eyy_approx)./2).^2 + exy_approx.^2);
    nanField = geff > 2.0 | geff == 0 | isnan(geff) | isnan(dudx_data(:));
    
    % --- Plot 1: Full-Map Displacement Gradient Fields with Grain Boundaries ---
    if cfg.enable_gradient_maps
        cl = [-0.15 0.15];
        figGrad = figure('Visible', 'off');
        newMtexFigure('layout', [2, 2]);
        plot(newEBSD, squeeze(H_exp(1,1,:))); caxis(cl); title(sprintf('[%s] du/dx (H_{11})', prefix));
        hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        
        nextAxis;
        plot(newEBSD, squeeze(H_exp(1,2,:))); caxis(cl); title(sprintf('[%s] du/dy (H_{12})', prefix));
        hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        
        nextAxis;
        plot(newEBSD, squeeze(H_exp(2,1,:))); caxis(cl); title(sprintf('[%s] dv/dx (H_{21})', prefix));
        hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);

        nextAxis;
        plot(newEBSD, squeeze(H_exp(2,2,:))); caxis(cl); title(sprintf('[%s] dv/dy (H_{22})', prefix));
        hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        colormap viridis; mtexColorbar; hold off;
        exportScaledFigure(figGrad, fullfile(outDir, [prefix, '_Displacement_Gradient_Fields.png']), 'nodatestring');
        close(figGrad);
    end
    
    % --- Plot 2: Residual Norm Difference to Each Primary Variant ---
    NormDiff = zeros(nPoints, N);
    for i_v = 1:N
        DefDiff_v = H_exp - H_id2d(:,:,i_v);
        NormDiff(:, i_v) = vecnorm(reshape(DefDiff_v, 4, nPoints), 2, 1);
        NormDiff(nanField(:), i_v) = NaN;
    end
    
    if cfg.enable_residual_norm_diff
        figNorm = figure('Visible', 'off');
        newMtexFigure('layout', [1, N]);
        for i_v = 1:N
            nextAxis;
            plot(newEBSD, NormDiff(:, i_v)); caxis([0 0.8]);
            title(sprintf('Diff to PV%d', pVars(i_v))); hold on;
            quiver(parentGrain, pT_k1_trace(pVars(i_v)), 'color', 'w');
        end
        colormap viridis; mtexColorbar; hold off;
        exportScaledFigure(figNorm, fullfile(outDir, [prefix, '_NormDiff_to_PVs.png']), 'nodatestring');
        close(figNorm);
    end
    
    % --- Map-Wide Shape Mismatch & Scalar Projection Factor ---
    He_flat = reshape(H_exp, 4, nPoints);
    shapeMismatch = zeros(N, nPoints);
    magnitudeX = zeros(N, nPoints);
    for v = 1:N
        Hv = Hv_flat(:, v);
        x_proj = (Hv' * He_flat) ./ (Hv' * Hv);
        magnitudeX(v, :) = x_proj;
        Mismatch = He_flat - (Hv * x_proj);
        shapeMismatch(v, :) = vecnorm(Mismatch, 2, 1);
    end
    shapeMismatch(:, nanField) = NaN;
    magnitudeX(:, nanField) = NaN;
    
    % --- Plot 3: Best Fit Primary Variant from Displacement Gradient Difference ---
    if cfg.enable_best_fit_primary
        threshP = 0.2;
        if isfield(cfg, 'best_fit_primary_thresh'), threshP = cfg.best_fit_primary_thresh; end
        magCutoffP = 0.05;
        if isfield(cfg, 'best_fit_primary_mag_cutoff'), magCutoffP = cfg.best_fit_primary_mag_cutoff; end
        
        [fitPrimary, bestFitPrimary] = min(NormDiff, [], 2);
        linIdx_p = sub2ind(size(magnitudeX), bestFitPrimary(:)', 1:nPoints);
        bestFitMagX = magnitudeX(linIdx_p)';
        
        validPrimary = (fitPrimary < threshP) & (bestFitMagX > magCutoffP) & ~isnan(bestFitMagX);
        
        plotColorGrad = ones(length(bestFitPrimary), 3);
        plotColorGrad(validPrimary, :) = colors_p(bestFitPrimary(validPrimary), :);
        
        figBestFit = figure('Visible', 'off');
        plot(newEBSD, plotColorGrad); colormap(colors_p); caxis([1 N]); mtexColorbar; hold on;
        plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 2);
        if magCutoffP > 0
            title(sprintf('[%s] Best Fit Variant from Disp Gradient (Res < %.2f, ProjFactor > %.2f)', prefix, threshP, magCutoffP), 'Interpreter', 'none');
        else
            title(sprintf('[%s] Best Fit Variant from Disp Gradient (Res < %.2f)', prefix, threshP), 'Interpreter', 'none');
        end
        hold off;
        exportScaledFigure(figBestFit, fullfile(outDir, [prefix, '_BestFitPrimary_FromDispGradient.png']), 'nodatestring');
        close(figBestFit);
    end
    
    % --- Plot 4 & 5: Advanced Map-Wide Shape Mismatch & Scalar Projection Factor ---
    
    % Evaluate combined cost metric
    weight_s = 1.0; if isfield(cfg, 'weight_shape'), weight_s = cfg.weight_shape; end
    weight_m = 0.0; if isfield(cfg, 'weight_mag'),   weight_m = cfg.weight_mag;   end
    
    combinedCost = (weight_s .* shapeMismatch) + (weight_m .* abs(magnitudeX - 1.0));
    
    [minCost, bestVariantIdx] = min(combinedCost, [], 1);
    
    linearIndices = sub2ind(size(magnitudeX), bestVariantIdx, 1:nPoints);
    bestMagnitudeX = magnitudeX(linearIndices);
    minMismatch = shapeMismatch(linearIndices); % Ensure minMismatch holds the actual mismatch of the chosen variant, not the cost
    
    if cfg.enable_shape_mismatch
        magCutoff = 0.05; 
        shapeCutoff = 0.5;
        isValid = (bestMagnitudeX > magCutoff) & (minMismatch < shapeCutoff);
        validIndices = find(isValid & ~isnan(minMismatch));
        
        plotColorNoBackdrop = ones(nPoints, 3);
        plotColorNoBackdrop(validIndices, :) = colors_p(bestVariantIdx(validIndices), :);
        
        figShape = figure('Visible', 'off');
        plot(newEBSD, plotColorNoBackdrop); colormap(colors_p); caxis([1 N]); mtexColorbar; hold on;
        plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
        title(sprintf('[%s] Active Twinning Variants (Shape Mismatch < %.2f)', prefix, shapeCutoff), 'Interpreter', 'none'); hold off;
        exportScaledFigure(figShape, fullfile(outDir, [prefix, '_Shape_Mismatch_Variants.png']), 'nodatestring');
        close(figShape);

        % Continuous Shape Mismatch Map
        shape_cRange = [0, 0.5];
        plotMism = minMismatch;
        plotMism(isnan(plotMism)) = shape_cRange(1);
        
        figMism = figure('Visible', 'off');
        plot(newEBSD, plotMism); colormap viridis; setColorRange(shape_cRange); mtexColorbar; hold on;
        plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
        title(sprintf('[%s] Best Fit Shape Mismatch Map', prefix), 'Interpreter', 'none'); hold off;
        exportScaledFigure(figMism, fullfile(outDir, [prefix, '_Shape_Mismatch_Map.png']), 'nodatestring');
        close(figMism);
    end
    
    % --- Optional: Weight Series Generation ---
    if isfield(cfg, 'weight_series_enable') && cfg.weight_series_enable
        w_shape = 1.0;
        magCutoff_s = 0.05; 
        shapeCutoff_s = 0.5;
        
        for w_mag = cfg.weight_series_mag_weights
            tempCost = (w_shape .* shapeMismatch) + (w_mag .* abs(magnitudeX - 1.0));
            [~, tempBestVarIdx] = min(tempCost, [], 1);
            
            tempLinIdx = sub2ind(size(magnitudeX), tempBestVarIdx, 1:nPoints);
            tempMagX = magnitudeX(tempLinIdx);
            tempMis = shapeMismatch(tempLinIdx);
            
            tempValid = (tempMagX > magCutoff_s) & (tempMis < shapeCutoff_s);
            tempValidIdx = find(tempValid & ~isnan(tempMis));
            
            tempColorMap = ones(nPoints, 3);
            tempColorMap(tempValidIdx, :) = colors_p(tempBestVarIdx(tempValidIdx), :);
            
            figSeries = figure('Visible', 'off');
            plot(newEBSD, tempColorMap); colormap(colors_p); caxis([1 N]); mtexColorbar; hold on;
            plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
            title(sprintf('[%s] Active Variants (Shape W:%.1f, Mag W:%.1f)', prefix, w_shape, w_mag), 'Interpreter', 'none'); hold off;
            exportScaledFigure(figSeries, fullfile(outDir, sprintf('%s_Variants_Shape%.1f_Mag%.1f.jpg', prefix, w_shape, w_mag)), 'nodatestring');
            close(figSeries);
        end
    end
    
    % --- Plot 5: Scalar Projection Factor x ---
    if cfg.enable_projection_factor
        proj_cRange = [0.05, 1.25];
        plotProj = bestMagnitudeX;
        plotProj(isnan(plotProj)) = proj_cRange(1);
        
        figProj = figure('Visible', 'off');
        plot(newEBSD, plotProj); colormap viridis; setColorRange(proj_cRange); mtexColorbar; hold on;
        plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
        title(sprintf('[%s] Scalar Projection Factor x (1.0 = Theoretical Shear)', prefix), 'Interpreter', 'none'); hold off;
        exportScaledFigure(figProj, fullfile(outDir, [prefix, '_Scalar_Projection_Factor.png']), 'nodatestring');
        close(figProj);
    end
    
    % --- Plot 6 & 7: Per-Variant Breakdown of Shape Mismatch & Projection Factor ---
    if isfield(cfg, 'enable_all_variants_breakdown') && cfg.enable_all_variants_breakdown
        nCols = min(3, N);
        nRows = ceil(N / nCols);
        
        % 1. Shape Mismatch across all individual variants
        shape_cRange = [0, 0.5];
        figAllMism = figure('Name', 'Shape Mismatch All Variants', 'Visible', 'off');
        newMtexFigure('layout', [nRows, nCols]);
        for v = 1:N
            if v > 1, nextAxis; end
            plotMism_v = shapeMismatch(v, :);
            plotMism_v(isnan(plotMism_v)) = shape_cRange(1);
            plot(newEBSD, plotMism_v); hold on;
            plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
            title(sprintf('Var %d Shape Mismatch', v));
            colormap(gca, 'viridis'); setColorRange(shape_cRange); hold off;
        end
        mtexColorbar;
        sgtitle(sprintf('[%s] Shape Mismatch across All Variants (Lower = Better Fit)', prefix), 'Interpreter', 'none');
        exportScaledFigure(figAllMism, fullfile(outDir, [prefix, '_Shape_Mismatch_All_Variants.png']), 'nodatestring');
        close(figAllMism);
        
        % 2. Scalar Projection Factor x across all individual variants
        proj_cRange = [0.05, 1.25];
        figAllMag = figure('Name', 'Scalar Projection Factor All Variants', 'Visible', 'off');
        newMtexFigure('layout', [nRows, nCols]);
        for v = 1:N
            if v > 1, nextAxis; end
            plotMag_v = magnitudeX(v, :);
            plotMag_v(isnan(plotMag_v)) = proj_cRange(1);
            plot(newEBSD, plotMag_v); hold on;
            plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
            title(sprintf('Var %d Projection Factor x', v));
            colormap(gca, 'viridis'); setColorRange(proj_cRange); hold off;
        end
        mtexColorbar;
        sgtitle(sprintf('[%s] Scalar Projection Factor across All Variants (1.0 = Theo Shear)', prefix), 'Interpreter', 'none');
        exportScaledFigure(figAllMag, fullfile(outDir, [prefix, '_Scalar_Projection_All_Variants.png']), 'nodatestring');
        close(figAllMag);
    end
    
    % --- Plot 8: Direct comparison of Variant 3 vs Variant 6 ---
    if isfield(cfg, 'enable_var3_vs_var6_comparison') && cfg.enable_var3_vs_var6_comparison && N >= 6
        vA = 3;
        vB = 6;
        
        diff_shape = shapeMismatch(vB, :) - shapeMismatch(vA, :); % < 0 means vB (6) has lower (better) shape mismatch
        diff_mag   = magnitudeX(vB, :) - magnitudeX(vA, :);
        
        shape_cRange = [0, 0.5];
        proj_cRange  = [0.05, 1.25];

        % Row 1: Shape Mismatch
        plotMism_vA = shapeMismatch(vA, :);
        plotMism_vA(isnan(plotMism_vA)) = shape_cRange(1);
        figure; plot(newEBSD, plotMism_vA); colormap(gca, 'viridis'); setColorRange(shape_cRange); mtexColorbar;
        title(sprintf('Var %d Shape Mismatch', vA)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ShapeMismatch', '_Var', num2str(vA), '.jpg']), 'nodatestring'); close(gcf);   

        plotMism_vB = shapeMismatch(vB, :);
        plotMism_vB(isnan(plotMism_vB)) = shape_cRange(1);
        figure; plot(newEBSD, plotMism_vB); colormap(gca, 'viridis'); setColorRange(shape_cRange); mtexColorbar;
        title(sprintf('Var %d Shape Mismatch', vB)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ShapeMismatch', '_Var', num2str(vB), '.jpg']), 'nodatestring'); close(gcf); 

        figure; plot(newEBSD, diff_shape); colormap(gca, 'parula'); 
        maxDiffS = max(abs(prctile(diff_shape, 2)), abs(prctile(diff_shape, 98)));
        if isnan(maxDiffS) || maxDiffS == 0, maxDiffS = 0.1; end
        setColorRange([-maxDiffS, maxDiffS]); mtexColorbar;
        title(sprintf('Diff (Var %d - Var %d)', vB, vA)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ShapeMismatch', '_Var', num2str(vB), '-Var', num2str(vA),'.jpg']), 'nodatestring'); close(gcf);

        % Row 2: Projection Factor
        plotMag_vA = magnitudeX(vA, :);
        plotMag_vA(isnan(plotMag_vA)) = proj_cRange(1);
        figure; plot(newEBSD, plotMag_vA); colormap(gca, 'viridis'); setColorRange(proj_cRange); mtexColorbar;
        title(sprintf('Var %d Proj Factor', vA)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ProjectionFactor', '_Var', num2str(vA), '.jpg']), 'nodatestring'); close(gcf);
        
        plotMag_vB = magnitudeX(vB, :);
        plotMag_vB(isnan(plotMag_vB)) = proj_cRange(1);
        figure; plot(newEBSD, plotMag_vB); colormap(gca, 'viridis'); setColorRange(proj_cRange); mtexColorbar;
        title(sprintf('Var %d Proj Factor', vB)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ProjectionFactor', '_Var', num2str(vB), '.jpg']), 'nodatestring'); close(gcf);
        
        figure; plot(newEBSD, diff_mag); colormap(gca, 'parula');
        maxDiffM = max(abs(prctile(diff_mag, 2)), abs(prctile(diff_mag, 98)));
        if isnan(maxDiffM) || maxDiffM == 0, maxDiffM = 0.5; end
        setColorRange([-maxDiffM, maxDiffM]); mtexColorbar;
        title(sprintf('Diff (Var %d - Var %d)', vB, vA)); hold on; plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1);
        exportScaledFigure(gcf, fullfile(outDir, [prefix, '_ProjectionFactor', '_Var', num2str(vB), '-Var', num2str(vA), '.jpg']), 'nodatestring'); close(gcf);



        
        % Also plot a scatter plot of Shape Mismatch vs Magnitude for the two variants
        % figScatter = figure('Name', 'Var 3 vs Var 6 Scatter', 'Visible', 'off');
        % valid_pts = ~isnan(shapeMismatch(vA, :)) & ~isnan(shapeMismatch(vB, :));
        % mag_A = magnitudeX(vA, valid_pts); mis_A = shapeMismatch(vA, valid_pts);
        % mag_B = magnitudeX(vB, valid_pts); mis_B = shapeMismatch(vB, valid_pts);
        % 
        % scatter(mag_A, mis_A, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        % scatter(mag_B, mis_B, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.2);
        % 
        % xline(0.05, 'k--', 'Mag Cutoff (0.05)'); % Default hardcoded cutoff in code
        % yline(0.5, 'k--', 'Shape Cutoff (0.5)'); % Default hardcoded cutoff in code
        % 
        % xlabel('Scalar Projection Factor x (1.0 = Theo Shear)');
        % ylabel('Shape Mismatch');
        % legend(sprintf('Variant %d', vA), sprintf('Variant %d', vB), 'Location', 'best');
        % title(sprintf('[%s] Scatter: Var %d vs Var %d', prefix, vA, vB), 'Interpreter', 'none');
        % grid on;
        % exportScaledFigure(figScatter, fullfile(outDir, [prefix, '_Scatter_Var', num2str(vB), '_vs_Var', num2str(vA), '.png']), 'nodatestring');
        % close(figScatter);
    end
    
    % --- Tilt Sensitivity Series Analysis ---
    if cfg.tilt_series.enable
        % Automatically resolve axis name string for labeling if not explicitly set
        vec = cfg.tilt_series.vector;
        if isfield(cfg.tilt_series, 'axis_label') && ~isempty(cfg.tilt_series.axis_label)
            axisStr = cfg.tilt_series.axis_label;
        elseif norm(vec - xvector) < 1e-5
            axisStr = 'Xaxis';
        elseif norm(vec - yvector) < 1e-5
            axisStr = 'Yaxis';
        elseif norm(vec - zvector) < 1e-5
            axisStr = 'Zaxis';
        else
            axisStr = sprintf('Vec_%.2f_%.2f_%.2f', vec.x, vec.y, vec.z);
        end
        
        fprintf('  -> Running Tilt Sensitivity Series along [%s] (%d angles: from %d° to %d° in %d° steps)...\n', ...
                axisStr, length(cfg.tilt_series.range_deg), min(cfg.tilt_series.range_deg), max(cfg.tilt_series.range_deg), ...
                cfg.tilt_series.range_deg(2)-cfg.tilt_series.range_deg(1));
        
        tiltDir = fullfile(outDir, sprintf('Tilt_Sensitivity_%s', axisStr), prefix);
        if ~exist(tiltDir, 'dir'), mkdir(tiltDir); end
        
        originalParentOri = parentGrain.meanOrientation;
        reforientation = orientation.byEuler(0*degree, 0*degree, 0*degree, CS);
        
        for tilt_deg = cfg.tilt_series.range_deg
            fprintf('     [Tilt %s %+2d°] Computing theoretical twinning shear & matching variants...\n', axisStr, tilt_deg);
            
            rot_tilt = rotation.byAxisAngle(cfg.tilt_series.vector, tilt_deg * degree);
            parentOri_tilted = rot_tilt * originalParentOri;
            % Note: Do NOT project to fundamental region here! Applying symmetry operators rotates
            % local crystal axes and silently permutes the variant numbering in tS_sym!
            % parentOri_tilted = parentOri_tilted.project2FundamentalRegion(reforientation);
            
            H_tilted = parentOri_tilted * tS_sym.displacementGradient;
            H_id_tilted = H_tilted(pVars);
            H_id2d_tilted = H_id_tilted.M(1:2, 1:2, :);
            Hv_flat_tilted = reshape(H_id2d_tilted, 4, N);
            
            shapeMismatch_tilted = zeros(N, nPoints);
            magnitudeX_tilted    = zeros(N, nPoints);
            for v = 1:N
                Ev = Hv_flat_tilted(:, v);
                x_val = (Ev' * He_flat) ./ (Ev' * Ev);
                magnitudeX_tilted(v, :) = x_val;
                Mismatch = He_flat - (Ev * x_val);
                shapeMismatch_tilted(v, :) = vecnorm(Mismatch, 2, 1);
            end
            shapeMismatch_tilted(:, nanField) = NaN;
            magnitudeX_tilted(:, nanField) = NaN;
            
            [minMismatch_t, bestVarIdx_t] = min(shapeMismatch_tilted, [], 1);
            linIdx_t = sub2ind(size(magnitudeX_tilted), bestVarIdx_t, 1:nPoints);
            bestMagX_t = magnitudeX_tilted(linIdx_t);
            
            isValid_t = (bestMagX_t > cfg.tilt_series.magCutoff) & (minMismatch_t < cfg.tilt_series.shapeCutoff);
            validIdx_t = find(isValid_t & ~isnan(minMismatch_t));
            
            % Plot 1: Active Twinning Variants under Tilt
            if cfg.enable_shape_mismatch
                plotColorTilt = ones(nPoints, 3);
                plotColorTilt(validIdx_t, :) = colors_p(bestVarIdx_t(validIdx_t), :);
                
                figT1 = figure('Visible', 'off');
                plot(newEBSD, plotColorTilt); colormap(colors_p); caxis([1 N]); mtexColorbar; hold on;
                plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
                title(sprintf('[%s] Active Variants (%s Tilt %+d°, Mismatch < %.2f)', prefix, axisStr, tilt_deg, cfg.tilt_series.shapeCutoff), 'Interpreter', 'none'); hold off;
                exportScaledFigure(figT1, fullfile(tiltDir, sprintf('%s_%s_ShapeMismatch_Variants_Tilt_%+ddeg.png', prefix, axisStr, tilt_deg)), 'nodatestring');
                close(figT1);
            end
            
            % Plot 2: Scalar Projection Factor x under Tilt
            if cfg.enable_projection_factor
                proj_cRange = [0.05, 1.25];
                bestMagX_plot = bestMagX_t;
                bestMagX_plot(bestMagX_t > 1.5 | isnan(bestMagX_plot)) = proj_cRange(1);
                
                figT2 = figure('Visible', 'off');
                plot(newEBSD, bestMagX_plot); colormap viridis; setColorRange(proj_cRange); mtexColorbar; hold on;
                plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
                title(sprintf('[%s] Scalar Projection Factor x (%s Tilt %+d°)', prefix, axisStr, tilt_deg), 'Interpreter', 'none'); hold off;
                exportScaledFigure(figT2, fullfile(tiltDir, sprintf('%s_%s_Scalar_Projection_Factor_Tilt_%+ddeg.png', prefix, axisStr, tilt_deg)), 'nodatestring');
                close(figT2);
            end
            
            % Plot 3: Shape Mismatch under Tilt
            if cfg.enable_residual_norm_diff
                shape_cRange = [0, cfg.tilt_series.shapeCutoff];
                minMismatch_plot = minMismatch_t;
                minMismatch_plot(minMismatch_t > 1 | isnan(minMismatch_plot)) = shape_cRange(1);
                
                figT3 = figure('Visible', 'off');
                plot(newEBSD, minMismatch_plot); colormap viridis; setColorRange(shape_cRange); mtexColorbar; hold on;
                plot(grains.boundary('indexed', 'notIndexed'), 'lineColor', 'k', 'linewidth', 1.5);
                title(sprintf('[%s] Shape Mismatch Map (%s Tilt %+d°)', prefix, axisStr, tilt_deg), 'Interpreter', 'none'); hold off;
                exportScaledFigure(figT3, fullfile(tiltDir, sprintf('%s_%s_ShapeMismatch_Map_Tilt_%+ddeg.png', prefix, axisStr, tilt_deg)), 'nodatestring');
                close(figT3);
            end
        end
        fprintf('  -> Successfully saved all %s tilt sensitivity maps to: %s\n', axisStr, tiltDir);
    end

    % Save summary MAT file for this step
    save(fullfile(outDir, [prefix, '_Twin_Analysis_Results.mat']), 'bestVariantIdx', 'bestMagnitudeX', 'minMismatch', 'tS_sym', 'cfg', 'shapeMismatch', 'magnitudeX');
    fprintf('  -> Completed twin mechanics evaluation for increment %s!\n', prefix);
end

fprintf('\n========================================================================\n');
fprintf('ALL TWINNING VARIANT ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
