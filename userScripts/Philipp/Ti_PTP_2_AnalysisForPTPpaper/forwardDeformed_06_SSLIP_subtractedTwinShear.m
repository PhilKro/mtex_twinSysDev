%% forwardDeformed_06_SSLIP_subtractedTwinShear.m
% Applies the SSLIP analysis method to forward deformed DIC data on TKD datasets with multiple grains,
% globally subtracting the theoretical displacement gradient of a selected twin variant (v) formed in the parent grain.
% This decouples twinning shear from secondary dislocation glide (slip) plasticity occurring in the presence of twinning.
%
% Outputs results to: Processed_Output_Deformed_Forward/SSLIP_SubtractedTwinShear_Results

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% User Configuration for Twin Subtraction
% Specify which primary {10-12} tensile twin variants to globally subtract before SSLIP optimization.
% Set to a single variant ID (e.g., target_variants = [1]) or evaluate a series (e.g., target_variants = 1:6).
target_variants = [1]; 

%% 1. Load Registry & Configure Step Selection
[registry, defaults] = a00_Dataset_Registry();

% Filter for non-excluded deformed steps with valid TKD maps and DIC groups
defSteps = registry(~[registry.isReference] & ~[registry.exclude] & ...
                    ~cellfun(@isempty, {registry.tkd_fname}) & ~cellfun(@isempty, {registry.dic_group}));

if isempty(defSteps)
    fprintf('========================================================================\n');
    fprintf('NOTICE: No valid deformed step entries found in a00_Dataset_Registry().\n');
    fprintf('========================================================================\n');
    return;
end

inDir  = fullfile(defaults.dir, defaults.outDir_deformed);
outDir = fullfile(defaults.dir, defaults.outDir_deformed, 'SSLIP_SubtractedTwinShear_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING TWIN-SHEAR-SUBTRACTED SSLIP BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

%% 2. Iterate Over Each Deformed Step
for stepIdx = 1:numel(defSteps)
    defStep = defSteps(stepIdx);
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('Processing Step [%s] (%d/%d) - DIC Group: %s\n', defStep.stepName, stepIdx, numel(defSteps), defStep.dic_group);
    fprintf('------------------------------------------------------------------------\n');
    
    matPath = fullfile(inDir, [defStep.stepName, '_deformedEBSD_withForwardDIC.mat']);
    if ~exist(matPath, 'file')
        warning('Aligned dataset not found at:\n  %s\nRun forwardDeformed_03_IMPORTER first for Step [%s]. Skipping.', matPath, defStep.stepName);
        continue;
    end
    
    fprintf('Loading aligned dataset: %s...\n', matPath);
    dataLoaded = load(matPath);
    newEBSD = dataLoaded.newEBSD;
    
    %% 3. Retrieve or Reconstruct Grains & Identify Primary Parent Grain
    grainThresh = 5*degree; 
    minSize = 300; % Threshold in pixels to isolate significant grains
    alpha = 1;

    if isfield(dataLoaded, 'grains')
        grains = dataLoaded.grains;
        fprintf('  -> Successfully loaded precomputed grains directly from .mat dataset.\n');
    else
        fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
        thresh_ci = defaults.thresh_ci;
        newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
        [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
    end

    % Identify primary Parent Grain (matrix grain with largest surface area/pixels)
    [~, parentIdx] = max(grains.numPixel);
    parentGrain = grains(parentIdx);
    fprintf('  -> Primary Parent Matrix Identified: Grain %d (%d px, (%.1f, %.1f) um).\n', ...
            parentGrain.id, parentGrain.numPixel, parentGrain.centroid(1), parentGrain.centroid(2));

    %% 4. Calculate Theoretical Twin Displacement Gradient Tensors ({10-12} Twinning)
    tS_1012 = twinSystem.hexagonal_1012(newEBSD.CS);
    tS_sym  = tS_1012.symmetrise; % All 6 primary variants
    
    % Transform theoretical twinning displacement gradients into global specimen reference frame
    H_twin_global = parentGrain.meanOrientation * tS_sym.displacementGradient;

    %% 5. Retrieve or Calculate Global Raw Displacement Gradients (H_raw)
    all_fields = fields(newEBSD.prop);
    
    % Search for precomputed gradient components first
    hxx_idx = find(contains(all_fields, {'_Hxx', '_dudx'}, 'IgnoreCase', true), 1);
    hxy_idx = find(contains(all_fields, {'_Hxy', '_dudy'}, 'IgnoreCase', true), 1);
    hyx_idx = find(contains(all_fields, {'_Hyx', '_dvdx'}, 'IgnoreCase', true), 1);
    hyy_idx = find(contains(all_fields, {'_Hyy', '_dvdy'}, 'IgnoreCase', true), 1);
    
    ux_idx = find(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true), 1);
    uy_idx = find(contains(all_fields, '_Displacements_y_displacement', 'IgnoreCase', true), 1);

    if isempty(ux_idx) || isempty(uy_idx)
        warning('Could not find displacement fields in newEBSD.prop for step %s. Skipping.', defStep.stepName);
        continue;
    end
    ux_field = all_fields{ux_idx};
    uy_field = all_fields{uy_idx};

    if ~isempty(hxx_idx) && ~isempty(hxy_idx) && ~isempty(hyx_idx) && ~isempty(hyy_idx)
        fprintf('  -> Using pre-computed DIC displacement gradient fields directly from dataset.\n');
        Hxx_raw = double(newEBSD.prop.(all_fields{hxx_idx}));
        Hxy_raw = double(newEBSD.prop.(all_fields{hxy_idx}));
        Hyx_raw = double(newEBSD.prop.(all_fields{hyx_idx}));
        Hyy_raw = double(newEBSD.prop.(all_fields{hyy_idx}));
    else
        fprintf('  -> Computing displacement gradient field via Gaussian smoothed displacements...\n');
        ebsd_grid = newEBSD.gridify;
        grid_size = size(ebsd_grid);
        
        U_raw = double(ebsd_grid.prop.(ux_field)); U_raw(U_raw == 0) = NaN;
        V_raw = double(ebsd_grid.prop.(uy_field)); V_raw(V_raw == 0) = NaN;
        U_grid = reshape(U_raw * ebsd_grid.dPos, grid_size);
        V_grid = reshape(V_raw * ebsd_grid.dPos, grid_size);
        
        filt_opt.filt_std = 4;
        try
            sm_data = filterDisplacements(U_grid, V_grid, filt_opt);
            U_grid = sm_data.U; V_grid = sm_data.V;
        catch
            warning('filterDisplacements not found. Using raw displacement grid for gradient calculation.');
        end
        
        [Hxx_grid, Hxy_grid] = gradient(U_grid, ebsd_grid.dPos, ebsd_grid.dPos);
        [Hyx_grid, Hyy_grid] = gradient(V_grid, ebsd_grid.dPos, ebsd_grid.dPos);
        
        grid_ids  = ebsd_grid.id(:);
        valid_idx = ~isnan(grid_ids) & grid_ids > 0;
        [~, loc] = ismember(newEBSD.id, grid_ids(valid_idx));
        
        tmp = Hxx_grid(:); Hxx_raw = tmp(valid_idx); Hxx_raw = Hxx_raw(loc);
        tmp = Hxy_grid(:); Hxy_raw = tmp(valid_idx); Hxy_raw = Hxy_raw(loc);
        tmp = Hyx_grid(:); Hyx_raw = tmp(valid_idx); Hyx_raw = Hyx_raw(loc);
        tmp = Hyy_grid(:); Hyy_raw = tmp(valid_idx); Hyy_raw = Hyy_raw(loc);
    end

    %% 6. Define Slip Systems for Titanium (HCP)
    sS_basal = slipSystem.basal(newEBSD.CS);
    sS_prism = slipSystem.prismaticA(newEBSD.CS);
    sS = [sS_basal; sS_prism];
    sS = symmetrise(sS, 'antipodal');

    valid_grains = grains(grains.numPixel > minSize);
    grains_ids = valid_grains.id;
    fprintf('  -> Ready to run subtracted SSLIP across %d significant grain(s).\n', length(grains_ids));

    % --- Plot 1: Displacement Gradient Field BEFORE Subtraction ---
    all_vals_raw = [Hxx_raw(:); Hxy_raw(:); Hyx_raw(:); Hyy_raw(:)];
    all_vals_raw = all_vals_raw(~isnan(all_vals_raw));
    cmin_raw = 0.2 * min(all_vals_raw);
    cmax_raw = 0.2 * max(all_vals_raw);
    
    figRaw = figure('Name', 'Displacement Gradient Verification BEFORE Subtraction', 'Visible', 'off', 'Position', [150, 150, 1200, 900]);
    newMtexFigure('layout', [2, 2]);
    raw_vals   = {Hxx_raw, Hxy_raw, Hyx_raw, Hyy_raw};
    raw_titles = {'du/dx (H_{11} Raw)', 'du/dy (H_{12} Raw)', 'dv/dx (H_{21} Raw)', 'dv/dy (H_{22} Raw)'};
    raw_cb_lbl = {'du/dx', 'du/dy', 'dv/dx', 'dv/dy'};
    for k = 1:4
        if k > 1, nextAxis; end
        plot(newEBSD, raw_vals{k}); hold on;
        plot(grains.boundary, 'lineW', 1, 'color', 'k');
        quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
        quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
        title(raw_titles{k});
        colormap(gca, 'parula');
        mtexColorbar('title', raw_cb_lbl{k});
        caxis([cmin_raw cmax_raw]);
        hold off;
    end
    sgtitle(sprintf('[%s] Experimental Displacement Gradients BEFORE Subtraction', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figRaw, fullfile(outDir, sprintf('%s_BeforeSubtraction_Gradients.png', defStep.stepName)), 'nodatestring');
    close(figRaw);

    %% 7. Loop over Selected Twin Variants to Subtract and Run SSLIP
    for var_idx = 1:length(target_variants)
        v = target_variants(var_idx);
        if v > numel(tS_sym)
            warning('Variant index %d out of bounds (max %d). Skipping.', v, numel(tS_sym));
            continue;
        end
        
        fprintf('\n   =====================================================================\n');
        fprintf('   Applying Global Subtraction for Twin Variant %d (of %d)\n', v, numel(tS_sym));
        fprintf('   =====================================================================\n');
        
        % Extract 2D theoretical displacement gradient tensor for variant v
        H_tv = H_twin_global(v).M(1:2, 1:2);
        fprintf('     Theoretical Twin Gradient Tensor H_id(v=%d):\n', v);
        fprintf('       [ %8.4f,  %8.4f ]\n       [ %8.4f,  %8.4f ]\n', H_tv(1,1), H_tv(1,2), H_tv(2,1), H_tv(2,2));
        
        % Rigorous Multiplicative Decoupling: F_slip = F_total * (F_twin)^-1
        F_twin_2d = eye(2) + H_tv;
        F_twin_inv = inv(F_twin_2d);
        
        a11 = F_twin_inv(1,1); a12 = F_twin_inv(1,2);
        a21 = F_twin_inv(2,1); a22 = F_twin_inv(2,2);
        
        F_tot_11 = 1 + Hxx_raw;
        F_tot_12 = Hxy_raw;
        F_tot_21 = Hyx_raw;
        F_tot_22 = 1 + Hyy_raw;
        
        F_slip_11 = F_tot_11 .* a11 + F_tot_12 .* a21;
        F_slip_12 = F_tot_11 .* a12 + F_tot_12 .* a22;
        F_slip_21 = F_tot_21 .* a11 + F_tot_22 .* a21;
        F_slip_22 = F_tot_21 .* a12 + F_tot_22 .* a22;
        
        % Residual displacement gradient driving secondary slip: H_slip = F_slip - I
        newEBSD.prop.Hxx_sub = F_slip_11 - 1;
        newEBSD.prop.Hxy_sub = F_slip_12;
        newEBSD.prop.Hyx_sub = F_slip_21;
        newEBSD.prop.Hyy_sub = F_slip_22 - 1;
        
        % Calculate twinned parent orientation (parent orientation transformed by twin variant misorientation)
        if isprop(tS_sym(v), 'misorientation') || ismethod(tS_sym(v), 'misorientation') || isfield(tS_sym(v), 'misorientation')
            twinMis = tS_sym(v).misorientation;
        else
            twinMis = tS_sym(v);
        end
        parentOri_twinned = parentGrain.meanOrientation * twinMis;
        
        % --- Save Visual Check Plot of Subtracted Fields ---
        Hxx_s = newEBSD.prop.Hxx_sub; Hxy_s = newEBSD.prop.Hxy_sub;
        Hyx_s = newEBSD.prop.Hyx_sub; Hyy_s = newEBSD.prop.Hyy_sub;
        all_vals_sub = [Hxx_s(:); Hxy_s(:); Hyx_s(:); Hyy_s(:)];
        all_vals_sub = all_vals_sub(~isnan(all_vals_sub));
        cmin_sub = 0.2 * min(all_vals_sub);
        cmax_sub = 0.2 * max(all_vals_sub);
        
        figSub = figure('Name', 'Decoupled Secondary Slip Gradients AFTER Subtraction', 'Visible', 'off', 'Position', [150, 150, 1200, 900]);
        newMtexFigure('layout', [2, 2]);
        sub_fields = {Hxx_s, Hxy_s, Hyx_s, Hyy_s};
        sub_titles = {'H_{11} Decoupled Slip', 'H_{12} Decoupled Slip', 'H_{21} Decoupled Slip', 'H_{22} Decoupled Slip'};
        sub_cb_lbl = {'H_{11}', 'H_{12}', 'H_{21}', 'H_{22}'};
        for k = 1:4
            if k > 1, nextAxis; end
            plot(newEBSD, sub_fields{k}); hold on;
            plot(grains.boundary, 'lineW', 1, 'color', 'k');
            quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
            quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
            title(sprintf('%s (Var %d)', sub_titles{k}, v));
            colormap(gca, 'parula');
            mtexColorbar('title', sub_cb_lbl{k});
            caxis([cmin_sub cmax_sub]);
            hold off;
        end
        sgtitle(sprintf('[%s] Multiplicatively Decoupled Secondary Slip Gradients (Variant %d)', defStep.stepName, v), 'Interpreter', 'none');
        exportScaledFigure(figSub, fullfile(outDir, sprintf('%s_SubtractedGradients_Var%d.png', defStep.stepName, v)), 'nodatestring');
        close(figSub);

        % --- Run SSLIP across individual grains using Decoupled Gradients and Twinned Parent Orientation ---
        for g_i = 1:length(grains_ids)
            clear IDoptions
            close all
            grainId = grains_ids(g_i);

            fprintf('\n     [Variant %d | Grain %d (%d px)] Initializing slip system geometry in twinned lattice...\n', v, grainId, grains(grainId).numPixel);

            % Project slip systems using the twinned parent orientation instead of original parent orientation
            sSLocal = parentOri_twinned * sS;

            threshold_SF = 0.1;
            stressState = stressTensor.uniaxial(yvector);
            
            % Reconfigure slip systems to ensure positive Schmid factor under assumed loading
            sf_signs = sign(sSLocal.SchmidFactor(stressState));
            sf_signs(sf_signs == 0) = 1;
            sSLocal.n = sf_signs .* sSLocal.n;
            
            SF_crit = abs(sSLocal.SchmidFactor(stressState)) > threshold_SF;

            sSLocal_filt = sSLocal(SF_crit);
            sS_filt = sS(SF_crit);

            if isempty(sS_filt)
                warning('     No slip systems with SF > %.2f for Grain %d. Skipping grain.', threshold_SF, grainId);
                continue;
            end

            % Slice EBSD to grain and gridify (automatically slices Hxx_sub, etc.)
            ebsd_grain = newEBSD(grains(grainId));
            ebsd_grain = ebsd_grain.gridify;
            ebsd_grain.orientations = repmat(parentOri_twinned, size(ebsd_grain.orientations));

            U = ebsd_grain.prop.(ux_field) * newEBSD.dPos; 
            V = ebsd_grain.prop.(uy_field) * newEBSD.dPos;

            plotBaseName = sprintf('%s_SSLIP_Grain%d_Var%d_%s', defStep.stepName, grainId, v, defStep.dic_group);

            IDoptions.IDMethod = 1;
            IDoptions.filterSize = 4;
            IDoptions.coarsegrain = 4;
            IDoptions.threshResidual = 0.02;
            IDoptions.minEeff = 0.01;
            IDoptions.posConstr = 1;
            IDoptions.plotSSLIP = 1;
            IDoptions.layout = [];
            IDoptions.plotTraces = 1;
            IDoptions.stress = stressState;
            IDoptions.casename = plotBaseName;
            IDoptions.plotDefGrad = 1;
            IDoptions.saveFig = 1;
            IDoptions.plotname = fullfile(outDir, plotBaseName);
            IDoptions.sizeAdjust = 0.5;
            IDoptions.plotResidual = 1;
            IDoptions.logscale = 1;
            IDoptions.logmin = 0.01;
            try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
            IDoptions.comment = sprintf('Step_%s_Grain_%d_SubVar_%d', defStep.stepName, grainId, v);
            IDoptions.NoSs = 1:length(sS_filt);

            % Feed subtracted gradient fields into SSLIP optimizer!
            IDoptions.Hxx = ebsd_grain.prop.Hxx_sub;
            IDoptions.Hxy = ebsd_grain.prop.Hxy_sub;
            IDoptions.Hyx = ebsd_grain.prop.Hyx_sub;
            IDoptions.Hyy = ebsd_grain.prop.Hyy_sub;

            fprintf('     Running optimizer on %d candidate systems with subtracted shear...\n', length(sS_filt));

            try
                [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);

                % Save calculation MAT object
                resMatPath = fullfile(outDir, [plotBaseName, '_SSLIP_results.mat']);
                save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId', 'v', 'H_tv');
                fprintf('     -> Saved subtracted optimization results to: %s\n', resMatPath);

                % Save CSV table
                SF = abs(SchmidFactor(sSLocal_filt, IDoptions.stress));
                SF = round(SF*1000)/1000;

                T = table((1:length(sS_filt))', arrayfun(@(x) x.b.char, sS_filt, 'UniformOutput', false), ...
                          arrayfun(@(x) x.n.char, sS_filt, 'UniformOutput', false), SF', ...
                          'VariableNames', {'SlipSystem_ID', 'Burgers_Vector_b', 'Plane_Normal_n', 'Schmid_Factor'});

                csvPath = fullfile(outDir, [plotBaseName, '_SlipSystems.csv']);
                writetable(T, csvPath);
            catch ME
                warning('     SSLIP calculation failed for Grain %d (Variant %d) in step %s:\n        %s', grainId, v, defStep.stepName, ME.message);
            end
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL TWIN-SHEAR-SUBTRACTED SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
