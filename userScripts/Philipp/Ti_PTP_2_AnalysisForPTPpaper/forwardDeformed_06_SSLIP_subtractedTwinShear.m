%% forwardDeformed_06_SSLIP_subtractedTwinShear.m
% Applies the SSLIP analysis method to forward deformed DIC data on TKD datasets with multiple grains,
% globally subtracting the theoretical displacement gradient of a selected twin variant (v) 
% formed in the primary parent grain.
% This decouples twinning shear from secondary dislocation glide (slip) plasticity.
%
% Outputs results directly to: Processed_Output_Deformed_Forward/SSLIP_SubtractedTwinShear_Results.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\SSLIP'); initSSLIP();

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% ========================================================================
%  USER CONFIGURATION & PRESETS (Modify as desired)
% =========================================================================
cfg = struct();

% --- PRESET A: Subtracted Twin Shear SSLIP ---
cfg.analysis_name    = 'SubtractedTwinShear_Var6_rotationEnabled';
cfg.target_variants  = [6];         % Twin variants to subtract (e.g., [1, 2] or [6])
cfg.shear_scaling_factor = 1.0;     % Twinning shear multiplier (1.0 = theoretical shear)

cfg.slip_system_mode = 'all';       % 'all', 'a_only', 'ca_only', or 'basal_pyr'
cfg.threshold_SF     = 0.01;        % Minimum Schmid Factor threshold
cfg.stress_vector    = yvector;     % Uniaxial stress loading axis
cfg.restrict2GrainId = 0;           % 0 = all, [1,2,3] grainIds
cfg.stressSign = 1;                 %-1 compression, 1 extension
cfg.enableRotation   = 1;           % enable rotation fields

cfg.sslip_method     = 1;           % IDMethod (1 = multi-slip optimization)
cfg.posConstr        = 1;           % Enforce non-negative shear activity
cfg.threshResidual   = 0.02;        % Residual threshold for correlation rate calculation
cfg.threshResidualFraction = 0.50;  % Threshold residual as fraction of effective shear
cfg.minEeff          = 0.01;        % Minimum effective shear strain required for metric evaluation
cfg.minSize          = 300;         % Minimum grain size in pixels to analyze
cfg.logscale         = 0;           % 
cfg.bundleSimilar    = 1;           % 
cfg.outputReport = 1;               %

% Lightweight JPG export settings
cfg.jpg_dpi          = 150;
cfg.jpg_target_px    = 1200;
cfg.clim_percentile  = 99;
cfg.custom_clim      = [0 0.3];          % Custom color limits for slip activity fields (e.g., [0, 0.05]). Leave empty [] for percentile scaling.
cfg.plotBoundariesOverDisps = 0;

% =========================================================================

%% 1. Load Registry & Configure Step Selection
% [registry, defaults] = a00_Dataset_Registry_DifferentialDIC();
[registry, defaults] = a00_Dataset_Registry();

% Filter for non-excluded deformed steps with valid TKD maps and DIC groups
defSteps = registry(~[registry.isReference] & ~[registry.exclude] & ...
                    ~cellfun(@isempty, {registry.tkd_fname}) & ~cellfun(@isempty, {registry.dic_group}));

if isempty(defSteps)
    fprintf('========================================================================\n');
    fprintf('NOTICE: No valid deformed step entries found in registry.\n');
    fprintf('========================================================================\n');
    return;
end

inDir  = fullfile(defaults.dir, defaults.outDir_deformed);
outDir = fullfile(defaults.dir, defaults.outDir_deformed, ['SSLIP_SubtractedTwinShear_Results_', cfg.analysis_name]);
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING TWIN-SHEAR-SUBTRACTED SSLIP BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

%% 2. Iterate Over Each Deformed Step
for stepIdx = 2 %1:numel(defSteps)
    defStep = defSteps(stepIdx);
    prefix = defStep.stepName;
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('Processing Step [%s] (%d/%d) - DIC Group: %s\n', prefix, stepIdx, numel(defSteps), defStep.dic_group);
    fprintf('------------------------------------------------------------------------\n');
    
    matPath = fullfile(inDir, [prefix, '_deformedEBSD_withForwardDIC.mat']);
    if ~exist(matPath, 'file')
        warning('Aligned dataset not found at:\n  %s\nRun forwardDeformed_03_IMPORTER first for Step [%s]. Skipping.', matPath, prefix);
        continue;
    end
    
    fprintf('Loading aligned dataset: %s...\n', matPath);
    dataLoaded = load(matPath);
    newEBSD = dataLoaded.newEBSD;
    
    %% 3. Retrieve or Reconstruct Grains & Identify Primary Parent Grain
    grainThresh = 5*degree; 
    alpha = 1;

    if isfield(dataLoaded, 'grains')
        grains = dataLoaded.grains;
        fprintf('  -> Successfully loaded precomputed grains directly from .mat dataset.\n');
    else
        fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
        thresh_ci = defaults.thresh_ci;
        newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
        [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', cfg.minSize, 'angle', grainThresh, 'alpha', alpha);
    end

    % Identify primary Parent Grain (matrix grain with largest surface area/pixels)
    [~, parentIdx] = max(grains.numPixel);
    parentGrain = grains(parentIdx);
    fprintf('  -> Primary Parent Matrix Identified: Grain %d (%d px, (%.1f, %.1f) um).\n', ...
            parentGrain.id, parentGrain.numPixel, parentGrain.centroid.x, parentGrain.centroid.y);

    % Save reference Grain Map figure for this step
    figG = figure('Visible', 'off');
    plot(grains, grains.meanOrientation); hold on;
    text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
    title(sprintf('%s: Analyzed Grains (>%d px) | Parent: Grain %d', prefix, cfg.minSize, parentGrain.id), 'Interpreter', 'none');
    hold off;
    exportScaledFigure(figG, fullfile(outDir, sprintf('%s_GrainMap.jpg', prefix)), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
    close(figG);

    %% 4. Calculate Theoretical Twin Displacement Gradient Tensors ({10-12} Twinning)
    tS_1012 = twinSystem.hexagonal_1012(newEBSD.CS);
    tS_sym  = tS_1012.symmetrise; % All 6 primary variants
    
    % Transform theoretical twinning displacement gradients into global specimen reference frame
    % using the orientation of the primary parent grain.
    H_twin_global = parentGrain.meanOrientation * tS_sym.displacementGradient;

    %% 5. Locate the Forward DIC Displacement Fields
    all_fields = fields(newEBSD.prop);
    ux_idx = find(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true), 1);
    uy_idx = find(contains(all_fields, '_Displacements_y_displacement', 'IgnoreCase', true), 1);
    
    if isempty(ux_idx) || isempty(uy_idx)
        warning('Could not find x and y displacement fields in newEBSD.prop for step %s. Skipping.', prefix);
        continue;
    end
    ux_field = all_fields{ux_idx};
    uy_field = all_fields{uy_idx};
    
    hxx_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hxy_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    hyx_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hyy_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    
    if ~isempty(hxx_idx) && ~isempty(hxy_idx) && ~isempty(hyx_idx) && ~isempty(hyy_idx)
        Hxx_raw = double(newEBSD.prop.(all_fields{hxx_idx}));
        Hxy_raw = double(newEBSD.prop.(all_fields{hxy_idx}));
        Hyx_raw = double(newEBSD.prop.(all_fields{hyx_idx}));
        Hyy_raw = double(newEBSD.prop.(all_fields{hyy_idx}));
        hasGradients = true;
    else
        warning('No raw displacement gradients found in newEBSD.prop. Run importer first.');
        continue;
    end

    %% 6. Define Candidate Slip Systems based on Toggled Mode
    CS = newEBSD.CS;
    sS_basal_all = slipSystem.basal(CS);
    sS_prism_all = slipSystem.prismaticA(CS);
    sS_pyra_all  = slipSystem.pyramidalA(CS);
    sS_ca_all    = slipSystem.pyramidalCA(CS);

    switch lower(cfg.slip_system_mode)
        case {'all', 'everything'}
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all; sS_ca_all];
            mode_str = 'All Systems (Basal, Prism, Pyr <a>, Pyr <c+a>)';
        case {'a', 'a_only', 'only_a'}
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all];
            mode_str = 'Only <a> Type Systems (Basal, Prism, Pyr <a>)';
        case {'ca', 'c+a', 'ca_only', 'only_ca'}
            sS_raw = sS_ca_all;
            mode_str = 'Only <c+a> Type Systems (Pyr <c+a>)';
        case {'basal_pyramidal', 'basal_pyr', 'only_basal_pyr'}
            sS_raw = [sS_basal_all; sS_pyra_all; sS_ca_all];
            mode_str = 'Basal & Pyramidal Systems (Excluding Prismatic)';
        otherwise
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all; sS_ca_all];
            mode_str = 'Default: All Systems';
    end

    sS = symmetrise(sS_raw, 'antipodal');
    fprintf('  -> Slip System Mode [%s]: %d candidate systems initialized.\n', cfg.slip_system_mode, length(sS));

    %% 7. Evaluate Significant Grains
    valid_grains = grains(grains.numPixel > cfg.minSize);
    if isempty(valid_grains)
        warning('No grains larger than %d pixels in step %s. Skipping.', cfg.minSize, prefix);
        continue;
    end

    grains_ids = valid_grains.id; 
    if isfield(cfg, 'restrict2GrainId') && ~isequal(cfg.restrict2GrainId, 0) && ~isempty(cfg.restrict2GrainId) && all(cfg.restrict2GrainId > 0)
        grains_ids = intersect(grains_ids, cfg.restrict2GrainId);
        if isempty(grains_ids)
            warning('None of the restricted grain IDs [%s] were found among significant grains. Skipping.', num2str(cfg.restrict2GrainId));
            continue;
        end
        fprintf('  -> Restricting analysis to %d specified grain(s).\n', length(grains_ids));
    else
        fprintf('  -> Analyzing all %d significant grains...\n', length(grains_ids));
    end
    
    combinedMatPath = fullfile(outDir, sprintf('%s_%s_CombinedSSLIP.mat', prefix, cfg.slip_system_mode));
    if exist(combinedMatPath, 'file')
        fprintf('  -> Loading combined SSLIP results from: %s\n', combinedMatPath);
        load(combinedMatPath, 'stepResults');
    else
        stepResults = struct();
    end

    %% 8. Loop over Selected Twin Variants to Subtract and Run SSLIP
    for var_idx = 1:length(cfg.target_variants)
        v = cfg.target_variants(var_idx);
        if v > numel(tS_sym)
            warning('Variant index %d out of bounds (max %d). Skipping.', v, numel(tS_sym));
            continue;
        end
        
        fprintf('\n   =====================================================================\n');
        fprintf('   Applying Global Subtraction for Twin Variant %d (of %d)\n', v, numel(tS_sym));
        fprintf('   =====================================================================\n');
        
        % Extract 2D theoretical displacement gradient tensor for variant v
        H_tv = H_twin_global(v).M(1:2, 1:2);
        H_tv_scaled = cfg.shear_scaling_factor * H_tv;
        
        fprintf('     Scaled Theoretical Twin Gradient Tensor H_id(v=%d, Scale=%.2f):\n', v, cfg.shear_scaling_factor);
        fprintf('       [ %8.4f,  %8.4f ]\n       [ %8.4f,  %8.4f ]\n', H_tv_scaled(1,1), H_tv_scaled(1,2), H_tv_scaled(2,1), H_tv_scaled(2,2));
        
        % Rigorous Multiplicative Decoupling: F_slip = F_total * (F_twin)^-1
        F_twin_2d = eye(2) + H_tv_scaled;
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

        newEBSD.prop.Hyy_sub = F_slip_22 - 1;

        for g_i = 1:length(grains_ids)
            close all;
            grainId = grains_ids(g_i);

            fprintf('\n   [Variant %d | Grain %d (%d px)] Initializing slip geometry...\n', v, grainId, grains(grainId).numPixel);

            % IMPORTANT FIX: Project global slip systems into the experimental local grain lattice orientation.
            % We DO NOT overwrite the experimental orientation with the theoretical twinned parent orientation,
            % as requested by the user, ensuring evaluation in true grain coordinates.
            sSLocal = grains(grainId).meanOrientation * sS;

            stressState = cfg.stressSign * stressTensor.uniaxial(cfg.stress_vector);
            
            % Reconfigure slip systems to ensure positive Schmid factor
            sf_vals = sSLocal.SchmidFactor(stressState);
            sf_signs = sign(sf_vals);
            sf_signs(sf_signs == 0) = 1;
            sSLocal.b = sf_signs .* sSLocal.b;
            
            SF_crit = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
            sSLocal_filt = sSLocal(SF_crit);
            sS_filt = sS(SF_crit);

            if isempty(sS_filt)
                warning('   No slip systems with SF >= %.2f for Step [%s] Grain %d. Skipping grain.', cfg.threshold_SF, prefix, grainId);
                continue;
            end

            % Extract EBSD data slice specific to this grain and gridify
            ebsd_grain = newEBSD(grains(grainId));
            ebsd_grain = ebsd_grain.gridify;

            % Convert displacement field to physical units via dPos
            U = ebsd_grain.prop.(ux_field) * newEBSD.dPos; 
            V = ebsd_grain.prop.(uy_field) * newEBSD.dPos;

            plotBaseName = sprintf('%s_Grain%d_Var%d_%s', prefix, grainId, v, cfg.slip_system_mode);

            % Configure New SSLIP Architecture Parameters
            ssl_cfg = struct();
            ssl_cfg.preprocess = struct();
            ssl_cfg.preprocess.filterSize  = 4;
            ssl_cfg.preprocess.coarsegrain = 4;
            
            ssl_cfg.solver = struct();
            ssl_cfg.solver.IDMethod       = cfg.sslip_method;
            ssl_cfg.solver.threshResidual = cfg.threshResidual;
            if isfield(cfg, 'threshResidualFraction')
                ssl_cfg.solver.threshResidualFraction = cfg.threshResidualFraction;
            end
            ssl_cfg.solver.posConstr      = cfg.posConstr;
            ssl_cfg.solver.enableRotation = cfg.enableRotation;
            

            
            % Setup DeformationData contract
            DefData = struct();
            DefData.Hxx = ebsd_grain.prop.Hxx_sub;
            DefData.Hxy = ebsd_grain.prop.Hxy_sub;
            DefData.Hyx = ebsd_grain.prop.Hyx_sub;
            DefData.Hyy = ebsd_grain.prop.Hyy_sub;

            variantGrainKey = sprintf('var_%d_grain_%d', v, grainId);

            try
                if false; %isfield(stepResults, variantGrainKey)
                    fprintf('   Loading existing SSLIP optimization results for Var %d Grain %d from combined mat...\n', v, grainId);
                    ebsdID       = stepResults.(variantGrainKey).ebsdID;
                    sSLocal_filt = stepResults.(variantGrainKey).sSLocal_filt;
                    optOut       = stepResults.(variantGrainKey).optOut;
                else
                    fprintf('   Running SSLIP optimization (%d candidate systems) with subtracted shear...\n', length(sS_filt));
                    [ebsdID, prepData, optOut] = SSLIP(ebsd_grain, DefData, sSLocal_filt, ssl_cfg);
                    
                    % optOut is returned for plotting backward compatibility
                    optOut.plotname = fullfile(outDir, plotBaseName);
                    optOut.stress = stressState;
                    optOut.clim_percentile = cfg.clim_percentile;
                    optOut.posConstr = cfg.posConstr;
                        
                    % Save deformation gradient summary figure using new modular function
                    if isfield(cfg, 'plotDefGrad') && cfg.plotDefGrad
                        plotSSLIP_DeformationFields(ebsdID, optOut);
                        defGrad_fig = gcf;
                        exportScaledFigure(defGrad_fig, fullfile(outDir, [plotBaseName, '_DefGrad.jpg']), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                        close(defGrad_fig);
                    end
    
                    % Save into the combined struct immediately after solving
                    stepResults.(variantGrainKey).ebsdID = ebsdID;
                    stepResults.(variantGrainKey).sSLocal_filt = sSLocal_filt;
                    stepResults.(variantGrainKey).optOut = optOut;
                    save(combinedMatPath, 'stepResults', 'cfg', '-v7.3');
                end
                
                % --- NEW MODULAR PLOTTING SECTION ---
                close all;
                optOut.grains = grains;
                optOut.plotDominant = 1;
                optOut.num_top_systems = 3;
                optOut.saveFig = 0;
                optOut.outputReport = cfg.outputReport;
                optOut.ori = grains(grainId).meanOrientation;
                optOut.custom_clim = cfg.custom_clim;
                
                if cfg.outputReport
                    optOut.outputReport = fullfile(outDir, [plotBaseName, '_SlipSystems.csv']);
                end

                % Call Slip Activities
                plotSSLIP_SlipActivities(ebsdID, sSLocal_filt, optOut)
                exportScaledFigure(gcf, [optOut.plotname, '_SlipActivities.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                close(gcf);
                
                optOut.bundleSimilar = cfg.bundleSimilar;
                optOut.plotGrainContext = 1;
                
                % Call Dominant Activities Plotter
                plotSSLIP_DominantActivities(ebsdID, sSLocal_filt, optOut);
                exportScaledFigure(gcf, [optOut.plotname, '_SlipActivities_Dominant.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                close(gcf);
                
                % Call Residual Plotter with context
                plotSSLIP_Residual(ebsdID, optOut);
                exportScaledFigure(gcf, [optOut.plotname, '_Residual.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                close(gcf);
                
                % Call Rotation Plotter if rotation was calculated
                if isfield(ebsdID.prop, 'rotationIDcor')
                    plotSSLIP_Rotation(ebsdID, optOut);
                    exportScaledFigure(gcf, [optOut.plotname, '_Rotation.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    close(gcf);
                end
                fprintf('   -> Complete!\n');
            catch ME
                fprintf(2, '   [ERROR] SSLIP processing failed for Grain %d in step %s:\n%s\n', grainId, prefix, getReport(ME));
            end
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL TWIN-SHEAR-SUBTRACTED SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
