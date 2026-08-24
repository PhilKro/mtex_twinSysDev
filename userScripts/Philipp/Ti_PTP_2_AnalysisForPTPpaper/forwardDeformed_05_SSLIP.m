%% forwardDeformed_05_SSLIP.m
% Performs multi-slip SSLIP optimization on all Deformed Forward datasets 
% (from forwardDeformed_03_IMPORTER).
% Features: Toggleable slip system filtering, Combined Error Index, lightweight JPG export,
% and comprehensive Markdown reporting.
% Outputs results directly to: Processed_Output_Deformed_Forward/SSLIP_Results.

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

% --- PRESET A: All Slip Systems with Low SF Threshold ---
cfg.analysis_name    = 'Standard_AllSlipSystems_LowSF_rotationEnabled';
cfg.slip_system_mode = 'all';       % 'all', 'a_only', 'ca_only', or 'basal_pyr'
cfg.threshold_SF     = 0.01;        % Minimum Schmid Factor threshold
cfg.stress_vector    = yvector;     % Uniaxial stress loading axis
cfg.restrict2GrainId = 0;           % 0 = all, [1,2,3] grainIds
cfg.stressSign       = 1;           % -1 compression, 1 extension
cfg.enableRotation   = 1;           % enable rotation fields

cfg.sslip_method     = 1;           % IDMethod (1 = multi-slip optimization)
cfg.posConstr        = 1;           % Enforce non-negative shear activity
cfg.threshResidual   = 0.02;        % Residual threshold for correlation rate calculation
cfg.threshResidualFraction = 0.50;  % Threshold residual as fraction of effective shear
cfg.minEeff          = 0.01;        % Minimum effective shear strain required for metric evaluation
cfg.minSize          = 300;         % Minimum grain size in pixels to analyze
cfg.logscale = 0; 
cfg.bundleSimilar = 1;


% Lightweight JPG export settings
cfg.jpg_dpi          = 150;
cfg.jpg_target_px    = 1200;
cfg.clim_percentile  = 99;
cfg.custom_clim      = [0 0.125];          % Custom color limits for slip activity fields (e.g., [0, 0.05]). Leave empty [] for percentile scaling.
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
outDir = fullfile(defaults.dir, defaults.outDir_deformed, ['SSLIP_Results_', cfg.analysis_name]);
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING DEFORMED FORWARD SSLIP BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

%% 2. Iterate Over Each Deformed Step
for stepIdx = 1:numel(defSteps)
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
    
    %% 3. Retrieve or Reconstruct Grains
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

    % Save reference Grain Map figure for this step
    figG = figure('Visible', 'off');
    plot(grains, grains.meanOrientation); hold on;
    text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
    title(sprintf('%s: Analyzed Grains (>%d px)', prefix, cfg.minSize), 'Interpreter', 'none');
    hold off;
    exportScaledFigure(figG, fullfile(outDir, sprintf('%s_GrainMap.jpg', prefix)), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
    close(figG);

    %% 4. Locate the Forward DIC Displacement Fields
    all_fields = fields(newEBSD.prop);
    % Displacement logic removed as per user preference for direct gradients
    
    hxx_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hxy_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    hyx_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hyy_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    
    if isempty(hxx_idx) || isempty(hxy_idx) || isempty(hyx_idx) || isempty(hyy_idx)
        warning('Could not find DIC gradients (Hxx, Hxy, Hyx, Hyy) for step %s. Skipping.', prefix);
        continue;
    end

    %% 5. Define Candidate Slip Systems based on Toggled Mode
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

    %% 6. Evaluate Significant Grains
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

    for g_i = 1:length(grains_ids)
        clear ssl_cfg;
        close all;
        grainId = grains_ids(g_i);

        fprintf('\n   [Grain %d (%d px)] Initializing slip geometry...\n', grainId, grains(grainId).numPixel);

        % Project global slip systems into local grain lattice orientation
        sSLocal = grains(grainId).meanOrientation * sS;

        stressState = cfg.stressSign * stressTensor.uniaxial(cfg.stress_vector);
        
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



        plotBaseName = sprintf('%s_Grain%d_%s', prefix, grainId, cfg.slip_system_mode);

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
        ssl_cfg.solver.stress         = stressState;
        
        % Setup DeformationData contract
        DefData = struct();
        DefData.Hxx = ebsd_grain.prop.(all_fields{hxx_idx});
        DefData.Hxy = ebsd_grain.prop.(all_fields{hxy_idx});
        DefData.Hyx = ebsd_grain.prop.(all_fields{hyx_idx});
        DefData.Hyy = ebsd_grain.prop.(all_fields{hyy_idx});
        
        grainKey = sprintf('grain_%d', grainId);

        try
            if false % isfield(stepResults, grainKey)
                fprintf('   Loading existing SSLIP optimization results for Grain %d from combined mat...\n', grainId);
                ebsdID       = stepResults.(grainKey).ebsdID;
                sSLocal_filt = stepResults.(grainKey).sSLocal_filt;
                optOut       = stepResults.(grainKey).optOut;
            else
                fprintf('   Running SSLIP optimizer for Grain %d...\n', grainId);
                [ebsdID, ~, optOut, sSLocal_filt] = SSLIP(ebsd_grain, DefData, sSLocal_filt, ssl_cfg);
                
                % optOut is returned for plotting backward compatibility
                optOut.plotname = fullfile(outDir, plotBaseName);
                optOut.stress = stressState;
                optOut.clim_percentile = cfg.clim_percentile;
                optOut.custom_clim = cfg.custom_clim;
                optOut.posConstr = cfg.posConstr;
                
                if isfield(cfg, 'plotDefGrad') && cfg.plotDefGrad
                    plotSSLIP_DeformationFields(ebsdID, optOut);
                    defGrad_fig = gcf;
                    exportScaledFigure(defGrad_fig, fullfile(outDir, [plotBaseName, '_DefGrad.jpg']), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    close(defGrad_fig);
                end
    
                % Save into the combined struct immediately after solving
                stepResults.(grainKey).ebsdID = ebsdID;
                stepResults.(grainKey).sSLocal_filt = sSLocal_filt;
                stepResults.(grainKey).optOut = optOut;
                save(combinedMatPath, 'stepResults', 'cfg', '-v7.3');
            end
            
            % --- NEW MODULAR PLOTTING SECTION ---
            close all;
            optOut.grains = grains;
            optOut.plotDominant = 1;
            optOut.num_top_systems = 3;
            optOut.saveFig = 0;
            optOut.ori = grains(grainId).meanOrientation;
            optOut.outputReport = fullfile(outDir, [plotBaseName, '_SlipSystems.csv']);
            
            % Call Slip Activities
            plotSSLIP_SlipActivities(ebsdID, sSLocal_filt, optOut)
            exportScaledFigure(gcf, [optOut.plotname, '_SlipActivities.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
            close(gcf);

            optOut.bundleSimilar = cfg.bundleSimilar;
            optOut.plotGrainContext = 1;

            % Call Dominant Activities Plotter
            plotSSLIP_DominantActivities(ebsdID, sSLocal_filt, optOut);
            exportScaledFigure(gcf, [optOut.plotname, '_SlipActivitiesDominant.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
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
            

            close all;
            fprintf('   -> Complete!\n');
        catch ME
            fprintf(2, '   [ERROR] SSLIP processing failed for Grain %d in step %s:\n%s\n', grainId, prefix, getReport(ME));
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL DEFORMED FORWARD SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
