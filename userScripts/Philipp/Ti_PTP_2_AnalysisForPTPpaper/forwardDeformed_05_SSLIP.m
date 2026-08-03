%% forwardDeformed_05_SSLIP.m
% Performs SSLIP analysis on all Deformed Forward datasets (from forwardDeformed_03_IMPORTER).
% Automatically iterates over all non-excluded deformed steps possessing both a TKD scan and a DIC group,
% evaluating slip system activities across all significant reconstructed grains.
% Outputs results directly to: Processed_Output_Deformed_Forward/SSLIP_Results.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

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
outDir = fullfile(defaults.dir, defaults.outDir_deformed, 'SSLIP_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING DEFORMED FORWARD SSLIP BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
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
    
    %% 3. Retrieve or Reconstruct Grains
    grainThresh = 5*degree; 
    minSize = 300; % Threshold in pixels to isolate significant grains
    alpha = 1;

    if isfield(dataLoaded, 'grains')
        grains = dataLoaded.grains;
        fprintf('  -> Successfully loaded precomputed, FZ-projected grains directly from .mat dataset.\n');
    else
        fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
        thresh_ci = defaults.thresh_ci;
        newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
        [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
    end

    % Save reference Grain Map figure for this step
    figG = figure('Visible', 'off');
    plot(grains, grains.meanOrientation); hold on;
    text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
    title(sprintf('%s: Reconstructed Grains (>%d px)', defStep.stepName, minSize), 'Interpreter', 'none');
    hold off;
    exportScaledFigure(figG, fullfile(outDir, [defStep.stepName, '_GrainMap_.png']), 'nodatestring');
    close(figG);

    %% 4. Locate the Single Forward DIC Displacement Fields
    all_fields = fields(newEBSD.prop);
    ux_idx = find(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true), 1);
    uy_idx = find(contains(all_fields, '_Displacements_y_displacement', 'IgnoreCase', true), 1);
    
    if isempty(ux_idx) || isempty(uy_idx)
        warning('Could not find x and y displacement fields in newEBSD.prop for step %s. Skipping.', defStep.stepName);
        continue;
    end
    ux_field = all_fields{ux_idx};
    uy_field = all_fields{uy_idx};
    fprintf('  -> Bound DIC Displacement fields: %s & %s\n', ux_field, uy_field);

    %% 5. Define Slip Systems for Titanium (HCP)
    sS_basal = slipSystem.basal(newEBSD.CS);
    sS_prism = slipSystem.prismaticA(newEBSD.CS);
    % Optional: Uncomment below to activate Pyramidal systems
    % sS_pyrA  = slipSystem.pyramidalA(newEBSD.CS);
    % sS_pyrCA = slipSystem.pyramidalCA(newEBSD.CS);

    sS = [sS_basal; sS_prism];
    sS = symmetrise(sS, 'antipodal');

    %% 6. Filter & Evaluate Grains
    valid_grains = grains(grains.numPixel > minSize);
    if isempty(valid_grains)
        warning('No grains larger than %d pixels in step %s. Skipping.', minSize, defStep.stepName);
        continue;
    end

    grains_ids = valid_grains.id;
    fprintf('  -> Found %d significant grain(s) for independent SSLIP fitting.\n', length(grains_ids));

    for g_i = 1:length(grains_ids)
        clear IDoptions
        close all
        grainId = grains_ids(g_i);

        fprintf('\n   [Grain %d (%d px)] Initializing local slip geometry...\n', grainId, grains(grainId).numPixel);

        % Project global slip systems into local grain lattice orientation
        sSLocal = grains(grainId).meanOrientation * sS;

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
            warning('   No slip systems with SF > %.2f for Step [%s] Grain %d. Skipping grain.', threshold_SF, defStep.stepName, grainId);
            continue;
        end

        % Extract EBSD data slice specific to this grain and gridify
        ebsd_grain = newEBSD(grains(grainId));
        ebsd_grain = ebsd_grain.gridify;

        % Convert displacement field to physical units via dPos
        U = ebsd_grain.prop.(ux_field) * newEBSD.dPos; 
        V = ebsd_grain.prop.(uy_field) * newEBSD.dPos;

        idName = sprintf('%s_Grain%d', defStep.stepName, grainId);
        plotBaseName = sprintf('%s_SSLIP_Grain%d_%s', defStep.stepName, grainId, defStep.dic_group);

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
        IDoptions.comment = sprintf('Step_%s_Grain_%d_DIC_%s', defStep.stepName, grainId, defStep.dic_group);
        IDoptions.NoSs = 1:length(sS_filt);

        % Detect and pass pre-computed Lagrangian displacement gradients if present in dataset
        propFields = fields(ebsd_grain.prop);
        hxx_idx = find(contains(propFields, {'_Hxx', '_dudx'}, 'IgnoreCase', true), 1);
        hxy_idx = find(contains(propFields, {'_Hxy', '_dudy'}, 'IgnoreCase', true), 1);
        hyx_idx = find(contains(propFields, {'_Hyx', '_dvdx'}, 'IgnoreCase', true), 1);
        hyy_idx = find(contains(propFields, {'_Hyy', '_dvdy'}, 'IgnoreCase', true), 1);

        if ~isempty(hxx_idx) && ~isempty(hxy_idx) && ~isempty(hyx_idx) && ~isempty(hyy_idx)
            fprintf('   -> Precomputed gradient fields found! Enabling Lagrangian SSLIP bypass.\n');
            IDoptions.Hxx = ebsd_grain.prop.(propFields{hxx_idx});
            IDoptions.Hxy = ebsd_grain.prop.(propFields{hxy_idx});
            IDoptions.Hyx = ebsd_grain.prop.(propFields{hyx_idx});
            IDoptions.Hyy = ebsd_grain.prop.(propFields{hyy_idx});
        else
            warning('   Could not find precomputed gradient fields. Using standard gradient calculation in SSLIP method that may contain artifacts due to forward interpolation of displacement field.');
        end

        

        fprintf('   Running optimization solver via SSLIP on %d candidate systems...\n', length(sS_filt));

        % Call SSLIP calculation core
        try
            [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);

            % Save calculation MAT object
            resMatPath = fullfile(outDir, [plotBaseName, '_SSLIP_results.mat']);
            save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId');
            fprintf('   -> Saved optimization results to: %s\n', resMatPath);

            % Build and save summary CSV report for local Schmid factors
            SF = abs(SchmidFactor(sSLocal_filt, IDoptions.stress));
            SF = round(SF*1000)/1000;

            T = table((1:length(sS_filt))', arrayfun(@(x) x.b.char, sS_filt, 'UniformOutput', false), ...
                      arrayfun(@(x) x.n.char, sS_filt, 'UniformOutput', false), SF', ...
                      'VariableNames', {'SlipSystem_ID', 'Burgers_Vector_b', 'Plane_Normal_n', 'Schmid_Factor'});

            csvPath = fullfile(outDir, [plotBaseName, '_SlipSystems.csv']);
            writetable(T, csvPath);
            fprintf('   -> Saved slip system table to: %s\n', csvPath);
        catch ME
            warning('   SSLIP calculation failed for Grain %d in step %s:\n      %s', grainId, defStep.stepName, ME.message);
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL DEFORMED FORWARD SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
