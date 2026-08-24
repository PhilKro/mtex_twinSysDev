%% undeformedReference_05_SSLIP.m
% Performs SSLIP analysis on the Undeformed Reference EBSD scan mapped across all DIC loading steps (from undeformedReference_03_IMPORTER).
% Evaluates slip system activity across multiple DIC increments on the undeformed reference grid.
% Outputs results directly to: Processed_Output_Undeformed/SSLIP_Results.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Load Registry & Locate Undeformed Reference Dataset
[registry, defaults] = a00_Dataset_Registry();
inDir  = fullfile(defaults.dir, defaults.outDir_undeformed);
outDir = fullfile(defaults.dir, defaults.outDir_undeformed, 'SSLIP_Results');
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
options = dataLoaded.options;
[~, matname_noext, ~] = fileparts(matname);

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

figG = figure('Visible', 'off');
plot(grains, grains.meanOrientation); hold on;
text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
title('Undeformed Reference: Reconstructed Grains', 'Interpreter', 'none'); hold off;
exportScaledFigure(figG, fullfile(outDir, 'Undeformed_Reference_GrainMap_.png'), 'nodatestring');
close(figG);

%% 3. Extract DIC Step Prefixes
all_fields = fields(newEBSD.prop);
displacement_fields_logical = contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true);
displacement_fields = all_fields(displacement_fields_logical);

prefix_list = cellfun(@(x) x(1:find(x=='_', 1)-1), displacement_fields, 'UniformOutput', false);
prefix_list = unique(prefix_list);
fprintf('Found %d DIC loading increments mapped to reference grid.\n', length(prefix_list));

%% 4. Define Slip Systems for Titanium (HCP)
sS_basal = slipSystem.basal(newEBSD.CS);
sS_prism = slipSystem.prismaticA(newEBSD.CS);
sS = [sS_basal; sS_prism];
sS = symmetrise(sS, 'antipodal');

%% 5. Filter & Evaluate Grains across Increments
valid_grains = grains(grains.numPixel > minSize); 
if isempty(valid_grains)
    error('No valid grains (>%d px) found for analysis.', minSize);
end

grains_ids = valid_grains.id; 
fprintf('Analyzing %d significant grains across all DIC increments...\n', length(grains_ids));

for i_i = 1:length(prefix_list)
    prefix = prefix_list{i_i};
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('Evaluating Reference SSLIP for DIC Increment [%s] (%d/%d)\n', prefix, i_i, length(prefix_list));
    fprintf('------------------------------------------------------------------------\n');
    
    for g_i = 1:length(grains_ids)
        grainId = grains_ids(g_i);
        
        fprintf('   [Grain %d (%d px)] Initializing slip geometry...\n', grainId, grains(grainId).numPixel);
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
            warning('   No slip systems with SF > %.2f for Grain %d. Skipping.', threshold_SF, grainId);
            continue;
        end
        
        ebsd_grain = newEBSD(grains(grainId));
        ebsd_grain = ebsd_grain.gridify;
        
        ux_fname = [prefix, '_Displacements_x_displacement'];
        uy_fname = [prefix, '_Displacements_y_displacement'];
        
        if ~isfield(ebsd_grain.prop, ux_fname) || ~isfield(ebsd_grain.prop, uy_fname)
            warning('   Missing displacement data %s for grain %d. Skipping.', ux_fname, grainId);
            continue;
        end
        
        U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos; 
        V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
        
        plotBaseName = sprintf('%s_SSLIP_Grain%d', prefix, grainId);
        
        clear cfg
        close all
        
        % -------------------------------------------------------------
        % Configure SSLIP run (Modern Architecture)
        % -------------------------------------------------------------
        cfg = struct();
        cfg.preprocessor.filterSize = 4;
        cfg.preprocessor.coarsegrain = 4;
        
        cfg.solver.IDMethod = 1; % Combined ID
        cfg.solver.threshResidual = 0.02;
        cfg.solver.minEeff = 0.01;
        cfg.solver.posConstr = 1; 
        cfg.solver.NoSs = 1:length(sS_filt);
        

        
        cfg.plot.stress = stressState;
        cfg.plot.casename = plotBaseName;
        cfg.plot.plotDefGrad = 1;
        cfg.plot.saveFig = 1;
        cfg.plot.plotname = fullfile(outDir, plotBaseName);
        cfg.plot.sizeAdjust = 0.5;
        cfg.plot.plotResidual = 1;
        cfg.plot.logscale = 1;
        cfg.plot.logmin = 0.01;
        try, cfg.plot.cmap = viridis(512); catch, cfg.plot.cmap = parula(256); end
        cfg.plot.comment = sprintf('RefGrain_%d_Inc_%s', grainId, prefix);
        cfg.plot.plotTraces = 1;
        
        fprintf('   Running SSLIP optimization (Vectorized)...\n');
        try
            [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, cfg);
            
            % External plotting decoupling
            plotSSLIP(ebsdID, sSLocal_filt, optOut);
            
            resMatPath = fullfile(outDir, [plotBaseName, '_SSLIP_results.mat']);
            save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId');
            
            SF = abs(SchmidFactor(sSLocal_filt, cfg.plot.stress));
            SF = round(SF*1000)/1000;
            
            T = table((1:length(sS_filt))', arrayfun(@(x) x.b.char, sS_filt, 'UniformOutput', false), ...
                      arrayfun(@(x) x.n.char, sS_filt, 'UniformOutput', false), SF', ...
                      'VariableNames', {'SlipSystem_ID', 'Burgers_Vector_b', 'Plane_Normal_n', 'Schmid_Factor'});
            
            writetable(T, fullfile(outDir, [plotBaseName, '_SlipSystems.csv']));
            fprintf('   -> Complete! Saved figures and CSV table to %s\n', outDir);
        catch ME
            warning('   SSLIP calculation failed for Grain %d in increment %s:\n      %s', grainId, prefix, ME.message);
        end
    end
end
fprintf('\n========================================================================\n');
fprintf('ALL UNDEFORMED REFERENCE SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
