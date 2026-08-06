% =========================================================================
% undeformedReference_06b_Replot_SSLIP.m
% Fast Interactive Re-Plotting Tool for Completed SSLIP Analyses
% =========================================================================
% This lightweight utility allows you to load previously calculated SSLIP
% optimization results (*_SSLIP_results.mat) and experiment with color limits,
% colormaps, scalebars, and visual formatting in SECONDS without re-running
% the multi-minute DIC loading, grain decoupling, or optimization solvers.
% =========================================================================

clear all; close all; clc;

%% 1. TARGET ANALYSIS FOLDER & RESULTS
% Specify the exact folder where your SSLIP results were saved:
results_dir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2\Processed_Output_Undeformed\SSLIP_DecoupledTwinShear_Standard_MultiSlip_AllSystems_noPosConstr';

% Find all saved SSLIP result files in this directory
matFiles = dir(fullfile(results_dir, '*_SSLIP_results.mat'));
if isempty(matFiles)
    error('No *_SSLIP_results.mat files found in: %s\nPlease check the folder path.', results_dir);
end
fprintf('Found %d saved SSLIP optimization dataset(s) ready for fast re-plotting.\n', length(matFiles));

%% 2. CUSTOMIZE PLOTTING & COLOR SCALE OPTIONS
% Adjust any visual parameters here without modifying the core optimization:
replot_cfg.use_custom_maxE   = false;  % Set true to lock a specific maximum color limit (e.g. 0.05)
replot_cfg.maxE              = 0.05;   % Fixed upper bound when use_custom_maxE is true
replot_cfg.clim_percentile   = 99;     % Automatic upper limit percentile when use_custom_maxE is false (default: 99)

replot_cfg.logscale          = 0;      % 0 = Linear scale (Recommended), 1 = Logarithmic scale
replot_cfg.plotTraces        = 1;      % 1 = Overlay slip plane traces and direction arrows, 0 = Off
replot_cfg.plotResidual      = 1;      % 1 = Generate residual error map, 0 = Skip residual
replot_cfg.saveFig           = 1;      % 1 = Overwrite saved PNG figures in the folder, 0 = Display only

%% 3. EXECUTION LOOP 
for i = 1:length(matFiles)
    filePath = fullfile(results_dir, matFiles(i).name);
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('[%d/%d] Loading pre-computed results: %s\n', i, length(matFiles), matFiles(i).name);
    
    % Load saved variables: 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId', etc.
    data = load(filePath);
    
    if ~isfield(data, 'ebsdID') || ~isfield(data, 'optOut') || ~isfield(data, 'sSLocal_filt')
        warning('File %s does not contain complete SSLIP results structure. Skipping.', matFiles(i).name);
        continue;
    end
    
    % Merge customized options into optOut
    optOut = data.optOut;
    if replot_cfg.use_custom_maxE
        optOut.maxE = replot_cfg.maxE;
    else
        if isfield(optOut, 'maxE'), optOut = rmfield(optOut, 'maxE'); end
    end
    optOut.clim_percentile = replot_cfg.clim_percentile;
    optOut.logscale        = replot_cfg.logscale;
    optOut.plotTraces      = replot_cfg.plotTraces;
    optOut.plotResidual    = replot_cfg.plotResidual;
    optOut.saveFig         = replot_cfg.saveFig;
    
    % Re-run plotting directly
    fprintf('   -> Generating plots with updated color scaling...\n');
    tic;
    plotSSLIP(data.ebsdID.prop.slipIDcor, data.ebsdID.prop.residualEeff, data.ebsdID, data.sSLocal_filt, optOut);
    fprintf('   -> Complete in %.2f seconds!\n', toc);
end

fprintf('\n========================================================================\n');
fprintf('ALL FAST RE-PLOTTING TASKS COMPLETE!\n');
fprintf('========================================================================\n');
