%% a03_forwardDeformed_Launch_Interactive_PrettyPlotter.m
% =========================================================================
% FORWARD DEFORMED TKD-DIC INTERACTIVE VISUALIZATION ORCHESTRATOR
% =========================================================================
% This script loads fully warped and consolidated forward-deformed EBSD/DIC
% domain structures from MAT archives (produced by forwardDeformed_03_IMPORTER.m
% and stored in Processed_Output_Deformed_Forward) and streams them directly 
% into the dual-mode Interactive Visualization Studio (Interactive_PrettyPlotter).
%
% Workflow:
%   1. Discovery: Scans Processed_Output_Deformed_Forward for precompiled
%      deformed EBSD domain MAT archives (*_deformedEBSD_withForwardDIC.mat).
%   2. Batch Selection: Presents an interactive selection dialog (listdlg)
%      enabling experimentalists to evaluate individual loading steps or run
%      semi-automated batch processing loops across multiple deformation increments.
%   3. Ingestion & Validation: Loads newEBSD (warped EBSD grid augmented with
%      DaVis DIC strains, displacements, and SEM rasters) and associated grain 
%      topologies without requiring computationally intensive re-alignment.
%   4. Studio Launch: Invokes Interactive_PrettyPlotter in serial progression
%      while preserving UI customization overrides between consecutive increments.
%
% Author: Philipp / Antigravity
% Date:   2026-08-06
% =========================================================================

clear variables; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\Philipp\Ti_PTP_2_AnalysisForPTPpaper');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\NanoMech_Aignment_Matlab\src');

% Establish uniform MTEX Spatial Conventions
plotx2east;
plotzIntoPlane;
try
    plottingConvention.default.east = xvector;
    plottingConvention.default.outOfScreen = -zvector;
catch
end

%% 1. Resolve Input Output Directories & Scan Deformed MAT Archives
try
    % [~, defaults] = a00_Dataset_Registry();
    [~, defaults] = a00_Dataset_Registry_DifferentialDIC();
    baseDir = defaults.dir;
    inputDir = fullfile(baseDir, defaults.outDir_deformed);
catch
    baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2';
    inputDir = fullfile(baseDir, 'Processed_Output_Deformed_Forward');
end

if ~exist(inputDir, 'dir')
    error('[Directory Error] Deformed data folder not found:\n  %s\nPlease run forwardDeformed_03_IMPORTER.m first.', inputDir);
end

exportDir = fullfile(inputDir, 'PrettyPlots_Interactive_Deformed');
if ~exist(exportDir, 'dir'), mkdir(exportDir); end

fprintf('========================================================================\n');
% FORWARD DEFORMED EBSD-DIC STUDIO LAUNCHER
fprintf('Scanning input directory:\n  %s\n', inputDir);

% Discover valid deformed EBSD/DIC MAT archives while ignoring alignment and config helpers
matFiles = dir(fullfile(inputDir, '*_deformedEBSD_withForwardDIC.mat'));
if isempty(matFiles)
    % Try fallback search for general combined EBSD mats if specific naming varies
    matFiles = dir(fullfile(inputDir, '*.mat'));
    matFiles = matFiles(~contains({matFiles.name}, 'alignment_data', 'IgnoreCase', true));
    matFiles = matFiles(~contains({matFiles.name}, 'specimen_mask', 'IgnoreCase', true));
    matFiles = matFiles(~contains({matFiles.name}, '_Params.mat', 'IgnoreCase', true));
end

if isempty(matFiles)
    error('No processed deformed MAT datasets found in directory:\n  %s\nEnsure forwardDeformed_03_IMPORTER.m has successfully generated output files.', inputDir);
end

file_names = {};
file_paths = {};
for i = 1:length(matFiles)
    fname = matFiles(i).name;
    % Extract human-readable deformation increment name
    dispName = regexprep(fname, '_deformedEBSD_withForwardDIC\.mat$', '');
    dispName = regexprep(dispName, '_combinedEBSD\.mat$', '');
    dispName = regexprep(dispName, '\.mat$', '');
    
    file_names{end+1} = dispName; %#ok<AGROW>
    file_paths{end+1} = fullfile(matFiles(i).folder, fname); %#ok<AGROW>
end

fprintf('Discovered %d consolidated deformed TKD-DIC datasets.\n\n', length(file_names));

%% 2. Interactive Step Selection (Single or Batch Mode)
[sel_indices, is_ok] = listdlg('ListString', file_names, ...
                               'SelectionMode', 'multiple', ...
                               'PromptString', 'Select Deformed Steps (Use Ctrl/Shift for Batch Loop):', ...
                               'Name', 'Deformed TKD-DIC Studio Selector', ...
                               'ListSize', [450, 320]);

if ~is_ok || isempty(sel_indices)
    fprintf('No deformation increments selected. Defaulting to first step [%s].\n', file_names{1});
    sel_indices = 1;
end

% Establish base GUI options configuration
opts = struct();
opts.exportDir          = exportDir;
opts.defaultBgField     = 'SE_raw';         % Prioritize SEM image or LapVar as background
opts.defaultOverlayMode = 'Scalar Property'; % Default to scalar property/strain mode for deformed datasets
opts.defaultIpfKey      = 'ipfColorKey';
opts.defaultIpfDir      = 'Y';
opts.batchMode          = (length(sel_indices) > 1);

%% 3. Execute Studio Presentation Loop
for idx_step = 1:length(sel_indices)
    sel_idx = sel_indices(idx_step);
    stepName = file_names{sel_idx};
    mat_path = file_paths{sel_idx};
    
    fprintf('\n========================================================================\n');
    fprintf('LOADING DEFORMED STEP [%d/%d]: [%s]\n', idx_step, length(sel_indices), stepName);
    fprintf('  File: %s\n', mat_path);
    
    cData = load(mat_path);
    
    % Robust EBSD grid structure extraction
    if isfield(cData, 'newEBSD')
        ebsd_grid = cData.newEBSD;
    elseif isfield(cData, 'ebsd_grid')
        ebsd_grid = cData.ebsd_grid;
    elseif isfield(cData, 'ebsd')
        ebsd_grid = cData.ebsd;
    else
        warning('[Skip] Neither newEBSD nor ebsd_grid found in file: %s', mat_path);
        continue;
    end
    
    % Ensure map grid architecture is maintained
    try, ebsd_grid = ebsd_grid.gridify; catch, end
    
    % Extract or evaluate grain topology
    if isfield(cData, 'grains') && ~isempty(cData.grains)
        grains = cData.grains;
    else
        fprintf('  -> Calculating grain boundary topology on-the-fly (>300 px, 5°)...\n');
        [grains, ~] = calcGrains(ebsd_grid, 'minPixel', 300, 'angle', 5*degree);
    end
   
    % Calculate Effective Shear Strain (geff) prior to studio launch

    props_list = fieldnames(ebsd_grid.prop);

    idx_eyy = find(contains(props_list, 'AxialStrain_yy', 'IgnoreCase', true), 1);
    idx_exx = find(contains(props_list, 'TransverseStrain_xx', 'IgnoreCase', true), 1);
    idx_exy = find(contains(props_list, 'ShearStrain_xy', 'IgnoreCase', true), 1);
    
    if isempty(idx_eyy) || isempty(idx_exx) || isempty(idx_exy)
        fprintf('  -> Strains not found in ebsd.props please check that they are imported correctly. No effective shear strain calculated\n');
    else
        geff = [];
        try
            Exx_interp = double(ebsd_grid.prop.(props_list{idx_exx})(:));
            Eyy_interp = double(ebsd_grid.prop.(props_list{idx_eyy})(:));
            Exy_interp = double(ebsd_grid.prop.(props_list{idx_exy})(:));
            geff = sqrt(((Exx_interp(:) - Eyy_interp(:))/2).^2 + ((Exy_interp(:) - Exy_interp(:))/2).^2);
            fprintf('     + Computed geff via DaVis gradient fields fallback.\n');
        catch ME
            warning(ME.identifier, 'Effective shear strain calculation failed.', ME.message);         
        end
        if ~isempty(geff) && numel(geff) == numel(ebsd_grid)
            ebsd_grid.prop.geff = reshape(geff, size(ebsd_grid));
            fprintf('     + Attached ebsd_grid.prop.geff for visualization studio.\n');
        end
    end
    
    
    % Verify fallback background layer availability
    available_props = fields(ebsd_grid.prop);
    if any(strcmpi(available_props, 'SE_raw'))
        opts.defaultBgField = 'SE_raw';
    elseif any(strcmpi(available_props, 'LapVar'))
        opts.defaultBgField = 'LapVar';
    elseif any(strcmpi(available_props, 'BS'))
        opts.defaultBgField = 'BS';
    else
        opts.defaultBgField = 'ci';
    end
    
    % Auto-select interesting default scalar overlay property if available
    if any(contains(available_props, 'geff', 'IgnoreCase', true))
        opts.defaultPropOverlay = 'geff';
    elseif any(contains(available_props, 'Shear_strain', 'IgnoreCase', true))
        idx_p = find(contains(available_props, 'Shear_strain', 'IgnoreCase', true), 1);
        opts.defaultPropOverlay = available_props{idx_p};
    elseif any(contains(available_props, 'du_dx', 'IgnoreCase', true))
        idx_p = find(contains(available_props, 'du_dx', 'IgnoreCase', true), 1);
        opts.defaultPropOverlay = available_props{idx_p};
    elseif any(contains(available_props, 'displacement', 'IgnoreCase', true))
        idx_p = find(contains(available_props, 'displacement', 'IgnoreCase', true), 1);
        opts.defaultPropOverlay = available_props{idx_p};
    end
    
    % change scan unit to nm
    ebsd_grid.scanUnit = 'nm';
    %% 4. Launch Interactive Visualization Studio
    fprintf('------------------------------------------------------------------------\n');
    fprintf('LAUNCHING INTERACTIVE STUDIO [%d/%d]: %s\n', idx_step, length(sel_indices), stepName);
    
    opts.datasetName = stepName;
    opts.currentStep = idx_step;
    opts.totalSteps  = length(sel_indices);
    
    % Launch interactive studio and retrieve user configuration overrides for subsequent increments
    [status, opts] = Interactive_PrettyPlotter(ebsd_grid, grains, opts);
    
    if strcmpi(status, 'exit')
        fprintf('\n[Batch Aborted] User exited batch processing at dataset %d (%s).\n', idx_step, stepName);
        break;
    elseif strcmpi(status, 'skip')
        fprintf('\n[Batch Skipped] Advancing directly to next dataset without saving.\n');
        continue;
    else
        fprintf('\n[Batch Step Complete] Advancing to step %d/%d while preserving interactive settings.\n', min(idx_step+1, length(sel_indices)), length(sel_indices));
    end
end

fprintf('========================================================================\n');
fprintf('ALL SELECTED FORWARD DEFORMATION INCREMENTS PROCESSED SUCCESSFULLY!\n');
fprintf('========================================================================\n');
