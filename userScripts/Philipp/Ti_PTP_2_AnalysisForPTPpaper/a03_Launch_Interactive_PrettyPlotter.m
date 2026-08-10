%% a03_Launch_Interactive_PrettyPlotter.m
% =========================================================================
% DATA ENGINEERING & ETL ORCHESTRATOR FOR INTERACTIVE PRETTY PLOTTER
% =========================================================================
% Following modern software design principles, this script handles all 
% Data Discovery, ETL (Extract, Transform, Load), Mask Application, and MAT
% Caching operations independently from the presentation UI layer.
%
% Workflow:
%   1. Discovery: Scans Ti_cAxis_PTP_2 for TKD .ang files while adhering to
%      exclusion rules defined in a00_Dataset_Registry().
%   2. User Selection: Presents an interactive dataset selector (listdlg / CLI).
%   3. Caching & Transformation: Loads precompiled combined MAT archives if
%      available. Otherwise, executes the full ingestion pipeline:
%      - Merges supplementary H5OINA property fields (_h5oina_props.mat).
%      - Applies standard TKD sample coordinate rotation (-70° + -7° tilt).
%      - Filters CI < 0.10, gridifies, and crops to specimen masks (_specimen_mask.mat).
%      - Pre-calculates grain boundaries (>300 px) and archives to disk.
%   4. Functional Launch: Invokes Interactive_PrettyPlotter(ebsd, grains, opts)
%      and captures the interactive studio controller handles.
%
% Author: Philipp / Antigravity
% Date:   2026-08-05
% =========================================================================

clear variables; close all; clc;
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

% MTEX Spatial Conventions
plotx2east;
plotzIntoPlane;

%% 1. Configure Directories & Locate Target Datasets
baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2';
if ~exist(baseDir, 'dir')
    try, [~, defs] = a00_Dataset_Registry(); baseDir = defs.dir; catch, baseDir = pwd; end
end

maskDir   = fullfile(baseDir, 'Masks');
cacheDir  = fullfile(baseDir, 'Combined_EBSD_MATs');
exportDir = fullfile(baseDir, 'PrettyPlots_Interactive');

if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
if ~exist(exportDir, 'dir'), mkdir(exportDir); end

fprintf('========================================================================\n');
fprintf('TKD DATASET DISCOVERY & ETL PIPELINE\n');
fprintf('Scanning: %s\n', baseDir);

angFiles = dir(fullfile(baseDir, '**', '*.ang'));
if isempty(angFiles)
    error('No .ang files discovered in directory %s.', baseDir);
end

% Filter against dataset registry exclusions
try
    [registry, ~] = a00_Dataset_Registry();
    excl_fnames = {registry([registry.exclude] & ~cellfun(@isempty, {registry.tkd_fname})).tkd_fname};
    clean_excluded = cellfun(@(x) strrep(x, '.ang', ''), excl_fnames, 'UniformOutput', false);
catch
    clean_excluded = {};
end

file_names = {};
file_paths = {};
for i = 1:length(angFiles)
    [~, bName, ~] = fileparts(angFiles(i).name);
    if ~any(strcmpi(clean_excluded, bName))
        file_names{end+1} = bName; %#ok<AGROW>
        file_paths{end+1} = fullfile(angFiles(i).folder, angFiles(i).name); %#ok<AGROW>
    end
end

if isempty(file_names)
    error('All discovered datasets are marked as excluded in the registry.');
end
fprintf('Found %d candidate TKD datasets.\n\n', length(file_names));

%%% 2. Interactive Dataset Selector (Single or Multi-Select Batching)
% Display list box dialog for quick specimen selection (Use Ctrl/Shift for batch mode)
[sel_indices, is_ok] = listdlg('ListString', file_names, ...
                               'SelectionMode', 'multiple', ...
                               'PromptString', 'Select TKD specimens (Select multiple / Ctrl+A for Batch Mode):', ...
                               'Name', 'EBSD ETL & Batch Selector', ...
                               'ListSize', [450, 350]);
if ~is_ok || isempty(sel_indices)
    fprintf('No dataset selected. Defaulting to first discovered dataset [%s].\n', file_names{1});
    sel_indices = 1;
end

% Define base persistent visual configuration across batch runs
opts = struct();
opts.exportDir      = exportDir;
opts.defaultBgField = 'BS';        % Fallback to 'ci' if BS is absent
opts.defaultIpfKey  = 'ipfColorKey';
opts.defaultIpfDir  = 'Y';
opts.batchMode      = (length(sel_indices) > 1);
opts.mask_spec      = true;        % Default enable specimen masking if available
opts.crop_spec      = true;        % Default enable specimen map cropping if available

%% 3. Execute ETL Pipeline & Semi-Automatic Batch Presentation Loop
for idx_step = 1:length(sel_indices)
    sel_idx = sel_indices(idx_step);
    baseName = file_names{sel_idx};
    ang_path = file_paths{sel_idx};
    mat_cache_path = fullfile(cacheDir, [baseName, '_combinedEBSD.mat']);
    
    fprintf('\n========================================================================\n');
    fprintf('PROCESSING DATASET [%d/%d]: [%s]\n', idx_step, length(sel_indices), baseName);
    
    if exist(mat_cache_path, 'file')
        fprintf('[Cache Hit] Loading pre-compiled combined EBSD domain object from MAT...\n');
        fprintf('  File: %s\n', mat_cache_path);
        cData     = load(mat_cache_path);
        ebsd      = cData.ebsd;
        ebsd_grid = cData.ebsd_grid;
        grains    = cData.grains;
        
        % Ensure specimen_location is cleanly attached even if loading from an older MAT archive
        mask_file = fullfile(maskDir, sprintf('%s_specimen_mask.mat', baseName));
        if ~isfield(ebsd_grid.prop, 'specimen_location') && exist(mask_file, 'file')
            m_data = load(mask_file, 'specimenlocation');
            if isfield(m_data, 'specimenlocation')
                spec_loc = reshape(m_data.specimenlocation, size(ebsd_grid));
                if isequal(size(spec_loc), size(ebsd_grid))
                    ebsd_grid.prop.specimen_location = logical(spec_loc);
                    save(mat_cache_path, 'ebsd', 'ebsd_grid', 'grains', 'baseName', '-v7.3');
                end
            end
        end
    else
        fprintf('[Cache Miss] Running ingestion transformation pipeline from raw .ang...\n');
        CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
        ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
        
        % Merge supplementary H5OINA scalar properties
        prop_file = fullfile(maskDir, sprintf('%s_h5oina_props.mat', baseName));
        if exist(prop_file, 'file')
            fprintf('  -> Merging H5OINA diagnostic layers (%s_h5oina_props.mat)...\n', baseName);
            props = load(prop_file);
            p_flds = fields(props);
            for p_i = 1:length(p_flds)
                fn = p_flds{p_i}; val = props.(fn);
                if isnumeric(val) && numel(val) == length(ebsd)
                    ebsd.prop.(fn) = val(:);
                    fprintf('     + Attached property field: ebsd.prop.%s\n', fn);
                end
            end
        end
        
        % Apply standard TKD sample geometry rotations (-70° + -7° tilt around X, 180° around Z)
        tkd_tilt = -7 * degree;
        ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
        ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
        
        % Basic quality thresholding
        % ebsd(ebsd.ci < 0.10) = 'notIndexed';
        ebsd_grid = ebsd.gridify;    
    
        % Integrate Specimen Masking (Non-destructive logical property attachment)
        mask_file = fullfile(maskDir, sprintf('%s_specimen_mask.mat', baseName));
        if exist(mask_file, 'file')
            fprintf('  -> Loading specimen mask boundary (%s_specimen_mask.mat)...\n', baseName);
            m_data = load(mask_file, 'specimenlocation');
            if isfield(m_data, 'specimenlocation')
                spec_loc = m_data.specimenlocation;
                spec_loc = reshape(spec_loc, size(ebsd_grid));
                if isequal(size(spec_loc), size(ebsd_grid))
                    ebsd_grid.prop.specimen_location = logical(spec_loc);
                    fprintf('     + Attached non-destructive ebsd_grid.prop.specimen_location for GUI toggles.\n');
                elseif numel(spec_loc) == length(ebsd)
                    error('Your specimen location is acting up please check.')
                end
            end
        else
            fprintf('  -> No specimen mask found. Retaining full spatial footprint.\n');
        end
    
        % Precompute Grain Topology (>300 px)
        fprintf('  -> Computing grain boundaries (>300 px, 5°)...\n');
        [grains, ~] = calcGrains(ebsd_grid, 'minPixel', 300, 'angle', 5*degree);
        
        % Archive fully consolidated domain dataset to MAT cache
        save(mat_cache_path, 'ebsd', 'ebsd_grid', 'grains', 'baseName', '-v7.3');
        fprintf('  -> Archived combined domain dataset to disk:\n     %s\n', mat_cache_path);
    end

    % change scan unit to nm
    ebsd_grid.scanUnit = 'nm';
    %% 4. Launch Functional Presentation Studio with Persistent State
    fprintf('------------------------------------------------------------------------\n');
    fprintf('LAUNCHING INTERACTIVE STUDIO [%d/%d]: %s\n', idx_step, length(sel_indices), baseName);
    
    opts.datasetName = baseName;
    opts.currentStep = idx_step;
    opts.totalSteps  = length(sel_indices);
    
    [status, opts] = Interactive_PrettyPlotter(ebsd_grid, grains, opts);
    
    if strcmpi(status, 'exit')
        fprintf('\n[Batch Aborted] User terminated batch processing at dataset %d (%s).\n', idx_step, baseName);
        break;
    elseif strcmpi(status, 'skip')
        fprintf('\n[Batch Skipped] Advancing directly to next dataset without saving.\n');
        continue;
    else
        fprintf('\n[Batch Complete] Advancing to dataset %d/%d while preserving UI settings.\n', min(idx_step+1, length(sel_indices)), length(sel_indices));
    end
end
fprintf('========================================================================\n');
fprintf('ALL SELECTED DATASETS PROCESSED SUCCESSFULLY!\n');
fprintf('========================================================================\n');
