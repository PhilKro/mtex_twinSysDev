%% forwardDeformed_03_IMPORTER.m
% Two-Phase Deformed Forward Workflow:
%   Phase 1: Upfront Batch Interactive Point Selection (User Interactive)
%   Phase 2: Unattended Heavy Computation & Eulerian Forward Warping (25x speedup via Delaunay caching)
% Outputs directly to: Processed_Output_Deformed_Forward.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\NanoMech_Aignment_Matlab\src');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% Load Master Configuration & Dataset Registry
% [registry, defaults] = a00_Dataset_Registry();
[registry, defaults] = a00_Dataset_Registry_DifferentialDIC();

% Filter for deformed state entries that have a valid associated TKD dataset and are not excluded
defSteps = registry(~[registry.isReference] & ~[registry.exclude] & ~cellfun(@isempty, {registry.tkd_fname}));
if isempty(defSteps)
    fprintf('========================================================================\n');
    fprintf('NOTICE: No valid deformed step entries found in a00_Dataset_Registry().\n');
    fprintf('========================================================================\n');
    return;
end

%% ========================================================================
%  PHASE 1: UPFRONT BATCH INTERACTIVE ALIGNMENT POINT SELECTION
% =========================================================================
fprintf('\n========================================================================\n');
fprintf('--- PHASE 1: BATCH INTERACTIVE ALIGNMENT POINT SELECTION (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

for stepIdx = 1:numel(defSteps)
    defStep = defSteps(stepIdx);
    outPath = fullfile(defaults.dir, defaults.outDir_deformed);
    if ~exist(outPath, 'dir'), mkdir(outPath); end
    
    alignSavePath = fullfile(outPath, [defStep.stepName, '_alignment_data_deformed.mat']);
    
    % Determine if alignment points are already saved or provided in registry
    hasSavedFile = exist(alignSavePath, 'file') || ...
                   (~isempty(defStep.sel_points_fname) && exist(defStep.sel_points_fname, 'file')) || ...
                   (~isempty(defStep.sel_points_fname) && exist(fullfile(defaults.dir, defStep.sel_points_fname), 'file'));
    
    if hasSavedFile
        fprintf('Step [%s] (%d/%d): Alignment points already exist. Skipping interactive clicking.\n', defStep.stepName, stepIdx, numel(defSteps));
        continue;
    end
    
    fprintf('\nStep [%s] (%d/%d): No alignment points found. Loading minimal data for selection window...\n', defStep.stepName, stepIdx, numel(defSteps));
    
    ang_path = fullfile(defaults.dir, [defStep.tkd_fname, '.ang']);
    mask_path = fullfile(defaults.dir, defStep.mask_fname);
    SEImageFullPath = fullfile(defaults.SEImagefolder, defStep.se_fname);
    
    if ~exist(ang_path, 'file')
        warning('EBSD file not found at %s. Cannot select points for %s.', ang_path, defStep.stepName);
        continue;
    end
    
    CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
    ebsd_tmp = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame','setting 2');


    % Load optional extra properties from .h5oina if extracted
    prop_file = strrep(mask_path, '_specimen_mask.mat', '_h5oina_props.mat');
    if exist(prop_file, 'file')
        props = load(prop_file);
        if isfield(props, 'BS') && length(props.BS) == length(ebsd_tmp), ebsd_tmp.prop.BS = props.BS; end
        if isfield(props, 'BC') && length(props.BC) == length(ebsd_tmp), ebsd_tmp.prop.BC = props.BC; end
        if isfield(props, 'FSE') && length(props.FSE) == length(ebsd_tmp), ebsd_tmp.prop.FSE = props.FSE; end
        if isfield(props, 'BSE') && length(props.BSE) == length(ebsd_tmp), ebsd_tmp.prop.BSE = props.BSE; end
        if isfield(props, 'MSE') && length(props.MSE) == length(ebsd_tmp), ebsd_tmp.prop.MSE = props.MSE; end
        if isfield(props, 'GradMag') && length(props.GradMag) == length(ebsd_tmp), ebsd_tmp.prop.GradMag = props.GradMag; end
        if isfield(props, 'LapVar') && length(props.LapVar) == length(ebsd_tmp), ebsd_tmp.prop.LapVar = props.LapVar; end
    end


    ebsd_tmp = rotate(ebsd_tmp, rotation.byAxisAngle(xvector, -70*degree + defaults.tkd_tilt), 'keepXY');
    ebsd_tmp(ebsd_tmp.ci < defaults.thresh_ci) = 'notIndexed';
    ebsd_tmp = rotate(ebsd_tmp, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    ebsd_tmp = ebsd_tmp.gridify;
    
    if exist(mask_path, 'file')
        load(mask_path, 'specimenlocation');
        try, ebsd_tmp(~specimenlocation) = 'notIndexed'; catch, end
    end    
    
    rawData_tmp = struct();
    rawData_tmp.SE_raw.fname = SEImageFullPath;
    rawData_tmp.SE_raw.psize = 1;
    rawData_tmp.SE_raw.dtype = 'SEM';
    imInfo = imfinfo(SEImageFullPath);
    rawData_tmp.SE_raw.imregion = [1 1 imInfo.Width-1 imInfo.Height-1];
    rawData_tmp.EBSD0.fname = defStep.tkd_fname;
    rawData_tmp.EBSD0.data = ebsd_tmp;
    rawData_tmp.EBSD0.dtype = 'EBSD';
    rawData_tmp.dataNames = {'SE_raw', 'EBSD0'};
    rawData_tmp = loadRawData(rawData_tmp);
    
    options_tmp = struct();
    options_tmp.Ndof = 6;
    options_tmp.Nsel = 4;
    options_tmp.case = 'EBSD';
    options_tmp.selectionIPF = 1;
    options_tmp.EBSD_selection = 'LapVar';
    options_tmp.colormap = 'gray';
    options_tmp.sel_points_fname = 0; % Force interactive selection
    options_tmp.in_situIncrementName = [defStep.stepName, '_', datestr(now,'yyyymmdd_HHMMSS')];
    
    alignment_data_deformed = pointSelectionForAlignment(rawData_tmp.SE_raw, rawData_tmp.EBSD0, options_tmp);
    alignment_data_EBSD = alignment_data_deformed; % Duplicate for backward compatibility
    save(alignSavePath, 'alignment_data_deformed', 'alignment_data_EBSD');
    fprintf('  -> Successfully collected & saved alignment points for %s at:\n     %s\n', defStep.stepName, alignSavePath);
    
    close all; % Clean up alignment figures before moving to next step
end

fprintf('\n========================================================================\n');
fprintf('ALL INTERACTIVE ALIGNMENTS COMPLETE! You can leave the computer now.\n');
fprintf('--- PHASE 2: UNATTENDED HEAVY COMPUTATION & FORWARD WARPING ---\n');
fprintf('========================================================================\n');


%% ========================================================================
%  PHASE 2: UNATTENDED HEAVY COMPUTATION & EULERIAN FORWARD WARPING
% =========================================================================
for stepIdx = 1:numel(defSteps)
    defStep = defSteps(stepIdx);
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('[PHASE 2] Executing Forward Warping for Step: %s (%d/%d)\n', defStep.stepName, stepIdx, numel(defSteps));
    fprintf('------------------------------------------------------------------------\n');
    
    input = struct();
    input.fname             = defStep.tkd_fname;
    input.SEImagefname      = defStep.se_fname;
    input.mask_fname        = defStep.mask_fname;
    input.dir               = defaults.dir;
    input.SEImagefolder     = defaults.SEImagefolder;
    input.DIC_filepath      = defaults.DIC_filepath;
    input.dic_group         = defStep.dic_group;
    input.outPath           = fullfile(defaults.dir, defaults.outDir_deformed);

    input.SEImageFullPath = fullfile(input.SEImagefolder, input.SEImagefname);
    input.ang_path        = fullfile(input.dir, [input.fname, '.ang']);
    input.mask_path       = fullfile(input.dir, input.mask_fname);

    options = struct();
    options.dateString = datestr(now,'yyyymmdd_HHMMSS');
    options.in_situIncrementName = [defStep.stepName, '_', options.dateString];

    outname = fullfile(input.outPath, sprintf('%s_deformed_ioConfig.json', options.in_situIncrementName));
    jsonText = jsonencode(input);
    fid = fopen(outname,'w'); if fid > 0, fprintf(fid, '%s', jsonText); fclose(fid); end

    %% EBSD Loading & Properties
    if ~exist(input.ang_path, 'file')
        warning('EBSD scan not found at %s. Skipping.', input.ang_path);
        continue;
    end
    CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
    ebsd = EBSD.load(input.ang_path, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

    prop_file = fullfile(input.dir, 'Masks', sprintf('%s_h5oina_props.mat', input.fname));
    if exist(prop_file, 'file')
        try
            props = load(prop_file);
            if isfield(props, 'BS') && length(props.BS) == length(ebsd), ebsd.prop.BS = props.BS; end
            if isfield(props, 'BC') && length(props.BC) == length(ebsd), ebsd.prop.BC = props.BC; end
            if isfield(props, 'FSE') && length(props.FSE) == length(ebsd), ebsd.prop.FSE = props.FSE; end
            if isfield(props, 'BSE') && length(props.BSE) == length(ebsd), ebsd.prop.BSE = props.BSE; end
            if isfield(props, 'MSE') && length(props.MSE) == length(ebsd), ebsd.prop.MSE = props.MSE; end
            if isfield(props, 'LapVar') && length(props.LapVar) == length(ebsd), ebsd.prop.LapVar = props.LapVar; end
            if isfield(props, 'GradMag') && length(props.GradMag) == length(ebsd), ebsd.prop.GradMag = props.GradMag; end
        catch, end
    end

    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + defaults.tkd_tilt), 'keepXY');
    ebsd(ebsd.ci < defaults.thresh_ci) = 'notIndexed';
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    refOrientation = orientation.byEuler(defaults.refOrientationEulerAngles, ebsd.CS);
    ebsd.orientations = ebsd.orientations.project2FundamentalRegion(refOrientation);
    ebsd = ebsd.gridify;

    if exist(input.mask_path, 'file')
        load(input.mask_path, 'specimenlocation');
        try, ebsd(~specimenlocation) = 'notIndexed'; catch, end
    end

    %% Load Alignment & Warp EBSD to SE Grid
    rawData = struct();
    rawData.SE_raw.fname = input.SEImageFullPath;
    rawData.SE_raw.psize = 1;
    rawData.SE_raw.dtype = 'SEM';
    imInfo = imfinfo(input.SEImageFullPath);
    rawData.SE_raw.imregion = [1 1 imInfo.Width-1 imInfo.Height-1];
    rawData.EBSD0.fname = input.fname;
    rawData.EBSD0.data = ebsd;
    rawData.EBSD0.dtype = 'EBSD';
    rawData.dataNames = {'SE_raw', 'EBSD0'};
    rawData = loadRawData(rawData);

    alignSavePath = fullfile(input.outPath, [defStep.stepName, '_alignment_data_deformed.mat']);
    if exist(alignSavePath, 'file')
        loaddata = load(alignSavePath);
        if isfield(loaddata, 'alignment_data_deformed'), alignment_data = loaddata.alignment_data_deformed;
        else, alignment_data = loaddata.alignment_data_EBSD; end
    elseif ~isempty(defStep.sel_points_fname) && exist(defStep.sel_points_fname, 'file')
        loaddata = load(defStep.sel_points_fname);
        if isfield(loaddata, 'alignment_data_deformed'), alignment_data = loaddata.alignment_data_deformed;
        else, alignment_data = loaddata.alignment_data_EBSD; end
    elseif ~isempty(defStep.sel_points_fname) && exist(fullfile(defaults.dir, defStep.sel_points_fname), 'file')
        loaddata = load(fullfile(defaults.dir, defStep.sel_points_fname));
        if isfield(loaddata, 'alignment_data_deformed'), alignment_data = loaddata.alignment_data_deformed;
        else, alignment_data = loaddata.alignment_data_EBSD; end
    else
        error('Missing alignment file for %s in unattended Phase 2.', defStep.stepName);
    end

    psize_new = rawData.SE_raw.psize;
    [Xnew, Ynew] = meshgrid(0:psize_new:max(rawData.SE_raw.X(:)), 0:psize_new:max(rawData.SE_raw.Y(:)));
    newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data);
    newEBSD.prop.originalId = newEBSD.prop.oldId;
    originalEBSD = ebsd;
    originalEBSD.prop.originalId = ebsd.id;
    newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');
    newEBSD = newEBSD.gridify;

    %% DIC Data Loading & Extraction (Unstructured Eulerian Forward Interpolation)
    fprintf('Reading HDF5 DIC dataset group [%s]...\n', defStep.dic_group);
    dicFilePath = input.DIC_filepath;
    try, grpInfo = h5info(dicFilePath, ['/' defStep.dic_group]); catch, continue; end
    
    allPaths = collectH5Paths(grpInfo);
    
    xcoordPath = findDatasetPath(allPaths, 'Coordinates/x');
    if isempty(xcoordPath), xcoordPath = ['/' defStep.dic_group '/Coordinates/x']; end
    
    ycoordPath = findDatasetPath(allPaths, 'Coordinates/y');
    if isempty(ycoordPath), ycoordPath = ['/' defStep.dic_group '/Coordinates/y']; end
    
    maskPath = findDatasetPath(allPaths, 'Coordinates/Mask');
    if isempty(maskPath), maskPath = ['/' defStep.dic_group '/Coordinates/Mask']; end
    
    uxPath = findDatasetPath(allPaths, 'Displacements/x-displacement');
    if isempty(uxPath), uxPath = ['/' defStep.dic_group '/Displacements/x-displacement']; end
    
    uyPath = findDatasetPath(allPaths, 'Displacements/y-displacement');
    if isempty(uyPath), uyPath = ['/' defStep.dic_group '/Displacements/y-displacement']; end

    try
        X_dic_raw = double(h5read(dicFilePath, xcoordPath));
        Y_dic_raw = double(h5read(dicFilePath, ycoordPath));
        u_dic_raw = double(h5read(dicFilePath, uxPath));
        v_dic_raw = double(h5read(dicFilePath, uyPath));
    catch ME
        warning('Could not read core DaVis coordinates/displacements. Skipping DIC import. (%s)', ME.message);
        continue;
    end
    
    try
        Mask_raw = double(h5read(dicFilePath, maskPath));
    catch
        Mask_raw = ones(size(X_dic_raw));
    end
    
    % --- Determine true spatial offsets from SE_raw mapping ---
    % SE_raw crop coordinates (1-based index mapped to 0-based physical origin)
    x0 = rawData.SE_raw.imregion(1) - 1; 
    y0 = rawData.SE_raw.imregion(2) - 1; 
    
    % Extract physical scaling
    psize = rawData.SE_raw.psize;
    if isempty(psize) || isnan(psize) || psize == 0, psize = 1; end
    
    % --- Valid Eulerian Forward Cloud Construction ---
    % 1. Must be physically tracked by DaVis (Mask > 0)
    % 2. Must have valid displacements (not NaN and not exactly 0.0)
    valid_mask = ~isnan(u_dic_raw) & ~isnan(v_dic_raw) & ...
                 (u_dic_raw ~= 0.0) & (v_dic_raw ~= 0.0) & ...
                 (Mask_raw > 0);
    
    % Map to physical EBSD reference frame (Standard Forward Warping: X + u, Y + v)
    X_def_phys = (X_dic_raw(valid_mask) + u_dic_raw(valid_mask) - x0) * psize;
    Y_def_phys = (Y_dic_raw(valid_mask) + v_dic_raw(valid_mask) - y0) * psize;
    
    % --- Pre-Interpolation Verification Scatter Plot ---
    fprintf('  -> Saving pre-interpolation scatter plot for raw DIC subset...\n');
    fig_raw = figure('Visible', 'off');
    
    grid_sz = size(newEBSD);
    x_ext = [min(newEBSD.x(:)), max(newEBSD.x(:))];
    y_ext = [min(newEBSD.y(:)), max(newEBSD.y(:))];
    SE_img = reshape(newEBSD.prop.SE_raw, grid_sz);
    
    % Plot the background SE image in grayscale
    ax1 = axes('Parent', fig_raw);
    imagesc(ax1, x_ext, y_ext, SE_img); 
    colormap(ax1, 'gray');
    axis(ax1, 'image'); axis(ax1, 'off');
    
    % Overlay scatter on a second transparent axis so colormaps don't conflict
    ax2 = axes('Parent', fig_raw, 'Position', ax1.Position);
    v_valid = v_dic_raw(valid_mask);
    scatter(ax2, X_def_phys, Y_def_phys, 10, v_valid, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
    axis(ax2, 'image'); axis(ax2, 'off');
    xlim(ax2, x_ext); ylim(ax2, y_ext);
    set(ax2, 'Color', 'none', 'YDir', ax1.YDir, 'XDir', ax1.XDir); % Sync axis directions!
    colormap(ax2, 'parula');
    
    title(ax2, sprintf('Step %s: Raw Unstructured DIC Points (Pre-Interpolation)', defStep.stepName), 'Interpreter', 'none');
    
    if ~isempty(v_valid)
        v_min = prctile(v_valid, 2);
        v_max = prctile(v_valid, 98);
        if v_max == v_min, v_max = v_min + 1; end
        caxis(ax2, [v_min, v_max]);
    end
    cb = colorbar(ax2); cb.Title.String = 'Raw Displacement y';
    exportScaledFigure(fig_raw, fullfile(input.outPath, [defStep.stepName, '_RawDIC_PreInterpolation_.png']), 'nodatestring');
    close(fig_raw);
    
    % --- Initialize Single Caching scatteredInterpolant ---
    fprintf('  -> Constructing 2D Delaunay Triangulation on Deformed Frame (%d valid tracked points)...\n', sum(valid_mask(:)));
    tic;
    % Using 'natural' inside for smooth, artifact-free Voronoi neighbor interpolation (fixes barcode effect)
    % and crucially 'none' for extrapolation to prevent edge smearing artifacts!
    F_dic = scatteredInterpolant(X_def_phys, Y_def_phys, ones(size(X_def_phys)), 'natural', 'none');
    t_del = toc;
    % --- Compute Concave Hull to Prevent Extrapolation Across Notches ---
    fprintf('  -> Computing Alpha Shape (Concave Hull) to remove convex hull artifacts...\n');
    tic;
    % Estimate DIC grid point spacing in physical units to set a tight Alpha radius
    step_x = nanmedian(diff(X_dic_raw(1,:))) * psize;
    if step_x == 0, step_x = nanmedian(diff(X_dic_raw(:,1))) * psize; end
    step_y = nanmedian(diff(Y_dic_raw(:,1))) * psize;
    if step_y == 0, step_y = nanmedian(diff(Y_dic_raw(1,:))) * psize; end
    nominal_spacing = max(abs(step_x), abs(step_y));
    if nominal_spacing == 0 || isnan(nominal_spacing), nominal_spacing = 10; end
    
    % A radius of ~2.5x the grid spacing safely wraps the true boundaries 
    % without poking holes inside the domain, but tightly crops the hourglass notches.
    ashape = alphaShape(X_def_phys, Y_def_phys, nominal_spacing * 2.5);
    in_hull = inShape(ashape, newEBSD.x, newEBSD.y);
    t_shape = toc;
    fprintf('  -> Concave hull mask generated in %.2f s (Grid spacing: %.1f um). Interpolating %d data fields...\n', t_shape, nominal_spacing, numel(allPaths));
    
    % --- Interpolate All DIC Data Fields ---
    tic;
    baseName = defStep.dic_group;
    rawData.DIC = struct(); 
    rawData.DIC.names = {};
    
    for k = 1:numel(allPaths)
        fullPath = allPaths{k};
        try, raw_field = double(h5read(dicFilePath, fullPath)); catch, continue; end
        
        % Check if field matches spatial dimensions
        if numel(raw_field) ~= numel(X_dic_raw)
            continue; 
        end
        
        % Construct valid struct name
        rel = strrep(fullPath, ['/' defStep.dic_group '/'], '');
        dsSuffix = strjoin(strsplit(rel, '/'), '_');
        shortName = matlab.lang.makeValidName([baseName '_' dsSuffix]);
        
        % Set all non-valid points (outside mask or exactly 0.0) to NaN explicitly
        raw_field(~valid_mask) = NaN;
        
        % Extract valid points and rapidly interpolate onto EBSD grid
        F_dic.Values = raw_field(valid_mask);
        interp_data = F_dic(newEBSD.x, newEBSD.y);
        
        % Destroy the false convex hull interpolation across the notches!
        interp_data(~in_hull) = NaN; 
        newEBSD.prop.(shortName) = interp_data;
        
        % Save to rawData to maintain backwards compatibility
        rawData.DIC.(shortName) = raw_field; 
        rawData.DIC.names{end+1} = shortName;
    end
    t_interp = toc;
    fprintf('  -> Completed Eulerian forward spatial interpolation in %.2f s\n', t_interp);

    %% Check & Apply Coordinate System Corrections via 2D Grid Dominance Analysis
    coord_x_fields = find(contains(fieldnames(newEBSD.prop), 'Coordinates_x', 'IgnoreCase', true));
    coord_y_fields = find(contains(fieldnames(newEBSD.prop), 'Coordinates_y', 'IgnoreCase', true));
    all_prop_names = fieldnames(newEBSD.prop);

    if ~isempty(coord_x_fields) && ~isempty(coord_y_fields)
        dic_x_field = all_prop_names{coord_x_fields(1)};
        dic_y_field = all_prop_names{coord_y_fields(1)};
        
        grid_sz = size(newEBSD);
        DIC_x = reshape(newEBSD.prop.(dic_x_field), grid_sz);
        DIC_y = reshape(newEBSD.prop.(dic_y_field), grid_sz);
        ebsd_x = reshape(newEBSD.x, grid_sz);
        ebsd_y = reshape(newEBSD.y, grid_sz);
        
        tmp = abs(diff(DIC_x, 1, 2)); var_DICx_horz = nanmedian(tmp(:));
        tmp = abs(diff(DIC_x, 1, 1)); var_DICx_vert = nanmedian(tmp(:));
        
        if var_DICx_vert > var_DICx_horz
            warning('DIC_x varies vertically rather than horizontally! HDF5 read transposition check should have handled this.');
        else
            fprintf('  -> Grid Dominance check confirms DaVis X and Y axes match EBSD X and Y spatial matrix orientation.\n');
        end
        
        tmp = diff(ebsd_x, 1, 2); dir_EBSD_x = sign(nanmedian(tmp(:)));
        tmp = diff(DIC_x,  1, 2); dir_DIC_x  = sign(nanmedian(tmp(:)));
        tmp = diff(ebsd_y, 1, 1); dir_EBSD_y = sign(nanmedian(tmp(:)));
        tmp = diff(DIC_y,  1, 1); dir_DIC_y  = sign(nanmedian(tmp(:)));
        
        s_x = 1; s_y = 1;
        if dir_DIC_x ~= dir_EBSD_x, s_x = -1; end
        if dir_DIC_y ~= dir_EBSD_y, s_y = -1; end
        
        if s_x == -1 || s_y == -1
            warning('Coordinate directions oppose between DIC and EBSD! Applying sign corrections (s_x = %d, s_y = %d)...', s_x, s_y);
            for f_idx = 1:numel(all_prop_names)
                fn = all_prop_names{f_idx};
                if contains(fn, 'x_displacement', 'IgnoreCase', true)
                    newEBSD.prop.(fn) = s_x .* newEBSD.prop.(fn);
                elseif contains(fn, 'y_displacement', 'IgnoreCase', true)
                    newEBSD.prop.(fn) = s_y .* newEBSD.prop.(fn);
                elseif (contains(fn, 'du', 'IgnoreCase', true) && contains(fn, 'dy', 'IgnoreCase', true)) || (contains(fn, 'dv', 'IgnoreCase', true) && contains(fn, 'dx', 'IgnoreCase', true))
                    newEBSD.prop.(fn) = (s_x * s_y) .* newEBSD.prop.(fn);
                elseif contains(fn, 'Shear_strain', 'IgnoreCase', true) || contains(fn, 'xy', 'IgnoreCase', true)
                    newEBSD.prop.(fn) = (s_x * s_y) .* newEBSD.prop.(fn);
                end
            end
            fprintf('  -> Sign corrections applied completely.\n');
        else
            fprintf('  -> Coordinate directions increase identically (s_x = +1, s_y = +1). No negative sign flips needed.\n');
        end
    end

    %% Grain Calculation & Fundamental Region (FZ) Projection
    grainThresh = 5*degree; 
    minSize = 300; % Threshold in pixels for significant grains
    alpha = 1;

    [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

    % replaced by a unified FZ in the Dataset Registry
    % fprintf('  -> Projecting intra-grain orientations to Fundamental Zone centered on grain mean orientation...\n');
    % valid_grains = grains(grains.numPixel > 0);
    % for g_idx = 1:length(valid_grains)
    %     currGrainID = valid_grains(g_idx).id;
    %     currGrainMeanOri = valid_grains(g_idx).meanOrientation;
    % 
    %     mask_grain = (newEBSD.grainId == currGrainID);
    %     if any(mask_grain(:))
    %         newEBSD(mask_grain).orientations = newEBSD(mask_grain).orientations.project2FundamentalRegion(currGrainMeanOri);
    %     end
    % end
    % 
    % % Re-run calcGrains post-projection
    % [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

    %% Verification Plotting & Saving
    SE_gray = mat2gray(newEBSD.prop.SE_raw);
    SE_rgb = repmat(SE_gray, 1, 1, 3);
    
    fig1 = figure('Visible', 'off');
    plot(newEBSD, SE_rgb); axis image off; hold on;
    plot(grains.boundary, 'lineColor', 'red', 'linewidth', 1.5);
    title(sprintf('Step %s: FZ-Projected EBSD Grains over Deformed SE', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(fig1, fullfile(input.outPath, [defStep.stepName, '_DeformedEBSD_SE_.png']), 'nodatestring');
    close(fig1);

    if ~isempty(uyPath)
        % Identify the sanitized struct field name for y_displacement
        rel_y = strrep(uyPath, ['/' defStep.dic_group '/'], '');
        dsSuffix_y = strjoin(strsplit(rel_y, '/'), '_');
        uy_prop_name = matlab.lang.makeValidName([defStep.dic_group '_' dsSuffix_y]);
        
        if isfield(newEBSD.prop, uy_prop_name)
            % Extract scalar data and compute percentiles
            uy_data = newEBSD.prop.(uy_prop_name)(:);
            valid_uy = uy_data(~isnan(uy_data));
            
            if ~isempty(valid_uy)
                uy_min = prctile(valid_uy, 2);
                uy_max = prctile(valid_uy, 98);
                if uy_max == uy_min, uy_max = uy_min + 1; end
                
                % Normalize scalar data to [0, 1]
                uy_norm = (uy_data - uy_min) ./ (uy_max - uy_min);
                uy_norm = max(0, min(1, uy_norm));
                
                % Convert to RGB using parula colormap
                cmap = parula(256);
                idx = round(1 + 255 * uy_norm);
                idx(isnan(idx)) = 1; % Temporary fallback for NaNs
                uy_rgb_flat = reshape(ind2rgb(idx, cmap), [], 3);
                
                % Define alpha mask (0.55 where valid, 0 where NaN)
                alpha_mask = 0.55 * double(~isnan(uy_data));
                alpha_3d = repmat(alpha_mask, [1, 3]);
                
                % Agglomerate mathematically
                SE_rgb_flat = reshape(SE_rgb, [], 3);
                final_rgb = (1 - alpha_3d) .* SE_rgb_flat + alpha_3d .* uy_rgb_flat;
                
                fig2 = figure('Visible', 'off');
                plot(newEBSD, final_rgb); axis image off; hold on;
                
                % Fake the colorbar scaling for the legend
                colormap('parula'); caxis([uy_min, uy_max]);
                title(sprintf('Step %s: Deformed EBSD with Forward-Mapped Displacement Y', defStep.stepName), 'Interpreter', 'none');
                mtexColorbar('title', 'Forward Displacement y (Eulerian)');
                
                exportScaledFigure(fig2, fullfile(input.outPath, [defStep.stepName, '_ForwardDIC_Overlay_.png']), 'nodatestring');
                close(fig2);
            end
        end
    end

    %% 4-Panel Coordinate System Visual Check
    coord_x_fields = find(contains(fieldnames(newEBSD.prop), 'Coordinates_x', 'IgnoreCase', true));
    coord_y_fields = find(contains(fieldnames(newEBSD.prop), 'Coordinates_y', 'IgnoreCase', true));
    all_prop = fieldnames(newEBSD.prop);
    if ~isempty(coord_x_fields) && ~isempty(coord_y_fields)
        dic_x_f = all_prop{coord_x_fields(1)};
        dic_y_f = all_prop{coord_y_fields(1)};
        newEBSD.prop.DICX = newEBSD.prop.(dic_x_f);
        newEBSD.prop.DICY = newEBSD.prop.(dic_y_f);
        
        grid_sz = size(newEBSD);
        x_grid = reshape(newEBSD.x, grid_sz);
        y_grid = reshape(newEBSD.y, grid_sz);
        x_ext = [nanmean(x_grid(:, 1)), nanmean(x_grid(:, end))];
        y_ext = [nanmean(y_grid(1, :)), nanmean(y_grid(end, :))];
        
        fig_check = figure('Visible', 'off', 'Position', [100, 100, 1200, 900]);
        subplot(2, 2, 1); imagesc(x_ext, y_ext, x_grid); title('MTEX ebsd.x (\mum)'); axis image off; colormap(gca, 'parula'); cb1 = colorbar; cb1.Title.String = 'X (\mum)';
        subplot(2, 2, 2); imagesc(x_ext, y_ext, y_grid); title('MTEX ebsd.y (\mum)'); axis image off; colormap(gca, 'parula'); cb2 = colorbar; cb2.Title.String = 'Y (\mum)';
        subplot(2, 2, 3); imagesc(x_ext, y_ext, reshape(newEBSD.prop.DICX, grid_sz)); title(['DaVis DICX: ' strrep(dic_x_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb3 = colorbar; cb3.Title.String = 'DIC X';
        subplot(2, 2, 4); imagesc(x_ext, y_ext, reshape(newEBSD.prop.DICY, grid_sz)); title(['DaVis DICY: ' strrep(dic_y_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb4 = colorbar; cb4.Title.String = 'DIC Y';
        
        for sp_i = 1:4
            subplot(2, 2, sp_i); hold on;
            quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
            quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
            hold off;
        end
        
        checkOutName = fullfile(input.outPath, [defStep.stepName, '_CoordSystem_4PanelCheck_.png']);
        exportScaledFigure(fig_check, checkOutName, 'nodatestring');
        close(fig_check);
        fprintf('  -> Saved 4-Panel Coordinate System Check figure to: %s\n', checkOutName);
    end

    %% 4-Panel Displacement Gradient System Visual Check
    all_prop = fieldnames(newEBSD.prop);
    dudx_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
    dudy_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));
    dvdx_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
    dvdy_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));

    if ~isempty(dudx_idx) && ~isempty(dudy_idx) && ~isempty(dvdx_idx) && ~isempty(dvdy_idx)
        dudx_f = all_prop{dudx_idx(end)};
        dudy_f = all_prop{dudy_idx(end)};
        dvdx_f = all_prop{dvdx_idx(end)};
        dvdy_f = all_prop{dvdy_idx(end)};
        
        grid_sz = size(newEBSD);
        x_grid = reshape(newEBSD.x, grid_sz);
        y_grid = reshape(newEBSD.y, grid_sz);
        x_ext = [nanmean(x_grid(:, 1)), nanmean(x_grid(:, end))];
        y_ext = [nanmean(y_grid(1, :)), nanmean(y_grid(end, :))];
        
        fig_grad = figure('Visible', 'off', 'Position', [150, 150, 1200, 900]);
        subplot(2, 2, 1); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudx_f), grid_sz)); title(['du/dx: ' strrep(dudx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb1 = colorbar; cb1.Title.String = 'du/dx';
        subplot(2, 2, 2); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudy_f), grid_sz)); title(['du/dy: ' strrep(dudy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb2 = colorbar; cb2.Title.String = 'du/dy';
        subplot(2, 2, 3); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdx_f), grid_sz)); title(['dv/dx: ' strrep(dvdx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb3 = colorbar; cb3.Title.String = 'dv/dx';
        subplot(2, 2, 4); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdy_f), grid_sz)); title(['dv/dy: ' strrep(dvdy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb4 = colorbar; cb4.Title.String = 'dv/dy';
        
        for sp_i = 1:4
            subplot(2, 2, sp_i); hold on;
            quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
            quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
            hold off;
        end
        
        gradOutName = fullfile(input.outPath, [defStep.stepName, '_DispGradients_4PanelCheck_.png']);
        exportScaledFigure(fig_grad, gradOutName, 'nodatestring');
        close(fig_grad);
        fprintf('  -> Saved 4-Panel Displacement Gradient Check figure to: %s\n', gradOutName);
    end
    finalOutName = fullfile(input.outPath, [defStep.stepName, '_deformedEBSD_withForwardDIC.mat']);
    save(finalOutName, 'newEBSD', 'originalEBSD', 'grains', 'options', 'defStep');
    fprintf('  -> Saved complete FZ-projected dataset & grain object to: %s\n', finalOutName);
end
fprintf('\n========================================================================\n');
fprintf('ALL DEFORMED FORWARD DATASETS COMPUTED AND SAVED SUCCESSFULLY!\n');
fprintf('========================================================================\n');

%% Local Support Functions
function allPaths = collectH5Paths(infoStruct)
    allPaths = {};
    if isfield(infoStruct, 'Datasets')
        for i = 1:numel(infoStruct.Datasets)
            p = [infoStruct.Name, '/', infoStruct.Datasets(i).Name];
            p = strrep(p, '//', '/');
            allPaths{end+1} = p;
        end
    end
    if isfield(infoStruct, 'Groups')
        for i = 1:numel(infoStruct.Groups)
            grp = infoStruct.Groups(i);
            allPaths = [allPaths, collectH5Paths(grp)]; %#ok<AGROW>
        end
    end
end

function fullPath = findDatasetPath(allPaths, requiredRelPath)
    req = strsplit(requiredRelPath, '/');
    reqTail = req(max(end-1,1):end);
    fullPath = '';
    for i = 1:numel(allPaths)
        p = allPaths{i};
        parts = strsplit(p, '/');
        if length(parts) >= length(reqTail)
            if isequal(parts(end-numel(reqTail)+1:end), reqTail)
                fullPath = p;
                return;
            end
        end
    end
end
