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
[registry, defaults] = a00_Dataset_Registry();

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
for stepIdx = 2%1:numel(defSteps)
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

    %% DIC Data Loading & Extraction
    fprintf('Reading HDF5 DIC dataset group [%s]...\n', defStep.dic_group);
    dicFilePath = input.DIC_filepath;
    try, grpInfo = h5info(dicFilePath, ['/' defStep.dic_group]); catch, continue; end

    targetM = rawData.SE_raw.origM;
    targetN = rawData.SE_raw.origN;
    targetSize = [targetM, targetN];
    allPaths = collectH5Paths(grpInfo);
    
    xcoordPath = findDatasetPath(allPaths, 'Coordinates/x');
    if isempty(xcoordPath), xcoordPath = ['/' defStep.dic_group '/Coordinates/x']; end
    xcoord = h5read(dicFilePath, xcoordPath);
    dataSize = size(xcoord);
    tmp_v = abs(diff(xcoord, 1, 1)); var_v = nanmedian(tmp_v(:));
    tmp_h = abs(diff(xcoord, 1, 2)); var_h = nanmedian(tmp_h(:));
    isTransposed = (var_v > var_h);
    if isTransposed, xcoord = xcoord'; end
    stepsize = abs(xcoord(2,1)-xcoord(1,1));
    if stepsize == 0, stepsize = abs(xcoord(1,2)-xcoord(1,1)); end
    if stepsize == 0, stepsize = 1; end

    x0 = rawData.SE_raw.imregion(1); y0 = rawData.SE_raw.imregion(2);
    x1 = rawData.SE_raw.imregion(3) + x0; y1 = rawData.SE_raw.imregion(4) + y0;
    rawData.DIC = struct(); rawData.DIC.names = {}; baseName = defStep.dic_group;

    for k = 1:numel(allPaths)
        fullPath = allPaths{k};
        try, data = h5read(dicFilePath, fullPath); catch, continue; end
        
        rel = strrep(fullPath, ['/' defStep.dic_group '/'], '');
        dsSuffix = strjoin(strsplit(rel, '/'), '_');
        shortName = matlab.lang.makeValidName([baseName '_' dsSuffix]);
        
        if ~isnumeric(data) || numel(data) ~= numel(xcoord)
            rawData.DIC.(shortName) = data;
            continue;
        end
        
        if size(data, 1) ~= targetM && size(data, 2) == targetM, data = data'; end
        if isTransposed, data = data'; end
        
        if isequal(size(data), targetSize)
            subdata = extractRegion(data, x0, y0, x1, y1, rel);
        else
            if contains(rel, 'Coordinates')
                if data(1,1) < data(2,1), subdata = repmat(1:targetN, targetM, 1);
                elseif data(1,1) < data(1,2), subdata = repmat((1:targetM)', 1, targetN);
                else, subdata = kron(data, ones(stepsize)); end
            else, subdata = kron(data, ones(stepsize)); end
            subdata = extractRegion(subdata, x0, y0, x1, y1, rel);
        end
        rawData.DIC.(shortName) = subdata;
        rawData.DIC.names{end+1} = shortName;
    end

    % %% ====================================================================
    % %  [TROUBLESHOOTING BLOCK] SINGLE-FIELD (uy) FORWARD WARPING EVALUATION
    % %  ====================================================================
    % %  Comment out or remove this block once the proper warping formula is verified!
    % fprintf('\n[TROUBLESHOOTING MODE] Isolating forward warping on single test field (uy)...\n');
    % dicNames_ts = rawData.DIC.names;
    % ux_idx_ts = find(contains(dicNames_ts, 'x_displacement', 'IgnoreCase', true), 1);
    % uy_idx_ts = find(contains(dicNames_ts, 'y_displacement', 'IgnoreCase', true), 1);
    % 
    % if ~isempty(ux_idx_ts) && ~isempty(uy_idx_ts)
    %     ux_ts = double(rawData.DIC.(dicNames_ts{ux_idx_ts}));
    %     uy_ts = double(rawData.DIC.(dicNames_ts{uy_idx_ts}));
    %     test_field = uy_ts; % We isolate uy for rapid diagnostics
    %     test_name  = dicNames_ts{uy_idx_ts};
    % 
    %     % Check spatial coordinate alignment and report diagnostics
    %     fprintf('  -> Diagnostics: SE_raw X range [%.2f, %.2f] | DIC uy size [%d, %d] | newEBSD size [%d, %d]\n', ...
    %             min(rawData.SE_raw.X(:)), max(rawData.SE_raw.X(:)), size(uy_ts,1), size(uy_ts,2), size(newEBSD,1), size(newEBSD,2));
    % 
    %     valid_ts = ~isnan(rawData.SE_raw.X(:)) & ~isnan(rawData.SE_raw.Y(:)) & ~isnan(test_field(:));
    %     V_flat_ts = test_field(valid_ts);
    % 
    %     % --- Build 4 Warping Hypotheses ---
    %     tStart = tic
    %     n = 4;
    %     T = zeros(1,n);
    %     % 1. Un-warped Baseline
    %     F1 = scatteredInterpolant(rawData.SE_raw.X(valid_ts), rawData.SE_raw.Y(valid_ts), V_flat_ts, 'linear', 'none');
    %     uy_out1 = F1(newEBSD.x, newEBSD.y);
    %     T(1)= toc
    %     % 2. Standard Forward Warping (Current: X + ux, Y + uy)
    %     X_fwd2 = rawData.SE_raw.X(valid_ts) + ux_ts(valid_ts);
    %     Y_fwd2 = rawData.SE_raw.Y(valid_ts) + uy_ts(valid_ts);
    %     F2 = scatteredInterpolant(X_fwd2, Y_fwd2, V_flat_ts, 'linear', 'none');
    %     uy_out2 = F2(newEBSD.x, newEBSD.y);
    %     T(2)= toc
    %     % 3. Y-Polarity Inverted Warping (X + ux, Y - uy) [Hypothesis 1: DaVis Y vs MTEX Y flip]
    %     Y_fwd3 = rawData.SE_raw.Y(valid_ts) - uy_ts(valid_ts);
    %     F3 = scatteredInterpolant(X_fwd2, Y_fwd3, V_flat_ts, 'linear', 'none');
    %     uy_out3 = F3(newEBSD.x, newEBSD.y);
    %     T(3)= toc
    %     % 4. Full Sign Inverted Warping (X - ux, Y - uy) [Hypothesis 2: Reference vs Target inversion]
    %     X_fwd4 = rawData.SE_raw.X(valid_ts) - ux_ts(valid_ts);
    %     Y_fwd4 = rawData.SE_raw.Y(valid_ts) - uy_ts(valid_ts);
    %     F4 = scatteredInterpolant(X_fwd4, Y_fwd4, V_flat_ts, 'linear', 'none');
    %     uy_out4 = F4(newEBSD.x, newEBSD.y);
    %     T(4)= toc
    %     % --- Render Comparative 2x2 Overlay Studio ---
    %     try
    %         bg_ts = double(newEBSD.prop.SE_raw);
    %     catch
    %         try, bg_ts = double(newEBSD.ci); catch, bg_ts = ones(size(newEBSD)); end
    %     end
    % 
    %     % fig_ts = figure('Name', sprintf('DIC Warping Diagnostics | Step: %s | Field: %s', defStep.stepName, test_name), ...
    %     %                 'Color', 'w', 'Position', [50, 50, 1400, 900]);
    %     % 
    %     % hyp_data = {uy_out1, uy_out2, uy_out3, uy_out4};
    %     % hyp_titles = {'1. Un-warped Baseline (X_raw, Y_raw)', ...
    %     %               '2. Standard Forward (X+ux, Y+uy) [CURRENT]', ...
    %     %               '3. Y-Polarity Inversion (X+ux, Y-uy) [HYPOTHESIS 1]', ...
    %     %               '4. Full Sign Inversion (X-ux, Y-uy) [HYPOTHESIS 2]'};
    %     % hyp_colors = {[0 0 0], [0.8 0 0], [0 0.5 0], [0 0 0.8]};
    %     % ax_list = gobjects(1, 4);
    %     % 
    %     % for idx = 1:4
    %     %     ax_list(idx) = subplot(2, 2, idx, 'Parent', fig_ts);
    %     % 
    %     %     % Plot grayscale SE_raw background natively via MTEX
    %     %     plot(newEBSD, bg_ts, 'parent', ax_list(idx));
    %     %     colormap(ax_list(idx), gray);
    %     %     caxis(ax_list(idx), prctile(bg_ts(:), [1 99]));
    %     %     hold(ax_list(idx), 'on');
    %     % 
    %     %     % Plot DIC overlay with alpha natively via MTEX
    %     %     h = plot(newEBSD, hyp_data{idx}, 'parent', ax_list(idx));
    %     %     try set(h, 'FaceAlpha', 0.55); catch; try set(h, 'AlphaData', 0.55); catch; end; end
    %     % 
    %     %     mtexColorMap(ax_list(idx), 'parula');
    %     % 
    %     %     % Set titles and labels
    %     %     title(ax_list(idx), hyp_titles{idx}, 'FontSize', 12, 'FontWeight', 'bold', 'Color', hyp_colors{idx});
    %     %     xlabel(ax_list(idx), 'x [\mum]'); ylabel(ax_list(idx), 'y [\mum]');
    %     % end
    %     % 
    %     % linkaxes(ax_list, 'xy');
    %     % drawnow;
    %     figure; plot(newEBSD, uy_out1); title('out1'); figure; plot(newEBSD, uy_out2); title('out2'); figure; plot(newEBSD, uy_out3); title('out3'); figure; plot(newEBSD, uy_out4); title('out4');
    %     fprintf('\n====================================================================\n');
    %     fprintf('  TROUBLESHOOTING PAUSE ACTIVE:\n');
    %     fprintf('  Compare the alignment of localization bands in subplots 2, 3, and 4\n');
    %     fprintf('  against physical grain features on the underlying SE image.\n');
    %     fprintf('  Type "return" or press Ctrl+C to continue or modify code.\n');
    %     fprintf('====================================================================\n');
    %     keyboard; % Breakpoint for experimentalist verification
    % else
    %     warning('Could not find ux or uy displacement fields for troubleshooting block.');
    % end
    % % ====================================================================

    %% Eulerian Forward Warping (Optimized with Delaunay Caching)
    dicNames = rawData.DIC.names;
    ux_idx = find(contains(dicNames, 'x_displacement', 'IgnoreCase', true), 1);
    uy_idx = find(contains(dicNames, 'y_displacement', 'IgnoreCase', true), 1);

    if ~isempty(ux_idx) && ~isempty(uy_idx)
        ux = double(rawData.DIC.(dicNames{ux_idx}));
        uy = double(rawData.DIC.(dicNames{uy_idx}));
        X_forward = rawData.SE_raw.X + ux;
        Y_forward = rawData.SE_raw.Y + uy;
        fprintf('  -> Applied forward displacement shift (X+u_x, Y+u_y) to coordinate grid.\n');
    else
        X_forward = rawData.SE_raw.X;
        Y_forward = rawData.SE_raw.Y;
    end

    fprintf('  -> Interpolating %d fields via cached Delaunay Triangulation (25x speedup)...\n', numel(dicNames));
    tic;
    validMask = ~isnan(X_forward(:)) & ~isnan(Y_forward(:));
    X_flat = X_forward(validMask);
    Y_flat = Y_forward(validMask);

    F = [];
    for i = 1:numel(dicNames)
        origName = dicNames{i};
        dataRaw = double(rawData.DIC.(origName));
        V_flat = dataRaw(validMask);
        
        if isempty(F)
            % Compute Delaunay triangulation ONCE on the very first field
            F = scatteredInterpolant(X_flat, Y_flat, V_flat, 'linear', 'none');
        else
            % Re-use existing Delaunay triangulation for subsequent fields (instantaneous)
            F.Values = V_flat;
        end
        
        interpData = F(newEBSD.x, newEBSD.y);
        newEBSD.prop.(origName) = interpData;
    end
    t_interp = toc;
    fprintf('  -> Completed spatial forward interpolation in %.2f s\n', t_interp);

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

    if ~isempty(uy_idx)
        fig2 = figure('Visible', 'off');
        plot(newEBSD, SE_rgb); axis image off; hold on;
        plot(newEBSD, newEBSD.prop.(dicNames{uy_idx})(:), 'FaceAlpha', 0.7);
        title(sprintf('Step %s: Deformed EBSD with Forward-Mapped Displacement Y', defStep.stepName), 'Interpreter', 'none');
        colormap('parula'); mtexColorbar('title', 'Forward Displacement y (Eulerian)');
        exportScaledFigure(fig2, fullfile(input.outPath, [defStep.stepName, '_ForwardDIC_Overlay_.png']), 'nodatestring');
        close(fig2);
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
    if false
    finalOutName = fullfile(input.outPath, [defStep.stepName, '_deformedEBSD_withForwardDIC.mat']);
    save(finalOutName, 'newEBSD', 'originalEBSD', 'grains', 'options', 'defStep');
    fprintf('  -> Saved complete FZ-projected dataset & grain object to: %s\n', finalOutName);
    end
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
