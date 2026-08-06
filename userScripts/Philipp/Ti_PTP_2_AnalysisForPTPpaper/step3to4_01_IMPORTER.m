%% step3to4_01_IMPORTER.m
% Importer & alignment script specifically correlating TKD Step 03 (BW158)
% with the standalone Step 3 to 4 DIC analysis file.
% Outputs directly to a dedicated subfolder: Processed_Output_Step3to4

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\NanoMech_Aignment_Matlab\src');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Load Master Configuration & Locate Step 03 TKD Dataset
[registry, defaults] = a00_Dataset_Registry();

% Locate Step_03_TKD (insitu_new Specimen 1 03 Map Data 7_BW158) in registry
stepIdx = find(strcmpi({registry.stepName}, 'Step_03_TKD') | ...
               contains({registry.tkd_fname}, '03 Map Data 7_BW158', 'IgnoreCase', true), 1, 'first');

if isempty(stepIdx)
    error('Could not find Step_03_TKD (BW158 Map Data 7) in a00_Dataset_Registry.');
end
stepInfo = registry(stepIdx);

fprintf('--- Running Step 3 to 4 Importer on TKD Scan: %s ---\n', stepInfo.tkd_fname);

%% 2. Setup Input Parameters, Override DIC Path & Dedicated Output Folder
input.fname             = stepInfo.tkd_fname;
input.SEImagefname      = stepInfo.se_fname;
input.mask_fname        = stepInfo.mask_fname;
input.dir               = defaults.dir;
input.SEImagefolder     = defaults.SEImagefolder;
input.DIC_filepath      = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\DIC exported data\FromImages\Processing_subset=19_step=6_relativeToFirst_median3\OnlyStep3to4_correspondsTo03and03tkdtifs0001.hdf5';

% Route outputs specifically to the Step 3 to 4 subfolder
input.outPath = fullfile(defaults.dir, 'Processed_Output_Step3to4');
if ~exist(input.outPath, 'dir'), mkdir(input.outPath); end

input.SEImageFullPath = fullfile(input.SEImagefolder, input.SEImagefname);
input.ang_path        = fullfile(input.dir, [input.fname, '.ang']);
input.mask_path       = fullfile(input.dir, input.mask_fname);

options.dateString = datestr(now,'yyyymmdd_HHMMSS');
options.in_situIncrementName = [input.fname, '_Step3to4_', options.dateString];

outname = fullfile(input.outPath, sprintf('%s_ioConfig.json', options.in_situIncrementName));
jsonText = jsonencode(input);
fid = fopen(outname,'w');
if fid > 0
    fprintf(fid, '%s', jsonText);
    fclose(fid);
end

%% 3. EBSD Data Loading & Filtering
if ~exist(input.ang_path, 'file')
    error('EBSD file not found: %s', input.ang_path);
end

CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
ebsd = EBSD.load(input.ang_path, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

% Load optional extra properties from .h5oina if extracted
prop_file = fullfile(input.dir, 'Masks', sprintf('%s_h5oina_props.mat', input.fname));
if ~exist(prop_file, 'file')
    prop_file = strrep(input.mask_path, '_specimen_mask.mat', '_h5oina_props.mat');
end
if exist(prop_file, 'file')
    fprintf('Loading extra properties from %s...\n', prop_file);
    try
        props = load(prop_file);
        if isfield(props, 'BS') && length(props.BS) == length(ebsd), ebsd.prop.BS = props.BS; end
        if isfield(props, 'BC') && length(props.BC) == length(ebsd), ebsd.prop.BC = props.BC; end
        if isfield(props, 'FSE') && length(props.FSE) == length(ebsd), ebsd.prop.FSE = props.FSE; end
        if isfield(props, 'BSE') && length(props.BSE) == length(ebsd), ebsd.prop.BSE = props.BSE; end
        if isfield(props, 'MSE') && length(props.MSE) == length(ebsd), ebsd.prop.MSE = props.MSE; end
        if isfield(props, 'LapVar') && length(props.LapVar) == length(ebsd), ebsd.prop.LapVar = props.LapVar; end
        if isfield(props, 'GradMag') && length(props.GradMag) == length(ebsd), ebsd.prop.GradMag = props.GradMag; end
        fprintf('  -> Successfully loaded optional H5OINA properties.\n');
    catch
        warning('Failed to load extra properties from %s', prop_file);
    end
end

ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + defaults.tkd_tilt), 'keepXY');
ebsd(ebsd.ci < defaults.thresh_ci) = 'notIndexed';
ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
refOrientation = orientation.byEuler(defaults.refOrientationEulerAngles, ebsd.CS);
ebsd.orientations = ebsd.orientations.project2FundamentalRegion(refOrientation);
ebsd = ebsd.gridify;

% Apply Specimen Mask
if exist(input.mask_path, 'file')
    fprintf('Applying specimen mask: %s\n', input.mask_path);
    load(input.mask_path, 'specimenlocation');
    try
        ebsd(~specimenlocation) = 'notIndexed';
    catch
        warning('Mask dimension mismatch, applying directly to 1D index.');
        ebsd(~specimenlocation) = 'notIndexed';
    end
else
    warning('No specimen mask found at %s', input.mask_path);
end

%% 4. Initial Quality Plot (IPF Y)
ipf = ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection = yvector;
figure('Visible', 'off');
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'), colors);
title('Step 03 (BW158) TKD - IPF Y');
exportScaledFigure(gcf, fullfile(input.outPath, [input.fname, '_IPF_Y.png']), 'nodatestring');

%% 5. Align TKD to SE Image
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

REF = 'SE_raw';
DEF = 'EBSD0';
options.Ndof = 6;
options.Nsel = 4;
options.case = 'EBSD';
options.selectionIPF = 1;
options.EBSD_selection = 'LapVar';
options.colormap = 'gray';

% Check multiple candidates for existing alignment points to prevent redundant interactive clicking
candPaths = {
    fullfile(input.outPath, 'Step_03_TKD_alignment_data_deformed.mat'), ...
    fullfile(defaults.dir, defaults.outDir_deformed, 'Step_03_TKD_alignment_data_deformed.mat'), ...
    fullfile(defaults.dir, stepInfo.sel_points_fname), ...
    stepInfo.sel_points_fname ...
};

loaded_alignment = false;
for p = 1:numel(candPaths)
    if ~isempty(candPaths{p}) && exist(candPaths{p}, 'file')
        fprintf('Found existing alignment data at: %s\nLoading alignment data...\n', candPaths{p});
        tmp_al = load(candPaths{p});
        if isfield(tmp_al, 'alignment_data_EBSD')
            alignment_data_EBSD = tmp_al.alignment_data_EBSD;
        elseif isfield(tmp_al, 'alignment_data_deformed')
            alignment_data_EBSD = tmp_al.alignment_data_deformed;
        end
        loaded_alignment = true;
        break;
    end
end

if ~loaded_alignment
    fprintf('No prior alignment points found. Starting interactive point selection...\n');
    options.sel_points_fname = 0;
    alignment_data_EBSD = pointSelectionForAlignment(rawData.(REF), rawData.(DEF), options);
    alignment_data_deformed = alignment_data_EBSD;
    save(fullfile(input.outPath, 'Step_03_TKD_alignment_data_deformed.mat'), 'alignment_data_deformed', 'alignment_data_EBSD');
    fprintf('Saved newly collected alignment points to: %s\n', fullfile(input.outPath, 'Step_03_TKD_alignment_data_deformed.mat'));
end

psize_new = rawData.SE_raw.psize;
[Xnew, Ynew] = meshgrid(0:psize_new:max(rawData.SE_raw.X(:)), 0:psize_new:max(rawData.SE_raw.Y(:)));
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data_EBSD);
newEBSD.prop.originalId = newEBSD.prop.oldId;
originalEBSD = ebsd;
originalEBSD.prop.originalId = ebsd.id;
newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');
newEBSD = newEBSD.gridify;

%% 6. Verification Plot: Aligned EBSD over SE Image
newGrains = calcGrains(newEBSD, 'angle', 5*degree);
figure('Visible', 'off');
plot(newEBSD, newEBSD.prop.SE_raw, 'micronbar', 'off'); hold on;
plot(newGrains.boundary, 'lineColor', 'red', 'linewidth', 1.5);
title('Aligned Step 03 EBSD boundaries over SE image'); hold off;
exportScaledFigure(gcf, fullfile(input.outPath, [input.fname, '_AlignedEBSD_SE.png']), 'nodatestring');

%% 7. DIC Data Loading & Import onto Aligned EBSD Grid
dicFilePath = input.DIC_filepath;
if ~exist(dicFilePath, 'file')
    error('DIC HDF5 file not found at: %s', dicFilePath);
end
info = h5info(dicFilePath);
groups = info.Groups;

targetM = rawData.SE_raw.origM;
targetN = rawData.SE_raw.origN;
targetSize = [targetM, targetN];

requiredDatasetPaths = {
    'Coordinates/x',
    'Coordinates/y',
    'Displacements/x-displacement',
    'Displacements/y-displacement',
    'Gradients/du÷dx - RC',
    'Gradients/du÷dy - RC',
    'Gradients/dv÷dx - RC',
    'Gradients/dv÷dy - RC',
    'Strains/Axial strain Ɛyy - RC',
    'Strains/Shear strain Ɛxy - RC',
    'Strains/Transverse strain Ɛxx - RC'
};

rawData.DIC = struct();
rawData.DIC.names = {};

firstGrpPath = groups(1).Name;
allPaths = collectH5Paths(info);
resolvedPaths = struct();
for i = 1:numel(requiredDatasetPaths)
    rel = requiredDatasetPaths{i};
    full_p = findDatasetPath(allPaths, rel);
    if isempty(full_p)
        full_p = [firstGrpPath, '/', rel];
    end
    resolvedPaths.(matlab.lang.makeValidName(rel)) = full_p;
end

xcoord = h5read(dicFilePath, resolvedPaths.Coordinates_x);
ycoord = h5read(dicFilePath, resolvedPaths.Coordinates_y);
dataSize = size(xcoord);
tmp_v = abs(diff(xcoord, 1, 1)); var_v = nanmedian(tmp_v(:));
tmp_h = abs(diff(xcoord, 1, 2)); var_h = nanmedian(tmp_h(:));
isTransposed = (var_v > var_h);
if isTransposed
    xcoord = xcoord';
    ycoord = ycoord';
end

stepsize = abs(xcoord(2,1)-xcoord(1,1));
if stepsize == 0, stepsize = abs(xcoord(1,2)-xcoord(1,1)); end
if stepsize == 0, stepsize = 1; end

x0 = rawData.SE_raw.imregion(1);
y0 = rawData.SE_raw.imregion(2);
x1 = rawData.SE_raw.imregion(3) + x0;
y1 = rawData.SE_raw.imregion(4) + y0;

for g = 1:numel(groups)
    grpName = groups(g).Name;
    baseName = strrep(grpName, '/', '');
    fprintf('  -> Loading DIC Group from Step 3 to 4 HDF5: %s\n', baseName);
    
    for k = 1:numel(requiredDatasetPaths)
        rel = requiredDatasetPaths{k};
        fullPath = [grpName, '/', rel];
        
        try
            data = h5read(dicFilePath, fullPath);
        catch
            warning('Could not read dataset: %s in %s', fullPath, dicFilePath);
            continue;
        end
        
        if size(data, 1) ~= targetM && size(data, 2) == targetM
            data = data';
        end
        if isTransposed, data = data'; end
        
        if isequal(size(data), targetSize)
            subdata = extractRegion(data, x0, y0, x1, y1, rel);
        else
            if contains(rel, 'Coordinates')
                if data(1,1) < data(2,1)
                    subdata = repmat(1:targetN, targetM, 1);
                elseif data(1,1) < data(1,2)
                    subdata = repmat((1:targetM)', 1, targetN);
                else
                    subdata = kron(data, ones(stepsize));
                end
                subdata = extractRegion(subdata, x0, y0, x1, y1, rel);
            else
                subdata = kron(data, ones(stepsize));
                subdata = extractRegion(subdata, x0, y0, x1, y1, rel);
            end
        end
        
        dsSuffix = strjoin(strsplit(rel, '/'), '_');
        shortName = matlab.lang.makeValidName([baseName '_' dsSuffix]);
        rawData.DIC.(shortName) = subdata;
        rawData.DIC.names{end+1} = shortName;
    end
end

% Interpolate each DIC dataset directly onto aligned EBSD grid
dicNames = rawData.DIC.names;
if ~isempty(dicNames)
    for i = 1:numel(dicNames)
        origName = dicNames{i};
        dicField.X = rawData.SE_raw.X;
        dicField.Y = rawData.SE_raw.Y;
        dicField.data = double(rawData.DIC.(origName));
        newEBSD = warpIncludeDataInMtex(newEBSD, dicField, origName);
    end
end

%% 8. Check & Apply Coordinate System Corrections via 2D Grid Dominance Analysis
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
        warning('DIC_x varies vertically rather than horizontally!');
    else
        fprintf('Grid Dominance check confirms DaVis X and Y axes match EBSD X and Y spatial matrix orientation.\n');
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
        fprintf('Coordinate directions increase identically (s_x = +1, s_y = +1). No negative sign flips needed.\n');
    end
end

%% 9. Grain Calculation & FZ Projection
grainThresh = 5*degree; 
minSize = 300; 
alpha = 1;

[grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

%% 10. 4-Panel Coordinate & Displacement Gradient System Visual Checks
if ~isempty(coord_x_fields) && ~isempty(coord_y_fields)
    dic_x_f = all_prop_names{coord_x_fields(1)};
    dic_y_f = all_prop_names{coord_y_fields(1)};
    newEBSD.prop.DICX = newEBSD.prop.(dic_x_f);
    newEBSD.prop.DICY = newEBSD.prop.(dic_y_f);
    
    grid_sz = size(newEBSD);
    x_grid = reshape(newEBSD.x, grid_sz);
    y_grid = reshape(newEBSD.y, grid_sz);
    x_ext = [nanmean(x_grid(:, 1)), nanmean(x_grid(:, end))];
    y_ext = [nanmean(y_grid(1, :)), nanmean(y_grid(end, :))];
    
    fig_check = figure('Name', 'Coordinate System Verification', 'Position', [100, 100, 1200, 900], 'Visible', 'off');
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
    
    checkOutName = fullfile(input.outPath, [input.fname, '_CoordSystem_4PanelCheck.png']);
    exportScaledFigure(fig_check, checkOutName, 'nodatestring');
    close(fig_check);
end

% Displacement gradients check
all_prop = fieldnames(newEBSD.prop);
dudx_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
dudy_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));
dvdx_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
dvdy_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));

if ~isempty(dudx_idx) && ~isempty(dudy_idx) && ~isempty(dvdx_idx) && ~isempty(dvdy_idx)
    dudx_f = all_prop{dudx_idx(1)};
    dudy_f = all_prop{dudy_idx(1)};
    dvdx_f = all_prop{dvdx_idx(1)};
    dvdy_f = all_prop{dvdy_idx(1)};
    
    grid_sz = size(newEBSD);
    x_grid = reshape(newEBSD.x, grid_sz);
    y_grid = reshape(newEBSD.y, grid_sz);
    x_ext = [nanmean(x_grid(:, 1)), nanmean(x_grid(:, end))];
    y_ext = [nanmean(y_grid(1, :)), nanmean(y_grid(end, :))];

    dudx_data = newEBSD.prop.(dudx_f);
    dudy_data = newEBSD.prop.(dudy_f);
    dvdx_data = newEBSD.prop.(dvdx_f);
    dvdy_data = newEBSD.prop.(dvdy_f);
    
    all_vals = [dudx_data(:); dudy_data(:); dvdx_data(:); dvdy_data(:)];
    all_vals = all_vals(~isnan(all_vals));
    cmin = 0.2 * min(all_vals);
    cmax = 0.2 * max(all_vals);
    
    fig_grad = figure('Name', 'Displacement Gradient Verification', 'Position', [150, 150, 1200, 900], 'Visible', 'off');
    subplot(2, 2, 1); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudx_f), grid_sz)); title(['du/dx: ' strrep(dudx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); colorbar; caxis([cmin cmax]);
    subplot(2, 2, 2); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudy_f), grid_sz)); title(['du/dy: ' strrep(dudy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); colorbar; caxis([cmin cmax]);
    subplot(2, 2, 3); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdx_f), grid_sz)); title(['dv/dx: ' strrep(dvdx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); colorbar; caxis([cmin cmax]);
    subplot(2, 2, 4); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdy_f), grid_sz)); title(['dv/dy: ' strrep(dvdy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); colorbar; caxis([cmin cmax]);
    for sp_i = 1:4
        subplot(2, 2, sp_i); hold on;
        quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
        quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
        hold off;
    end
    
    gradOutName = fullfile(input.outPath, [input.fname, '_DispGradients_4PanelCheck.png']);
    exportScaledFigure(fig_grad, gradOutName, 'nodatestring');
    close(fig_grad);
end

%% 11. Save Aligned Dataset
finalOutName = fullfile(input.outPath, 'step3to4_EBSD_withDIC.mat');
save(finalOutName, 'newEBSD', 'originalEBSD', 'grains', 'options');
fprintf('\n========================================================================\n');
fprintf('STEP 3 to 4 IMPORT & ALIGNMENT COMPLETE!\n');
fprintf('Saved fully aligned dataset & FZ grains to:\n  %s\n', finalOutName);
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
