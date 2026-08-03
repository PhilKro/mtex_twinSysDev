%% undeformedReference_03_IMPORTER.m
% Aligns the Undeformed Reference TKD/EBSD dataset to the Reference SE image,
% and imports ALL DIC HDF5 datasets directly onto this undeformed grid without forward extrapolation.
% Outputs to a dedicated Undeformed directory: Processed_Output_Undeformed.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\NanoMech_Aignment_Matlab\src');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% Load Master Configuration & Dataset Registry
[registry, defaults] = a00_Dataset_Registry();

% Locate the reference dataset in the registry (ignoring any marked as excluded)
refIdx = find([registry.isReference] & ~[registry.exclude], 1, 'first');
if isempty(refIdx)
    error('No non-excluded reference step found in a00_Dataset_Registry (where isReference == true and exclude == false).');
end
refStep = registry(refIdx);
fprintf('--- Running Undeformed Reference Importer on Step: %s ---\n', refStep.stepName);

%% Setup Input Parameters & Paths
input.fname             = refStep.tkd_fname;
input.SEImagefname      = refStep.se_fname;
input.mask_fname        = refStep.mask_fname;
input.sel_points_fname  = refStep.sel_points_fname;
input.dir               = defaults.dir;
input.SEImagefolder     = defaults.SEImagefolder;
input.DIC_filepath      = defaults.DIC_filepath;

% Route outputs specifically to the undeformed folder
input.outPath = fullfile(defaults.dir, defaults.outDir_undeformed);
if ~exist(input.outPath, 'dir'), mkdir(input.outPath); end

input.SEImageFullPath = fullfile(input.SEImagefolder, input.SEImagefname);
input.ang_path        = fullfile(input.dir, [input.fname, '.ang']);
input.mask_path       = fullfile(input.dir, input.mask_fname);

options.dateString = datestr(now,'yyyymmdd_HHMMSS');
options.in_situIncrementName = [input.fname,'_',options.dateString];

outname = fullfile(input.outPath, sprintf('%s_undeformed_ioConfig.json', options.in_situIncrementName));
jsonText = jsonencode(input);
fid = fopen(outname,'w');
if fid > 0
    fprintf(fid, '%s', jsonText);
    fclose(fid);
end

%% EBSD Data Loading
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
ebsd = EBSD.load(input.ang_path, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

% Load optional extra properties from .h5oina if extracted
prop_file = fullfile(input.dir, 'Masks', sprintf('%s_h5oina_props.mat', input.fname));
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

% Apply Mask
if exist(input.mask_path, 'file')
    fprintf('Applying specimen mask...\n');
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

%% Initial Quality Plots
ipf = ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection = yvector;
figure;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'), colors);
title('Undeformed Reference IPF Y');
exportScaledFigure(gcf, fullfile(input.outPath, [input.fname, '_', options.dateString, '_ebsd_.png']), 'nodatestring');

%% Align Reference TKD to Reference SE Image
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

if ~isempty(input.sel_points_fname)
    % Resolve path if relative
    if ~exist(input.sel_points_fname, 'file') && exist(fullfile(input.dir, input.sel_points_fname), 'file')
        options.sel_points_fname = fullfile(input.dir, input.sel_points_fname);
    else
        options.sel_points_fname = input.sel_points_fname;
    end
end

alignment_data_EBSD = pointSelectionForAlignment(rawData.(REF), rawData.(DEF), options);
save(fullfile(input.outPath, [options.in_situIncrementName, '_', options.case, '_alignment_data_EBSD']), 'alignment_data_EBSD');

psize_new = rawData.SE_raw.psize;
[Xnew, Ynew] = meshgrid(0:psize_new:max(rawData.SE_raw.X(:)), 0:psize_new:max(rawData.SE_raw.Y(:)));
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data_EBSD);
newEBSD.prop.originalId = newEBSD.prop.oldId;
originalEBSD = ebsd;
originalEBSD.prop.originalId = ebsd.id;
newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');
newEBSD = newEBSD.gridify;

%% Verification Plotting
newGrains = calcGrains(newEBSD, 'angle', 5*degree);
figure;
plot(newEBSD, newEBSD.prop.SE_raw, 'micronbar', 'off');
hold on;
plot(newGrains.boundary, 'lineColor', 'red', 'linewidth', 1.5);
title('Aligned Reference EBSD boundaries over Undeformed SE image');
exportScaledFigure(gcf, fullfile(input.outPath, [input.fname, '_', options.dateString, '_AlignedEBSD_SE_.png']), 'nodatestring');

%% DIC Data Loading & Import onto Undeformed Reference Grid
dicFilePath = input.DIC_filepath;
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
    full = findDatasetPath(allPaths, rel);
    if isempty(full)
        full = [firstGrpPath, '/', rel];
    end
    resolvedPaths.(matlab.lang.makeValidName(rel)) = full;
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
    
    % Skip if this group is marked as excluded in the dataset registry (e.g. B00010)
    if any([registry.exclude] & strcmpi({registry.dic_group}, baseName))
        fprintf('  -> Skipping excluded DIC group: %s (as configured in registry)\n', baseName);
        continue;
    end
    
    for k = 1:numel(requiredDatasetPaths)
        rel = requiredDatasetPaths{k};
        fullPath = [grpName, '/', rel];
        
        try
            data = h5read(dicFilePath, fullPath);
        catch
            warning('Could not read %s', fullPath);
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

% Interpolate each DIC dataset directly onto undeformed reference TKD grid (without forward extrapolation)
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

%% Check & Apply Coordinate System Corrections via 2D Grid Dominance Analysis
% Determine whether DaVis DIC coordinates (/Coordinates/x and /Coordinates/y) are swapped (X <-> Y)
% or inverted in sign relative to the physical EBSD map coordinates (ebsd.x, ebsd.y).
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
    
    % Check dominant rate of change along rows (dim 1, vertical) vs columns (dim 2, horizontal)
    % ebsd.x varies along dim 2; ebsd.y varies along dim 1.
    tmp = abs(diff(DIC_x, 1, 2)); var_DICx_horz = nanmedian(tmp(:));
    tmp = abs(diff(DIC_x, 1, 1)); var_DICx_vert = nanmedian(tmp(:));
    
    if var_DICx_vert > var_DICx_horz
        warning('DIC_x varies vertically rather than horizontally! HDF5 read transposition check should have handled this.');
    else
        fprintf('Grid Dominance check confirms DaVis X and Y axes match EBSD X and Y spatial matrix orientation.\n');
    end
    
    % Evaluate direction signs after ensuring proper axis assignment
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

%% Grain Calculation & Fundamental Region (FZ) Projection
grainThresh = 5*degree; 
minSize = 300; 
alpha = 1;

[grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

% replaced by a unified FZ in the Dataset Registry
% fprintf('Projecting intra-grain orientations to Fundamental Zone (FZ) centered on grain mean orientation...\n');
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
% [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

%% Final Visualization & Saving
figure;
SE_gray = mat2gray(newEBSD.prop.SE_raw);
SE_rgb = repmat(SE_gray, 1, 1, 3);
plot(newEBSD, SE_rgb); 
axis image off; hold on;
last_y_disp = dicNames{end}; 
maxStrainTest = newEBSD.prop.(last_y_disp);
plot(newEBSD, maxStrainTest(:), 'FaceAlpha', 0.7);
title('Undeformed Reference TKD with Final Step Displacement Y');
colormap('parula');
mtexColorbar('title', 'Displacement y (Lagrangian)');

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
    
    fig_check = figure('Name', 'Coordinate System Verification', 'Position', [100, 100, 1200, 900]);
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
    
    checkOutName = fullfile(input.outPath, [input.fname, '_', options.dateString, '_CoordSystem_4PanelCheck_.png']);
    exportScaledFigure(fig_check, checkOutName, 'nodatestring');
    fprintf('Saved 4-Panel Coordinate System Check figure to: %s\n', checkOutName);
end

%% 4-Panel Displacement Gradient System Visual Check
plot_inc = 7

all_prop = fieldnames(newEBSD.prop);
dudx_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
dudy_idx = find(contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));
dvdx_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true) & ~contains(all_prop, 'dy', 'IgnoreCase', true));
dvdy_idx = find(contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true) & ~contains(all_prop, 'dx', 'IgnoreCase', true));

if ~isempty(dudx_idx) && ~isempty(dudy_idx) && ~isempty(dvdx_idx) && ~isempty(dvdy_idx)
    try
        dudx_f = all_prop{dudx_idx(plot_inc)};
        dudy_f = all_prop{dudy_idx(plot_inc)};
        dvdx_f = all_prop{dvdx_idx(plot_inc)};
        dvdy_f = all_prop{dvdy_idx(plot_inc)};
    catch
        dudx_f = all_prop{dudx_idx(end)};
        dudy_f = all_prop{dudy_idx(end)};
        dvdx_f = all_prop{dvdx_idx(end)};
        dvdy_f = all_prop{dvdy_idx(end)};
    end
    
    grid_sz = size(newEBSD);
    x_grid = reshape(newEBSD.x, grid_sz);
    y_grid = reshape(newEBSD.y, grid_sz);
    x_ext = [nanmean(x_grid(:, 1)), nanmean(x_grid(:, end))];
    y_ext = [nanmean(y_grid(1, :)), nanmean(y_grid(end, :))];

    % Collect all gradient data
    dudx_data = newEBSD.prop.(dudx_f);
    dudy_data = newEBSD.prop.(dudy_f);
    dvdx_data = newEBSD.prop.(dvdx_f);
    dvdy_data = newEBSD.prop.(dvdy_f);
    
    all_vals = [dudx_data(:); dudy_data(:); dvdx_data(:); dvdy_data(:)];
    
    % Ignore NaNs
    all_vals = all_vals(~isnan(all_vals));
    
    % Global limits
    cmin = 0.2*min(all_vals);
    cmax = 0.2*max(all_vals);
    
    fig_grad = figure('Name', 'Displacement Gradient Verification', 'Position', [150, 150, 1200, 900]);
    subplot(2, 2, 1); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudx_f), grid_sz)); title(['du/dx: ' strrep(dudx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb1 = colorbar; cb1.Title.String = 'du/dx';
    caxis([cmin cmax])
    subplot(2, 2, 2); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dudy_f), grid_sz)); title(['du/dy: ' strrep(dudy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb2 = colorbar; cb2.Title.String = 'du/dy';
    caxis([cmin cmax])
    subplot(2, 2, 3); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdx_f), grid_sz)); title(['dv/dx: ' strrep(dvdx_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb3 = colorbar; cb3.Title.String = 'dv/dx';
    caxis([cmin cmax])
    subplot(2, 2, 4); imagesc(x_ext, y_ext, reshape(newEBSD.prop.(dvdy_f), grid_sz)); title(['dv/dy: ' strrep(dvdy_f, '_', '\_')]); axis image off; colormap(gca, 'parula'); cb4 = colorbar; cb4.Title.String = 'dv/dy';
    caxis([cmin cmax])
    for sp_i = 1:4
        subplot(2, 2, sp_i); hold on;
        quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
        quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
        hold off;
    end
    
    gradOutName = fullfile(input.outPath, [input.fname, '_', options.dateString, '_DispGradients_4PanelCheck_.png']);
    exportScaledFigure(fig_grad, gradOutName, 'nodatestring');
    fprintf('Saved 4-Panel Displacement Gradient Check figure to: %s\n', gradOutName);
end

finalOutName = fullfile(input.outPath, [options.in_situIncrementName, '_undeformedEBSD_withAllDIC.mat']);
save(finalOutName, 'newEBSD', 'originalEBSD', 'grains', 'options');
fprintf('Saved Undeformed Reference dataset & FZ grains to: %s\n', finalOutName);

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
