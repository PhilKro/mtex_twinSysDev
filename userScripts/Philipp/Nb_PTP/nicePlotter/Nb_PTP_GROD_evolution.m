%% Nb_PTP_GROD_evolution.m
% Script to calculate and plot Grain Reference Orientation Deviation (GROD)
% across multiple increments. It uses a fixed reference orientation from the
% first increment for all subsequent calculations to ensure consistency and
% track evolution from the initial state.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Parameters
if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir = fullfile(baseDir, 'GROD_Analysis');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% --- Analysis Parameters ---
grainThresh_ref = 20 * degree;  % Stricter threshold for establishing stable reference grains
grainThresh_inc = 20 * degree; % Liberal threshold for defining GROD calculation area
minSize = 10;                  % Minimum grain size in pixels
tkd_tilt = -10 * degree;
thresh_ci = 0.13;
CS = crystalSymmetry.load('Nb-Niobium.cif');

%% 3. Find and Sort Data Files
angFiles = dir(fullfile(baseDir, '*263.ang'));
% Sort files by name to ensure they are processed in temporal order
[~, sortOrder] = sort({angFiles.name});
angFiles = angFiles(sortOrder);

%% 4. Establish Reference Orientations from the First Increment
if isempty(angFiles)
    error('No .ang files found in the specified directory: %s', baseDir);
end

fprintf('--- Establishing reference orientations from: %s ---\n', angFiles(1).name);
[~, baseName_ref, ~] = fileparts(angFiles(1).name);
% Load the first file to define the reference state
ebsd_ref = EBSD.load(fullfile(baseDir, angFiles(1).name), CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');

% Standard corrections
ebsd_ref = rotate(ebsd_ref, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
ebsd_ref(ebsd_ref.ci < thresh_ci) = 'notIndexed';
ebsd_ref = rotate(ebsd_ref, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');

% Apply specimen mask if it exists
mask_file_ref = fullfile(baseDir, sprintf('%s_specimen_mask.mat', baseName_ref));
if exist(mask_file_ref, 'file')
    fprintf('Applying specimen mask to reference data...\n');
    load(mask_file_ref, 'specimenlocation');
    ebsd_ref(~specimenlocation) = 'notIndexed';
else
    warning('No specimen mask found for %s. Proceeding with unmasked data.', baseName_ref);
end

% Calculate reference grains with a stricter threshold for stability
[grains_ref, ~] = calcGrains(ebsd_ref, 'minPixel', minSize, 'angle', grainThresh_ref);
ref_orientations = grains_ref.meanOrientation;
fprintf('Found %d reference grains to track.\n\n', length(grains_ref));

%% 5. Loop Through All Increments to Calculate and Plot GROD
for i = 1:length(angFiles)
    currentFile = angFiles(i).name;
    ang_path = fullfile(baseDir, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    fprintf('--- Processing %d/%d: %s ---\n', i, length(angFiles), currentFile);
    
    % Load current EBSD data and apply corrections
    ebsd_inc = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd_inc = rotate(ebsd_inc, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd_inc(ebsd_inc.ci < thresh_ci) = 'notIndexed';
    ebsd_inc = rotate(ebsd_inc, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Apply specimen mask if it exists
    mask_file_inc = fullfile(baseDir, sprintf('%s_specimen_mask.mat', baseName));
    if exist(mask_file_inc, 'file')
        fprintf('  -> Applying specimen mask...\n');
        load(mask_file_inc, 'specimenlocation');
        ebsd_inc(~specimenlocation) = 'notIndexed';
    else
        warning('  -> No specimen mask found for %s. Proceeding with unmasked data.', baseName);
    end

    % Calculate grains for the current increment using the liberal threshold
    [grains_inc, ebsd_inc.grainId] = calcGrains(ebsd_inc, 'minPixel', minSize, 'angle', grainThresh_inc);
    
    % --- Map current grains to the fixed reference orientations ---
    % Create a list of reference orientations indexed by the current grainId
    oriRefList = orientation.nan(max(grains_inc.id), 1, CS);
    for k = 1:length(grains_inc)
        g_inc = grains_inc(k);
        if g_inc.phase == -1
            continue
        end
        % Find the best matching reference orientation by minimum misorientation angle
        [~, best_match_idx] = min(angle(g_inc.meanOrientation, ref_orientations));
        oriRefList(g_inc.id) = ref_orientations(best_match_idx);
    end
    
    % --- Calculate GROD using the fixed reference orientations ---
    grod = ebsd_inc.calcGROD(grains_inc, oriRefList);
    
    % --- Plotting ---
    
    % Plot 1: GROD Angle Map
    fig1 = figure('Name', [baseName ' GROD Angle'], 'Visible', 'off');
    plot(ebsd_inc, grod.angle./degree, 'micronbar', 'off');
    hold on; plot(grains_inc.boundary, 'lineWidth', 1.5); hold off;
    % mtexColorMap LaboTeX;
    mtexColorMap viridis
    mtexColorbar('title', 'GROD Angle (deg)');
    title(sprintf('GROD Angle - %s', baseName), 'Interpreter', 'none');
    exportScaledFigure(fig1, fullfile(outDir, [baseName, '_GROD_Angle_viridis.png']), 'noDateString');
    
    % Plot 2: GROD Axis Map (Crystal Coordinates)
    fig2 = figure('Name', [baseName ' GROD Axis Map'], 'Visible', 'off');
    colorKey = HSVDirectionKey(CS, 'antipodal');
    alpha_val = min(grod.angle./degree/10, 1); % Scale transparency by angle
    plot(ebsd_inc, colorKey.direction2color(grod.axis), 'micronbar', 'off', 'FaceAlpha', alpha_val);
    hold on; plot(grains_inc.boundary, 'lineWidth', 1.5); hold off;
    title(sprintf('GROD Axis (Crystal Coords) - %s', baseName), 'Interpreter', 'none');
    exportScaledFigure(fig2, fullfile(outDir, [baseName, '_GROD_AxisMap.png']), 'noDateString');
    
    % Plot 3: GROD Axis Pole Figure
    fig3 = figure('Name', [baseName ' GROD Axis PF'], 'Visible', 'off');
    
    % Plot the weighted density contour map as the background
    plot(grod.axis, 'weights', grod.angle, 'contourf', 'fundamentalRegion', 'antipodal', 'halfwidth', 10*degree, 'DisplayName', 'Weighted Density');
    mtexColorbar('title', 'Weighted M.U.D.');
    hold on;
    
    % Plot 4: KAM map
    kam = ebsd_inc.KAM / degree;

    % lets plot it
    figure;
    plot(ebsd_inc,kam,'micronbar','off')
    hold on
    plot(grains_inc.boundary,'lineWidth',1.5)
    setColorRange([0,5])
    mtexColorbar
    mtexColorMap LaboTeX
    hold off
    exportScaledFigure(fig2, fullfile(outDir, [baseName, 'KAM_Map.png']), 'noDateString');

    % --- Add scatter points for axes with the highest angles ---
    angle_thresh_80 = prctile(grod.angle, 80);

    % Select and plot axes for top 20%
    axes_top_20 = grod.axis(grod.angle >= angle_thresh_80);
    if ~isempty(axes_top_20)
        plot(axes_top_20,'contourf','fundamentalRegion','antipodal','figSize','small', 'DisplayName', 'Top 20% Angle Axes');
    end
   
    title(sprintf('GROD Axis of top 20% of misorientation angles - %s', baseName), 'Interpreter', 'none');
    legend('Location', 'best');
    exportScaledFigure(fig3, fullfile(outDir, [baseName, '_GROD_AxisPF_top20.png']), 'noDateString');
    
    close(fig1, fig2, fig3);
    fprintf('  -> Saved GROD plots for %s\n', baseName);
end

fprintf('\nFinished processing all files.\n');