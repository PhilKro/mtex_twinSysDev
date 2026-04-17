%% Initialization & MTEX Preferences
clear; close all; clc;

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

if exist('master_pwd', 'var')
    baseDir = master_pwd;
else
    baseDir = pwd;
end

% Crystal symmetry for Rhenium
CS = crystalSymmetry('6/mmm', [2.761 2.761 4.458], 'x||a', 'mineral', 'Re');

% Common settings
thresh_ci = 0.12;
grainThresh = 5 * degree;
minPixelSize = 5;

%% 1. Load TKD Data (Reference EBSD)
fprintf('Loading TKD Data...\n');
inOut_dir = fullfile(baseDir, 'onedriveData','Philipp','Re_PTP','Nanomegas_angs');
tkd_file = fullfile(baseDir, 'onedriveData', 'Philipp', 'Re_PTP', 'DIC_Aligned_Data', 'data_Re_7deg_PTP3', 'EBSD_angs', 'Re_7deg_insitu Specimen 1 8 Map Data 11_BW158.ang');
tkd_tilt = -11 * degree; % Typical TKD tilt, adjust if needed

ebsd_tkd = EBSD.load(tkd_file, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');

% Perform standard TKD corrections for tilt and misaligned axes
ebsd_tkd = rotate(ebsd_tkd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
ebsd_tkd(ebsd_tkd.ci < thresh_ci) = 'notIndexed';
ebsd_tkd = rotate(ebsd_tkd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');


%% 2. Load Nanomegas Data (ASTAR / TEM)
fprintf('Loading Nanomegas Data...\n');
nm_dir = inOut_dir;

% Toggle between the zoomed in and zoomed out files:
% nm_fname = 'Nanomegas_300x250px_3nm_0c6deg_15pre_10kx (0-0 300-250)_full.ang'; % Zoomed out
nm_fname = 'Nanomegas_200x150px_1nm_0c6deg_15pre_10kx (0-0 200-150).ang'; % Zoomed in

nm_path = fullfile(nm_dir, nm_fname);

% Nanomegas specific settings (Backside Viewing Adjustment)
% TKD viewing direction is +z, Nanomegas viewing direction is -z.
% To align them, the sample must be "flipped". A point inversion of everything 
% corresponds to inverting both X and Y spatially (180 around Z).
% A physical flip is usually a 180 deg rotation around Y (left-right) or X (up-down).
% Toggle the appropriate one that aligns your specific maps.
point_inversion = true;  % 180 deg around Z (inverts X and Y)
flip_left_right = false; % 180 deg around Y (flips X and Z)
flip_up_down    = false; % 180 deg around X (flips Y and Z)

% Additional fine-tuning rotations
xrot = 0 * degree;
yrot = 0 * degree;
zrot = 175 * degree;

% Note: ASTAR/Nanomegas usually maps well with 'setting 4' in MTEX
ebsd_nm = EBSD.load(nm_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 4');
% ebsd_nm = loadEBSD_ang(nm_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 4');
ebsd_nm(ebsd_nm.ci < thresh_ci) = 'notIndexed';

% Perform backside correction
if point_inversion
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(zvector, 180*degree));
end
if flip_left_right
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(yvector, 180*degree));
end
if flip_up_down
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(xvector, 180*degree));
end

% Perform additional corrections for TEM tilt and misaligned axes
if xrot ~= 0
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(xvector, xrot));
end
if yrot ~= 0
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(yvector, yrot));
end
if zrot ~= 0
    ebsd_nm = rotate(ebsd_nm, rotation.byAxisAngle(zvector, zrot));
end



%% Project Orientations to Fundamental Region
fprintf('Projecting Orientations to Fundamental Region...\n');
reforientation = orientation.byEuler(0*degree,0*degree,0*degree, CS);
ebsd_tkd.orientations = ebsd_tkd.orientations.project2FundamentalRegion(reforientation);
ebsd_nm.orientations = ebsd_nm.orientations.project2FundamentalRegion(reforientation);

ebsd_tkd = ebsd_tkd.gridify;
ebsd_nm = ebsd_nm.gridify;
%% 3. Grain Reconstruction
% Extract base filenames for output
[~, tkd_file_base, ~] = fileparts(tkd_file);
[~, nm_file_base, ~] = fileparts(nm_path);
fprintf('Reconstructing Grains...\n');

% TKD Grains
[grains_tkd, ebsd_tkd.grainId] = calcGrains(ebsd_tkd, 'angle', grainThresh, 'minPixel', minPixelSize);
grains_tkd = smooth(grains_tkd, 5);

% Nanomegas Grains
[grains_nm, ebsd_nm.grainId] = calcGrains(ebsd_nm, 'angle', grainThresh, 'minPixel', minPixelSize, 'alpha', 1);
grains_nm = smooth(grains_nm, 5);

% Nanomegas Garbage indexing pseudosymmetry
nanoMegas_pseudoSymmetry = orientation.byAxisAngle(Miller(-1,1,0,0, grains_nm.CS, 'UVTW'), 65.33*degree);
[ebsd_nm,grains_nm,numChanged] = cleanUpPseudoSym(ebsd_nm,grains_nm,nanoMegas_pseudoSymmetry);
numChanged

%% 4. Plotting Section
fprintf('Plotting Maps and Crystal Shapes...\n');

ipfKey = ipfColorKey(CS);
ipfKey.inversePoleFigureDirection = zvector; % IPF Y map

% Create shared crystal shape geometry
cS1 = crystalShape.hex(CS);

% --- TKD Plots ---
fig_tkd = figure('Name', 'TKD EBSD Analysis', 'Position', [100 100 1200 500]);
newMtexFigure('layout', [1, 2]);

nextAxis(1, 1);
plot(ebsd_tkd, ipfKey.orientation2color(ebsd_tkd.orientations), 'micronbar', 'on');
hold on;
plot(grains_tkd.boundary, 'linewidth', 1.5, 'linecolor', 'k');
title('TKD Map - IPF Y');

nextAxis(1, 2);
plot(grains_tkd, ipfKey.orientation2color(grains_tkd.meanOrientation));
hold on;
plot(grains_tkd.boundary, 'linewidth', 1.5, 'linecolor', 'k');
plot(grains_tkd, 0.4 * cS1, 'FaceAlpha', 0.8, 'linewidth', 1.5);
title('TKD - Grains & Crystal Shapes');
exportScaledFigure(fig_tkd, fullfile(inOut_dir, ['TKD_', tkd_file_base, '_IPFY_GrainsCrystalShapes.jpg']), 'nodatestring');

% --- Nanomegas Plots ---
fig_nm = figure('Name', 'Nanomegas ASTAR Analysis', 'Position', [150 150 1200 500]);
newMtexFigure('layout', [1, 2]);

nextAxis(1, 1);
plot(ebsd_nm, ipfKey.orientation2color(ebsd_nm.orientations), 'micronbar', 'on');
hold on;
plot(grains_nm.boundary, 'linewidth', 1.5, 'linecolor', 'k');
title('Nanomegas Map - IPF Y');

nextAxis(1, 2);
plot(grains_nm, ipfKey.orientation2color(grains_nm.meanOrientation));
hold on;
plot(grains_nm.boundary, 'linewidth', 1.5, 'linecolor', 'k');
% Adjust scaling parameter based on nm field of view
plot(grains_nm, 0.4 * cS1, 'FaceAlpha', 0.8, 'linewidth', 1.5);
exportScaledFigure(fig_nm, fullfile(inOut_dir, ['Nanomegas_', nm_file_base, '_IPFY_GrainsCrystalShapes.jpg']), 'nodatestring');
title('Nanomegas - Grains & Crystal Shapes');

%% 5. Identify Parent Grains & Twin Systems
disp('Identifying parent grains and defining twin systems...');

% TKD Parent
[~, parentId_tkd] = max(grains_tkd.numPixel);
parentGrain_tkd = grains_tkd(parentId_tkd);
parentGrains_tkd = grains_tkd(angle(parentGrain_tkd.meanOrientation, grains_tkd.meanOrientation) < 8*degree);

% Nanomegas Parent
[~, parentId_nm] = max(grains_nm.numPixel);
parentGrain_nm = grains_nm(parentId_nm);
parentGrains_nm = grains_nm(angle(parentGrain_nm.meanOrientation, grains_nm.meanOrientation) < 8*degree);

% Twin System definition ({11-21} for Re)
tS_1121 = twinSystem.hexagonal_1121_0001(CS);
tS_sym = tS_1121.symmetrise;

%% 6. Variant Determination (TKD)
disp('Calculating variants for TKD data...');
primaryTwin_tkd = parentGrain_tkd.meanOrientation * tS_sym;
secondTwin_tkd = primaryTwin_tkd * tS_sym;

primary_twin_oris_tkd = primaryTwin_tkd.orientation;
secondary_twin_oris_tkd = secondTwin_tkd.orientation;

tolerance = 8 * degree;
isPrimary_tkd = false(length(grains_tkd), 1);
isSecondary_tkd = false(length(grains_tkd), 1);
primaryVariantId_tkd = zeros(length(grains_tkd), 1);
secondaryVariantId_tkd = zeros(length(grains_tkd), 2);

for i = 1:length(grains_tkd)
    if any(grains_tkd(i).id == parentGrains_tkd.id), continue; end
    [minAngle, bestFit] = min(angle(grains_tkd(i).meanOrientation, primary_twin_oris_tkd));
    if minAngle < tolerance
        isPrimary_tkd(i) = true;
        primaryVariantId_tkd(i) = primaryTwin_tkd(bestFit).variantId;
    end
end

for i = 1:length(grains_tkd)
    if any(grains_tkd(i).id == parentGrains_tkd.id) || isPrimary_tkd(i), continue; end
    [minAngle, bestFit] = min(angle(grains_tkd(i).meanOrientation, secondary_twin_oris_tkd));
    if minAngle < tolerance
        isSecondary_tkd(i) = true;
        secondaryVariantId_tkd(i,:) = [secondTwin_tkd(bestFit).parent.variantId, secondTwin_tkd(bestFit).variantId];
    end
end

%% 7. Variant Determination (Nanomegas)
disp('Calculating variants for Nanomegas data...');
primaryTwin_nm = parentGrain_nm.meanOrientation * tS_sym;
secondTwin_nm = primaryTwin_nm * tS_sym;

primary_twin_oris_nm = primaryTwin_nm.orientation;
secondary_twin_oris_nm = secondTwin_nm.orientation;

isPrimary_nm = false(length(grains_nm), 1);
isSecondary_nm = false(length(grains_nm), 1);
primaryVariantId_nm = zeros(length(grains_nm), 1);
secondaryVariantId_nm = zeros(length(grains_nm), 2);

for i = 1:length(grains_nm)
    if any(grains_nm(i).id == parentGrains_nm.id), continue; end
    [minAngle, bestFit] = min(angle(grains_nm(i).meanOrientation, primary_twin_oris_nm));
    if minAngle < tolerance
        isPrimary_nm(i) = true;
        primaryVariantId_nm(i) = primaryTwin_nm(bestFit).variantId;
    end
end

for i = 1:length(grains_nm)
    if any(grains_nm(i).id == parentGrains_nm.id) || isPrimary_nm(i), continue; end
    [minAngle, bestFit] = min(angle(grains_nm(i).meanOrientation, secondary_twin_oris_nm));
    if minAngle < tolerance
        isSecondary_nm(i) = true;
        secondaryVariantId_nm(i,:) = [secondTwin_nm(bestFit).parent.variantId, secondTwin_nm(bestFit).variantId];
    end
end

%% 8. Plotting Twin Maps
fprintf('Plotting Twin Variant Maps...\n');
fig_twins = figure('Name', 'Twin Variant Comparison', 'Position', [200 200 1200 500]);
newMtexFigure('layout', [1, 2]);

% Colormap for primary variants
varColors = lines(length(tS_sym));

% --- TKD Twin Map ---
nextAxis(1, 1);
plot(grains_tkd, grains_tkd.meanOrientation, 'facealpha', 0.3); hold on;
plot(parentGrains_tkd, 'facecolor', [0.8 0.8 0.8]);
text(parentGrains_tkd, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

primaryGrains_tkd = grains_tkd(isPrimary_tkd);
varIds_tkd = primaryVariantId_tkd(isPrimary_tkd);
uniqueVars_tkd = unique(varIds_tkd);
for v_idx = 1:length(uniqueVars_tkd)
    variant = uniqueVars_tkd(v_idx);
    if variant == 0, continue; end
    g_var = primaryGrains_tkd(varIds_tkd == variant);
    hold on;
    plot(g_var, 'facecolor', varColors(variant,:));
    text(g_var, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
end

secondaryGrains_tkd = grains_tkd(isSecondary_tkd);
if ~isempty(secondaryGrains_tkd)
    hold on
    plot(secondaryGrains_tkd, 'facecolor', 'm');
    secIds = secondaryVariantId_tkd(isSecondary_tkd, :);
    labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    text(secondaryGrains_tkd, labels, 'HorizontalAlignment', 'center');
end
hold on
plot(grains_tkd.boundary, 'linewidth', 1.5, 'linecolor', 'k');
title('TKD - Twin Variants');
exportScaledFigure(fig_twins, fullfile(inOut_dir, ['TwinVariants_TKD_', tkd_file_base, '.jpg']), 'nodatestring');
hold off;

% --- Nanomegas Twin Map ---
nextAxis(1, 2);
plot(grains_nm, grains_nm.meanOrientation, 'facealpha', 0.3); hold on;
plot(parentGrains_nm, 'facecolor', [0.8 0.8 0.8]);
text(parentGrains_nm, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

primaryGrains_nm = grains_nm(isPrimary_nm);
varIds_nm = primaryVariantId_nm(isPrimary_nm);
uniqueVars_nm = unique(varIds_nm);
for v_idx = 1:length(uniqueVars_nm)
    variant = uniqueVars_nm(v_idx);
    if variant == 0, continue; end
    g_var = primaryGrains_nm(varIds_nm == variant);
    hold on
    plot(g_var, 'facecolor', varColors(variant,:));
    text(g_var, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
end

secondaryGrains_nm = grains_nm(isSecondary_nm);
if ~isempty(secondaryGrains_nm)
    hold on
    plot(secondaryGrains_nm, 'facecolor', 'm');
    secIds = secondaryVariantId_nm(isSecondary_nm, :);
    labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    text(secondaryGrains_nm, labels, 'HorizontalAlignment', 'center');
end
hold on
plot(grains_nm.boundary, 'linewidth', 1.5, 'linecolor', 'k');
title('Nanomegas - Twin Variants');
exportScaledFigure(fig_twins, fullfile(inOut_dir, ['TwinVariants_Nanomegas_', nm_file_base, '.jpg']),'target_px', 4000, 'nodatestring');
hold off;
