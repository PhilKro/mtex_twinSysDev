%% Twin Analysis of Rhenium EBSD Data
% This script performs a twin analysis on a Calcite EBSD dataset.

% Clear workspace and close figures
clear all
close all

% Set MTEX preferences for axis directions
plotx2east
plotzIntoPlane

%% 1. Load and Pre-process Data
disp('Loading and pre-processing EBSD data...');

% CM-39
fname = "C:\Users\phkr\MTEX_workEnv\data\LuizMoralez\CM-39_clean.ang";

% Load EBSD data
ebsd = EBSD.load(fname, 'convertEuler2SpatialReferenceFrame','setting 2');

% Clean data based on a confidence index (CI) threshold
ebsd(ebsd.ci < 0.15) = 'notIndexed';

% Gridify the data for faster processing and visualization
ebsd = ebsd.gridify;

% Rectangular region selection
region = [320,138,300,300];
ebsd = ebsd(inpolygon(ebsd, region));

% some initial test plot
ipf=ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection=yvector;
colors = ipf.orientation2color(ebsd.orientations);


%figure;
plot(ebsd,colors)
% exportScaledFigure(gcf, [fname, 'ebsd.png'])
%% 2. Grain Reconstruction
disp('Reconstructing grains...');

% Define parameters for grain reconstruction
grainThresh = 5*degree; 
minPixel = 5; % minimum number of pixels per grain
alpha = 1; % alpha parameter for filling during reconstruction

% Calculate grains from the indexed EBSD data
[grains, ebsd.grainId] = calcGrains(ebsd, 'angle', grainThresh, 'minPixel', minPixel, 'alpha', alpha);

% Smooth grain boundaries
grains = smooth(grains, 5);

% Restrict to Calcite phase to avoid multi-phase errors later
grains = grains('Calcite');

% Plot reconstructed grains
figure;
grainColors = ipf.orientation2color(grains.meanOrientation);
plot(grains, grainColors);
title('Reconstructed Grains in Calcite Sample');
% exportScaledFigure(gcf, [fname, 'grains.png'])

% Plot grain boundary angles
figure;
plot(grains); hold on; plot(grains.boundary, angle(grains.boundary.misorientation)/degree); mtexColorbar
%% 3. Define twin system


CS = ebsd.CS;

possible_tS_1018 = twinSystem.calculateTheoreticalTwins(CS, Miller(1,0,-1,8, CS), 1:8)

twin_E_lowShear = possible_tS_1018(1); % q = 5

possible_tS_1014 = twinSystem.calculateTheoreticalTwins(CS, Miller(1,0,-1,4, CS), 1:8)

twin_rp_lowShear = possible_tS_1014(1);

possible_tS_1012 = twinSystem.calculateTheoreticalTwins(CS, Miller(1,0,-1,2, CS), 1:8)

twin_fn_lowShear = possible_tS_1012(2);

%% 4. Boundary Identification
gB = grains.boundary('Calcite', 'Calcite');

% Identify boundaries for each of the 3 low shear twin types
is_E_twin = isTwinning(gB, -twin_E_lowShear.parentTwinMisorientation, 8*degree); % negative because of definition of the tS misorientation as improper rotation and the non-centrosymmetric point group of Calcite
is_rp_twin = isTwinning(gB, -twin_rp_lowShear.parentTwinMisorientation, 8*degree);
is_fn_twin = isTwinning(gB, -twin_fn_lowShear.parentTwinMisorientation, 8*degree);

% Plot grain boundaries colored by twin type
figure;
plot(grains, grainColors, 'micronbar', 'off');
hold on;
plot(gB.subSet(is_E_twin), 'linecolor', 'r', 'linewidth', 2, 'displayName', 'E (10-18)');
plot(gB.subSet(is_rp_twin), 'linecolor', 'g', 'linewidth', 2, 'displayName', 'r'' (10-14)');
plot(gB.subSet(is_fn_twin), 'linecolor', 'b', 'linewidth', 2, 'displayName', 'f- (10-12)');
hold off;
title('Grains and Twin Boundaries');
legend('Location', 'best');