%% Twin Analysis of Rhenium EBSD Data
% This script performs a twin analysis on a Rhenium EBSD dataset. It
% identifies a parent grain, calculates primary and secondary twin
% orientations based on the {11-21} twinning system, and visualizes the
% parent, primary, and secondary twins in the grain map.

% Clear workspace and close figures
clear all
close all

% Set MTEX preferences for axis directions
plotx2east
plotzIntoPlane

%% 1. Load and Pre-process Data
disp('Loading and pre-processing EBSD data...');

% Define file name and crystal symmetry for Rhenium
% Please ensure the .ang file is in the MATLAB path.
fname = 'userScripts/TwinClassPresentation/pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% Load EBSD data
ebsd = EBSD.load(fname, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

% Clean data based on a confidence index (CI) threshold
ebsd(ebsd.ci < 0.15) = 'notIndexed';

tkd_tilt = -20 * degree;
% Perform corrections for tkd tilt and misaligned axes
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

% Gridify the data for faster processing and visualization
ebsd = ebsd.gridify;

% some initial test plot
ipf=ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection=yvector;
colors = ipf.orientation2color(ebsd.orientations);

figure;
plot(ebsd,colors)
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

% Plot reconstructed grains
figure;
plot(grains, grains.meanOrientation);
title('Reconstructed Grains in Rhenium Sample');

%% 3. Identify Parent Grain and Define Twin Systems
disp('Identifying parent grain and defining twin systems...');

% Identify the largest grain as the parent grain for the analysis
[~, parentId] = max(grains.numPixel);
parentGrain = grains(parentId);
parentGrains = grains(angle(parentGrain.meanOrientation, grains.meanOrientation) < 8*degree);

% Create a twinSystem object for the {11-21} twinning in Rhenium
tS_1121 = twinSystem.hexagonal_1121_0001(CS);

% Symmetrise to get all unique twin variants
tS_sym = tS_1121.symmetrise;

% Generate the first layer of twins from the parent grain
grain_sym_twin = parentGrain * tS_sym;

% Generate the second layer of twins from the first layer
analysis_twin = grain_sym_twin * tS_sym;

%% 4. Twin Identification and Analysis
disp('Identifying primary and secondary twins...');

% Get orientations of the primary (first layer) twins
% implement sth like this
% primary_twins = unique(analysis_twin.parent);
% primary_twin_oris = primary_twins.orientation;
% this is the shortcut
primary_twin_oris = grain_sym_twin.orientation;

% Get orientations of the secondary (second layer) twins
secondary_twin_oris = analysis_twin.orientation;

% Define tolerance for orientation matching
tolerance = 8 * degree;

% Initialize logical arrays to classify grains
isPrimary = false(length(grains), 1);
isSecondary = false(length(grains), 1);
primaryVariantId = zeros(length(grains), 1);
secondaryVariantId = zeros(length(grains), 2);

% --- Identify Primary Twins ---
% Iterate through all grains to find primary twins
for i = 1:length(grains)
    % Skip the parent grain
    if grains(i).id == parentGrain.id, continue; end
    
    % Find the best matching primary twin variant for the current grain
    [minAngle, bestFit] = min(angle(grains(i).meanOrientation, primary_twin_oris));
    
    % If the match is within tolerance, classify as a primary twin
    if minAngle < tolerance
        isPrimary(i) = true;
        primaryVariantId(i) = grain_sym_twin(bestFit).variantId;
    end
end

% --- Identify Secondary Twins ---
% Iterate through all grains to find secondary twins among the remaining ones
for i = 1:length(grains)
    % Skip parent and primary grains
    if any(grains(i).id == parentGrains.id) || isPrimary(i), continue; end
    
    % Check if the grain orientation matches any secondary twin orientation
    [minAngle, bestFit] = min(angle(grains(i).meanOrientation, secondary_twin_oris));
    
    % If the match is within tolerance, classify as a secondary twin
    if minAngle < tolerance
        isSecondary(i) = true;
        secondaryVariantId(i,:) = [analysis_twin(bestFit).parent.variantId, analysis_twin(bestFit).variantId];
    end
end

%% 5. Visualization
disp('Visualizing the twin analysis results...');

% Create a new figure for the final plot
figure;
% Plot the base grain map with some transparency
plot(grains, grains.meanOrientation, 'facealpha', 0.3)
hold on


% Define colormap for primary twin variants
colors = lines(length(primary_twin_oris));

% Plot the parent grain in a grey tone
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

% Plot the identified primary twins with variant-specific colors
primaryGrains = grains(isPrimary);
variantIdsOfPrimaryGrains = primaryVariantId(isPrimary);
uniqueVariants = unique(variantIdsOfPrimaryGrains);

for v_idx = 1:length(uniqueVariants)
    variant = uniqueVariants(v_idx);
    if variant == 0, continue; end % Skip if variant ID is 0
    
    % Select grains belonging to the current variant
    grainsOfVariant = primaryGrains(variantIdsOfPrimaryGrains == variant);
    
    if ~isempty(grainsOfVariant)
        hold on
        plot(grainsOfVariant, 'facecolor', colors(variant,:));
        text(grainsOfVariant, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
    end
end

% Plot the identified secondary twins in pink
secondaryGrains = grains(isSecondary);
if ~isempty(secondaryGrains)
    hold on
    plot(secondaryGrains, 'facecolor', 'm'); % 'm' for magenta/pink
    secIds = secondaryVariantId(isSecondary, :);
    labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    text(secondaryGrains, labels, 'HorizontalAlignment', 'center');
end

hold off
title('Primary and Secondary Twin Analysis of Rhenium');
legend off % Hide legend for clarity

disp('Analysis complete.');
