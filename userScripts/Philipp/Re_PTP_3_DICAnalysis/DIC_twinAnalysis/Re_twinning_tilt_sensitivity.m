%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tilt Sensitivity Analysis for Twinning Identification               %%%
%%%                                                                     %%%
%%% This script investigates the effect of an assumed sample tilt       %%%
%%% correction error on the identification of twin variants based on    %%%
%%% experimental strain data. It iterates through a series of tilt      %%%
%%% corrections, recalculates the theoretical twin shear for each tilt, %%%
%%% and compares it against the DIC strain field.                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% 1. Setup and Configuration
% MTEX preferences
plotx2east;
plotzIntoPlane;

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% --- Analysis Parameters ---
% Define the range of tilt corrections to test (in degrees)
% This will be applied as an additional rotation around the x-axis.
tilt_series = -10:5:10; 

% Primary variants to consider in the analysis
pVars = 1:6;

% Cutoff criteria for identifying a valid twin pixel
magCutoff = 0.22; 
shapeCutoff = 0.5; 

% --- Data and Paths ---
% inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3_Subset19_Step6_relative2First_firstOrderSF';
inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';

fnames = {fullfile(inputPath,'inc_0_newEBSD__withStrain_05-08-2025_10-27'),...
    fullfile(inputPath,'inc_4_newEBSD__withStrain_05-08-2025_23-58'),...
    fullfile(inputPath,'inc_6_newEBSD__withStrain_05-12-2025_15-07'),...
    fullfile(inputPath,'inc_7_newEBSD__withStrain_05-12-2025_15-16'),...
    fullfile(inputPath,'inc_8_newEBSD__withStrain_05-12-2025_15-23')};

% fnames = {fullfile(inputPath,'inc_4_newEBSD__withStrain_04-15-2026_16-48'),...
%     fullfile(inputPath,'inc_6_newEBSD__withStrain_04-15-2026_17-01'),...
%     fullfile(inputPath,'inc_7_newEBSD__withStrain_04-15-2026_17-06')};

N=length(fnames);

for i = 1:N
    loaded_data = load([fnames{i}, '_MTEX6.1.mat']);
    % loaded_data = load(fnames{i});
    ebsd{i}   = loaded_data.newEBSD;
    reforientation = orientation.byEuler(0*degree,0*degree,0*degree, ebsd{i}.CS);
    ebsd{i}.orientations = ebsd{i}.orientations.project2FundamentalRegion(reforientation);

    grains{i} = loaded_data.grains;
    area{i}   = loaded_data.area;
end
% ---------------------------------

CS = ebsd{1}.CS;
clear data ebsd0


% Define output directory for plots
outputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\MTEX_Analysis_output\Tilt_Sensitivity';
if ~exist(outputPath, 'dir'), mkdir(outputPath); end

% Select which increment to analyze
inc = 3;
comment = '_filt3';

%% 2. Initial Data Preparation (Done Once)
fprintf('--- Preparing Initial Data for Increment %d ---\n', inc);

CS = grains{inc}.CS;

% Ensure notIndexed pixels are handled correctly
ebsd{inc}(ebsd{inc}.orientations.isnan) = 'notIndexed';

% Identify the parent grain (assuming it's the largest one)
[~, parentId] = max(grains{inc}.numPixel);
parentGrain = grains{inc}(parentId);
originalParentOri = parentGrain.meanOrientation;
fprintf('Original parent grain ID: %d\n', parentId);

% Define the {11-21} twinning system for Rhenium
tS_1121 = twinSystem.hexagonal_1121_0001(CS);
tS_sym = tS_1121.symmetrise;

% Prepare experimental strain tensor E_exp from DIC data
% This is done only once as the experimental data does not change.
try
    nPoints = numel(ebsd{inc}.prop.exx2);
    E_exp = zeros(2, 2, nPoints);
    E_exp(1, 1, :) = ebsd{inc}.prop.exx2(:);
    E_exp(2, 1, :) = ebsd{inc}.prop.exy2(:);
    E_exp(1, 2, :) = ebsd{inc}.prop.exy2(:);
    E_exp(2, 2, :) = ebsd{inc}.prop.eyy2(:);
    Ee_flat = reshape(E_exp, 4, nPoints);
    fprintf('Using exx2/exy2/eyy2 for experimental strain.\n');
catch
    warning('exx2/exy2/eyy2 not found. Falling back to exx/exy/eyy.');
    nPoints = numel(ebsd{inc}.prop.exx);
    E_exp = zeros(2, 2, nPoints);
    E_exp(1, 1, :) = ebsd{inc}.prop.exx(:);
    E_exp(2, 1, :) = ebsd{inc}.prop.exy(:);
    E_exp(1, 2, :) = ebsd{inc}.prop.exy(:);
    E_exp(2, 2, :) = ebsd{inc}.prop.eyy(:);
    Ee_flat = reshape(E_exp, 4, nPoints);
end

% Define a colormap for the variants
colors = lines(numel(pVars));

%% 3. Tilt Sensitivity Loop
fprintf('\n--- Starting Tilt Sensitivity Analysis ---\n');

for tilt_correction_deg = tilt_series
    fprintf('Analyzing with tilt correction: %+d degrees\n', tilt_correction_deg);
    
    % --- a. Apply Tilt Correction to Parent Orientation ---
    tilt_rad = tilt_correction_deg * degree;
    rot = rotation.byAxisAngle(xvector, tilt_rad);
    
    % Rotate the original parent orientation
    parentOri_tilted = rot * originalParentOri;
    parentOri_tilted= parentOri_tilted.project2FundamentalRegion(reforientation);
    % --- b. Recalculate Theoretical Twin Tensors ---
    % Generate twin orientations from the tilted parent
    primaryTwin_tilted = parentOri_tilted * tS_sym;

    % Get theoretical displacement gradients and Green-Lagrange strain
    H_tilted = primaryTwin_tilted.displacementGradient;
    F_id_tilted = H_tilted(pVars) + tensor.eye(size(H_tilted(pVars)));
    F_id_tilted = deformationGradientTensor(F_id_tilted);
    E_id_tilted = F_id_tilted.GreenStrain;
    
    % Get the 2D projection of the theoretical strain
    E_id2d_tilted = E_id_tilted.M(1:2, 1:2, :);
    Ev_flat_tilted = reshape(E_id2d_tilted, 4, numel(pVars));

    % --- c. Perform Shape and Magnitude Analysis ---
    shapeMismatch = zeros(numel(pVars), nPoints);
    magnitudeX = zeros(numel(pVars), nPoints);

    for v = 1:numel(pVars)
        Ev = Ev_flat_tilted(:, v);
        
        % Calculate scalar projection factor 'x'
        x = (Ev' * Ee_flat) ./ (Ev' * Ev);
        magnitudeX(v, :) = x;
        
        % Calculate shape mismatch
        Proj = Ev * x;
        Mismatch = Ee_flat - Proj;
        shapeMismatch(v, :) = vecnorm(Mismatch, 2, 1);
    end

    % Find best fit variant based on minimum shape mismatch
    [minMismatch, bestVariantIdx] = min(shapeMismatch, [], 1);
    
    % Get the magnitude 'x' for the best-fit variant
    linearIndices = sub2ind(size(magnitudeX), bestVariantIdx, 1:nPoints);
    bestMagnitudeX = magnitudeX(linearIndices);

    % Identify valid pixels based on cutoff criteria
    isValid = (bestMagnitudeX > magCutoff) & (minMismatch < shapeCutoff);
    validIndices = find(isValid);
    nanMask = isnan(minMismatch);
    validIndices = validIndices(~nanMask(validIndices));

    % --- d. Plot and Save Results for Current Tilt ---
    
    % Plot 1: Best Fit Variant Map
    plotColor = ones(nPoints, 3); % White background
    plotColor(validIndices, :) = colors(bestVariantIdx(validIndices), :);
    
    fig1 = figure('Visible', 'off');
    plot(ebsd{inc}, plotColor);
    title(sprintf('Best Fit Variant (Tilt: %+d°)', tilt_correction_deg));
    colormap(colors);
    caxis([1 numel(pVars)]);
    mtexColorbar;
    
    saveFigure(fullfile(outputPath, sprintf('inc%d_BestFitMap_tilt_%+ddeg%s.png', inc, tilt_correction_deg, comment)));
    close(fig1);

    % Plot 2: Scalar Projection Factor 'x' Map
    fig2 = figure('Visible', 'off');
    plot(ebsd{inc}, bestMagnitudeX);
    title(sprintf('Scalar Projection Factor ''x'' (Tilt: %+d°)', tilt_correction_deg));
    colormap turbo;
    setColorRange([0.2, 1.35]);
    mtexColorbar;
    
    saveFigure(fullfile(outputPath, sprintf('inc%d_ProjectionFactor_tilt_%+ddeg%s.png', inc, tilt_correction_deg, comment)));
    close(fig2);
    
    % Plot 3: Shape Mismatch Map
    fig3 = figure('Visible', 'off');
    plot(ebsd{inc}, minMismatch);
    title(sprintf('Shape Mismatch of Best Fit (Tilt: %+d°)', tilt_correction_deg));
    colormap(viridis);
    setColorRange([0, shapeCutoff]);
    mtexColorbar;
    
    saveFigure(fullfile(outputPath, sprintf('inc%d_ShapeMismatch_tilt_%+ddeg%s.png', inc, tilt_correction_deg, comment)));
    close(fig3);
end

fprintf('\n--- Tilt sensitivity analysis complete. Plots saved to: %s ---\n', outputPath);