%% =========================================================================
% standalone_01_GND_Analysis_BW158.m
% Standalone Geometrically Necessary Dislocation (GND) Analysis for HCP Ti
% =========================================================================
% Target file: insitu_new Specimen 1 05a 30kV Map Data 12_BW158.ang
%
% Key Features & HCP Mechanics:
%   1. Loads TKD EBSD scan, corrects orientation alignment, and applies masks.
%   2. Computes local crystal curvature using half-quadratic smoothing.
%   3. Defines all primary HCP dislocation systems (Basal <a>, Prism <a>,
%      Pyr <a>, and Pyr <c+a> 1st/2nd order).
%   4. CRITICAL: Explicitly incorporates correct Burgers vector lengths (|b|).
%      Unlike BCC/FCC where all primary |b| are identical, HCP <c+a> systems
%      have |b| ~ 5.53 Å compared to <a> systems with |b| ~ 2.95 Å.
%      Line energies are scaled proportionally to |b|^2 (screws) and
%      |b|^2 / (1 - nu) (edges) to prevent artificial over-assignment of <c+a>.
%   5. Generates publication-ready figures for Total GNDs, Edge/Screw ratio,
%      and individual dislocation families.
% =========================================================================

clear all; close all; clc;

%% 1. MTEX & PLOTTING CONVENTIONS
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. CONFIGURATION & DIRECTORIES
% Input Data Paths
baseDir   = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2';
ang_file  = fullfile(baseDir, 'insitu_new Specimen 1 05a 30kV Map Data 12_BW158.ang');
maskDir   = fullfile(baseDir, 'Masks');

% Output Path
outDir    = fullfile(baseDir, 'Processed_Output_GND_BW158');
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output directory: %s\n', outDir);
end

% TKD Alignment & EBSD Cleaning Parameters
cfg.tkd_tilt     = -7 * degree;  % TKD tilt correction (e.g., -7 to -11 degrees)
cfg.thresh_ci    = 0.10;         % Minimum Confidence Index (CI) threshold
cfg.minSize      = 50;           % Minimum grain size in pixels
cfg.grainThresh  = 5 * degree;   % Grain boundary threshold angle
cfg.smooth_alpha = 3.0;          % Half-quadratic filter smoothing parameter
cfg.nu           = 0.30;         % Poisson's ratio for Titanium
cfg.force_recalc = false;        % Set true to bypass saved .mat and recalculate GNDs
cfg.caxis_limits = [14.0 16.5];  % Log10 limits (m^-2) for plotting density maps
cfg.FZ_eulerAngles = [4.9181,0.0672,1.3421]; % FZ center such that the GNDs and SSLIP results in the twin line up (extracted from parentOri_twinned of undeformedReference_06_SSLIP_subtractedTwinShear.m)
cfg.dislocation_mode = 'ca_only'; % Options: 'all' (all <a> & <c+a>), 'a_only', 'ca_only', 'basal_prism', 'basal_pyr'

[~, baseName, ~] = fileparts(ang_file);
mode_tag     = ['_', lower(cfg.dislocation_mode)];
gnd_mat_path = fullfile(outDir, sprintf('GND_Results_%s%s.mat', baseName, mode_tag));

%% 3. LOAD OR CALCULATE GNDs
if exist(gnd_mat_path, 'file') && ~cfg.force_recalc
    fprintf('Loading existing GND calculation results from:\n  %s\n', gnd_mat_path);
    load(gnd_mat_path, 'ebsd', 'ebsdSGND', 'grains', 'dS', 'dSRot', 'rho', 'factor', 'kappa', 'cfg');
else
    fprintf('========================================================================\n');
    fprintf('STARTING NEW GND ANALYSIS FOR: %s\n', baseName);
    fprintf('========================================================================\n');
    
    if ~exist(ang_file, 'file')
        error('Input .ang file not found at:\n  %s\nPlease verify the file path.', ang_file);
    end
    
    % --- 3.1 Load EBSD and apply crystal symmetry ---
    % Titanium symmetry with room-temperature lattice parameters (c/a ~ 1.586)
    CS = crystalSymmetry('6/mmm', [2.95 2.95 4.68], 'x||a', 'mineral', 'Ti');
    fprintf('Loading TKD EBSD map...\n');
    ebsd = EBSD.load(ang_file, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % --- 3.2 Align axes and apply TKD tilt correction ---
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + cfg.tkd_tilt), 'keepXY');
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % --- 3.3 Apply confidence index (CI) filtering ---
    ebsd(ebsd.ci < cfg.thresh_ci) = 'notIndexed';

    fz_CenterOri = orientation.byEuler(cfg.FZ_eulerAngles, CS);
    ebsd.orientations = ebsd.orientations.project2FundamentalRegion(fz_CenterOri);
    
    % --- 3.4 Check and apply specimen mask (if available) ---
    % mask_file = fullfile(maskDir, sprintf('%s_specimen_mask.mat', baseName));
    % if exist(mask_file, 'file')
    %     fprintf('  -> Applying interactive specimen mask to remove vacuum artifacts...\n');
    %     maskData = load(mask_file, 'specimenlocation');
    %     if isfield(maskData, 'specimenlocation')
    %         if length(maskData.specimenlocation) == length(ebsd)
    %             ebsd.prop.specimenlocation = maskData.specimenlocation;
    %             ebsd(~ebsd.prop.specimenlocation) = 'notIndexed';
    %         else
    %             warning('Mask dimension mismatch. Proceeding unmasked.');
    %         end
    %     end
    % else
    %     fprintf('  -> No specimen mask found in %s. Proceeding with unmasked data.\n', maskDir);
    % end
    
    % --- 3.5 Reconstruct Grains ---
    fprintf('Reconstructing grains (Threshold: %.1f deg, Min Size: %d px)...\n', cfg.grainThresh/degree, cfg.minSize);
    [grains, ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', cfg.grainThresh, 'minPixel', cfg.minSize,'boundary','tight');
    
    % --- 3.6 Curvature Calculation via Half-Quadratic Smoothing ---
    fprintf('Applying half-quadratic orientation smoothing (alpha = %.1f)...\n', cfg.smooth_alpha);
    F = halfQuadraticFilter;
    F.alpha = cfg.smooth_alpha;
    ebsdSGND = gridify(smooth(ebsd('indexed'), F));
    
    fprintf('Computing curvature tensor field (kappa)...\n');
    kappa = ebsdSGND.curvature;
    
    % --- 3.7 Define HCP Dislocation Systems & Burgers Vectors ---
    fprintf('Constructing HCP dislocation systems & assigning elastic line energies...\n');
    sSBasal   = symmetrise(slipSystem.basal(CS), 'antipodal');
    sSPrism   = symmetrise(slipSystem.prismaticA(CS), 'antipodal');
    sSPyraA   = symmetrise(slipSystem.pyramidalA(CS), 'antipodal');
    sSPyraCA  = symmetrise(slipSystem.pyramidalCA(CS), 'antipodal');
    sSPyraCA2 = symmetrise(slipSystem.pyramidal2CA(CS), 'antipodal');
    
    % Combine systems according to configured user toggle
    switch lower(cfg.dislocation_mode)
        case 'a_only'
            sSelected = [sSBasal; sSPrism; sSPyraA];
            fprintf('  -> Dislocation Mode [a_only]: Evaluating exclusively <a> type systems (Basal, Prism, Pyr <a>)\n');
        case 'ca_only'
            sSelected = [sSPyraCA; sSPyraCA2];
            fprintf('  -> Dislocation Mode [ca_only]: Evaluating exclusively <c+a> type systems (Pyr 1st & 2nd <c+a>)\n');
        case 'basal_prism'
            sSelected = [sSBasal; sSPrism];
            fprintf('  -> Dislocation Mode [basal_prism]: Evaluating exclusively Basal <a> and Prismatic <a> systems\n');
        case 'basal_pyr'
            sSelected = [sSBasal; sSPyraA; sSPyraCA; sSPyraCA2];
            fprintf('  -> Dislocation Mode [basal_pyr]: Evaluating Basal <a>, Pyramidal <a>, and Pyramidal <c+a> systems\n');
        case 'all'
            sSelected = [sSBasal; sSPrism; sSPyraA; sSPyraCA; sSPyraCA2];
            fprintf('  -> Dislocation Mode [all]: Evaluating all <a> and <c+a> systems\n');
        otherwise
            warning('Unknown dislocation_mode "%s". Defaulting to "all".', cfg.dislocation_mode);
            sSelected = [sSBasal; sSPrism; sSPyraA; sSPyraCA; sSPyraCA2];
    end
    dS = dislocationSystem(sSelected);
    
    % --- 3.8 CRITICAL BURGERS VECTOR & LINE ENERGY WEIGHTING ---
    % In HCP Ti, <a> dislocations have |b_a| = 2.95 Å, while <c+a> dislocations have
    % |b_c+a| = sqrt(a^2 + c^2) = 5.53 Å.
    % Since elastic strain energy per unit length scales as E ~ G * |b|^2 (for screws)
    % and E ~ G * |b|^2 / (1 - nu) (for edges), we MUST assign dS.u explicitly using
    % the exact Burgers vector norms to penalize high-energy <c+a> dislocations!
    
    dS(dS.isScrew).u = (dS(dS.isScrew).b.norm).^2;
    dS(dS.isEdge).u  = (dS(dS.isEdge).b.norm).^2 / (1 - cfg.nu);
    
    fprintf('  -> Verified %d unique dislocation systems in mode [%s]:\n', length(dS), cfg.dislocation_mode);
    a_idx = eq(dS.b, Miller(1,1,-2,0,CS,'UVTW'));
    if any(a_idx)
        fprintf('     * <a> type systems   : |b| ~ %.2f Å | Energy weights: %.2f (screw) / %.2f (edge)\n', ...
                mean(dS(a_idx).b.norm), mean(dS(a_idx).b.norm)^2, mean(dS(a_idx).b.norm)^2 / (1 - cfg.nu));
    end
    ca_idx = eq(dS.b, Miller(1,1,-2,3,CS,'UVTW'));
    if any(ca_idx)
        fprintf('     * <c+a> type systems : |b| ~ %.2f Å | Energy weights: %.2f (screw) / %.2f (edge)\n', ...
                mean(dS(ca_idx).b.norm), mean(dS(ca_idx).b.norm)^2, mean(dS(ca_idx).b.norm)^2 / (1 - cfg.nu));
    end
    
    % Rotate dislocation systems into specimen coordinates at each map pixel
    dSRot = ebsdSGND.orientations * dS;
    
    % --- 3.9 Fit Dislocation Density Tensor ---
    fprintf('Solving optimization problem (fitDislocationSystems) across EBSD map...\n');
    tic;
    [rho, factor] = fitDislocationSystems(kappa, dSRot);
    fprintf('  -> Optimization complete in %.2f seconds!\n', toc);
    
    % Save calculation to disk for fast replotting
    save(gnd_mat_path, 'ebsd', 'ebsdSGND', 'grains', 'dS', 'dSRot', 'rho', 'factor', 'kappa', 'cfg');
    fprintf('Saved complete GND dataset to:\n  %s\n', gnd_mat_path);
end

%% 4. COMPREHENSIVE VISUALIZATION & FIGURE EXPORT
fprintf('\n========================================================================\n');
fprintf('GENERATING PUBLICATION-READY GND MAPS\n');
fprintf('========================================================================\n');

% Calculate total density and energy fields
ebsdSGND.prop.totalGNDdensity       = factor * sum(abs(rho), 2, 'omitnan');
ebsdSGND.prop.totalDislocationEnergy = factor * sum(abs(rho .* dSRot.u), 2, 'omitnan');

% --- FIGURE 1: TOTAL GND DENSITY MAP ---
fprintf('Plotting Figure 1: Total GND Density...\n');
figTotal = figure;
plot(ebsdSGND, log10(ebsdSGND.prop.totalGNDdensity), 'micronbar', 'off');
hold on;
plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w');
hold off;
mtexColorMap parula;
caxis(cfg.caxis_limits);
cb = mtexColorbar('title', 'log_{10}(\rho_{GND} [m^{-2}])');
title(sprintf('Total GND Density (%s) - %s', cfg.dislocation_mode, baseName), 'Interpreter', 'none', 'FontSize', 14);
drawNow(gcm, 'figSize', 'large');

exportScaledFigure(figTotal, fullfile(outDir, [baseName, mode_tag, '_Total_GND_Density.png']), 'nodatestring');
close(figTotal);

% --- FIGURE 2: TOTAL EDGE vs. TOTAL SCREW GND DENSITY ---
fprintf('Plotting Figure 2: Total Edge vs. Total Screw Density...\n');
figES = figure;
newMtexFigure('layout', [1, 2]);

GND_E = factor * sum(abs(rho(:, dS.isEdge)), 2, 'omitnan');
GND_S = factor * sum(abs(rho(:, dS.isScrew)), 2, 'omitnan');

nextAxis(1, 1);
plot(ebsdSGND, log10(GND_E), 'micronbar', 'off');
hold on; plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w'); hold off;
mtexColorMap parula; caxis(cfg.caxis_limits);
title('Total Edge Dislocation Density (\rho_{Edge})', 'FontSize', 12);

nextAxis(1, 2);
plot(ebsdSGND, log10(GND_S), 'micronbar', 'off');
hold on; plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w'); hold off;
mtexColorMap parula; caxis(cfg.caxis_limits);
title('Total Screw Dislocation Density (\rho_{Screw})', 'FontSize', 12);

mtexColorbar('title', 'log_{10}(\rho_{GND} [m^{-2}])');
drawNow(gcm, 'figSize', 'large');

exportScaledFigure(gcf, fullfile(outDir, [baseName, mode_tag, '_Edge_vs_Screw_GND.png']), 'nodatestring');
close(gcf);

% --- FIGURE 3: BREAKDOWN BY HCP DISLOCATION FAMILIES ---
fprintf('Plotting Figure 3: Dislocation Family Breakdown...\n');
figFamilies = figure;
newMtexFigure('layout', [3, 3]);

% Helper to compute density for a mask of dislocation systems
calcFamily = @(mask) factor * sum(abs(rho(:, mask)), 2, 'omitnan');

% Family Definitions & Names for HCP Ti
CS_m = ebsdSGND.CS;
f_masks = { ...
    dS.isScrew & dS.b.norm < 4.0, ...                                                          % 1. <a> Screw
    dS.isEdge  & dS.b.norm < 4.0 & angle(dS.l, Miller(1,-1,0,0,CS_m,'UVTW')) < 1e-3, ...      % 2. Basal <a> Edge
    dS.isEdge  & dS.b.norm < 4.0 & angle(dS.l, Miller(0,0,0,1,CS_m,'UVTW')) < 1e-3, ...       % 3. Prism <a> Edge
    dS.isEdge  & dS.b.norm < 4.0 & angle(dS.l, Miller(1,-1,0,2,CS_m,'UVTW')) < 1e-3, ...      % 4. Pyr <a> Edge
    dS.isScrew & dS.b.norm > 4.0, ...                                                          % 5. <c+a> Screw
    dS.isEdge  & dS.b.norm > 4.0 & angle(dS.l, Miller(1,-1,0,0,CS_m,'UVTW')) < 1e-3, ...      % 6. Pyr2 <c+a> Edge
    dS.isEdge  & dS.b.norm > 4.0 & angle(dS.l, Miller(1,-1,0,0,CS_m,'UVTW')) >= 1e-3 ...       % 7. Pyr1 <c+a> Edge
};

f_titles = { ...
    '<a> Screw (\rho_{<a>, screw})', ...
    'Basal <a> Edge (\rho_{basal, <a>})', ...
    'Prismatic <a> Edge (\rho_{prism, <a>})', ...
    'Pyramidal <a> Edge (\rho_{pyr, <a>})', ...
    '<c+a> Screw (\rho_{<c+a>, screw})', ...
    'Pyr2 <c+a> Edge (\rho_{pyr2, <c+a>})', ...
    'Pyr1 <c+a> Edge (\rho_{pyr1, <c+a>})' ...
};

for idx = 1:length(f_masks)
    row_i = ceil(idx / 3);
    col_i = mod(idx - 1, 3) + 1;
    nextAxis(row_i, col_i);
    
    fam_density = calcFamily(f_masks{idx});
    plot(ebsdSGND, log10(fam_density), 'micronbar', 'off');
    hold on; plot(grains.boundary, 'linewidth', 1.0, 'linecolor', 'w'); hold off;
    mtexColorMap parula; caxis(cfg.caxis_limits);
    title(f_titles{idx}, 'FontSize', 10);
end

% Plot Total GND in position (3,2) for reference
nextAxis(3, 2);
plot(ebsdSGND, log10(ebsdSGND.prop.totalGNDdensity), 'micronbar', 'off');
hold on; plot(grains.boundary, 'linewidth', 1.0, 'linecolor', 'w'); hold off;
mtexColorMap parula; caxis(cfg.caxis_limits);
title(sprintf('Total GND (%s)', cfg.dislocation_mode), 'FontSize', 10, 'FontWeight', 'bold');

mtexColorbar('title', 'log_{10}(\rho_{GND} [m^{-2}])');
drawNow(gcm, 'figSize', 'large');

exportScaledFigure(gcf, fullfile(outDir, [baseName, mode_tag, '_HCP_Dislocation_Families.png']), 'nodatestring');
close(gcf);

% --- FIGURE 4: CURVATURE TENSOR COMPONENTS (kappa) ---
fprintf('Plotting Figure 4: Curvature Tensor Components (kappa)...\n');
figCurv = figure;
newMtexFigure('nrows', 3, 'ncols', 3);

for i = 1:3
    for j = 1:3
        nextAxis(i, j);
        plot(ebsdSGND, kappa{i, j}, 'micronBar', 'off');
        hold on; plot(grains.boundary, 'linewidth', 1.0, 'linecolor', 'k'); hold off;
        title(sprintf('\\kappa_{%d%d}', i, j), 'FontSize', 11);
    end
end
setColorRange([-0.15, 0.15]);
mtexColorMap parula;
mtexColorbar('title', 'Curvature [rad/\mu m]');
drawNow(gcm, 'figSize', 'large');

exportScaledFigure(gcf, fullfile(outDir, [baseName, mode_tag, '_Curvature_Tensor_kappa.png']), 'nodatestring');
close(gcf);

fprintf('========================================================================\n');
fprintf('GND ANALYSIS COMPLETE! All results and figures saved in:\n  %s\n', outDir);
fprintf('========================================================================\n');
