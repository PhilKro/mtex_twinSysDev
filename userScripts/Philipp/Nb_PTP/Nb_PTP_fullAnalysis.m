%% Nb_PTP_fullAnalysis.m
% This script compares the outputs of the SSLIP determination of dislocation 
% slip (based on undeformed reference TKD + DIC) with GND calculations from 
% deformed TKD scans at identical increments.
%
% It loads previously processed data if available to save computational time.

clear
close all

%% Initialization & Mtex Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Directories and Input Setup
baseDir = fullfile('onedriveData', 'Philipp', 'Nb_PTP');
outDir  = fullfile(baseDir, 'Processed_Output');
if ~exist(outDir, 'dir')
    warning('directory %s does not exist', outDir)
end

% Provide the name of the MAT file containing the fully merged EBSD & DIC 
% data from the Undeformed script
undeformed_dic_mat = fullfile(outDir, 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');

% Define increments (mapping TKD .ang files to concurrent DIC .hdf5 prefixes)
increments = struct();

increments(1).name = 'Increment 0';
increments(1).ang_file = 'insitu Specimen 1 0 Map Data 1_BW263.ang';
increments(1).dic_prefix = 'B0001';

increments(2).name = 'Increment 1';
increments(2).ang_file = 'insitu Specimen 1 1 Map Data 2_BW263.ang';
increments(2).dic_prefix = 'B0003';

increments(3).name = 'Increment 2';
increments(3).ang_file = 'insitu Specimen 1 2 Map Data 3_BW263.ang';
increments(3).dic_prefix = 'B0005';

increments(4).name = 'Increment 3';
increments(4).ang_file = 'insitu Specimen 1 3 Map Data 4_BW263.ang';
increments(4).dic_prefix = 'B0007';

increments(5).name = 'Increment 4';
increments(5).ang_file = 'insitu Specimen 1 4 Map Data 5_BW263.ang';
increments(5).dic_prefix = 'B0009';

increments(6).name = 'Increment 5';
increments(6).ang_file = 'insitu Specimen 1 5 Map Data 6_BW263.ang';
increments(6).dic_prefix = 'B0011';

increments(7).name = 'Increment 6';
increments(7).ang_file = 'insitu Specimen 1 7 Map Data 7_BW263.ang';
increments(7).dic_prefix = 'B0014';

% Choose which increment to analyze right now (1 through 7)
inc_idx = 4;
currentInc = increments(inc_idx);
fprintf('\n========================================\n');
fprintf('Analyzing %s (DIC Prefix: %s)\n', currentInc.name, currentInc.dic_prefix);
fprintf('========================================\n');

%% Universal Parameters
tkd_tilt = -10 * degree;
thresh_ci = 0.13;
grainThresh = 5 * degree;
minSize = 1000;
alpha = 1;

%% -------------------------------------------------------------------------
% PART 2: SSLIP Calculation on Undeformed Data + DIC
% -------------------------------------------------------------------------
sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', currentInc.dic_prefix));

if exist(sslip_mat_path, 'file')
    fprintf('Loading existing SSLIP data from %s...\n', sslip_mat_path);
    load(sslip_mat_path, 'ebsdID', 'sSLocalRef', 'optOut');
else
    fprintf('Calculating SSLIP data for %s...\n', currentInc.name);
    
    if ~exist(undeformed_dic_mat, 'file')
        error('Undeformed EBSD with DIC data not found at %s.', undeformed_dic_mat);
    end
    
    % Load aligned base data
    undefData = load(undeformed_dic_mat, 'newEBSD');
    ebsdUndef = undefData.newEBSD;
    
    % Clean and reconstruct grains
    ebsdUndef(ebsdUndef.ci < thresh_ci) = 'notIndexed';
    [grainsUndef, ebsdUndef.grainId] = calcGrains(ebsdUndef, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
    
    % Setup blank container for the fully combined SSLIP field
    ebsdID = ebsdUndef.gridify; 
    sS_BCC = symmetrise(slipSystem.bcc(ebsdUndef.CS), 'antipodal');
    ebsdID.prop.slipIDcor = zeros(length(sS_BCC), length(ebsdID));
    
    % We assume the 2 largest grains make up the bicrystal. We perform 
    % SSLIP sequentially on each grain so we can rotate the slip systems.
    [~, sortedIdx] = sort(grainsUndef.numPixel, 'descend');
    bicrystalGrains = grainsUndef(sortedIdx(1:min(2, length(grainsUndef))));
    
    sSLocalRef = sS_BCC; % Backup structure for reference
    
    for g_i = 1:length(bicrystalGrains)
        gId = bicrystalGrains(g_i).id;
        fprintf('  -> Running SSLIP for Grain ID %d...\n', gId);
        
        ebsd_g = ebsdUndef(grainsUndef(gId)).gridify;
        
        % Extract Displacements (Fixing previous bug: U is x, V is y)
        U = ebsd_g.prop.([currentInc.dic_prefix, '_Displacements_x_displacement']) * ebsdUndef.dPos;
        V = ebsd_g.prop.([currentInc.dic_prefix, '_Displacements_y_displacement']) * ebsdUndef.dPos;
        
        sSLocal_g = bicrystalGrains(g_i).meanOrientation * sS_BCC;
        
        % SSLIP Options
        IDoptions = struct();
        IDoptions.IDMethod = 1;
        IDoptions.filterSize = 4;
        IDoptions.coarsegrain = 4;
        IDoptions.threshResidual = 0.02;
        IDoptions.minEeff = 0.01;
        IDoptions.posConstr = 1;
        IDoptions.plotSSLIP = 1;
        IDoptions.saveFig = 0;
        IDoptions.plotDefGrad = 0;
        IDoptions.plotResidual = 0;
        IDoptions.NoSs = 1:length(sS_BCC);
        IDoptions.stress = stressTensor.uniaxial(yvector);      
        
        % Perform SSLIP locally
        [ebsdID_g, optOut] = SSLIP(ebsd_g, U, V, sSLocal_g, IDoptions);
        
        % Stitch local map back into the full ebsd grid using XY coordinates
        warning('The intersect between the ebsdId and ebsdID_g coordinates does not exist. The positions of x and y might be recalculated somewhere in SSLIP?')
        [~, idx_all, idx_g] = intersect([ebsdID.x(:), ebsdID.y(:)], [ebsdID_g.x(:), ebsdID_g.y(:)], 'rows');
        ebsdID.prop.slipIDcor(:, idx_all) = ebsdID_g.prop.slipIDcor(:, idx_g);
    end
    
    save(sslip_mat_path, 'ebsdID', 'sSLocalRef', 'optOut');
end

%% -------------------------------------------------------------------------
% PART 3: GND Calculation on Pristine Deformed TKD Data
% -------------------------------------------------------------------------
gnd_mat_path = fullfile(outDir, sprintf('GND_Data_%s.mat', currentInc.dic_prefix));

if exist(gnd_mat_path, 'file')
    fprintf('Loading existing GND data from %s...\n', gnd_mat_path);
    load(gnd_mat_path, 'ebsdDeformed', 'rho', 'factor', 'dSRot', 'dS', 'grainsDef');
else
    fprintf('Calculating GND data for %s...\n', currentInc.name);
    
    ang_path = fullfile(baseDir, currentInc.ang_file);
    if ~exist(ang_path, 'file')
        error('Deformed EBSD file not found at %s', ang_path);
    end
    
    CS = crystalSymmetry.load('Nb-Niobium.cif');
    ebsdDeformed = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsdDeformed = rotate(ebsdDeformed, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsdDeformed(ebsdDeformed.ci < thresh_ci) = 'notIndexed';
    ebsdDeformed = rotate(ebsdDeformed, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Project to fundamental reference orientation
    targetFile = fullfile(baseDir, 'refOri.mat');
    if exist(targetFile, 'file')
        refData = load(targetFile, 'refOrientation');
        ebsdDeformed.orientations = ebsdDeformed.orientations.project2FundamentalRegion(refData.refOrientation);
    else
        error('Reference orientation file not found at %s. Please define and save a reference orientation first.', targetFile);
    end
    
    % Clean and reconstruct grains
    [grainsDef, ebsdDeformed.grainId] = calcGrains(ebsdDeformed, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);

    % GND Processing (Filtered to remove noise artifacts)
    fprintf('  -> Computing curvature and fitting dislocations...\n');
    F = halfQuadraticFilter;
    F.alpha=3;
    ebsdSGND = gridify(smooth(ebsdDeformed('indexed'), F));
    kappa = ebsdSGND.curvature;
    
    dS = dislocationSystem.bcc(ebsdDeformed.CS);
    nu = 0.3;
    dS(dS.isEdge).u = 1 / (1 - nu);
    dS(dS.isScrew).u = 1;
    
    dSRot = ebsdSGND.orientations * dS;
    [rho, factor] = fitDislocationSystems(kappa, dSRot);
    
    ebsdDeformed = ebsdSGND; % Propagate smoothed orientations for plotting
    save(gnd_mat_path, 'ebsdDeformed', 'rho', 'factor', 'dSRot', 'dS', 'grainsDef');
end

%% -------------------------------------------------------------------------
% PART 4: Visualizing the comparison
% -------------------------------------------------------------------------
fprintf('Plotting comparison side-by-side...\n');

% 1. Find the top SSLIP system globally across the map
sslip_activity_sum = sum(abs(ebsdID.prop.slipIDcor), 2, 'omitnan');
[~, max_sslip_idx] = max(sslip_activity_sum);
sS_max = sSLocalRef(max_sslip_idx); 
sslip_data = abs(ebsdID.prop.slipIDcor(max_sslip_idx, :));

% 2. Find the top GND system
gnd_density_sum = sum(abs(rho), 1, 'omitnan');
[~, max_gnd_idx] = max(gnd_density_sum);
dS_max = dS(max_gnd_idx);
gnd_data = factor * abs(rho(:, max_gnd_idx));

figH = figure;
newMtexFigure('layout', [1, 2]);

% Left Panel: SSLIP Output
nextAxis;
[grainsUndefPlot, ~] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
plot(ebsdID, sslip_data, 'micronbar', 'off');
hold on;
plot(grainsUndefPlot.boundary, 'linewidth', 1.5, 'linecolor', 'w');
mtexColorMap viridis;
caxis([0, prctile(sslip_data(sslip_data > 0), 98)]); % Trim top 2% of extreme artifacts
mtexColorbar('title', 'Slip Activity (\gamma)');
title(sprintf('Top SSLIP System:\nb=[%s] n=(%s)', char(sS_max.b.round), char(sS_max.n.round)));

% Right Panel: GND Output
nextAxis;
plot(ebsdDeformed, log10(gnd_data), 'micronbar', 'off');
hold on;
plot(grainsDef.boundary, 'linewidth', 1.5, 'linecolor', 'w');
mtexColorMap parula;
caxis([14.5, 16.5]); % Limits mapped for typical log10(m^-2) densities
mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
typeStr = {'Edge', 'Screw'};
title(sprintf('Top GND System:\n%s, b=[%s] l=[%s]', typeStr{2 - dS_max.isEdge}, char(dS_max.b.round), char(dS_max.l.round)));

drawNow(gcm, 'figSize', 'large');

figPath = fullfile(outDir, sprintf('Comparison_SSLIP_vs_GND_%s.png', currentInc.dic_prefix));
exportScaledFigure(gcf, figPath, 'nodatestring');
fprintf('Saved comparison figure to %s\n', figPath);