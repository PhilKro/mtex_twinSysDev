%% Nb_PTP_fullAnalysis.m
% This script compares the outputs of the SSLIP determination of dislocation 
% slip (based on undeformed reference TKD + DIC) with GND calculations from 
% deformed TKD scans at identical increments.
%
% It loads previously processed data if available to save computational time.

if ~exist('is_master_run', 'var')
    clear
end
close all

%% Initialization & Mtex Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Directories and Input Setup
% Use pwd to create an absolute path and bypass OneDrive relative path bugs
if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    warning('No master_pwd exists taking native one.')
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir  = fullfile(baseDir, 'Processed_Output');
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created missing directory: %s\n', outDir);
end

% Provide the name of the MAT file containing the fully merged EBSD & DIC 
% data from the Undeformed script
if ~exist('undeformed_dic_mat', 'var')
    undeformed_dic_mat = fullfile(outDir, 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');
end

if ~exist('is_master_run', 'var')
    % Skip block if this is called by the master runner script 
    % Define increments (mapping TKD .ang files to concurrent DIC .hdf5 prefixes)
    increments = struct();
    
    increments(1).name = 'Increment 0';
    increments(1).ang_file = 'insitu Specimen 1 0 Map Data 1_BW263.ang';
    increments(1).dic_prefix = 'B0001';
    
    increments(2).name = 'Increment 1';
    increments(2).ang_file = '';
    increments(2).dic_prefix = 'B0002';
    
    increments(3).name = 'Increment 2';
    increments(3).ang_file = 'insitu Specimen 1 1 Map Data 2_BW263.ang';
    increments(3).dic_prefix = 'B0003';
    
    increments(4).name = 'Increment 3';
    increments(4).ang_file = '';
    increments(4).dic_prefix = 'B0004';
    
    increments(5).name = 'Increment 4';
    increments(5).ang_file = 'insitu Specimen 1 2 Map Data 3_BW263.ang';
    increments(5).dic_prefix = 'B0005';
    
    increments(6).name = 'Increment 5';
    increments(6).ang_file = '';
    increments(6).dic_prefix = 'B0006';
    
    increments(7).name = 'Increment 6';
    increments(7).ang_file = 'insitu Specimen 1 3 Map Data 4_BW263.ang';
    increments(7).dic_prefix = 'B0007';
    
    increments(8).name = 'Increment 7';
    increments(8).ang_file = '';
    increments(8).dic_prefix = 'B0008';
    
    increments(9).name = 'Increment 8';
    increments(9).ang_file = 'insitu Specimen 1 4 Map Data 5_BW263.ang';
    increments(9).dic_prefix = 'B0009';
    
    increments(10).name = 'Increment 9';
    increments(10).ang_file = '';
    increments(10).dic_prefix = 'B0010';
    
    increments(11).name = 'Increment 10';
    increments(11).ang_file = 'insitu Specimen 1 5 Map Data 6_BW263.ang';
    increments(11).dic_prefix = 'B0011';
    
    increments(12).name = 'Increment 11';
    increments(12).ang_file = '';
    increments(12).dic_prefix = 'B0012';
    
    increments(13).name = 'Increment 12';
    increments(13).ang_file = '';
    increments(13).dic_prefix = 'B0013';
    
    increments(14).name = 'Increment 13';
    increments(14).ang_file = 'insitu Specimen 1 7 Map Data 7_BW263.ang';
    increments(14).dic_prefix = 'B0014';
end

% Choose which increment to analyze right now (1 through 14)
if ~exist('is_master_run', 'var')
    inc_idx = 2;
    currentInc = increments(inc_idx);
end

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

% if exist(sslip_mat_path, 'file')
try 
    % error()
    load(sslip_mat_path, 'ebsdID', 'sSLocalRef', 'optOut');
    fprintf('Loading existing SSLIP data from %s...\n', sslip_mat_path);
%else
catch
    fprintf('Calculating SSLIP data for %s...\n', currentInc.name);
    
    if ~exist(undeformed_dic_mat, 'file') && ~exist('is_master_run', 'var') 
        error('Undeformed EBSD with DIC data not found at %s.', undeformed_dic_mat);
    end
    
    % Load aligned base data
    if ~exist('is_master_run', 'var')
        undefData = load(undeformed_dic_mat, 'newEBSD');
        ebsdUndef = undefData.newEBSD;
    else
        ebsdUndef = newEBSD;
    end
    
    % Clean and reconstruct grains
    ebsdUndef(ebsdUndef.ci < thresh_ci) = 'notIndexed';
    [grainsUndef, ebsdUndef.grainId] = calcGrains(ebsdUndef, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
    
    % Setup blank container for the fully combined SSLIP field
    ebsdID = ebsdUndef.gridify; 
    sS_BCC = symmetrise(slipSystem.bcc(ebsdUndef.CS), 'antipodal');
    ebsdID.prop.slipIDcor = zeros(length(sS_BCC), length(ebsdID));
    ebsdID.prop.residualEeff = nan(length(ebsdID), 1); % MTEX prefers column vectors for properties
    
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
        
        % Stitch local map back into the full ebsd grid using nearest neighbor 
        % mapping, since SSLIP's coarse-graining changes the grid coordinates
        % and truncates edges.
        in_grain_idx = find(ebsdID.grainId == gId);
        xq = ebsdID.x(in_grain_idx);
        yq = ebsdID.y(in_grain_idx);
        
        xc = ebsdID_g.x(:);
        yc = ebsdID_g.y(:);
        
        idx_nearest = knnsearch([xc, yc], [xq, yq]);
        ebsdID.prop.slipIDcor(:, in_grain_idx) = ebsdID_g.prop.slipIDcor(:, idx_nearest);
        
        resEeff_g = ebsdID_g.prop.residualEeff(:); % Force column vector to match dimensionality
        ebsdID.prop.residualEeff(in_grain_idx) = resEeff_g(idx_nearest);
    end
    
    save(sslip_mat_path, 'ebsdID', 'sSLocalRef', 'optOut');
end

%% -------------------------------------------------------------------------
% PART 3: GND Calculation on Pristine Deformed TKD Data
% -------------------------------------------------------------------------
if ~isempty(currentInc.ang_file)
    gnd_mat_path = fullfile(outDir, sprintf('GND_Data_%s.mat', currentInc.dic_prefix));
    
    % if exist(gnd_mat_path, 'file')
    try
        % error()
        load(gnd_mat_path, 'ebsdDeformed', 'rho', 'factor', 'dSRot', 'dS', 'grainsDef');
        fprintf('Loading existing GND data from %s...\n', gnd_mat_path);
    % else
    catch
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
else
    fprintf('No ANG file provided for %s. Skipping GND calculation.\n', currentInc.name);
end

%% -------------------------------------------------------------------------
% PART 4: Visualizing the comparison
% -------------------------------------------------------------------------
fprintf('Plotting comparison side-by-side...\n');

% Define the layout shape for each group (e.g., [2, 2] means a 2x2 block 
% of SSLIPs next to a 2x2 block of GNDs, creating a total 2x4 field).
plot_layout = [2, 2];
num_top = prod(plot_layout);

% 1. Find the top SSLIP systems globally across the map
sslip_activity_sum = sum(abs(ebsdID.prop.slipIDcor), 2, 'omitnan');
[~, sorted_sslip_idx] = sort(sslip_activity_sum, 'descend');
top_sslip_idx = sorted_sslip_idx(1:num_top);

% Determine global color limit for the top SSLIP systems to unify their colorbars
sslip_max_vals = zeros(num_top, 1);
for i = 1:num_top
    ssdata = abs(ebsdID.prop.slipIDcor(top_sslip_idx(i), :));
    sslip_max_vals(i) = prctile(ssdata(ssdata > 0), 98);
end
global_sslip_max = max(sslip_max_vals);
if isnan(global_sslip_max) || global_sslip_max <= 0, global_sslip_max = 1e-4; end

% 2. Find the top GND systems
if ~isempty(currentInc.ang_file)
    gnd_density_sum = sum(abs(rho), 1, 'omitnan');
    [~, sorted_gnd_idx] = sort(gnd_density_sum, 'descend');
    top_gnd_idx = sorted_gnd_idx(1:num_top);
end

% Figure 1: Top SSLIP systems
figSSLIP = figure;
newMtexFigure('layout', plot_layout);

[grainsUndefPlot, ~] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

for i = 1:num_top
    row_idx = ceil(i / plot_layout(2));
    col_idx = mod(i - 1, plot_layout(2)) + 1;
    
    nextAxis(row_idx, col_idx);
    sS_max = sSLocalRef(top_sslip_idx(i)); 
    sslip_data = abs(ebsdID.prop.slipIDcor(top_sslip_idx(i), :));
    
    plot(ebsdID, sslip_data, 'micronbar', 'off');
    hold on;
    plot(grainsUndefPlot.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    
    % Plot slip trace and slip direction rotated for each local grain
    sSLocal_grains = grainsUndefPlot.meanOrientation * sS_max;
    quiver(grainsUndefPlot, sSLocal_grains.trace, 'color', 'r', 'linewidth', 2);
    quiver(grainsUndefPlot, sSLocal_grains.b.normalize, 'color', 'w', 'linewidth', 2);
    
    mtexColorMap viridis;
    
    title(sprintf('Rank %d SSLIP System:\nb=[%s] n=(%s)', i, char(sS_max.b.round), char(sS_max.n.round)));
end

setColorRange([0, global_sslip_max]);
mtexColorbar('title', 'Slip Activity (\gamma)');
drawNow(gcm, 'figSize', 'large');

figPathSSLIP = fullfile(outDir, sprintf('Top_SSLIP_%s.jpg', currentInc.dic_prefix));
exportScaledFigure(gcf, figPathSSLIP, 'target_px', 3500, 'noDateString');
fprintf('Saved SSLIP figure to %s\n', figPathSSLIP);

% Figure 2: Top GND systems
if ~isempty(currentInc.ang_file)
    figGND = figure;
    newMtexFigure('layout', plot_layout);
    
    for i = 1:num_top
        row_idx = ceil(i / plot_layout(2));
        col_idx = mod(i - 1, plot_layout(2)) + 1;
        
        nextAxis(row_idx, col_idx);
        dS_max = dS(top_gnd_idx(i));
        gnd_data = factor * abs(rho(:, top_gnd_idx(i)));
        
        plot(ebsdDeformed, log10(gnd_data), 'micronbar', 'off');
        hold on;
        plot(grainsDef.boundary, 'linewidth', 1.5, 'linecolor', 'w');
        mtexColorMap parula;
        
        typeStr = {'Edge', 'Screw'};
        title(sprintf('Rank %d GND System:\n%s, b=[%s] l=[%s]', i, typeStr{2 - dS_max.isEdge}, char(dS_max.b.round), char(dS_max.l.round)));
    end
    
    setColorRange([14.5, 16.5]); % Limits mapped for typical log10(m^-2) densities
    mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
    drawNow(gcm, 'figSize', 'large');
    
    figPathGND = fullfile(outDir, sprintf('Top_GND_%s.jpg', currentInc.dic_prefix));
    exportScaledFigure(gcf, figPathGND, 'target_px', 3500, 'noDateString');
    fprintf('Saved GND figure to %s\n', figPathGND);
end