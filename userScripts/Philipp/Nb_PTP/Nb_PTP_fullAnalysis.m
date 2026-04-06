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
    inc_idx = 5;
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
    % just throw an error if you want to run the analysis despite an
    % existing .mat file.
    error()
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
    sS_BCC_plain = slipSystem.bcc(ebsdUndef.CS);
    
    %restrict bcc slip systems
    sSRestriction = 1:3;
    sS_BCC_select = sS_BCC_plain(sSRestriction);
    if length(sS_BCC_plain) > length(sS_BCC_select)
        warning('Number of bcc slip systems has been reduced to %f', length(sS_BCC_select))
    end

    sS_BCC = symmetrise(sS_BCC_select, 'antipodal');
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
        
        % assuming uniaxial tension in x-direction, reconfigure slip systems (by changing signs of b and n) to be
        % "positive". i.e., slip amplitudes will be positive under uniaxial
        % tension and Schmid factors will be positive
        loadDir = yvector;
        sSLocal_g.b = sign(cos(angle(loadDir,sSLocal_g.b))) .* sSLocal_g.b;
        sSLocal_g.n = sign(sSLocal_g.SchmidFactor(loadDir)) .* sSLocal_g.n;

        % SSLIP Options
        IDoptions = struct();
        IDoptions.IDMethod = 1;
        IDoptions.filterSize = 4;
        IDoptions.coarsegrain = 4;
        IDoptions.threshResidual = 0.02;
        IDoptions.minEeff = 0.01;

        % this pre-filter by loading direction, only taking the dir that
        % has positive SF
        IDoptions.posConstr = 0;

        IDoptions.plotSSLIP = 1;
        IDoptions.saveFig = 0;
        IDoptions.plotDefGrad = 0;
        IDoptions.plotResidual = 0;
        IDoptions.NoSs = 1:length(sS_BCC);
        IDoptions.stress = stressTensor.uniaxial(loadDir);      
        
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
        
        % Attempt to load and apply the specimen mask
        [ang_dir, ang_name, ~] = fileparts(ang_path);
        mask_file = fullfile(ang_dir, sprintf('%s_specimen_mask.mat', ang_name));

        if exist(mask_file, 'file')
            fprintf('Applying specimen mask to remove vacuum artifacts...\n');
            load(mask_file, 'specimenlocation');
            ebsdDeformed.prop.specimenlocation = specimenlocation;
            % Filter the EBSD map completely before grain reconstruction
            ebsdDeformed(ebsd.prop.specimenlocation) = 'notIndexed';
        else
            warning('No specimen mask found for %s. Proceeding with unmasked data.', currentInc.name);
        end

        % Now proceed to calcGrains(ebsd, ...)


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
plot_layout = [2, 3];
num_top = prod(plot_layout);

% --- MODIFIED SSLIP PLOTTING ---
% This section is modified to plot top systems PER GRAIN, not globally.

% 1. Reconstruct grains and identify left/right bicrystal grains
[grainsUndefPlot, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
[~, sorted_grain_idx] = sort(grainsUndefPlot.numPixel, 'descend');
top2_idx = sorted_grain_idx(1:min(2, length(grainsUndefPlot)));
[~, left_to_right] = sort(grainsUndefPlot(top2_idx).centroid.x);
bicrystal_idx = top2_idx(left_to_right);

left_g_id = grainsUndefPlot(bicrystal_idx(1)).id;
right_g_id = grainsUndefPlot(bicrystal_idx(2)).id;

mask_g1 = (ebsdID.grainId == left_g_id);
mask_g2 = (ebsdID.grainId == right_g_id);

% 2. Find top active slip systems for each grain individually
slip_mags_g1 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g1(:))), 2, 'omitnan');
slip_mags_g2 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g2(:))), 2, 'omitnan');

[~, sorted_g1_idx] = sort(slip_mags_g1, 'descend');
[~, sorted_g2_idx] = sort(slip_mags_g2, 'descend');

top_sslip_g1_idx = sorted_g1_idx(1:num_top);
top_sslip_g2_idx = sorted_g2_idx(1:num_top);

% 3. Calculate Schmid Factors and Ranks for all systems in each grain
stress = stressTensor.uniaxial(yvector); % As defined in the SSLIP calculation part
ori_g1 = grainsUndefPlot(left_g_id).meanOrientation;
ori_g2 = grainsUndefPlot(right_g_id).meanOrientation;

sSLocal_g1 = ori_g1 * sSLocalRef;
sSLocal_g2 = ori_g2 * sSLocalRef;

SF_g1 = abs(sSLocal_g1.SchmidFactor(stress));
SF_g2 = abs(sSLocal_g2.SchmidFactor(stress));

[~, SF_g1_rank_order] = sort(SF_g1, 'descend');
[~, SF_g2_rank_order] = sort(SF_g2, 'descend');

% Find max SF for delta calculation
max_sf_g1 = max(SF_g1);
max_sf_g2 = max(SF_g2);

% 4. Determine a unified color limit for the plots
all_slip_data = [];
for r = 1:num_top
    d1 = abs(ebsdID.prop.slipIDcor(top_sslip_g1_idx(r), mask_g1(:)));
    d2 = abs(ebsdID.prop.slipIDcor(top_sslip_g2_idx(r), mask_g2(:)));
    all_slip_data = [all_slip_data; d1(:); d2(:)];
end
global_sslip_max = prctile(all_slip_data(all_slip_data > 0), 98);
if isempty(global_sslip_max) || isnan(global_sslip_max) || global_sslip_max <= 0
    global_sslip_max = 1e-4;
end

% 5. Find the top GND systems (unchanged)
if ~isempty(currentInc.ang_file)
    gnd_density_sum = sum(abs(rho), 1, 'omitnan');
    [~, sorted_gnd_idx] = sort(gnd_density_sum, 'descend');
    top_gnd_idx = sorted_gnd_idx(1:num_top);
end

% Figure 1: Top SSLIP systems per grain
figSSLIP = figure;
newMtexFigure('layout', plot_layout);
cmap_viridis = viridis(256);
light_grey = [0.8 0.8 0.8];

for i = 1:num_top
    row_idx = ceil(i / plot_layout(2));
    col_idx = mod(i - 1, plot_layout(2)) + 1;
    nextAxis(row_idx, col_idx);
    
    % Get the i-th ranked active system for each grain
    sys_g1_idx = top_sslip_g1_idx(i);
    sys_g2_idx = top_sslip_g2_idx(i);
    
    % Create a combined plot map
    plot_data = nan(size(ebsdID));
    plot_data(mask_g1) = abs(ebsdID.prop.slipIDcor(sys_g1_idx, mask_g1));
    plot_data(mask_g2) = abs(ebsdID.prop.slipIDcor(sys_g2_idx, mask_g2));
    
    % --- New plotting logic with grey background ---
    % 1. Create base image (light grey for bicrystal, white for rest)
    final_rgb = ones(length(ebsdID), 3); % Start with white
    bicrystal_mask = ~isnan(ebsdID.prop.residualEeff);
    final_rgb(bicrystal_mask, :) = repmat(light_grey, sum(bicrystal_mask), 1);

    % 2. Normalize slip data to get colormap indices
    norm_data = plot_data / global_sslip_max;
    norm_data(isnan(norm_data)) = 0;
    norm_data(norm_data > 1) = 1;
    norm_data(norm_data < 0) = 0;
    color_indices = round(norm_data * 255) + 1;

    % 3. Find pixels with activity and apply viridis color
    is_overlay_pixel = (plot_data > 0);
    overlay_colors = cmap_viridis(color_indices(is_overlay_pixel), :);
    final_rgb(is_overlay_pixel, :) = overlay_colors;

    % 4. Plot the composite RGB image
    plot(ebsdID, final_rgb, 'micronbar', 'off');
    % --- End of new plotting logic ---

    hold on;
    plot(grainsUndefPlot.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    % Get slip systems in local grain coordinates for plotting
    sS_g1_local = grainsUndefPlot(left_g_id).meanOrientation * sSLocalRef(sys_g1_idx);
    sS_g2_local = grainsUndefPlot(right_g_id).meanOrientation * sSLocalRef(sys_g2_idx);
    
    % Plot quivers for each grain's specific system
    quiver(grainsUndefPlot(left_g_id), sS_g1_local.trace, 'color', 'r', 'linewidth', 2);
    b1 = sS_g1_local.b.normalize;
    b1.antipodal = true;
    quiver(grainsUndefPlot(left_g_id), b1, 'color', 'w', 'linewidth', 2);
    quiver(grainsUndefPlot(right_g_id), sS_g2_local.trace, 'color', 'r', 'linewidth', 2);
    b2 = sS_g2_local.b.normalize;
    b2.antipodal = true;
    quiver(grainsUndefPlot(right_g_id), b2, 'color', 'w', 'linewidth', 2);
    
    % Build the detailed title string
    sS_g1_ref = sSLocalRef(sys_g1_idx);
    sS_g2_ref = sSLocalRef(sys_g2_idx);
    
    sf_val_g1 = SF_g1(sys_g1_idx);
    sf_val_g2 = SF_g2(sys_g2_idx);
    
    sf_rank_g1 = find(SF_g1_rank_order == sys_g1_idx, 1);
    sf_rank_g2 = find(SF_g2_rank_order == sys_g2_idx, 1);
    
    delta_sf_g1 = max_sf_g1 - sf_val_g1;
    delta_sf_g2 = max_sf_g2 - sf_val_g2;
    
    mag_g1 = slip_mags_g1(sys_g1_idx);
    mag_g2 = slip_mags_g2(sys_g2_idx);
    
    title_str1 = sprintf('Left (G%d) Sys %d: b=[%s] n=(%s), Mag=%.4f, SF=%.2f (Rank %d, \\DeltaSF=%.2f)', ...
        left_g_id, sys_g1_idx, char(sS_g1_ref.b.round), char(sS_g1_ref.n.round), mag_g1, sf_val_g1, sf_rank_g1, delta_sf_g1);
    title_str2 = sprintf('Right (G%d) Sys %d: b=[%s] n=(%s), Mag=%.4f, SF=%.2f (Rank %d, \\DeltaSF=%.2f)', ...
        right_g_id, sys_g2_idx, char(sS_g2_ref.b.round), char(sS_g2_ref.n.round), mag_g2, sf_val_g2, sf_rank_g2, delta_sf_g2);
    
    title(sprintf('SSLIP Activity Rank %d\n%s\n%s', i, title_str1, title_str2), 'FontSize', 8);
end

% Manually add a colorbar as we are plotting RGB data
cax = axes('Parent', figSSLIP, 'Position', [0.92 0.15 0.02 0.7]);
colormap(cax, cmap_viridis);
caxis(cax, [0, global_sslip_max]);
cb = colorbar(cax);
cb.Label.String = 'Slip Activity (\gamma)';
set(cax, 'Visible', 'off');

drawNow(gcm, 'figSize', 'large');

figPathSSLIP = fullfile(outDir, sprintf('Top_SSLIP_perGrain_%s.jpg', currentInc.dic_prefix));
exportScaledFigure(gcf, figPathSSLIP, 'target_px', 3500, 'noDateString');
fprintf('Saved SSLIP per-grain figure to %s\n', figPathSSLIP);

% Figure 2: Top GND systems
if ~isempty(currentInc.ang_file)
    figGND = figure;
    newMtexFigure('layout', [2, 1]);
    
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

% Figure 3: Screw and Edge GNDs
if ~isempty(currentInc.ang_file)
    fig_ES_GND = figure;
    newMtexFigure('layout', [1, 2]);
    
    % one plot with the accumulative density of screw dislocations dS.isScrew
    % the other with the accumulative density of edge dislocations
    % dS.isEdge

    GND_E = factor * sum(abs(rho(:, dS.isEdge)), 2, 'omitnan');
    GND_S = factor * sum(abs(rho(:, dS.isScrew)), 2, 'omitnan');

    nextAxis(1, 1);
    plot(ebsdDeformed, log10(GND_E), 'micronbar', 'off');
    hold on;
    plot(grainsDef.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    mtexColorMap parula;
    title('Total Edge GND Density');

    nextAxis(1, 2);
    plot(ebsdDeformed, log10(GND_S), 'micronbar', 'off');
    hold on;
    plot(grainsDef.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    mtexColorMap parula;
    title('Total Screw GND Density');
    
    setColorRange([14.5, 16.5]); % Limits mapped for typical log10(m^-2) densities
    mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
    drawNow(gcm, 'figSize', 'large');
    
    figPathGND_ES = fullfile(outDir, sprintf('Edge_Screw_GND_%s.jpg', currentInc.dic_prefix));
    exportScaledFigure(gcf, figPathGND_ES, 'target_px', 3500, 'noDateString');
    fprintf('Saved Edge and Screw GND figure to %s\n', figPathGND_ES);
end