% N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\DIC exported data\Processing_subset=21_step=5%% Re_PTP3_prettyIPF_loop.m
% Script to loop over Re_PTP_3 .ang files and generate pretty IPF Y maps
% overlaid on a tunable greyscale Confidence Index (CI) background.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Tuning Parameters
baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
outDir = fullfile(baseDir, 'PrettyPlots');

if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output directory: %s\n', outDir);
end

% --- "PHOTOSHOP CURVES" TUNING PARAMETERS ---
% 1. Black Point (0.0 to 1.0): CI values below this will be pure black.
blackPoint = 0.05; 

% 2. White Point (0.0 to 1.0): CI values above this will be pure white.
whitePoint = 0.98; 

% 3. Gamma (Midtones): 
%    < 1.0 makes the mid-range lighter (good if the map is too dark)
%    > 1.0 makes the mid-range darker (good if the map looks washed out)
%    = 1.0 applies no gamma curve.
gammaVal = 1.5; 

% 4. IPF Overlay Opacity (0.0 to 1.0)
ipf_alpha = 0.55;
% --------------------------------------------

%% 3. Load Data & Loop over .ang files
% Define Crystal Symmetry for Rhenium
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

tkd_tilt = -7 * degree; % Adjust this depending on your exact TKD geometry setup if needed
thresh_ci = 0.1;

% Find all .ang files matching the insitu Re sequence (searches recursively)
angFiles = dir(fullfile(baseDir, '**', 'Re_7deg_insitu Specimen*.ang'));

if isempty(angFiles)
    error('No .ang files found matching "Re_7deg_insitu Specimen*.ang" in %s', baseDir);
end

for i = 1:length(angFiles)
    currentFile = angFiles(i).name;
    ang_path = fullfile(angFiles(i).folder, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    fprintf('Processing %d/%d: %s\n', i, length(angFiles), currentFile);
    
    % Load EBSD Data
    ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % Apply standard rotations for TKD alignment
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Base filtering
    ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
    

    ebsd = ebsd.gridify;
    % --- Calculate Tuned CI Greyscale (Curves Adjustment) ---
    ci_raw = ebsd.ci;
    ci_min = min(ci_raw(:));
    ci_max = max(ci_raw(:));
    
    ci_norm = (ci_raw - ci_min) ./ (ci_max - ci_min);
    ci_tuned = (ci_norm - blackPoint) ./ (whitePoint - blackPoint);
    ci_tuned(ci_tuned < 0) = 0; % clip darks
    ci_tuned(ci_tuned > 1) = 1; % clip highlights
    ci_tuned = ci_tuned .^ gammaVal; % Apply gamma
    
    cibw = [ci_tuned(:), ci_tuned(:), ci_tuned(:)];
    cibw(~ebsd.isIndexed, :) = 0; % Mask out the background (black)


    % --- Apply Specimen Mask ---
    % Look for a specimen mask .mat file containing 'specimenlocation'
    mask_file = fullfile(angFiles(i).folder, sprintf('%s_specimen_mask.mat', baseName));
    
    if exist(mask_file, 'file')
        
        fprintf('  -> Applying specimen mask...\n');
        load(mask_file, 'specimenlocation');
        ebsd(~specimenlocation) = 'notIndexed';
        
        % --- Crop to the bounding box of the specimen location ---
        valid_ebsd = ebsd(specimenlocation);
        if length(valid_ebsd) > 0
            xmin = min(valid_ebsd.x); xmax = max(valid_ebsd.x);
            ymin = min(valid_ebsd.y); ymax = max(valid_ebsd.y);
            ebsd = ebsd(inpolygon(ebsd, [xmin, ymin, xmax-xmin, ymax-ymin]));
        end
    else
        warning('  -> No specimen mask found at %s. Proceeding with unmasked data.', mask_file);
    end
    
    % --- Plotting IPF Y ---
    ipf = ipfColorKey(ebsd.CS);
    ipf.inversePoleFigureDirection = yvector;
    
    fig = figure('Name', currentFile, 'Visible', 'off');
    plot(ebsd, cibw, 'micronbar', 'off');
    hold on;

    colors = ipf.orientation2color(ebsd('indexed').orientations);
    plot(ebsd('indexed'), colors, 'FaceAlpha', ipf_alpha);
    
    title(['IPF Y - ', baseName], 'Interpreter', 'none');
    
    % --- Export Figure ---
    out_name = fullfile(outDir, sprintf('%s_Y_IPF_TunedCI_noMicron.jpg', baseName)); %sprintf('%s_Y_IPF_TunedCI.jpg', baseName)
    exportScaledFigure(gcf, out_name, 'noDateString');
    close(fig);
    
    fprintf('  -> Saved plot to %s\n', out_name);
    
    % --- Twin Network Analysis ---
    % Project orientations to Fundamental Region
    reforientation = orientation.byEuler(0*degree, 0*degree, 0*degree, CS);
    ebsd('indexed').orientations = ebsd('indexed').orientations.project2FundamentalRegion(reforientation);
    
    % Reconstruct Grains
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 5*degree, 'minPixel', 5);
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 5*degree, 'minPixel', 5);
    grains = smooth(grains, 1);
    
    % Extract only internal boundaries between indexed grains
    gB = grains.boundary(CS.mineral, CS.mineral);

    % Define {11-21} Twin System for Re
    tS_1121 = twinSystem.hexagonal_1121(CS);

    % Initialize twinNetwork with maxGen = 1
    network = twinNetwork(grains, gB, tS_1121, ebsd('indexed'), 'thresh', 8*degree, 'maxGen', 1);

    % Step 1: Building Initial Topology
    disp('    Step 1: Building Initial Topology...');
    network.buildInitialGraph('maxGen', 1);
    if ~isempty(network.Graph) && numedges(network.Graph) > 0
        network.plotGraph(network.Graph);
        figG = gcf;
        title(sprintf('[%s] Step 1: Initial Raw Topology', baseName), 'Interpreter', 'none');
        exportScaledFigure(figG, fullfile(outDir, sprintf('%s_01_Raw_Topology.png', baseName)), 'noDateString');
        close(figG);
    end

    % Step 2: Extracting Spanning Tree
    disp('    Step 2: Extracting Spanning Tree...');
    network.extractSpanningTree();
    if ~isempty(network.mstGraph) && numedges(network.mstGraph) > 0
        network.plotGraph(network.mstGraph);
        figMST = gcf;
        title(sprintf('[%s] Step 2: Minimum Spanning Tree', baseName), 'Interpreter', 'none');
        exportScaledFigure(figMST, fullfile(outDir, sprintf('%s_02_Minimum_Spanning_Tree.png', baseName)), 'noDateString');
        close(figMST);
    end

    % Step 3: Resolving Node States and Variants
    disp('    Step 3: Resolving Node States and Variants...');
    network.resolveNodeStates();
    if ~isempty(network.mstGraph) && numedges(network.mstGraph) > 0
        network.plotGraph(network.mstGraph);
        figRes = gcf;
        title(sprintf('[%s] Step 3: MST after FZ Propagation', baseName), 'Interpreter', 'none');
        exportScaledFigure(figRes, fullfile(outDir, sprintf('%s_03_Resolved_Node_States.png', baseName)), 'noDateString');
        close(figRes);
    end

    % Step 4: Pruning and Rewiring Outliers
    disp('    Step 4: Pruning and Rewiring Outliers...');
    network.pruneAndRewireOutliers('outlierMode', 'theoretical', 'maxErrThreshold', 25*degree, 'debug', true, 'rejectGenIncrease', true);
    if ~isempty(network.mstGraph) && numedges(network.mstGraph) > 0
        network.plotGraph(network.mstGraph);
        figPrune = gcf;
        title(sprintf('[%s] Step 4: Pruned Graph Topology', baseName), 'Interpreter', 'none');
        exportScaledFigure(figPrune, fullfile(outDir, sprintf('%s_04_Pruned_Rewired.png', baseName)), 'noDateString');
        close(figPrune);
    end

    % Step 5: Reducing to Quotient Graph
    disp('    Step 5: Reducing to Quotient Graph...');
    network.reduceToQuotient();
    if ~isempty(network.QuotientGraph) && numedges(network.QuotientGraph) > 0
        network.plotGraph(network.QuotientGraph, 'noGrainsBG', true);
        figQ = gcf;
        title(sprintf('[%s] Step 5: Final Quotient Graph', baseName), 'Interpreter', 'none');
        exportScaledFigure(figQ, fullfile(outDir, sprintf('%s_05_Quotient_Graph.png', baseName)), 'noDateString');
        close(figQ);
    end

    % Step 6: Plotting grains with twin generation assignment (with text labels)
    disp('    Step 6: Plotting grains with twin generation assignment...');
    network.plot('twinsOnly', false, 'noText', false);
    figGen = gcf;
    title(sprintf('[%s] Twin Generations & Network Map', baseName), 'Interpreter', 'none');
    out_name_twins = fullfile(outDir, sprintf('%s_TwinVariants_twinNetwork.jpg', baseName));
    exportScaledFigure(figGen, out_name_twins, 'noDateString');
    close(figGen);
    fprintf('  -> Saved twin network plot to %s\n', out_name_twins);

    % Step 7: Plotting grains with twin generation assignment (without text labels)
    disp('    Step 7: Plotting grains with twin generation assignment without text...');
    network.plot('twinsOnly', false, 'noText', true, 'lineWidth', 2);
    figGenNoText = gcf;
    title(sprintf('[%s] Twin Generations & Network Map (No Text)', baseName), 'Interpreter', 'none');
    out_name_notext = fullfile(outDir, sprintf('%s_TwinVariants_twinNetwork_noText.jpg', baseName));
    exportScaledFigure(figGenNoText, out_name_notext, 'noDateString');
    close(figGenNoText);
    fprintf('  -> Saved twin network (no text) plot to %s\n', out_name_notext);

    % % --- (Legacy) Twin Variant Identification ---
    % % Project orientations to Fundamental Region
    % reforientation = orientation.byEuler(0*degree, 0*degree, 0*degree, CS);
    % ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
    % 
    % % Reconstruct Grains
    % [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 5*degree, 'minPixel', 5);
    % [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 5*degree, 'minPixel', 5);
    % % grains = smooth(grains, 5);
    % 
    % % Identify Parent Grain & Twin Systems
    % [~, parentId] = max(grains.numPixel);
    % parentGrain = grains(parentId);
    % parentGrains = grains(angle(parentGrain.meanOrientation, grains.meanOrientation) < 8*degree);
    % 
    % tS_1121 = twinSystem.hexagonal_1121(CS);
    % tS_sym = tS_1121.symmetrise;
    % 
    % % Calculate primary and secondary variants
    % primaryTwin = parentGrain.meanOrientation * tS_sym;
    % secondTwin = primaryTwin * tS_sym;
    % 
    % primary_twin_oris = primaryTwin.orientation;
    % secondary_twin_oris = secondTwin.orientation;
    % 
    % tolerance = 13 * degree;
    % isPrimary = false(length(grains), 1);
    % isSecondary = false(length(grains), 1);
    % primaryVariantId = zeros(length(grains), 1);
    % secondaryVariantId = zeros(length(grains), 2);
    % 
    % for j = 1:length(grains)
    %     if any(grains(j).id == parentGrains.id), continue; end
    %     [minAngle, bestFit] = min(angle(grains(j).meanOrientation, primary_twin_oris));
    %     if minAngle < tolerance
    %         isPrimary(j) = true;
    %         primaryVariantId(j) = primaryTwin(bestFit).variantId;
    %     end
    % end
    % 
    % for j = 1:length(grains)
    %     if any(grains(j).id == parentGrains.id) || isPrimary(j), continue; end
    %     [minAngle, bestFit] = min(angle(grains(j).meanOrientation, secondary_twin_oris));
    %     if minAngle < tolerance
    %         isSecondary(j) = true;
    %         secondaryVariantId(j,:) = [secondTwin(bestFit).parent.variantId, secondTwin(bestFit).variantId];
    %     end
    % end
    % 
    % % --- Plotting Twin Map ---
    % variant_alpha = 1;
    % fig_twins = figure('Name', ['Twin Variants - ' currentFile], 'Visible', 'off');
    % % plot(ebsd, cibw, 'micronbar', 'on');
    % % hold on;
    % 
    % plot(parentGrains, 'facecolor', [0.8 0.8 0.8], 'facealpha', variant_alpha);
    % 
    % varColors = lines(length(tS_sym));
    % 
    % primaryGrains = grains(isPrimary);
    % varIds = primaryVariantId(isPrimary);
    % uniqueVars = unique(varIds);
    % for v_idx = 1:length(uniqueVars)
    %     variant = uniqueVars(v_idx);
    %     if variant == 0, continue; end
    %     g_var = primaryGrains(varIds == variant);
    %     hold on;
    %     plot(g_var, 'facecolor', varColors(variant,:), 'facealpha', variant_alpha);
    % end
    % 
    % secondaryGrains = grains(isSecondary);
    % if ~isempty(secondaryGrains)
    %     hold on
    %     plot(secondaryGrains, 'facecolor', 'w', 'facealpha', variant_alpha);
    %     secIds = secondaryVariantId(isSecondary, :);
    %     %labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    %     %text(secondaryGrains, labels, 'HorizontalAlignment', 'center');
    % end
    % 
    % hold on;
    % plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'k');
    % 
    % title(['Twin Variants - ', baseName], 'Interpreter', 'none');
    % 
    % % --- Export Twin Map Figure ---
    % out_name_twins = fullfile(outDir, sprintf('%s_TwinVariants_TunedCI.jpg', baseName));
    % exportScaledFigure(fig_twins, out_name_twins, 'noDateString');
    % close(fig_twins);
    % 
    % fprintf('  -> Saved twin variants plot to %s\n', out_name_twins);
end

fprintf('\nFinished generating pretty IPF plots for Re_PTP_3!\n');
