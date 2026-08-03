%% forwardDeformed_04_TwinNetworkAnalysis.m
% Performs TwinNetwork2 graph-based topology and generation analysis on Deformed Forward TKD datasets.
% Integrates crystallographic orientation identification with {10-12} deformation twinning variant analysis in HCP Titanium.
% Operates directly on EBSD/TKD orientation maps (DIC strain data is optional and not required).
% Outputs results to: Processed_Output_Deformed_Forward/Twin_Analysis_Results.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\geometry');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Load Registry & Configure Step Selection
[registry, defaults] = a00_Dataset_Registry();

% Filter for non-excluded deformed steps with valid TKD maps (DIC presence optional)
defSteps = registry(~[registry.isReference] & ~[registry.exclude] & ~cellfun(@isempty, {registry.tkd_fname}));

if isempty(defSteps)
    fprintf('========================================================================\n');
    fprintf('NOTICE: No valid deformed TKD step entries found in a00_Dataset_Registry().\n');
    fprintf('========================================================================\n');
    return;
end

inDir  = fullfile(defaults.dir, defaults.outDir_deformed);
outDir = fullfile(defaults.dir, defaults.outDir_deformed, 'Twin_Analysis_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING DEFORMED FORWARD TWIN NETWORK BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

%% 2. Iterate Over Each Deformed Step
for stepIdx = 1:numel(defSteps)
    defStep = defSteps(stepIdx);
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('Processing Step [%s] (%d/%d) - TKD Map: %s\n', defStep.stepName, stepIdx, numel(defSteps), defStep.tkd_fname);
    fprintf('------------------------------------------------------------------------\n');
    
    ebsd = [];
    grains = [];
    
    % Attempt to load precomputed aligned EBSD object first
    matPath = fullfile(inDir, [defStep.stepName, '_deformedEBSD_withForwardDIC.mat']);
    if exist(matPath, 'file')
        fprintf('  -> Loading aligned EBSD dataset: %s...\n', matPath);
        dataLoaded = load(matPath);
        ebsd = dataLoaded.newEBSD;
        if isfield(dataLoaded, 'grains'), grains = dataLoaded.grains; end
    else
        fprintf('  -> Aligned MAT file not found. Loading TKD scan directly from .ang file...\n');
        ang_path = fullfile(defaults.dir, defStep.tkd_fname, [defStep.tkd_fname, '.ang']);
        if ~exist(ang_path, 'file')
            warning('EBSD scan not found at:\n  %s\nSkipping Step [%s].', ang_path, defStep.stepName);
            continue;
        end
        CS = crystalSymmetry('6/mmm', [2.95 2.95 4.68], 'x||a', 'mineral', 'Ti');
        ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
        ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + defaults.tkd_tilt), 'keepXY');
        ebsd(ebsd.ci < defaults.thresh_ci) = 'notIndexed';
        ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
        
        mask_path = fullfile(defaults.dir, defStep.mask_fname);
        if exist(mask_path, 'file')
            load(mask_path, 'specimenlocation');
            try, ebsd(~specimenlocation) = 'notIndexed'; catch, end
        end
        ebsd = ebsd.gridify;
    end
    
    %% 3. Retrieve or Reconstruct Grains & Extract Boundaries
    grainThresh = 5*degree; 
    minSize = 300; % Threshold in pixels to isolate significant grains
    alpha = 1;

    if ~isempty(grains) && isa(grains, 'grain2d')
        fprintf('  -> Successfully loaded precomputed grains directly from dataset.\n');
    else
        fprintf('  -> Calculating grains dynamically (thresh = %.1f deg, minSize = %d px)...\n', grainThresh/degree, minSize);
        ebsd(ebsd.ci < defaults.thresh_ci) = 'notIndexed';
        [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha);
    end
    
    CS = ebsd.CS;
    gB = grains.boundary(CS.mineral, CS.mineral);
    
    %% 4. Define Twin System & Initialize TwinNetwork2
    fprintf('  -> Initializing {10-12} tensile twin system and TwinNetwork2 engine...\n');
    tS_1012 = twinSystem.hexagonal_1012(CS);
    
    network = TwinNetwork2(grains, gB, tS_1012, ebsd, 'thresh', 5*degree);
    
    %% 5. Execute Twin Network Topology Pipeline
    disp('    Step 1: Building Initial Topology...');
    G = network.buildInitialGraph();
    network.plotGraph(G);
    figG = gcf;
    title(sprintf('[%s] Step 1: Initial Raw Topology', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figG, fullfile(outDir, [defStep.stepName, '_01_Raw_Topology.png']), 'nodatestring');
    close(figG);
    
    disp('    Step 2: Extracting Spanning Tree...');
    mstG = network.extractSpanningTree(G);
    network.plotGraph(mstG);
    figMST = gcf;
    title(sprintf('[%s] Step 2: Minimum Spanning Tree', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figMST, fullfile(outDir, [defStep.stepName, '_02_Minimum_Spanning_Tree.png']), 'nodatestring');
    close(figMST);
    
    disp('    Step 3: Resolving Node States and Variants...');
    mstG = network.resolveNodeStates(mstG);
    network.plotGraph(mstG);
    figRes = gcf;
    title(sprintf('[%s] Step 3: MST after FZ Propagation', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figRes, fullfile(outDir, [defStep.stepName, '_03_Resolved_Node_States.png']), 'nodatestring');
    close(figRes);
    
    % disp('    Step 4: Pruning and Rewiring Outliers...');
    % mstG = network.pruneAndRewireOutliers(mstG, 'debug', false);
    % network.plotGraph(mstG);
    % figPrune = gcf;
    % title(sprintf('[%s] Step 4: MST after Pruning and Rewiring', defStep.stepName), 'Interpreter', 'none');
    % exportScaledFigure(figPrune, fullfile(outDir, [defStep.stepName, '_04_Pruned_Rewired.png']), 'nodatestring');
    % close(figPrune);
    
    disp('    Step 5: Reducing to Quotient Graph...');
    Q = network.reduceToQuotient(mstG);
    network.plotGraph(Q);
    figQ = gcf;
    title(sprintf('[%s] Step 5: Final Quotient Graph', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figQ, fullfile(outDir, [defStep.stepName, '_05_Quotient_Graph.png']), 'nodatestring');
    close(figQ);
    
    disp('    Step 6: Plotting grains with twin generation assignment...');
    network.plot('twinsOnly', false);
    figGen = gcf;
    title(sprintf('[%s] Step 6: Twin Generations & Network Map', defStep.stepName), 'Interpreter', 'none');
    exportScaledFigure(figGen, fullfile(outDir, [defStep.stepName, '_06_Twin_Generations.png']), 'nodatestring');
    close(figGen);
    
    % Save summary MAT file for this step
    outMatPath = fullfile(outDir, [defStep.stepName, '_TwinNetwork2_Results.mat']);
    save(outMatPath, 'network', 'G', 'mstG', 'Q', 'tS_1012');
    fprintf('  -> Successfully completed & saved TwinNetwork2 evaluation for Step [%s] to:\n     %s\n', defStep.stepName, outMatPath);
    
    close all; % Cleanup any remaining open handles
end

fprintf('\n========================================================================\n');
fprintf('ALL DEFORMED TWIN NETWORK ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
