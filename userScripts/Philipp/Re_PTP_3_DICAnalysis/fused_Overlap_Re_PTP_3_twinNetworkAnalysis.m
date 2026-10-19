%% initialize mtex
clear
clc
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
dir = 'C:\Users\phkr\MTEX_workEnv\userScripts\Philipp\Re_PTP_3_DICAnalysis\';
file = 'Re 4a 25pct 158bw_fused_Re.ang';

path = fullfile(dir, file);
[~, baseName, ~] = fileparts(file);

% Define output directory for figures
outDir = fullfile(dir, 'output');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is usually safe (changed to -1 because ang is wierd)
thresh_ci = -1;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Re
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2'); %
% ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
% Read the raw .ang file directly to find points where phi1 is exactly 0.00000
fid = fopen(path, 'r');
raw_data = textscan(fid, '%f %*[^\n]', 'HeaderLines', 38);
fclose(fid);
raw_phi1 = raw_data{1};

if length(raw_phi1) == length(ebsd)
    ebsd(raw_phi1 == 0.0) = 'notIndexed';
else
    warning('Raw phi1 length does not match ebsd length. Skipping direct 0.00000 filter.');
end

ebsd.phase(ebsd.isnan) = -1;
% rot = rotation.byAxisAngle(zvector,180*degree);
% ebsd = rotate(ebsd,rot,'keepXY');
% ebsd = rotate(ebsd,rot);

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
% indebsd = ebsd('indexed');
% if length(indebsd) > 0
%     refebsd = indebsd(1);
%     reforientation = refebsd.orientations;
% else
%     disp('Warning: No indexed points found! Using a default reference orientation.');
%     % Create a rotation first to avoid MTEX 7 parsing bugs in orientation.byEuler
% rot = rotation.byEuler(101.3*degree,10.3*degree,260.7*degree);
%% initialize mtex
clear
clc
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
dir = 'C:\Users\phkr\MTEX_workEnv\userScripts\Philipp\Re_PTP_3_DICAnalysis\';
file = 'Re 4a 25pct 158bw_fused_Re.ang';

path = fullfile(dir, file);
[~, baseName, ~] = fileparts(file);

% Define output directory for figures
outDir = fullfile(dir, 'output');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is usually safe (changed to -1 because ang is wierd)
thresh_ci = -1;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Re
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2'); %
% ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
% Read the raw .ang file directly to find points where phi1 is exactly 0.00000
fid = fopen(path, 'r');
raw_data = textscan(fid, '%f %*[^\n]', 'HeaderLines', 38);
fclose(fid);
raw_phi1 = raw_data{1};

if length(raw_phi1) == length(ebsd)
    ebsd(raw_phi1 == 0.0) = 'notIndexed';
else
    warning('Raw phi1 length does not match ebsd length. Skipping direct 0.00000 filter.');
end

ebsd.phase(ebsd.isnan) = -1;

reforientation = orientation.byEuler(4.14*degree, 100.01*degree, 179.92*degree, CS);
ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
ebsd0 = ebsd;

% reforientation = orientation.byEuler(184.137*degree,79.9861*degree,180.076*degree, CS);
% original_baseName = baseName;
% sym_refs = reforientation.symmetrise;
% for iter_var = 1:length(sym_refs)
%     current_ref = sym_refs(iter_var);
%     [phi1, Phi, phi2] = Euler(current_ref);
%     fprintf('\n=======================================================\n');
%     fprintf('Running Variant %02d / %02d\n', iter_var, length(sym_refs));
%     fprintf('Reference Orientation Euler: %.2f, %.2f, %.2f (degrees)\n', phi1/degree, Phi/degree, phi2/degree);
%     fprintf('=======================================================\n\n');
% baseName = sprintf('%s_var%02d', original_baseName, iter_var);

% ebsd = ebsd0;
% ebsd.orientations = ebsd.orientations.project2FundamentalRegion(current_ref);


%% Crop EBSD Data
% Define your rectangular region: [xmin, ymin, width, height]
% uncomment the lines below to apply the crop
% region = [-6.5, -5.7, 0.8, 0.8]; 
% region = [-4.9, -6.2, 0.6, 0.8]; 
% region = [-4.0, -3.2, 2.7, 2.7]; %for pillar 48
% ebsd = ebsd(inpolygon(ebsd, region));
% ebsd0 = ebsd;


%% some initial test plot %%%

ipf=ipfColorKey(ebsd.CS);

figure;
ipf.inversePoleFigureDirection=xvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors)
title('IPF X')

ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
nextAxis
plot(ebsd('indexed'),colors)
title('IPF Y')

ipf.inversePoleFigureDirection = zvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
nextAxis
plot(ebsd('indexed'),colors)
title('IPF Z')

nextAxis
cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
cibw = cibw.xyz;
ebsd.prop.ciColor = cibw;
plot(ebsd,cibw)
colormap gray
exportScaledFigure(gcf, fullfile(dir, [baseName, '_EBSD_IPF.jpg']), 'target_px', 4000, 'noDateString')

% ipf with ci greyscale
% figure;
% plot(ebsd,cibw,'micronbar','off')
% hold on
% ipf.inversePoleFigureDirection = yvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
% title('IPF Y')
% exportScaledFigure(gcf, fullfile(dir, [baseName, '_IPF_CI_overlay.jpg']), 'target_px', 4000, 'noDateString')

%% grain reconstruction
% min angle used for grain reconstruction
grainThresh = 8*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 10;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 2.2;

[grains,ebsd.grainId,~] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha,'boundary','tight');
grains = grains.smooth(5);

% take only boundaries between indexed grains
gB = grains.boundary(CS.mineral,CS.mineral);

%%% test plot %%%

figure;
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(grains.meanOrientation);
% plot(ebsd,cibw,'micronbar','off')
% hold on
plot(grains,colors)%,'FaceAlpha',0.7)

% create and plot crystal shapes
cS1 = crystalShape.hex(ebsd.CS);

hold on
%plottablecS = grains.meanOrientation*cS1*0.5
plot(grains,0.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
exportScaledFigure(gcf, fullfile(dir, [baseName, '_grain_shapes.jpg']), 'target_px', 4000, 'noDateString')
% plot(grains(~isParentGrain),0.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% plot(grains(isParentGrain),0.2*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% % plot boundaries
figure;
plot(grains.boundary)
text(grains,grains.id)
exportScaledFigure(gcf, fullfile(dir, [baseName, '_grain_ids.jpg']), 'target_px', 4000, 'noDateString')

%% Graph-Based Twin Analysis (11-21 Twins)

% 1. Define the 11-21 twin system
tS_1121 = twinSystem.hexagonal_1121(CS);

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

%end % End of loop over symmetrical reference orientations
