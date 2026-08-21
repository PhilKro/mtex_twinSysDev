%% initialize mtex
clear 
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
dir = 'userScripts/Kamila';
file = 'Sample48 Site 3 Map Data 6_Map_Data_6_CircMask_EntireMap_NoGB_BW158.ang';
% dir = 'onedriveData/Kamila';
% file = 'pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';

path = fullfile(dir, file);
[~, baseName, ~] = fileparts(file);
% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.10;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Re
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
ebsd.phase(ebsd.isnan) = -1;
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');
ebsd = rotate(ebsd,rot);

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
indebsd = ebsd('indexed');
refebsd = indebsd(1);
reforientation = refebsd.orientations;
clear indebsd
reforientation = orientation.byEuler(101.3*degree,10.3*degree,260.7*degree);
% reforientation = orientation.byEuler(270.2*degree,79.7*degree,92.3*degree);

ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
ebsd0=ebsd;


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
figure;
plot(ebsd,cibw,'micronbar','off')
hold on
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
title('IPF Y')
exportScaledFigure(gcf, fullfile(dir, [baseName, '_IPF_CI_overlay.jpg']), 'target_px', 4000, 'noDateString')

%% grain reconstruction
% min angle used for grain reconstruction
grainThresh = 8*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 10;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 2.2;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha,'boundary','tight');

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

% 2. Initialize the TwinNetwork
% This class neatly encapsulates the graph-theoretical topology extraction,
% multi-region decomposition, and exact reference frame propagation.
network = TwinNetwork(grains, gB, tS_1121, ebsd, 'thresh', 8*degree);

% 3. Execute the topological pipeline step-by-step
disp('Step 1: Building Initial Topology...');
G = network.buildInitialGraph();
network.plotGraph(G);
title('Step 1: Initial Raw Topology');

disp('Step 2: Extracting Spanning Tree...');
mstG = network.extractSpanningTree(G);
network.plotGraph(mstG);
title('Step 2: Minimum Spanning Tree');

disp('Step 3: Resolving Node States and Variants...');
mstG = network.resolveNodeStates(mstG);
network.plotGraph(mstG);
title('Step 3: MST after FZ Propagation');

disp('Step 3.5: Plotting grains with twin generation assignment before pruning...');
network.plot('twinsOnly', false);
title('Step 3.5: Twin generations');

% disp('Step 4: Pruning and Rewiring Outliers...');
% mstG = network.pruneAndRewireOutliers(mstG, 'debug', true);
% network.plotGraph(mstG);
% title('Step 4: Pruned Graph Topology');

disp('Step 4: Pruning and Rewiring Outliers...');
mstG = network.pruneAndRewireOutliers(mstG, 'outlierMode', 'theoretical', 'maxErrThreshold', 25*degree, 'debug', true, 'rejectGenIncrease', true);
network.plotGraph(mstG);
title('Step 4: Pruned Graph Topology');

disp('Step 5: Reducing to Quotient Graph...');
Q = network.reduceToQuotient(mstG);
network.plotGraph(Q);
title('Step 5: Final Quotient Graph');

disp('Step 6: Plotting grains with twin generation assignment...');
network.plot('twinsOnly', false);
title('Step 6: Twin generations');
