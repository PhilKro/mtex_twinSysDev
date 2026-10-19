%% initialize mtex
clear 
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
% dir = 'N:\EBSD-Analysis-2\Tijmen\Re tensile 1\h5oina\';
% file = 'Re tensile 1 Specimen 1 Site 4 Map Data 4_Map_Data_4_BW203.ang';
dir = 'N:\EBSD-Analysis-2\Tijmen\Re HSR 1';
file = "Re HSR twin right region 2 Site 1 Map Data 7_Map_Data_7_BW123.ang";
% dir = 'onedriveData/Kamila';
% file = 'pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';

path = fullfile(dir, file);

% Create analysis directory
analysis_dir = fullfile(dir, 'analysis');
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end


% define tilt, only for TKD
tkd_tilt = 70 * degree;

% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.10;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Re
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

ebsd = rotate(ebsd,rot);

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
ebsd0=ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;

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


% ipf with ci greyscale
figure;
plot(ebsd,cibw,'micronbar','off')
hold on
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
title('IPF Y')

%% grain reconstruction
% min angle used for grain reconstruction
grainThresh = 8*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 9;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 2;

%%% perform reconstruction %%%
ebsd = ebsd0;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

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
% plot(grains(~isParentGrain),0.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% plot(grains(isParentGrain),0.2*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% % plot boundaries
plot(grains.boundary)

figure;
plot(grains,colors)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)
%% fill ebsd data
ebsd_filled = fill(ebsd, grains);
figure; plot(ebsd,ebsd.orientations); nextAxis; plot(ebsd_filled, ebsd_filled.orientations); 

%% close all

% Save existing figures before closing
figHandles = findobj('Type', 'figure');
for i = 1:length(figHandles)
    saveas(figHandles(i), fullfile(analysis_dir, sprintf('Figure_part1_%d.jpg', figHandles(i).Number)));
end
close all
figure;
plot(grains,colors)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)

%% Graph-Based Twin Analysis (11-21 Twins)

% --- Settings ---
% maxGenerationLimit = 2; % Enforces how deep the twinning can go (not in TwinNetwork just yet)

% 1. Define the 11-21 twin system
tS_1121 = twinSystem.hexagonal_1121(CS);
% tS_1121 = twinSystem.hexagonal_1122(CS);
% 2. Initialize the TwinNetwork
% This class neatly encapsulates the graph-theoretical topology extraction,
% multi-region decomposition, and exact reference frame propagation.
network = twinNetwork(grains, gB, tS_1121, ebsd, ...
    'thresh', 8*degree);

% 3. Execute the topological pipeline step-by-step
disp('Step 1: Building Initial Topology...');
network.buildInitialGraph();
network.plotGraph(network.Graph);
title('Step 1: Initial Raw Topology');

disp('Step 2: Extracting Spanning Tree...');
network.extractSpanningTree();
network.plotGraph(network.mstGraph);
title('Step 2: Minimum Spanning Tree');

disp('Step 3: Resolving Node States and Variants...');
network.resolveNodeStates();
network.plotGraph(network.mstGraph);
title('Step 3: MST after FZ Propagation');

% disp('Step 4: Pruning and Rewiring Outliers...');
% network.pruneAndRewireOutliers('debug', true);
% % network.pruneAndRewireOutliers('targetNode', 56, 'debug', true);
% network.plotGraph(network.mstGraph);
% title('Step 4: MST after Pruning and Rewiring');

disp('Step 5: Reducing to Quotient Graph...');
network.reduceToQuotient();
network.plotGraph(network.QuotientGraph, 'showBackground', false, 'layout', 'layered', 'scaleNodesByArea', true);
title('Step 5: Final Quotient Graph (Hierarchical, Area-scaled)');

disp('Step 6: Plotting grains with twin generation assignment...');
network.plot('twinsOnly', false);
title('Step 6: Twin generations');

disp('Step 7: Calculating and plotting Schmid Factors...');
% Calculate Schmid factors for uniaxial tension along the x-direction
% sigma = stressTensor.uniaxial(xvector);
sigma = stressTensor.uniaxial(yvector);
[sf, iw] = network.analyzeMechanics(sigma);

figure;
plot(network.grains, sf, 'micronbar', 'off');
colormap jet;
clim([0 0.5]);
mtexColorbar;
title('Twin Schmid Factor Map (Colorbar)');

figure;
plot(network.grains, sf, 'micronbar', 'off');
colormap jet;
clim([0 0.5]);
title('Twin Schmid Factor Map (Labels)');
hold on;
valid_idx = find(~isnan(sf));
if ~isempty(valid_idx)
    sf_valid = sf(valid_idx);
    x_cents = network.grains.centroid.x(valid_idx);
    y_cents = network.grains.centroid.y(valid_idx);
    labels = cellstr(num2str(sf_valid, '%.2f'));
    text(x_cents, y_cents, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold');
end
mtexColorbar
hold off;

% network.plotSigmaSections();


% Save remaining figures
figHandles = findobj('Type', 'figure');
for i = 1:length(figHandles)
    % Try to get figure title for meaningful name, fallback to number
    ax = get(figHandles(i), 'CurrentAxes');
    figName = sprintf('Figure_part2_%d', figHandles(i).Number);
    if ~isempty(ax)
        t = get(ax, 'Title');
        if ~isempty(t) && ~isempty(t.String)
            cleanTitle = strrep(t.String, ' ', '_');
            cleanTitle = strrep(cleanTitle, ':', '');
            figName = cleanTitle;
        end
    end
    saveas(figHandles(i), fullfile(analysis_dir, [figName '.jpg']));
    %saveas(figHandles(i), fullfile(analysis_dir, [figName '.fig'])); % Optional: save as fig too
end
