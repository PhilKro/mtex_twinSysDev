%% initialize mtex
clear 
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
% dir = 'userScripts/Kamila';
% file = 'Sample48 Site 3 Map Data 6_Map_Data_6_CircMask_EntireMap_NoGB_BW158.ang';
dir = 'onedriveData/Kamila';
file = 'pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';

path = fullfile(dir, file);

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
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
indebsd = ebsd('indexed');
refebsd = indebsd(1);
reforientation = refebsd.orientations;
clear indebsd
reforientation = orientation.byEuler(101.3*degree,10.3*degree,260.7*degree);
% reforientation = orientation.byEuler(270.2*degree,79.7*degree,92.3*degree);

ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
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
minSize = 40;
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

%% Graph-Based Twin Analysis (11-21 Twins)

% 1. Define the 11-21 twin system
tS_1121 = twinSystem.hexagonal_1121_0001(CS);
tS_sym = tS_1121.symmetrise; % Get all variants
twinMisoris = tS_sym.parentTwinMisorientation;
numVariants = length(twinMisoris);

% 2. Evaluate variants using absolute mean orientations to avoid FZ projection issues
% Extract all valid boundary segments
validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
gid1 = gB.grainId(validEdges, 1);
gid2 = gB.grainId(validEdges, 2);

ind1 = id2ind(grains, gid1);
ind2 = id2ind(grains, gid2);

O1 = grains(ind1).meanOrientation;
O2 = grains(ind2).meanOrientation;

% Vectorized angle calculation across all boundary pairs
twinMatch_valid = zeros(length(gid1), numVariants);
for v = 1:numVariants
    twinMatch_valid(:, v) = angle(O1 * twinMisoris(v), O2, 'noSym2');
end

% Map the results back to the original gB array size
twinMatch = inf(length(gB), numVariants);
twinMatch(validEdges, :) = twinMatch_valid;

[minAngle, variantId] = min(twinMatch, [], 2);

% 3. Identify twin boundaries based on a threshold angle
twinThresh = 5 * degree;
isTwinBoundary = minAngle < twinThresh;

% 4. Plot the grain map with different twin variants
figure;
plot(grains, grains.meanOrientation);
hold on
% Plot normal grain boundaries
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

% Define a colormap for the variants
cmap = lines(numVariants);

% Plot each variant with a different color
for v = 1:numVariants
    isCurrentVariant = isTwinBoundary & (variantId == v);
    if any(isCurrentVariant)
        plot(gB(isCurrentVariant), 'lineColor', cmap(v,:), 'lineWidth', 3, 'DisplayName', sprintf('Variant %d', v));
    end
end
title('11-21 Twin Variants')
hold off

%% 5. Build the Quotient Graph
% Extract the grain IDs for the boundaries that have a twinning relationship
twinPairs = gB.grainId(isTwinBoundary, :);
twinPairsVariants = variantId(isTwinBoundary);

% Filter out boundaries to unindexed regions (grainId == 0)
validEdges = twinPairs(:,1) > 0 & twinPairs(:,2) > 0;
twinPairs = twinPairs(validEdges, :);
twinPairsVariants = twinPairsVariants(validEdges);

% Find unique grain pairs that share a twin boundary
[uPairs, ~, ic] = unique(sort(twinPairs, 2), 'rows');
uVariants = zeros(size(uPairs, 1), 1);
for i = 1:size(uPairs, 1)
    uVariants(i) = mode(twinPairsVariants(ic == i));
end

% A. Trace Edge Variants to build Super Nodes (Orientational Regions)
activeGrains = unique(uPairs(:));
equivEdges = [];

for v_idx = 1:length(activeGrains)
    v = activeGrains(v_idx);
    
    % Find edges connecting to v
    edgesWithV = find(uPairs(:,1) == v | uPairs(:,2) == v);
    neighborsOfV = zeros(length(edgesWithV), 1);
    variantsOfV = uVariants(edgesWithV);
    
    for i = 1:length(edgesWithV)
        idx = edgesWithV(i);
        if uPairs(idx, 1) == v
            neighborsOfV(i) = uPairs(idx, 2);
        else
            neighborsOfV(i) = uPairs(idx, 1);
        end
    end
    
    % If two neighbors connect to v via the SAME variant, they are orientational equivalent
    for k = 1:numVariants
        neighborsK = neighborsOfV(variantsOfV == k);
        if length(neighborsK) > 1
            % Connect these equivalent neighbors
            for i = 1:length(neighborsK)-1
                equivEdges = [equivEdges; neighborsK(i), neighborsK(i+1)];
            end
        end
    end
end

% Include self loops to ensure all active grains are in the equivalence graph
if isempty(equivEdges)
    equivEdges = [activeGrains, activeGrains];
else
    equivEdges = [equivEdges; [activeGrains, activeGrains]];
end

% Ensure the graph can index up to the max grain ID
maxGid = max(grains.id);
equivEdges = [equivEdges; [maxGid, maxGid]];

% Construct equivalence graph and find connected components (Super Nodes)
G_equiv = graph(equivEdges(:,1), equivEdges(:,2));
bins = conncomp(G_equiv); % Array where bins(grainId) gives the Super Node ID

% B. Crown the Parent (Generation 0)
superNodes = unique(bins(activeGrains));
superNodeArea = zeros(max(superNodes), 1);

for i = 1:length(activeGrains)
    gid = activeGrains(i);
    sn = bins(gid);
    ind = id2ind(grains, gid);
    superNodeArea(sn) = superNodeArea(sn) + grains(ind).grainSize;
end

% The Super Node with maximum total area is Generation 0
[~, parentSuperNode] = max(superNodeArea);

% C. Build the Quotient Graph and compute Generations
qPairs = [bins(uPairs(:,1))', bins(uPairs(:,2))'];
% Remove self loops in Quotient Graph
isSelfLoop = qPairs(:,1) == qPairs(:,2);
qPairs(isSelfLoop, :) = [];

if isempty(qPairs)
    Q_G = graph();
else
    Q_G = graph(qPairs(:,1), qPairs(:,2));
end

% Ensure graph is large enough
if numnodes(Q_G) < max(bins)
    Q_G = addnode(Q_G, max(bins) - numnodes(Q_G));
end

% Calculate distances from Parent
dist_from_parent = distances(Q_G, parentSuperNode);

% Map generation back to active grains
grainGen = inf(maxGid, 1);
for i = 1:length(activeGrains)
    gid = activeGrains(i);
    sn = bins(gid);
    grainGen(gid) = dist_from_parent(sn);
end

% Assign generation as a property to ALL grains for plotting
genProp = inf(length(grains), 1);
activeInds = id2ind(grains, activeGrains);
genProp(activeInds) = grainGen(activeGrains);

% D. Plot Generations
figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.2); 
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

maxGenPlot = max(genProp(genProp < inf));
if ~isempty(maxGenPlot)
    cmapGen = lines(maxGenPlot + 1);
    for g = 0:maxGenPlot
        gInds = find(genProp == g);
        if ~isempty(gInds)
            if g == 0
                plot(grains(gInds), [0.3 0.3 0.3], 'DisplayName', 'Parent (Gen 0)');
            else
                plot(grains(gInds), cmapGen(g+1, :), 'DisplayName', sprintf('Generation %d', g));
            end
        end
    end
end
title('Twin Generations (Quotient Graph)')
legend('show', 'Location', 'northeastoutside');
hold off;

% E. Find Spawning Variants and Plot
spawnVariant = zeros(length(grains), 1);
for i = 1:length(activeGrains)
    gid = activeGrains(i);
    gInd = id2ind(grains, gid);
    gGen = grainGen(gid);
    
    if gGen > 0 && gGen < inf
        % Find neighbor with Gen == gGen - 1
        edgesWithV = find(uPairs(:,1) == gid | uPairs(:,2) == gid);
        for j = 1:length(edgesWithV)
            idx = edgesWithV(j);
            u = uPairs(idx, 1);
            if u == gid, u = uPairs(idx, 2); end
            
            if grainGen(u) == gGen - 1
                spawnVariant(gInd) = uVariants(idx);
                break;
            end
        end
    end
end

figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.2);
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

for v = 1:numVariants
    vInds = find(spawnVariant == v);
    if ~isempty(vInds)
        plot(grains(vInds), cmap(v, :), 'DisplayName', sprintf('Spawning Variant %d', v));
    end
end

parentInds = find(genProp == 0);
if ~isempty(parentInds)
    plot(grains(parentInds), [0.3 0.3 0.3], 'DisplayName', 'Parent Matrix');
end

title('Spawning Twin Variants')
legend('show', 'Location', 'northeastoutside');
hold off;
