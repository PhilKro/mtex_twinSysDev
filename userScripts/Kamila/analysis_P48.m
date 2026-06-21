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
    twinMatch_valid(:, v) = angle(O1 * twinMisoris(v), O2);
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
% 2. Identify Uncolored Twin Boundaries (Any Variant)
validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
gid1 = gB.grainId(validEdges, 1);
gid2 = gB.grainId(validEdges, 2);

ind1 = id2ind(grains, gid1);
ind2 = id2ind(grains, gid2);

O1 = grains(ind1).meanOrientation;
O2 = grains(ind2).meanOrientation;

% Calculate fit against any variant. MTEX angle naturally handles symmetry.
twinMatch_valid = zeros(length(gid1), numVariants);
for v = 1:numVariants
    % Dropping 'noSym2' allows true physical distance calculation across FZ boundaries
    twinMatch_valid(:, v) = angle(O1 * twinMisoris(v), O2);
end

% A boundary is a twin boundary if ANY variant matches within 5 degrees
twinThresh = 5 * degree;
[minAngle_valid, ~] = min(twinMatch_valid, [], 2);
isTwinBoundary_valid = minAngle_valid < twinThresh;

% Map back to original gB array size
isTwinBoundary = false(length(gB), 1);
isTwinBoundary(validEdges) = isTwinBoundary_valid;

% 3. Filtered Bipartite Projection (Topological Hierarchy)
% Extract unique grain pairs connected by a twin boundary
twinPairs = gB.grainId(isTwinBoundary, :);
[uPairs, ~, ~] = unique(sort(twinPairs, 2), 'rows');

activeGrains = unique(uPairs(:));
equivEdges = [];

% Trace paths of length 2: grains sharing a common twin neighbor
for v_idx = 1:length(activeGrains)
    w = activeGrains(v_idx);
    
    % Find all grains connected to w via a twin boundary
    edgesWithW = uPairs(:,1) == w | uPairs(:,2) == w;
    neighborsOfW = uPairs(edgesWithW, :);
    neighborsOfW = neighborsOfW(neighborsOfW ~= w); % Extract just the neighbor IDs
    
    % If w connects to multiple grains, check if they are the same orientation
    if length(neighborsOfW) > 1
        % Evaluate every pair of neighbors
        pairs = nchoosek(neighborsOfW, 2);
        for p = 1:size(pairs, 1)
            u = pairs(p, 1);
            v = pairs(p, 2);
            
            O_u = grains(id2ind(grains, u)).meanOrientation;
            O_v = grains(id2ind(grains, v)).meanOrientation;
            
            % If the physical misorientation is < 15 degrees, they belong to the same Super Node!
            if angle(O_u, O_v) < 15 * degree
                equivEdges = [equivEdges; u, v];
            end
        end
    end
end

% Include self loops to ensure all active grains exist in the equivalence graph
if isempty(equivEdges)
    equivEdges = [activeGrains, activeGrains];
else
    equivEdges = [equivEdges; [activeGrains, activeGrains]];
end

maxGid = max(grains.id);
equivEdges = [equivEdges; [maxGid, maxGid]];

% Construct equivalence graph and find Super Nodes
G_equiv = graph(equivEdges(:,1), equivEdges(:,2));
bins = conncomp(G_equiv); % Array where bins(grainId) gives the Super Node ID

% 4. Crown the Parent (Generation 0)
superNodes = unique(bins(activeGrains));
superNodeArea = zeros(max(superNodes), 1);

for i = 1:length(activeGrains)
    gid = activeGrains(i);
    sn = bins(gid);
    ind = id2ind(grains, gid);
    superNodeArea(sn) = superNodeArea(sn) + grains(ind).grainSize;
end

[~, parentSuperNode] = max(superNodeArea);

% 5. Build the Quotient Graph and compute Generations
qPairs = [bins(uPairs(:,1))', bins(uPairs(:,2))'];
% Remove self loops (edges between grains in the same super node)
isSelfLoop = qPairs(:,1) == qPairs(:,2);
qPairs(isSelfLoop, :) = [];

if isempty(qPairs)
    Q_G = graph();
else
    Q_G = graph(qPairs(:,1), qPairs(:,2));
end

if numnodes(Q_G) < max(bins)
    Q_G = addnode(Q_G, max(bins) - numnodes(Q_G));
end

% BFS distances from the Parent Super Node
dist_from_parent = distances(Q_G, parentSuperNode);

% Map generation back to all grains
grainGen = inf(maxGid, 1);
for i = 1:length(activeGrains)
    gid = activeGrains(i);
    sn = bins(gid);
    grainGen(gid) = dist_from_parent(sn);
end

genProp = inf(length(grains), 1);
activeInds = id2ind(grains, activeGrains);
genProp(activeInds) = grainGen(activeGrains);

% 6. Reference Frame Propagation (FZ-Agnostic Variant Assignment)
% Establish a single, global reference orientation from the Parent Matrix
parentGrainsInds = find(genProp == 0);
globalParentOri = mean(grains(parentGrainsInds).meanOrientation, 'robust');

spawnVariant = zeros(length(grains), 1);

% Compute variants for each generation sequentially
maxGenPlot = max(genProp(genProp < inf));

% Dictionary to store the ideal orientation of each grain for exact propagation
idealOri = repmat(globalParentOri, maxGid, 1);

for gen = 1:maxGenPlot
    genGids = grains(genProp == gen).id;
    
    for i = 1:length(genGids)
        child_gid = genGids(i);
        child_ind = id2ind(grains, child_gid);
        O_child = grains(child_ind).meanOrientation;
        
        % Find a parent grain of this child (generation - 1)
        edgesWithChild = find(uPairs(:,1) == child_gid | uPairs(:,2) == child_gid);
        for e = 1:length(edgesWithChild)
            idx = edgesWithChild(e);
            parent_gid = uPairs(idx, 1);
            if parent_gid == child_gid, parent_gid = uPairs(idx, 2); end
            
            if grainGen(parent_gid) == gen - 1
                % Evaluate variant using the parent's ideal theoretical orientation
                % This enforces a consistent FZ chain across the whole hierarchy!
                O_ideal_parent = idealOri(parent_gid);
                
                fit = zeros(1, numVariants);
                for v = 1:numVariants
                    fit(v) = angle(O_ideal_parent * twinMisoris(v), O_child);
                end
                
                [~, bestVar] = min(fit);
                spawnVariant(child_ind) = bestVar;
                
                % Store this child's ideal orientation for the next generation
                idealOri(child_gid) = O_ideal_parent * twinMisoris(bestVar);
                break; % Found the spawning variant
            end
        end
    end
end

% 7. Plot Generations
figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.2); 
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

if ~isempty(maxGenPlot)
    cmapGen = lines(maxGenPlot + 1);
    for g = 0:maxGenPlot
        gInds = find(genProp == g);
        if ~isempty(gInds)
            if g == 0
                hold on
                plot(grains(gInds), 'faceColor',[0.3 0.3 0.3], 'DisplayName', 'Parent (Gen 0)');
            else
                hold on
                plot(grains(gInds), 'faceColor', cmapGen(g+1, :), 'DisplayName', sprintf('Generation %d', g));
            end
        end
    end
end
title('Twin Generations (Filtered Bipartite Projection)')
legend('show', 'Location', 'northeastoutside');
hold off;

% 8. Plot Spawning Variants
figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.2);
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

cmapVar = lines(numVariants);
for v = 1:numVariants
    vInds = find(spawnVariant == v);
    if ~isempty(vInds)
        hold on
        plot(grains(vInds), 'faceColor',cmapVar(v, :), 'DisplayName', sprintf('Spawning Variant %d', v));
    end
end

if ~isempty(parentGrainsInds)
    hold on
    plot(grains(parentGrainsInds), 'faceColor',[0.3 0.3 0.3], 'DisplayName', 'Parent Matrix');
end

title('Globally Consistent Spawning Variants')
legend('show', 'Location', 'northeastoutside');
hold off;
