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

% Determine type 1 secondary twin indices
idx = nchoosek(1:length(tS_sym.rotAxis), 2);
valid = eq(tS_sym.rotAxis(idx(:,1)), tS_sym.rotAxis(idx(:,2)), 'noSymmetry', 'antipodal');
cozonalPairs = idx(valid, :);

% 2. Identify Twin Boundaries (Local Pixel Misorientation)
% Using MTEX's internal symmetry handling, we check the deviation 
% of the local boundary misorientation from the theoretical twin system.
twinThresh = 5 * degree;
isTwinBoundary = angle(gB.misorientation, tS_1121.parentTwinMisorientation) < twinThresh;

% Filter out boundaries to unindexed regions
validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
isTwinBoundary = isTwinBoundary & validEdges;

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

variantPath = cell(maxGid, 1);
parentGids = grains(parentGrainsInds).id;
for i = 1:length(parentGids)
    variantPath{parentGids(i)} = [];
end

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
                % Evaluate variant using the parent's actual orientation, BUT projected 
                % to its ideal reference frame to prevent FZ flipping!
                O_ideal_parent = idealOri(parent_gid);
                O_actual_parent = grains(id2ind(grains, parent_gid)).meanOrientation;
                O_parent_fixed = O_actual_parent.project2FundamentalRegion(O_ideal_parent);
                
                fit = zeros(1, numVariants);
                for v = 1:numVariants
                    fit(v) = angle(O_parent_fixed * twinMisoris(v), O_child);
                end
                
                [~, bestVar] = min(fit);
                spawnVariant(child_ind) = bestVar;
                
                % Store this child's ideal orientation for the next generation
                idealOri(child_gid) = O_parent_fixed * twinMisoris(bestVar);
                
                % Record the variant path
                parentPath = variantPath{parent_gid};
                variantPath{child_gid} = [parentPath, bestVar];
                break; % Found the spawning variant
            end
        end
    end
end

%% 7. Plot 1: Generations
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

% Plot 2: Twin Graph Overlay
figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.5);
hold on;
plot(grains.boundary, 'lineColor', 'k', 'lineWidth', 1);

% Draw edges between centroids for all twin boundaries
ind1 = id2ind(grains, uPairs(:,1));
ind2 = id2ind(grains, uPairs(:,2));

x1 = grains(ind1).centroid.x;
y1 = grains(ind1).centroid.y;
x2 = grains(ind2).centroid.x;
y2 = grains(ind2).centroid.y;

plot([x1'; x2'], [y1'; y2'], 'Color', 'r', 'LineWidth', 2);

% Plot nodes
scatter(grains(activeInds).centroid.x, grains(activeInds).centroid.y, 20, 'k', 'filled');
title('Twin Graph Structure Overlay');
hold off;

% Plot 3: Quotient Graph Topology
figure;
p = plot(Q_G, 'Layout', 'force', 'MarkerSize', 6);
highlight(p, parentSuperNode, 'NodeColor', 'r', 'MarkerSize', 10);
title('Quotient Graph Topology (Red = Parent Matrix)');

% Plot 4: Variant-Colored Boundaries
% Map boundary segments to variants
gB_twin = gB(isTwinBoundary);
bnd_gid1 = gB_twin.grainId(:,1);
bnd_gid2 = gB_twin.grainId(:,2);

bnd_variant = zeros(length(gB_twin), 1);
for i = 1:length(gB_twin)
    g1 = bnd_gid1(i);
    g2 = bnd_gid2(i);
    
    if grainGen(g1) > grainGen(g2)
        bnd_variant(i) = spawnVariant(id2ind(grains, g1));
    elseif grainGen(g2) > grainGen(g1)
        bnd_variant(i) = spawnVariant(id2ind(grains, g2));
    end
end

figure;
plot(grains, grains.meanOrientation, 'faceAlpha', 0.2);
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

cmapVar = lines(numVariants);
for v = 1:numVariants
    idx = (bnd_variant == v);
    if any(idx)
        hold on
        plot(gB_twin(idx), 'lineColor', cmapVar(v,:), 'lineWidth', 3, 'DisplayName', sprintf('Variant %d', v));
    end
end
title('Twin Boundaries Colored by Spawning Variant')
legend('show', 'Location', 'northeastoutside');
hold off;

% 11. Plot 4: Grains Filled by Twin Variant (Type 1 Cozonal Analysis)
figure;
plot(grains, 'faceColor', 'w');
hold on;
plot(grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);

cmapVar = lines(numVariants);
labels = cell(length(grains), 1);

parentIndsPlot = [];
primaryIndsPlot = cell(numVariants, 1);
type1IndsPlot = cell(numVariants, 1);
type2IndsPlot = [];

for i = 1:length(activeGrains)
    gid = activeGrains(i);
    gInd = id2ind(grains, gid);
    gen = grainGen(gid);
    path = variantPath{gid};
    
    if gen == 0
        parentIndsPlot = [parentIndsPlot; gInd];
        labels{gInd} = 'P';
    elseif gen == 1
        v1 = path(1);
        primaryIndsPlot{v1} = [primaryIndsPlot{v1}; gInd];
        labels{gInd} = sprintf('T_{%d}', v1);
    elseif gen == 2
        v1 = path(1);
        v2 = path(2);
        
        isCozonal = (v1 == v2);
        for k = 1:size(cozonalPairs, 1)
            if (v1 == cozonalPairs(k,1) && v2 == cozonalPairs(k,2)) || ...
               (v1 == cozonalPairs(k,2) && v2 == cozonalPairs(k,1))
                isCozonal = true;
                break;
            end
        end
        
        if isCozonal
            type1IndsPlot{v1} = [type1IndsPlot{v1}; gInd];
        else
            type2IndsPlot = [type2IndsPlot; gInd];
        end
        labels{gInd} = sprintf('T_{%d,%d}', v1, v2);
    elseif gen > 2
        type2IndsPlot = [type2IndsPlot; gInd];
        str = sprintf('%d,', path);
        labels{gInd} = sprintf('T_{%s}', str(1:end-1));
    end
end

if ~isempty(parentIndsPlot)
    hold on
    plot(grains(parentIndsPlot), 'faceColor', [0.8 0.8 0.8]);
end

for v = 1:numVariants
    if ~isempty(primaryIndsPlot{v})
        hold on
        plot(grains(primaryIndsPlot{v}), 'faceColor', cmapVar(v, :));
    end
end

for v = 1:numVariants
    if ~isempty(type1IndsPlot{v})
        hold on
        plot(grains(type1IndsPlot{v}), 'faceColor', cmapVar(v, :) * 0.5);
    end
end

if ~isempty(type2IndsPlot)
    hold on
    plot(grains(type2IndsPlot), 'faceColor', [1 0 1]);
end

activeInds = id2ind(grains, activeGrains);
text(grains(activeInds), labels(activeInds), 'FontSize', 8, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

title('Grains by Twin Variant (Type 1 Cozonal = Dark, Type 2 = Pink)');
hold off;
