clear
close all
mtexdata twins silent
CS = ebsd('indexed').CS;

% reconstruct grains
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',5*degree);

grains = smooth(grains,1);

% plot the grains
plot(grains,grains.meanOrientation)

%% Define {10-12} Twin System
% K1 = {10-12}, eta1 = <-1011>
cs = grains(1).CS;
tS = twinSystem.hexagonal_1012_1012(cs);
tS_sym = tS.symmetrise;

%% Group Twin Related Grains
pairs = grains.neighbors;

% misorientation between neighbors
mori = inv(grains(pairs(:,1)).meanOrientation) .* grains(pairs(:,2)).meanOrientation;

% check for twin relationship
isTwin = min(angle_outer(mori, [tS_sym.parentTwinMisorientation],'noSym2'), [], 2) < 5*degree;

% group grains
G = graph(pairs(isTwin,1), pairs(isTwin,2), [], length(grains));
componentId = conncomp(G);

% plot the groups
figure
plot(grains, ind2color(componentId, viridis(max(componentId))))
title('Twin Related Grain Clusters')

%% Identify Parent Grains
% Weights for criteria: [Size, Variants, InteractionWork]
weights = [1,1, 0]; 

% Stress tensor for Interaction Work (Tension along Y)
sigma = stressTensor.uniaxial(yvector);
% sigma = stressTensor.uniaxial(yvector);

% Initialize results
parentIds = [];
parentStats = [];

% Loop over each cluster
for i = 1:max(componentId)
    % Extract grains in the current cluster
    clusterMask = componentId == i;
    if sum(clusterMask) == 0, continue; end
    
    % Indices of grains in the full 'grains' list
    gIndices = find(clusterMask);
    cGrains = grains(gIndices);
    
    if length(cGrains) == 1
        parentIds = [parentIds; cGrains.id];
        continue;
    end
    
    % Initialize scores
    n = length(cGrains);
    s_size = zeros(n, 1);
    s_vars = zeros(n, 1);
    s_iw = zeros(n, 1);
    
    % 1. Size Criterion
    % Sum pixels of grains with similar orientation (< 8 deg)
    d = angle_outer(cGrains.meanOrientation, cGrains.meanOrientation);
    isSimilar = d < 8*degree;
    s_size = sum(isSimilar .* cGrains.numPixel', 2);
    
    % Prepare for 2 & 3
    % Filter neighbors to only those within the cluster
    inCluster = all(ismember(pairs, gIndices), 2);
    cPairs = pairs(inCluster, :);
    
    for k = 1:n
        currentInd = gIndices(k);
        currentOri = cGrains(k).meanOrientation;
        
        % Find neighbors of current grain within cluster
        nbMask = cPairs(:,1) == currentInd | cPairs(:,2) == currentInd;
        if ~any(nbMask), continue; end
        
        % Get neighbor indices
        nbRows = cPairs(nbMask, :);
        nbInds = unique(nbRows(:));
        nbInds(nbInds == currentInd) = [];
        
        % Calculate misorientations: inv(current) * neighbor
        nbOris = grains(nbInds).meanOrientation;
        mori = inv(currentOri) .* nbOris;
        
        % Determine twin variants
        ptm = tS_sym.parentTwinMisorientation;
        [minAng, vIdx] = min(angle_outer(mori, ptm,'noSym2'), [], 2);
        
        % Filter for valid twins
        isTwinRel = minAng < 5*degree;
        validVIdx = vIdx(isTwinRel);
        
        if isempty(validVIdx), continue; end
        
        % 2. Variant Criterion
        s_vars(k) = length(unique(validVIdx));
        
        % 3. Interaction Work Criterion
        % Activate twin systems: currentOri * tS_sym(variant)
        activeTS = currentOri * tS_sym(validVIdx);
        
        % Calculate IW
        iw = activeTS.interactionWork(sigma);
        
        % Mean IWs
        s_iw(k) = mean(iw);
    end
    
    % Normalize scores [0, 1]
    norm_s = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);
    
    totalScore = weights(1) * norm_s(s_size) + ...
                 weights(2) * norm_s(s_vars) + ...
                 weights(3) * norm_s(s_iw);
                 
    [~, bestIdx] = max(totalScore);
    
    % Store stats for the best grain
    parentStats = [parentStats; cGrains(bestIdx).id, s_size(bestIdx), s_vars(bestIdx), s_iw(bestIdx)];

    % Include all grains with < 8 deg misorientation to the identified parent
    isParent = angle(cGrains.meanOrientation, cGrains(bestIdx).meanOrientation) < 8*degree;
    parentIds = [parentIds; cGrains(isParent).id];
end

% Highlight Parents
parentGrains = grains(ismember(grains.id, parentIds));
hold on
plot(parentGrains, 'faceColor', 'r','faceAlpha', 0.5, 'lineWidth', 2)

% Annotate criteria values
if ~isempty(parentStats)
    [~, loc] = ismember(parentStats(:,1), grains.id);
    txt = arrayfun(@(i) sprintf('Sz:%d\nV:%d\nIW:%.1f', ...
        parentStats(i,2), parentStats(i,3), 100*parentStats(i,4)), ...
        1:size(parentStats,1), 'UniformOutput', false);
    text(grains(loc), txt, 'color','w','FontSize', 12, 'FontWeight', 'bold', 'Halo', true);
end

%% Extract and Plot Cluster for Grain 42
% Find the component ID of grain 42
clusterId = componentId(42);

% Extract all grains in this cluster
clusterGrains = grains(componentId == clusterId);

% Identify the representative parent grain (for orientation reference)
isRepParent = ismember(clusterGrains.id, parentStats(:,1));
if ~any(isRepParent)
    % Fallback if no parent identified in stats
    [~, idx] = max(clusterGrains.area);
    repParent = clusterGrains(idx);
else
    repParent = clusterGrains(isRepParent);
end

% Identify all parent grains (including those with similar orientation)
allParents = clusterGrains(ismember(clusterGrains.id, parentIds));
if isempty(allParents), allParents = repParent; end

% Plot the cluster with identified parent and twins
figure
plot(allParents, 'faceColor', [0.7 0.7 0.7])
hold on

% Plot twins with variant colors
colors = lines(length(tS_sym));
for i = 1:length(clusterGrains)
    g = clusterGrains(i);
    if ismember(g.id, allParents.id), continue; end
    
    % Determine the twin variant
    mori = inv(repParent.meanOrientation) * g.meanOrientation;
    [~, vIdx] = min(angle_outer(mori, tS_sym.parentTwinMisorientation,'noSym2'));
    hold on
    plot(g, 'faceColor', colors(vIdx,:))
    text(g, sprintf('T_{%d}', vIdx), 'HorizontalAlignment', 'center', 'FontSize', 8)
end

% plot(grains(42).boundary, 'lineWidth', 3, 'lineColor', 'r')
title(['Cluster #' num2str(clusterId) ' containing Grain 42'])

%% Calculate Net Deformation

% Initialize net displacement gradient tensor
netH = tensor(zeros(3,3)); 
totalArea = sum(clusterGrains.area);

% Sum up deformation from all twins
for i = 1:length(clusterGrains)
    g = clusterGrains(i);
    
    % Skip parent grains (no deformation relative to itself)
    if ismember(g.id, allParents.id), continue; end
    
    % Determine the twin variant
    mori = inv(repParent.meanOrientation) * g.meanOrientation;
    [~, vIdx] = min(angle_outer(mori, tS_sym.parentTwinMisorientation,'noSym2'));
    
    % Calculate displacement gradient for this variant in specimen coordinates
    grained_twin = repParent * tS_sym(vIdx);
    H = grained_twin.displacementGradient;
    
    % Add weighted contribution
    netH = netH + H * (g.area / totalArea);
end

disp('Net Displacement Gradient Tensor (Specimen Frame):');
disp(netH)

% Visualize the net deformation as a vector field
bnd = clusterGrains.boundary;
x_bnd = bnd.x; y_bnd = bnd.y;

% Create a grid over the cluster
xmin = min(x_bnd); xmax = max(x_bnd);
ymin = min(y_bnd); ymax = max(y_bnd);
step = min(xmax-xmin, ymax-ymin) / 15;
[xGrid, yGrid] = meshgrid(xmin:step:xmax, ymin:step:ymax);

% Filter points inside the cluster
k = boundary(x_bnd, y_bnd, 0.5);
in = inpolygon(xGrid, yGrid, x_bnd(k), y_bnd(k));
x = xGrid(in);
y = yGrid(in);

% Calculate displacement vectors relative to centroid
centroid = repParent.centroid;
r = [x - centroid.x, y - centroid.y, zeros(size(x))]'; % 3xN
dispVec = double(netH) * r; % 3xN

% Plot vector field
hold on
quiver(x, y, dispVec(1,:)', dispVec(2,:)', 'Color', 'k', 'LineWidth', 1.5)
% text(repParent, 'Net Deformation Field', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
