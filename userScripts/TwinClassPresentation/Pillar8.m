clear
close all

% Set MTEX preferences for axis directions
plotx2east
plotzIntoPlane

% Define file name and crystal symmetry for Titanium
fname = 'userScripts/TwinClassPresentation/2024-04-30_Kamila TKD Ti Pillar8 Site 1 Map Data 31_BW158 1.ang';
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

% Load EBSD data
ebsd = EBSD.load(fname, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

% Clean data based on a confidence index (CI) threshold
ebsd(ebsd.ci < 0.15) = 'notIndexed';

tkd_tilt = -20 * degree;
% Perform corrections for tkd tilt and misaligned axes
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

% Gridify the data for faster processing and visualization
ebsd = ebsd.gridify;

% reconstruct grains
[grains, ebsd.grainId] = calcGrains(ebsd,'angle',5*degree, 'minPixel', 5);

grains = smooth(grains,1);

% plot the grains
plot(grains,grains.meanOrientation)

%% Define {11-22} Twin System
% K1 = {11-22}, eta1 = <11-2-3>
cs = grains.CS;
tS = twinSystem.hexagonal_1122_1124(cs);
tS_sym = tS.symmetrise;

%% Group Twin Related Grains
pairs = grains.neighbors;

% misorientation between neighbors
mori = inv(grains(pairs(:,1)).meanOrientation) .* grains(pairs(:,2)).meanOrientation;

% check for twin relationship
isTwin = min(angle_outer(mori, [tS_sym.parentTwinMisorientation]), [], 2) < 5*degree;

% group grains
G = graph(pairs(isTwin,1), pairs(isTwin,2), [], length(grains));
componentId = conncomp(G);

% plot the groups
figure
plot(grains, ind2color(componentId, viridis(max(componentId))))
title('Twin Related Grain Clusters')

%% Identify Parent Grains
% Weights for criteria: [Size, Variants, InteractionWork]
weights = [0, 0, 1]; 

% Stress tensor for Interaction Work (Tension along Y)
% sigma = stressTensor.uniaxial(yvector);
sigma = - stressTensor.uniaxial(yvector);

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
        [minAng, vIdx] = min(angle_outer(mori, ptm), [], 2);
        
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

%% Twin Analysis and Visualization
disp('Starting twin analysis and visualization...');

% Initialize classification arrays
isParent = ismember(grains.id, parentIds);
isPrimary = false(length(grains), 1);
% Store [parent_cluster_id, variant_id]
primaryVariantInfo = zeros(length(grains), 2); 
tolerance = 8 * degree;

% --- Identify Primary Twins ---
disp('Identifying primary twins...');

% Identify clusters that have parent statistics
[~, loc] = ismember(parentStats(:,1), grains.id);
relevantClusters = unique(componentId(loc));

for i = relevantClusters(:).'
    
    % Find the representative parent grain ID for this cluster from parentStats
    clusterGrainIDs = grains(componentId == i).id;
    is_in_cluster = ismember(parentStats(:,1), clusterGrainIDs);
    
    rep_parent_id = parentStats(is_in_cluster, 1);
    parentOri = grains(rep_parent_id).meanOrientation;
    
    % Generate theoretical primary twin orientations
    primary_twins = parentOri * tS_sym;
    primary_twin_oris = primary_twins.orientation;
    
    % Get indices of non-parent grains in the current cluster
    clusterGrainInds = find(componentId == i & ~isParent.');
    
    % Iterate through grains in the cluster to find primary twins
    for k = 1:length(clusterGrainInds)
        gInd = clusterGrainInds(k);
        
        % Find the best matching primary twin variant
        [minAngle, bestFit] = min(angle(grains(gInd).meanOrientation, primary_twin_oris));
        
        % If the match is within tolerance, classify as a primary twin
        if minAngle < tolerance
            isPrimary(gInd) = true;
            primaryVariantInfo(gInd, :) = [i, primary_twins(bestFit).variantId];
        end
    end
end

%% Visualization
disp('Visualizing the twin analysis results...');

% Create a new figure for the final plot
figure;
hold on

% Plot unidentified grains in bright red
unidentifiedGrains = grains(~isParent & ~isPrimary);
if ~isempty(unidentifiedGrains)
    plot(unidentifiedGrains, 'facecolor', 'r', 'DisplayName', 'Unidentified');
end

% Plot parent grains in grey and label them
parentGrains = grains(isParent);
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
parentClusters = unique(componentId(isParent));
for i = 1:length(parentClusters)
    pID = parentClusters(i);
    if pID == 0, continue; end
    grainsOfParentCluster = grains(isParent' & componentId == pID);
    if ~isempty(grainsOfParentCluster)
        text(grainsOfParentCluster, ['P_{' num2str(pID) '}'], 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    end
end

% Plot primary twins with variant-specific colors and labels
colors = lines(length(tS_sym));
primaryInfo = primaryVariantInfo(isPrimary,:);
uniquePV = unique(primaryInfo, 'rows');

for k = 1:size(uniquePV, 1)
    pID = uniquePV(k, 1);
    vID = uniquePV(k, 2);
    if pID == 0 || vID == 0, continue; end
    
    grainsOfPV = grains(isPrimary' & primaryVariantInfo(:,1) == pID & primaryVariantInfo(:,2) == vID);
    
    if ~isempty(grainsOfPV)
        hold on
        plot(grainsOfPV, 'facecolor', colors(vID, :));
        label = sprintf('P_{%d}T_{%d}', pID, vID);
        text(grainsOfPV, label, 'HorizontalAlignment', 'center');
    end
end

hold off
title('Primary Twin Analysis of Pillar8');
legend off
