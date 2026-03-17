%%%%----
% THE FOLLOWING SCRIPT IS GARBAGE, THANK YOU FOR YOUR ATTENTION TO THIS
% MATTER
%%%%----

clear
close all

% Set MTEX preferences for axis directions
plotx2east
plotzIntoPlane

% Define file name and crystal symmetry for Titanium
% fname = 'onedriveData/Kamila/2024-04-30_Kamila TKD Ti Pillar8 Site 1 Map Data 31_BW158 1.ang';
% boundingbox = [0.75 0.45 6.1 12.8];

fname = 'onedriveData/Kamila/2024-04-30_Kamila TKD Ti Specimen 3 Site 2 Map Data 10_BW123.ang';
boundingbox = [0.3 0.45 6.3 9.8];

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
ebsd = ebsd(inpolygon(ebsd, boundingbox));


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
isTwin = min(angle_outer(mori, [tS_sym.parentTwinMisorientation]), [], 2) < 10*degree;

% group grains
G = graph(pairs(isTwin,1), pairs(isTwin,2), [], length(grains));
componentId = conncomp(G);

% plot the groups
figure
plot(grains, ind2color(componentId, viridis(max(componentId))))
title('Twin Related Grain Clusters')

%% Identify Parent Grains
% Weights for criteria: [Size, Variants, InteractionWork]
weights = [1, 0, 0]; 

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
isSecondary = false(length(grains), 1);

% Store [cluster_id, variant_id, refParentInd]
primaryVariantInfo = zeros(length(grains), 3); 
% Store [cluster_id, primary_var_id, secondary_var_id, refPrimInd]
secondaryVariantInfo = zeros(length(grains), 4); 

tolerance = 8 * degree;

% --- Identify Primary Twins ---
disp('Identifying primary twins...');

for i = 1:length(grains)
    if isParent(i), continue; end
    
    cID = componentId(i);
    if cID == 0, continue; end
    
    % Find parents in the same cluster
    clusterParents = find(isParent(:) & componentId(:) == cID);
    if isempty(clusterParents), continue; end
    
    % Find neighbors of this grain
    nbInds = pairs(pairs(:,1) == i, 2);
    nbInds = [nbInds; pairs(pairs(:,2) == i, 1)];
    
    parentNbs = intersect(nbInds, clusterParents);
    
    if ~isempty(parentNbs)
        % Use largest adjacent parent
        [~, maxIdx] = max(grains(parentNbs).numPixel);
        refParentInd = parentNbs(maxIdx);
    else
        % Fallback: use the largest parent in the cluster
        [~, maxIdx] = max(grains(clusterParents).numPixel);
        refParentInd = clusterParents(maxIdx);
    end
    
    refParentOri = grains(refParentInd).meanOrientation;
    
    % Check twin variant
    mori = inv(refParentOri) * grains(i).meanOrientation;
    [minAngle, bestFit] = min(angle_outer(mori, [tS_sym.parentTwinMisorientation]));
    
    if minAngle < tolerance
        isPrimary(i) = true;
        primaryVariantInfo(i, :) = [cID, bestFit, refParentInd];
    end
end

% --- Identify Secondary Twins ---
disp('Identifying secondary twins...');

for i = 1:length(grains)
    if isParent(i) || isPrimary(i), continue; end
    
    cID = componentId(i);
    if cID == 0, continue; end
    
    % Find primary twins in the same cluster
    clusterPrims = find(isPrimary(:) & componentId(:) == cID);
    if isempty(clusterPrims), continue; end
    
    % Find neighbors
    nbInds = pairs(pairs(:,1) == i, 2);
    nbInds = [nbInds; pairs(pairs(:,2) == i, 1)];
    
    primNbs = intersect(nbInds, clusterPrims);
    
    if ~isempty(primNbs)
        % Use largest adjacent primary twin
        [~, maxIdx] = max(grains(primNbs).numPixel);
        refPrimInd = primNbs(maxIdx);
    else
        % Fallback: use the largest primary twin in the cluster
        [~, maxIdx] = max(grains(clusterPrims).numPixel);
        refPrimInd = clusterPrims(maxIdx);
    end
    
    refPrimOri = grains(refPrimInd).meanOrientation;
    
    % Check twin variant
    mori = inv(refPrimOri) * grains(i).meanOrientation;
    [minAngle, bestFit] = min(angle_outer(mori, [tS_sym.parentTwinMisorientation]));
    
    if minAngle < tolerance
        isSecondary(i) = true;
        pVar = primaryVariantInfo(refPrimInd, 2);
        secondaryVariantInfo(i, :) = [cID, pVar, bestFit, refPrimInd];
    end
end

%% Visualization
disp('Visualizing the twin analysis results...');

% Create a new figure for the final plot
figure;
plot(grains, grains.meanOrientation, 'facealpha', 0.3)
hold on

% Plot unidentified grains in bright red
unidentifiedGrains = grains(~isParent & ~isPrimary & ~isSecondary);
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
    grainsOfParentCluster = grains(isParent(:) & componentId(:) == pID);
    if ~isempty(grainsOfParentCluster)
        text(grainsOfParentCluster, ['P_{' num2str(pID) '}'], 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    end
end

% Define colormap for primary twin variants
colors = lines(length(tS_sym));

% Plot primary twins with variant-specific colors and labels
primaryInfo = primaryVariantInfo(isPrimary,:);
uniquePV = unique(primaryInfo(:, 1:2), 'rows');

for k = 1:size(uniquePV, 1)
    pID = uniquePV(k, 1);
    vID = uniquePV(k, 2);
    if pID == 0 || vID == 0, continue; end
    
    grainsOfPV = grains(isPrimary(:) & primaryVariantInfo(:,1) == pID & primaryVariantInfo(:,2) == vID);
    
    if ~isempty(grainsOfPV)
        hold on
        plot(grainsOfPV, 'facecolor', colors(vID, :));
        label = sprintf('P_{%d}T_{%d}', pID, vID);
        text(grainsOfPV, label, 'HorizontalAlignment', 'center');
    end
end

% Plot secondary twins in pink
secondaryGrains = grains(isSecondary);
if ~isempty(secondaryGrains)
    hold on
    plot(secondaryGrains, 'facecolor', 'm');
    secInfo = secondaryVariantInfo(isSecondary, :);
    for k = 1:size(secInfo, 1)
        label = sprintf('T_{%d,%d}', secInfo(k,2), secInfo(k,3));
        text(secondaryGrains(k), label, 'HorizontalAlignment', 'center');
    end
end

hold off
title('Primary and Secondary Twin Analysis of Pillar8');
legend off

%% 6. Quiver K1 Plane Traces
disp('Plotting traces of K1 planes for primary and secondary twins...');

% Create a new figure for the traces
figure;
plot(grains, grains.meanOrientation, 'facealpha', 0.3)
hold on
plot(grains.boundary, 'lineWidth', 1.5);

% Plot parent grains
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
for i = 1:length(parentClusters)
    pID = parentClusters(i);
    if pID == 0, continue; end
    grainsOfParentCluster = grains(isParent(:) & componentId(:) == pID);
    if ~isempty(grainsOfParentCluster)
        text(grainsOfParentCluster, ['P_{' num2str(pID) '}'], 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    end
end

% --- Plot Primary Twin K1 Traces ---
primaryInds = find(isPrimary);
for i = 1:length(primaryInds)
    gInd = primaryInds(i);
    vID = primaryVariantInfo(gInd, 2);
    refParentInd = primaryVariantInfo(gInd, 3);
    
    if vID == 0 || refParentInd == 0, continue; end
    
    refParentOri = grains(refParentInd).meanOrientation;
    k1_crystal = tS_sym(vID).k1;
    
    % Transform to specimen coordinates using reference parent orientation
    k1_specimen = refParentOri * k1_crystal;
    k1_specimen.antipodal = 1;
    
    trace_dir = cross(k1_specimen, zvector);
    hold on
    quiver(grains(gInd), trace_dir, 'autoScale', 'off', 'color', colors(vID,:), 'LineWidth', 2);
end

% --- Plot Secondary Twin K1 Traces ---
secondaryInds = find(isSecondary);
for i = 1:length(secondaryInds)
    gInd = secondaryInds(i);
    s_vID = secondaryVariantInfo(gInd, 3);
    refPrimInd = secondaryVariantInfo(gInd, 4);
    
    if s_vID == 0 || refPrimInd == 0, continue; end
    
    refPrimOri = grains(refPrimInd).meanOrientation;
    k1_crystal = tS_sym(s_vID).k1;
    
    % Transform using the reference primary twin orientation
    k1_specimen = refPrimOri * k1_crystal;
    k1_specimen.antipodal = 1;
    
    trace_dir = cross(k1_specimen, zvector);
    hold on
    quiver(grains(gInd), trace_dir, 'autoScale', 'off', 'color', 'm', 'LineWidth', 2);
end

hold off
title('K1 Plane Traces of Identified Twins');
