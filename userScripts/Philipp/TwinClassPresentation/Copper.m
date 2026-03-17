%still needs some wort

clear
close all
mtexdata copper silent
CS = ebsd('indexed').CS;

% reconstruct grains
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',5*degree);

grains = smooth(grains,1);

% plot the grains
plot(grains,grains.meanOrientation)

%% Define {111} Twin System
% K1 = {111}, K2 = {111}
cs = grains(1).CS;
tS = twinSystem.fcc_111_111(cs);
tS_sym = tS.symmetrise;

%% Group Twin Related Grains
pairs = grains.neighbors;

% misorientation between neighbors
mori = inv(grains(pairs(:,1)).meanOrientation) .* grains(pairs(:,2)).meanOrientation;

% check for twin relationship
isTwin = min(angle_outer(mori, [-tS_sym.parentTwinMisorientation]), [], 2) < 5*degree;

% group grains
G = graph(pairs(isTwin,1), pairs(isTwin,2), [], length(grains));
componentId = conncomp(G);

% plot the groups
figure
plot(grains, ind2color(componentId, viridis(max(componentId))))
title('Twin Related Grain Clusters')

%% Identify Parent Grains
% Weights for criteria: [Size, Variants, InteractionWork]
weights = [1, 1, 1]; 

% Stress tensor for Interaction Work (Tension along Y)
sigma = stressTensor.uniaxial(yvector);
%sigma = - stressTensor.uniaxial(yvector);

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