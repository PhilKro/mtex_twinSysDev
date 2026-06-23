classdef TwinNetwork < handle
    % TwinNetwork - A class to represent and analyze graph-based twinning networks
    % It integrates with the twinSystem class to extract physical hierarchies
    % and deformation metrics from EBSD data.

    properties
        % Inputs
        grains      % grains object
        gB          % grain boundary object
        tS          % twinSystem object
        CS          % crystalSymmetry object
        
        % Settings
        settings    % Struct with configuration toggles
        
        % Twin System Attributes
        twinMisoris % Array of misorientations for all variants
        numVariants % Total number of variants
        doubleTwinMisoris % Misorientations for double twins
        doubleTwinPairsVariants % Variant pairs for double twins
        
        % Graph Topology
        activeGrains        % Array of grain IDs involved in twinning
        uPairs              % Single twin physical boundary pairs
        uDoublePairs        % Double twin physical boundary pairs
        virtualNodeBase     % Node ID base for virtual (ghost) nodes
        Q_G                 % Quotient Graph
        G_equiv             % Equivalence Graph
        bins                % Super Node assignments for each grain
        compBins            % Twin network connected component assignments
        uniqueComps         % Unique component IDs
        superNodeArea       % Area of each physical Super Node
        
        % Computed Attributes
        grainGen            % Generation depth of each grain
        spawnVariant        % The exact variant that spawned the grain
        variantPath         % Cell array of variant sequences (e.g. [1, 5])
        idealOri            % The ideal reference frame orientation
        parentGrainsInds    % Indices of Parent (Gen 0) grains
        parentSuperNodes    % Array of Parent Super Node IDs
        
        % Deformation Metrics
        defGradient         % Deformation gradient tensors for each grain
    end

    methods
        function obj = TwinNetwork(grains, gB, tS, varargin)
            % TwinNetwork Constructor
            obj.grains = grains;
            obj.gB = gB;
            obj.tS = tS;
            
            % Default Settings
            p = inputParser;
            addParameter(p, 'detectConsumedTwins', true, @islogical);
            addParameter(p, 'maxGenerationLimit', 2, @isnumeric);
            addParameter(p, 'twinThresh', 5 * degree, @isnumeric);
            parse(p, varargin{:});
            obj.settings = p.Results;
            
            % Setup twin system attributes
            obj.CS = tS.CS;
            tS_sym = tS.symmetrise;
            obj.twinMisoris = tS_sym.parentTwinMisorientation;
            obj.numVariants = length(obj.twinMisoris);
                       
            % Double Twin Variants setup
            obj.doubleTwinMisoris = [];
            obj.doubleTwinPairsVariants = [];
            if obj.settings.detectConsumedTwins
                for i = 1:obj.numVariants
                    for j = 1:obj.numVariants
                        misori = obj.twinMisoris(i) * obj.twinMisoris(j);
                        % A double twin must be physically distinct from the parent orientation.
                        % If it resolves to a low angle, it's just untwinning back to the parent (e.g. T4,4)!
                        if misori.angle > 15 * degree
                            obj.doubleTwinMisoris = [obj.doubleTwinMisoris; misori];
                            obj.doubleTwinPairsVariants = [obj.doubleTwinPairsVariants; i, j];
                        end
                    end
                end
            end
        end
        
        function obj = buildTopology(obj)
            % Identifies boundaries and builds the Equivalence and Quotient Graphs
            
            validEdges = obj.gB.grainId(:,1) > 0 & obj.gB.grainId(:,2) > 0;
            
            % 1. Primary Twins
            isTwinBoundary = any(angle(obj.gB.misorientation, obj.twinMisoris) < obj.settings.twinThresh, 2);
            isTwinBoundary = isTwinBoundary & validEdges;
            
            twinEdges = obj.gB.grainId(isTwinBoundary, :);
            obj.uPairs = unique(sort(twinEdges, 2), 'rows');
            
            % 2. Double Twins (Consumed Primary Twins)
            isDoubleTwinBoundary = false(length(obj.gB), 1);
            if obj.settings.detectConsumedTwins
                for k = 1:length(obj.doubleTwinMisoris)
                    match = angle(obj.gB.misorientation, obj.doubleTwinMisoris(k)) < obj.settings.twinThresh;
                    isDoubleTwinBoundary = isDoubleTwinBoundary | match;
                end
                isDoubleTwinBoundary = isDoubleTwinBoundary & validEdges & ~isTwinBoundary;
            end
            
            doubleTwinEdges = obj.gB.grainId(isDoubleTwinBoundary, :);
            obj.uDoublePairs = unique(sort(doubleTwinEdges, 2), 'rows');
            
            % 3. Equivalence Graph (Grains connected by low angle boundaries)
            obj.activeGrains = unique([obj.uPairs(:); obj.uDoublePairs(:)]);
            
            if isempty(obj.activeGrains)
                warning('No active twinning grains found.');
                return;
            end
            
            lowAngleThresh = 15 * degree;
            isEquivBoundary = obj.gB.misorientation.angle < lowAngleThresh & validEdges;
            equivBoundaryPairs = obj.gB.grainId(isEquivBoundary, :);
            
            equivEdges = [];
            for i = 1:size(equivBoundaryPairs, 1)
                if ismember(equivBoundaryPairs(i,1), obj.activeGrains) && ismember(equivBoundaryPairs(i,2), obj.activeGrains)
                    equivEdges = [equivEdges; sort(equivBoundaryPairs(i,:))];
                end
            end
            
            % Self loops to guarantee inclusion
            if isempty(equivEdges)
                equivEdges = [obj.activeGrains, obj.activeGrains];
            else
                equivEdges = [equivEdges; [obj.activeGrains, obj.activeGrains]];
            end
            
            maxGid = max(obj.grains.id);
            equivEdges = [equivEdges; [maxGid, maxGid]];
            
            obj.G_equiv = graph(equivEdges(:,1), equivEdges(:,2));
            obj.bins = conncomp(obj.G_equiv);
            
            % 4. Quotient Graph
            qPairs = [obj.bins(obj.uPairs(:,1))', obj.bins(obj.uPairs(:,2))'];
            if ~isempty(qPairs)
                isSelfLoop = qPairs(:,1) == qPairs(:,2);
                qPairs(isSelfLoop, :) = [];
            end
            
            obj.virtualNodeBase = max(obj.bins);
            numVirtualNodes = size(obj.uDoublePairs, 1);
            
            for i = 1:numVirtualNodes
                vNode = obj.virtualNodeBase + i;
                u = obj.uDoublePairs(i, 1);
                v = obj.uDoublePairs(i, 2);
                qPairs = [qPairs; obj.bins(u), vNode; vNode, obj.bins(v)];
            end
            
            if isempty(qPairs)
                obj.Q_G = graph();
            else
                obj.Q_G = graph(qPairs(:,1), qPairs(:,2));
            end
            
            maxNodeID = obj.virtualNodeBase + numVirtualNodes;
            if numnodes(obj.Q_G) < maxNodeID
                obj.Q_G = addnode(obj.Q_G, maxNodeID - numnodes(obj.Q_G));
            end
            
            % 5. Decompose Networks
            obj.compBins = conncomp(obj.Q_G);
            obj.uniqueComps = unique(obj.compBins);
            
            % Precompute super node areas
            obj.superNodeArea = zeros(obj.virtualNodeBase, 1);
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                sn = obj.bins(gid);
                ind = id2ind(obj.grains, gid);
                obj.superNodeArea(sn) = obj.superNodeArea(sn) + obj.grains(ind).grainSize;
            end
        end
        
        function obj = propagateReferences(obj)
            % Executes the BFS and 2-pass variant assignment
            if isempty(obj.activeGrains)
                return;
            end
            
            maxGid = max(obj.grains.id);
            obj.grainGen = inf(maxGid, 1);
            obj.spawnVariant = zeros(maxGid, 1);
            obj.variantPath = cell(maxGid, 1);
            obj.idealOri = orientation.nan(maxGid, 1, obj.CS);
            
            obj.parentGrainsInds = [];
            obj.parentSuperNodes = [];
            
            for c = 1:length(obj.uniqueComps)
                compID = obj.uniqueComps(c);
                nodesInComp = find(obj.compBins == compID);
                physicalNodesInComp = nodesInComp(nodesInComp <= obj.virtualNodeBase);
                
                if isempty(physicalNodesInComp)
                    continue;
                end
                
                activeGrainsInComp = obj.activeGrains(ismember(obj.bins(obj.activeGrains), physicalNodesInComp));
                if isempty(activeGrainsInComp)
                    continue;
                end
                
                [~, maxIdx] = max(obj.superNodeArea(physicalNodesInComp));
                mainParentSuperNode = physicalNodesInComp(maxIdx);
                
                % Grab the orientation of the main parent
                mainParentGids = obj.activeGrains(obj.bins(obj.activeGrains) == mainParentSuperNode);
                mainParentInds = id2ind(obj.grains, mainParentGids);
                globalParentOri = mean(obj.grains(mainParentInds).meanOrientation, 'robust');
                
                % Find ALL Super Nodes in this component that align with the parent orientation (isolated shards)
                compParentSuperNodes = [];
                for i = 1:length(physicalNodesInComp)
                    sn = physicalNodesInComp(i);
                    snGids = obj.activeGrains(obj.bins(obj.activeGrains) == sn);
                    snInds = id2ind(obj.grains, snGids);
                    snOri = mean(obj.grains(snInds).meanOrientation, 'robust');
                    
                    if angle(snOri, globalParentOri) < 15 * degree
                        compParentSuperNodes = [compParentSuperNodes; sn];
                    end
                end
                
                obj.parentSuperNodes = [obj.parentSuperNodes; compParentSuperNodes];
                
                % Multi-source BFS from all parent shards simultaneously
                dist_matrix = distances(obj.Q_G, compParentSuperNodes);
                if size(dist_matrix, 1) > 1
                    dist_from_parent = min(dist_matrix, [], 1);
                else
                    dist_from_parent = dist_matrix;
                end
                
                for i = 1:length(activeGrainsInComp)
                    gid = activeGrainsInComp(i);
                    g = dist_from_parent(obj.bins(gid));
                    if g > obj.settings.maxGenerationLimit
                        obj.grainGen(gid) = inf;
                    else
                        obj.grainGen(gid) = g;
                    end
                end
                
                compParentGids = activeGrainsInComp(obj.grainGen(activeGrainsInComp) == 0);
                compParentInds = id2ind(obj.grains, compParentGids);
                obj.parentGrainsInds = [obj.parentGrainsInds; compParentInds];
                
                globalParentOri = mean(obj.grains(compParentInds).meanOrientation, 'robust');
                obj.idealOri(compParentGids) = globalParentOri;
                for i = 1:length(compParentGids)
                    obj.variantPath{compParentGids(i)} = [];
                end
                
                compGens = obj.grainGen(activeGrainsInComp);
                maxGenComp = max(compGens(compGens < inf));
                if isempty(maxGenComp) || maxGenComp == 0
                    continue;
                end
                
                for gen = 1:maxGenComp
                    genGids = activeGrainsInComp(obj.grainGen(activeGrainsInComp) == gen);
                    successfulGids = [];
                    
                    % Pass 1: Direct Physical Edges
                    for i = 1:length(genGids)
                        child_gid = genGids(i);
                        child_ind = id2ind(obj.grains, child_gid);
                        O_child = obj.grains(child_ind).meanOrientation;
                        
                        foundSpawningVariant = false;
                        edgesWithChild = find(obj.uPairs(:,1) == child_gid | obj.uPairs(:,2) == child_gid);
                        for e = 1:length(edgesWithChild)
                            idx = edgesWithChild(e);
                            parent_gid = obj.uPairs(idx, 1);
                            if parent_gid == child_gid, parent_gid = obj.uPairs(idx, 2); end
                            
                            if obj.grainGen(parent_gid) == gen - 1 && ~isempty(obj.variantPath{parent_gid}) || (gen == 1 && obj.grainGen(parent_gid) == 0)
                                O_ideal_parent = obj.idealOri(parent_gid);
                                O_actual_parent = obj.grains(id2ind(obj.grains, parent_gid)).meanOrientation;
                                O_parent_fixed = O_actual_parent.project2FundamentalRegion(O_ideal_parent);
                                
                                fit = zeros(1, obj.numVariants);
                                for v = 1:obj.numVariants
                                    fit(v) = angle(O_parent_fixed * obj.twinMisoris(v), O_child);
                                end
                                
                                [~, bestVar] = min(fit);
                                obj.spawnVariant(child_ind) = bestVar;
                                obj.idealOri(child_gid) = O_parent_fixed * obj.twinMisoris(bestVar);
                                obj.variantPath{child_gid} = [obj.variantPath{parent_gid}, bestVar];
                                foundSpawningVariant = true;
                                break;
                            end
                        end
                        
                        if foundSpawningVariant
                            successfulGids = [successfulGids; child_gid];
                            continue;
                        end
                        
                        % Double Twins
                        edgesWithChild = find(obj.uDoublePairs(:,1) == child_gid | obj.uDoublePairs(:,2) == child_gid);
                        for e = 1:length(edgesWithChild)
                            idx = edgesWithChild(e);
                            parent_gid = obj.uDoublePairs(idx, 1);
                            if parent_gid == child_gid, parent_gid = obj.uDoublePairs(idx, 2); end
                            
                            if obj.grainGen(parent_gid) == gen - 2 && ~isempty(obj.variantPath{parent_gid}) || (gen == 2 && obj.grainGen(parent_gid) == 0)
                                O_ideal_parent = obj.idealOri(parent_gid);
                                O_actual_parent = obj.grains(id2ind(obj.grains, parent_gid)).meanOrientation;
                                O_parent_fixed = O_actual_parent.project2FundamentalRegion(O_ideal_parent);
                                
                                fit = zeros(size(obj.doubleTwinMisoris, 1), 1);
                                for v = 1:length(obj.doubleTwinMisoris)
                                    fit(v) = angle(O_parent_fixed * obj.doubleTwinMisoris(v), O_child);
                                end
                                
                                [~, bestDoubleIdx] = min(fit);
                                v1 = obj.doubleTwinPairsVariants(bestDoubleIdx, 1);
                                v2 = obj.doubleTwinPairsVariants(bestDoubleIdx, 2);
                                
                                obj.spawnVariant(child_ind) = v2;
                                obj.idealOri(child_gid) = O_parent_fixed * obj.doubleTwinMisoris(bestDoubleIdx);
                                obj.variantPath{child_gid} = [obj.variantPath{parent_gid}, v1, v2];
                                successfulGids = [successfulGids; child_gid];
                                break;
                            end
                        end
                    end
                    
                    % Pass 2: Broadcast
                    failedGids = setdiff(genGids, successfulGids);
                    for i = 1:length(failedGids)
                        child_gid = failedGids(i);
                        child_ind = id2ind(obj.grains, child_gid);
                        sn = obj.bins(child_gid);
                        
                        peer_gids = successfulGids(obj.bins(successfulGids) == sn);
                        if ~isempty(peer_gids)
                            peer_gid = peer_gids(1);
                            peer_ind = id2ind(obj.grains, peer_gid);
                            
                            obj.spawnVariant(child_ind) = obj.spawnVariant(peer_ind);
                            obj.idealOri(child_gid) = obj.idealOri(peer_gid);
                            obj.variantPath{child_gid} = obj.variantPath{peer_gid};
                        end
                    end
                end
            end
        end
        
        function obj = calculateDeformation(obj)
            % Compute total twin deformation gradient for each grain
            if isempty(obj.activeGrains)
                return;
            end
            
            tS_sym = obj.tS.symmetrise;
            obj.defGradient = cell(length(obj.grains), 1);
            
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                gInd = id2ind(obj.grains, gid);
                path = obj.variantPath{gid};
                
                % Parent orientation is required to transform into specimen coordinates
                if isempty(obj.idealOri(gid)) || isnan(obj.idealOri(gid).phi1)
                    pOri = obj.grains(gInd).meanOrientation;
                else
                    pOri = obj.idealOri(gid);
                end
                
                % Call the clean graph-friendly deformation cascader on the symmetrized array
                obj.defGradient{gInd} = tS_sym.calcDeformationSequence(path, pOri);
            end
        end
        
        function plotQuotientGraph(obj)
            if isempty(obj.Q_G) || numnodes(obj.Q_G) == 0
                return;
            end
            figure;
            p = plot(obj.Q_G, 'Layout', 'force', 'MarkerSize', 6);
            if ~isempty(obj.parentSuperNodes)
                highlight(p, obj.parentSuperNodes, 'NodeColor', 'r', 'MarkerSize', 10);
            end
            title('Quotient Graph Topology (Red = Parent Matrices)');
        end
        
        function plotGenerations(obj)
            genProp = inf(length(obj.grains), 1);
            activeInds = id2ind(obj.grains, obj.activeGrains);
            genProp(activeInds) = obj.grainGen(obj.activeGrains);
            
            compGens = genProp(genProp < inf);
            if isempty(compGens)
                maxGenPlot = 0;
            else
                maxGenPlot = max(compGens);
            end
            
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0.2); 
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            if maxGenPlot > 0
                cmapGen = lines(maxGenPlot + 1);
                for gen = 0:maxGenPlot
                    inds = activeInds(genProp(activeInds) == gen);
                    if ~isempty(inds)
                        hold on;
                        plot(obj.grains(inds), 'faceColor', cmapGen(gen + 1, :), 'faceAlpha', 0.8);
                    end
                end
            end
            title('Twinning Generations');
        end
        
        function plotTwinVariants(obj)
            % Plots Primary and Secondary variants with text labels
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0.2);
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            parentIndsPlot = [];
            type1IndsPlot = [];
            type2IndsPlot = [];
            labels = cell(length(obj.grains), 1);

            % Type 1 cozonal pairs (hexagonal specific or general if implemented)
            tS_sym = obj.tS.symmetrise;
            idx = nchoosek(1:length(tS_sym.rotAxis), 2);
            valid = eq(tS_sym.rotAxis(idx(:,1)), tS_sym.rotAxis(idx(:,2)), 'noSymmetry', 'antipodal');
            cozonalPairs = idx(valid, :);
            
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                gInd = id2ind(obj.grains, gid);
                gen = obj.grainGen(gid);
                path = obj.variantPath{gid};
                
                if gen == 0
                    parentIndsPlot = [parentIndsPlot; gInd];
                    labels{gInd} = 'P';
                elseif gen == 1
                    type1IndsPlot = [type1IndsPlot; gInd];
                    v1 = path(1);
                    labels{gInd} = sprintf('T_{%d}', v1);
                elseif gen == 2
                    v1 = path(1);
                    v2 = path(2);
                    isCozonal = any((cozonalPairs(:,1) == v1 & cozonalPairs(:,2) == v2) | ...
                                    (cozonalPairs(:,1) == v2 & cozonalPairs(:,2) == v1));
                    if isCozonal
                        type1IndsPlot = [type1IndsPlot; gInd];
                    else
                        type2IndsPlot = [type2IndsPlot; gInd];
                    end
                    labels{gInd} = sprintf('T_{%d,%d}', v1, v2);
                elseif gen <= obj.settings.maxGenerationLimit && gen > 2
                    type2IndsPlot = [type2IndsPlot; gInd];
                    str = sprintf('%d,', path);
                    labels{gInd} = sprintf('T_{%s}', str(1:end-1));
                end
            end
            
            if ~isempty(parentIndsPlot)
                hold on;
                plot(obj.grains(parentIndsPlot), 'faceColor', [0.7 0.7 0.7], 'faceAlpha', 0.8);
            end
            if ~isempty(type1IndsPlot)
                hold on;
                plot(obj.grains(type1IndsPlot), 'faceColor', [0 0.5 0.8], 'faceAlpha', 0.8);
            end
            if ~isempty(type2IndsPlot)
                hold on;
                plot(obj.grains(type2IndsPlot), 'faceColor', [1 0.4 0.8], 'faceAlpha', 0.8);
            end
            text(obj.grains, labels, 'FontSize', 8, 'BackgroundColor', 'w');
            title('Twinning Variants (Blue=Type1, Pink=Type2)');
        end
    end
end
