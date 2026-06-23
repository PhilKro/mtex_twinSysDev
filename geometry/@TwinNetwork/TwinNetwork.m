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
            
            topo = TwinNetwork.generateTopologyStruct(obj.gB, obj.grains, obj.settings, obj.twinMisoris, obj.doubleTwinMisoris);
            
            obj.uPairs = topo.uPairs;
            obj.uDoublePairs = topo.uDoublePairs;
            obj.activeGrains = topo.activeGrains;
            obj.G_equiv = topo.G_equiv;
            obj.bins = topo.bins;
            obj.virtualNodeBase = topo.virtualNodeBase;
            obj.Q_G = topo.Q_G;
            obj.compBins = topo.compBins;
            obj.uniqueComps = topo.uniqueComps;
            obj.superNodeArea = topo.superNodeArea;
            
            if isempty(obj.activeGrains)
                warning('No active twinning grains found.');
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
                        
                        [success, bestVar, O_ideal, newPath] = obj.resolveVariantMatch(...
                            child_gid, gen, O_child, obj.uPairs, obj.twinMisoris, (1:obj.numVariants)', 1);
                        
                        if success
                            obj.spawnVariant(child_ind) = bestVar;
                            obj.idealOri(child_gid) = O_ideal;
                            obj.variantPath{child_gid} = newPath;
                            successfulGids = [successfulGids; child_gid];
                            continue;
                        end
                        
                        % Double Twins
                        [success, bestVar, O_ideal, newPath] = obj.resolveVariantMatch(...
                            child_gid, gen, O_child, obj.uDoublePairs, obj.doubleTwinMisoris, obj.doubleTwinPairsVariants, 2);
                        
                        if success
                            obj.spawnVariant(child_ind) = bestVar;
                            obj.idealOri(child_gid) = O_ideal;
                            obj.variantPath{child_gid} = newPath;
                            successfulGids = [successfulGids; child_gid];
                            continue;
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
            
            % Plot generation layers (Parent = Gen 0, Primary = Gen 1...)
            
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0.2); 
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            genProp = inf(max(obj.grains.id), 1);
            genProp(obj.activeGrains) = obj.grainGen(obj.activeGrains);
            maxGenPlot = max(genProp(genProp < inf));
            
            if ~isempty(maxGenPlot)
                cmapGen = lines(maxGenPlot + 1);
                for g = 0:maxGenPlot
                    gInds = id2ind(obj.grains, find(genProp == g));
                    if ~isempty(gInds)
                        if g == 0
                            hold on
                            plot(obj.grains(gInds), 'faceColor', [0.3 0.3 0.3], 'DisplayName', 'Parent (Gen 0)');
                        else
                            hold on
                            plot(obj.grains(gInds), 'faceColor', cmapGen(g+1, :), 'DisplayName', sprintf('Generation %d', g));
                        end
                    end
                end
            end
            title('Twin Generations (Filtered Bipartite Projection)')
            legend('show', 'Location', 'northeastoutside');
            hold off;
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
    
    methods (Access = private)
        function [success, spawnVar, O_ideal_child, newPath] = resolveVariantMatch(obj, child_gid, gen, O_child, edgeList, misoriArray, variantMap, genStep)
            success = false;
            spawnVar = 0;
            O_ideal_child = orientation.empty;
            newPath = [];
            
            edgesWithChild = find(edgeList(:,1) == child_gid | edgeList(:,2) == child_gid);
            if isempty(edgesWithChild)
                return;
            end
            
            for e = 1:length(edgesWithChild)
                idx = edgesWithChild(e);
                parent_gid = edgeList(idx, 1);
                if parent_gid == child_gid
                    parent_gid = edgeList(idx, 2); 
                end
                
                % Guard clause: Validate parent generation based on step size
                isValidParent = (obj.grainGen(parent_gid) == gen - genStep && ~isempty(obj.variantPath{parent_gid})) || ...
                                (gen == genStep && obj.grainGen(parent_gid) == 0);
                if ~isValidParent
                    continue;
                end
                
                O_ideal_parent = obj.idealOri(parent_gid);
                O_actual_parent = obj.grains(id2ind(obj.grains, parent_gid)).meanOrientation;
                O_parent_fixed = O_actual_parent.project2FundamentalRegion(O_ideal_parent);
                
                numMisoris = length(misoriArray);
                fit = zeros(numMisoris, 1);
                for v = 1:numMisoris
                    fit(v) = angle(O_parent_fixed * misoriArray(v), O_child);
                end
                
                [~, bestIdx] = min(fit);
                
                vars = variantMap(bestIdx, :);
                spawnVar = vars(end);
                
                O_ideal_child = O_parent_fixed * misoriArray(bestIdx);
                newPath = [obj.variantPath{parent_gid}, vars];
                success = true;
                return;
            end
        end
    end
    
    methods (Static)
        function topo = generateTopologyStruct(gB, grains, settings, twinMisoris, doubleTwinMisoris)
            % Identifies boundaries and builds the Equivalence and Quotient Graphs
            topo = struct();
            validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
            
            % 1. Primary Twins
            isTwinBoundary = any(angle(gB.misorientation, twinMisoris) < settings.twinThresh, 2);
            isTwinBoundary = isTwinBoundary & validEdges;
            
            twinEdges = gB.grainId(isTwinBoundary, :);
            topo.uPairs = unique(sort(twinEdges, 2), 'rows');
            
            % 2. Double Twins (Consumed Primary Twins)
            isDoubleTwinBoundary = false(length(gB), 1);
            if settings.detectConsumedTwins
                for k = 1:length(doubleTwinMisoris)
                    match = angle(gB.misorientation, doubleTwinMisoris(k)) < settings.twinThresh;
                    isDoubleTwinBoundary = isDoubleTwinBoundary | match;
                end
                isDoubleTwinBoundary = isDoubleTwinBoundary & validEdges & ~isTwinBoundary;
            end
            
            doubleTwinEdges = gB.grainId(isDoubleTwinBoundary, :);
            topo.uDoublePairs = unique(sort(doubleTwinEdges, 2), 'rows');
            
            % 3. Equivalence Graph (Grains connected by low angle boundaries)
            topo.activeGrains = unique([topo.uPairs(:); topo.uDoublePairs(:)]);
            
            if isempty(topo.activeGrains)
                topo.G_equiv = graph();
                topo.bins = [];
                topo.Q_G = graph();
                topo.virtualNodeBase = 0;
                topo.compBins = [];
                topo.uniqueComps = [];
                topo.superNodeArea = [];
                return;
            end
            
            lowAngleThresh = 15 * degree;
            isEquivBoundary = gB.misorientation.angle < lowAngleThresh & validEdges;
            equivBoundaryPairs = gB.grainId(isEquivBoundary, :);
            
            equivEdges = [];
            for i = 1:size(equivBoundaryPairs, 1)
                if ismember(equivBoundaryPairs(i,1), topo.activeGrains) && ismember(equivBoundaryPairs(i,2), topo.activeGrains)
                    equivEdges = [equivEdges; sort(equivBoundaryPairs(i,:))];
                end
            end
            
            % Self loops to guarantee inclusion
            if isempty(equivEdges)
                equivEdges = [topo.activeGrains, topo.activeGrains];
            else
                equivEdges = [equivEdges; [topo.activeGrains, topo.activeGrains]];
            end
            
            maxGid = max(grains.id);
            equivEdges = [equivEdges; [maxGid, maxGid]];
            
            topo.G_equiv = graph(equivEdges(:,1), equivEdges(:,2));
            topo.bins = conncomp(topo.G_equiv);
            
            % 4. Quotient Graph
            qPairs = [topo.bins(topo.uPairs(:,1))', topo.bins(topo.uPairs(:,2))'];
            if ~isempty(qPairs)
                isSelfLoop = qPairs(:,1) == qPairs(:,2);
                qPairs(isSelfLoop, :) = [];
            end
            
            topo.virtualNodeBase = max(topo.bins);
            numVirtualNodes = size(topo.uDoublePairs, 1);
            
            for i = 1:numVirtualNodes
                vNode = topo.virtualNodeBase + i;
                u = topo.uDoublePairs(i, 1);
                v = topo.uDoublePairs(i, 2);
                qPairs = [qPairs; topo.bins(u), vNode; vNode, topo.bins(v)];
            end
            
            if isempty(qPairs)
                topo.Q_G = graph();
            else
                topo.Q_G = graph(qPairs(:,1), qPairs(:,2));
            end
            
            maxNodeID = topo.virtualNodeBase + numVirtualNodes;
            if numnodes(topo.Q_G) < maxNodeID
                topo.Q_G = addnode(topo.Q_G, maxNodeID - numnodes(topo.Q_G));
            end
            
            % 5. Decompose Networks
            topo.compBins = conncomp(topo.Q_G);
            topo.uniqueComps = unique(topo.compBins);
            
            % Precompute super node areas
            topo.superNodeArea = zeros(topo.virtualNodeBase, 1);
            for i = 1:length(topo.activeGrains)
                gid = topo.activeGrains(i);
                sn = topo.bins(gid);
                ind = id2ind(grains, gid);
                topo.superNodeArea(sn) = topo.superNodeArea(sn) + grains(ind).grainSize;
            end
        end
    end
end
