classdef TwinNetwork < handle
    % TwinNetwork - A class to represent and analyze graph-based twinning networks
    % It integrates with the twinSystem class to extract physical hierarchies
    % and deformation metrics from EBSD data.
    %
    % RECOMMENDED WORKFLOW:
    %   tN = TwinNetwork(grains, gB, twinSystems);
    %   tN = tN.analyze();
    %
    % The individual methods to build and analyze the graph (buildTopology, propagateReferences, calculateDeformation) 
    % remain available for step-by-step execution and custom workflows.

    properties
        % Inputs
        grains      % grains object
        gB          % grain boundary object
        twinSystems % Array or cell array of twinSystem objects
        CS          % crystalSymmetry object
        
        % Settings
        settings    % Struct with configuration toggles
        
        % Twin System Attributes
        twinRegistry % Struct array mapping variant sequence topologies
        
        % Graph Topology
        activeGrains        % Array of grain IDs involved in twinning
        physicalEdges       % Single unified array of all twinning boundaries (Nx2)
        edgeToRegistryMap   % containers.Map mapping uint64 bitshifted edges to array of registry indices
        virtualNodeBase     % Node ID base for virtual (ghost) nodes
        Q_G                 % Quotient Graph
        compBins            % Twin network connected component assignments
        uniqueComps         % Unique component IDs
        
        % Computed Attributes
        grainGen            % Generation depth of each grain
        spawnVariant        % The exact variant that spawned the grain
        variantPath         % Cell array of variant sequences (e.g. [1, 5])
        idealOri            % The ideal reference frame orientation
        parentGrainsInds    % Indices of Parent (Gen 0) grains
        
        % Deformation Metrics
        defGradient         % Deformation gradient tensors for each grain
    end

    methods
        function obj = TwinNetwork(grains, gB, twinSystems, varargin)
            % TwinNetwork Constructor
            obj.grains = grains;
            obj.gB = gB;
            if ~iscell(twinSystems)
                obj.twinSystems = {twinSystems};
            else
                obj.twinSystems = twinSystems;
            end
            
            % Default Settings
            p = inputParser;
            addParameter(p, 'detectDoubleTwins', false, @islogical); % Detects boundaries between parent and double twins of the passed twinSystems. Can also be done manually in the detectedTwinPaths (e.g. {[1], [2]} becomes {[1], [2], [1, 1], [2, 2]}).
            addParameter(p, 'maxGenerationLimit', 2, @isnumeric); % the depth to which twin generations are searched, 2 -> secondary twinning is resolved, 3 -> teritary twinning, etc.
            addParameter(p, 'twinThresh', 5 * degree, @isnumeric); % misorientation deviation from the twin oreinatation that is allowed for the boundary to be recognized as twin
            addParameter(p, 'detectedTwinPaths', {}, @iscell); % Specify which boundaries between parent and nested twins of the passed twinSystems are detected as edges of the graph. E.g. for two passed twinSystems the default is {[1], [2]}, can also be specified as {[1], [2], [1, 2], [2, 1], [1, 2, 1]} to detect different nested twin boudaries.            
            parse(p, varargin{:});
            obj.settings = p.Results;
            
            % Setup twin system attributes
            obj.CS = obj.twinSystems{1}.CS;
            obj.twinRegistry = struct('name', {}, 'genStep', {}, 'sysPath', {}, 'misoris', {}, 'variantMap', {}, 'edges', {});
            
            detectedPaths = obj.settings.detectedTwinPaths;
            if isempty(detectedPaths)
                % Auto-populate primary twins into the detected paths if no paths are provided
                for i = 1:length(obj.twinSystems)
                    detectedPaths{end+1} = [i]; % Primary twin
                end
                
                % Add Double Twins auto-generation if flag is set, matching old behavior
                if obj.settings.detectDoubleTwins
                    for i = 1:length(obj.twinSystems)
                        detectedPaths{end+1} = [i, i];
                    end
                end
            end
            
            % Populate twinRegistry based on detectedPaths
            for i = 1:length(detectedPaths)
                path = detectedPaths{i};
                genStep = length(path);
                
                % Recursively generate compounded misorientations
                [compMisoris, compVarMap] = TwinNetwork.computeCompoundTwins(obj.twinSystems, path);
                
                % Filter out identity mappings (low angle)
                validIdx = compMisoris.angle > 15 * degree;
                if ~any(validIdx)
                    continue;
                end
                
                regEntry.name = sprintf('Path [%s]', join(string(path), ','));
                regEntry.genStep = genStep;
                regEntry.sysPath = path;
                regEntry.misoris = compMisoris(validIdx);
                regEntry.variantMap = compVarMap(validIdx, :);
                regEntry.edges = [];
                
                obj.twinRegistry(end+1) = regEntry;
            end
        end
        
        function obj = analyze(obj)
            % Executes the standard full twin analysis workflow.
            % 
            % Workflow Steps:
            % 1. buildTopology(): Identifies boundaries, computes Equivalence Graphs, 
            %    and groups connected components. Yields graph structure and component bins.
            % 2. propagateReferences(): Runs BFS traversal to assign generations, 
            %    variant paths, and ideal orientations. Yields exact grain hierarchies.
            % 3. calculateDeformation(): Derives total deformation gradient tensors 
            %    based on the assigned variant paths. Yields continuous mechanic tensors.
            
            obj = obj.buildTopology();
            obj = obj.propagateReferences();
            obj = obj.calculateDeformation();
        end
        
        function obj = buildTopology(obj)
            % Identifies boundaries and builds the Equivalence and Quotient Graphs
            
            topo = TwinNetwork.generateTopologyStruct(obj.gB, obj.grains, obj.settings, obj.twinRegistry);
            
            obj.twinRegistry = topo.twinRegistry;
            obj.physicalEdges = topo.physicalEdges;
            obj.edgeToRegistryMap = topo.edgeToRegistryMap;
            obj.activeGrains = topo.activeGrains;
            obj.virtualNodeBase = topo.virtualNodeBase;
            obj.Q_G = topo.Q_G;
            obj.compBins = topo.compBins;
            obj.uniqueComps = topo.uniqueComps;
            
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
            idToIndMap = zeros(maxGid, 1);
            idToIndMap(obj.grains.id) = 1:length(obj.grains);
            allOris = obj.grains.meanOrientation;
            allSizes = obj.grains.grainSize;
            
            obj.grainGen = inf(maxGid, 1);
            obj.spawnVariant = zeros(maxGid, 1);
            obj.variantPath = cell(maxGid, 1);
            obj.idealOri = orientation.nan(maxGid, 1, obj.CS);
            
            parentGrainsIndsCell = cell(length(obj.uniqueComps), 1);
            
            for c = 1:length(obj.uniqueComps)
                compID = obj.uniqueComps(c);
                nodesInComp = find(obj.compBins == compID);
                
                % Since bins == grainId, physical nodes are exactly the grainIds
                activeGrainsInComp = nodesInComp(nodesInComp <= obj.virtualNodeBase);
                if isempty(activeGrainsInComp)
                    continue;
                end
                
                % Find the largest grain by grainSize
                compGrainsInds = idToIndMap(activeGrainsInComp);
                grainSizes = allSizes(compGrainsInds);
                [~, maxIdx] = max(grainSizes);
                
                mainParentGid = activeGrainsInComp(maxIdx);
                mainParentInd = idToIndMap(mainParentGid);
                globalParentOri = allOris(mainParentInd);
                
                % Find ALL grains in this component that align with the parent orientation
                compParentGids = [];
                for i = 1:length(activeGrainsInComp)
                    gid = activeGrainsInComp(i);
                    gInd = idToIndMap(gid);
                    gOri = allOris(gInd);
                    
                    if angle(gOri, globalParentOri) < 15 * degree
                        compParentGids = [compParentGids; gid];
                    end
                end
                
                % Multi-source BFS from all parent grains simultaneously
                dist_matrix = distances(obj.Q_G, compParentGids);
                if size(dist_matrix, 1) > 1
                    dist_from_parent = min(dist_matrix, [], 1);
                else
                    dist_from_parent = dist_matrix;
                end
                
                for i = 1:length(activeGrainsInComp)
                    gid = activeGrainsInComp(i);
                    g = dist_from_parent(gid);
                    if g > obj.settings.maxGenerationLimit
                        obj.grainGen(gid) = inf;
                    else
                        obj.grainGen(gid) = g;
                    end
                end
                
                compParentGids = activeGrainsInComp(obj.grainGen(activeGrainsInComp) == 0);
                compParentInds = idToIndMap(compParentGids);
                parentGrainsIndsCell{c} = compParentInds;
                
                globalParentOri = mean(allOris(compParentInds), 'robust');
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
                    
                    % Iterate through physical edges connected to this child
                    for i = 1:length(genGids)
                        child_gid = genGids(i);
                        child_ind = idToIndMap(child_gid);
                        O_child = allOris(child_ind);
                        
                        childEdges = obj.physicalEdges(obj.physicalEdges(:,1) == child_gid | obj.physicalEdges(:,2) == child_gid, :);
                        
                        for e = 1:size(childEdges, 1)
                            id1 = childEdges(e, 1);
                            id2 = childEdges(e, 2);
                            parent_gid = id1;
                            if parent_gid == child_gid
                                parent_gid = id2;
                            end
                            
                            key = bitshift(uint64(max(id1, id2)), 32) + uint64(min(id1, id2));
                            
                            if isKey(obj.edgeToRegistryMap, key)
                                regIndices = obj.edgeToRegistryMap(key);
                                
                                [success, bestVar, O_ideal, newPath] = obj.resolveVariantMatch(...
                                    child_gid, parent_gid, gen, O_child, regIndices, allOris, idToIndMap);
                                
                                if success
                                    obj.spawnVariant(child_ind) = bestVar;
                                    obj.idealOri(child_gid) = O_ideal;
                                    obj.variantPath{child_gid} = newPath;
                                    successfulGids = [successfulGids; child_gid];
                                    break;
                                end
                            end
                        end
                    end
                end
            end
            
            obj.parentGrainsInds = cell2mat(parentGrainsIndsCell);
        end
        
        function obj = calculateDeformation(obj)
            % Compute total twin deformation gradient for each grain
            if isempty(obj.activeGrains)
                return;
            end
            
            maxGid = max(obj.grains.id);
            idToIndMap = zeros(maxGid, 1);
            idToIndMap(obj.grains.id) = 1:length(obj.grains);
            allOris = obj.grains.meanOrientation;
            
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                gInd = idToIndMap(gid);
                path = obj.variantPath{gid};
                
                % Parent orientation is required to transform into specimen coordinates
                if isempty(obj.idealOri(gid)) || isnan(obj.idealOri(gid).phi1)
                    pOri = allOris(gInd);
                else
                    pOri = obj.idealOri(gid);
                end
                
                % Cascade deformation for multi-system paths
                F = eye(3);
                currentOri = pOri;
                for step = 1:size(path, 1)
                    sysIdx = path(step, 1);
                    varIdx = path(step, 2);
                    tS_sym = obj.twinSystems{sysIdx}.symmetrise;
                    
                    F_step = tS_sym.calcDeformationSequence(varIdx, currentOri);
                    F = F_step * F;
                    currentOri = currentOri * tS_sym.parentTwinMisorientation(varIdx);
                end
                
                obj.defGradient{gInd} = F;
            end
        end
        
        function plotQuotientGraph(obj)
            if isempty(obj.Q_G) || numnodes(obj.Q_G) == 0
                return;
            end
            figure;
            p = plot(obj.Q_G, 'Layout', 'force', 'MarkerSize', 6);
            
            if isempty(obj.parentGrainsInds)
                title('Quotient Graph Topology');
                return;
            end
            
            parentGids = obj.grains(obj.parentGrainsInds).id;
            highlight(p, parentGids, 'MarkerSize', 10, 'NodeColor', [0.3 0.3 0.3]);
            
            dist_matrix = distances(obj.Q_G, parentGids);
            if size(dist_matrix, 1) > 1
                nodeGen = min(dist_matrix, [], 1)';
            else
                nodeGen = dist_matrix';
            end
            
            maxGen = max(nodeGen(nodeGen < inf));
            if isempty(maxGen) || isnan(maxGen)
                maxGen = 0;
            end
            cmapGen = lines(maxGen + 1);
            
            for g = 0:maxGen
                % Physical nodes
                physNodes = find(nodeGen == g & (1:length(nodeGen))' <= obj.virtualNodeBase);
                if ~isempty(physNodes)
                    highlight(p, physNodes, 'NodeColor', cmapGen(g+1, :));
                end
                
                % Virtual nodes (faded by blending with white)
                virtNodes = find(nodeGen == g & (1:length(nodeGen))' > obj.virtualNodeBase);
                if ~isempty(virtNodes)
                    colorBlend = cmapGen(g+1, :) * 0.5 + [1 1 1] * 0.5;
                    highlight(p, virtNodes, 'NodeColor', colorBlend);
                end
            end
            
            title('Quotient Graph Topology (Node size = Parent, Colors = Generations, Faded = Virtual)');
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
            % Plots twin variants with text labels
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0);
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            labels = cell(length(obj.grains), 1);
            parentIndsPlot = [];
            twinIndsPlot = [];
            colorArray = zeros(length(obj.grains), 1);
            
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                gInd = id2ind(obj.grains, gid);
                gen = obj.grainGen(gid);
                path = obj.variantPath{gid};
                
                if gen == 0
                    parentIndsPlot = [parentIndsPlot; gInd];
                    labels{gInd} = 'P';
                elseif gen > 0 && ~isempty(path)
                    twinIndsPlot = [twinIndsPlot; gInd];
                    % Label as T_{v1, v2, ...} where v is the variant index
                    str = sprintf('%d,', path(:, 2));
                    labels{gInd} = sprintf('T_{%s}', str(1:end-1));
                    
                    % Give a color based on the final variant index
                    colorArray(gInd) = path(end, 2);
                end
            end
            
            if ~isempty(parentIndsPlot)
                hold on;
                plot(obj.grains(parentIndsPlot), 'faceColor', [0.7 0.7 0.7], 'faceAlpha', 0.8, 'DisplayName', 'P');
            end
            
            if ~isempty(twinIndsPlot)
                hold on;
                % Unique final variants
                uVars = unique(colorArray(twinIndsPlot));
                cmap = lines(max(uVars));
                for v = uVars'
                    vInds = twinIndsPlot(colorArray(twinIndsPlot) == v);
                    hold on
                    plot(obj.grains(vInds), 'faceColor', cmap(v, :), 'faceAlpha', 0.8);
                end
            end
            text(obj.grains, labels, 'FontSize', 8, 'BackgroundColor', 'w');
            title('Twinning Variants (Colored by Final Spawned Variant)');
        end
    end
    
    methods (Access = private)
        function [success, spawnVar, O_ideal_child, newPath] = resolveVariantMatch(obj, child_gid, parent_gid, gen, O_child, registryIndices, allOris, idToIndMap)
            success = false;
            spawnVar = 0;
            O_ideal_child = orientation.empty;
            newPath = [];
            
            for r = registryIndices
                regEntry = obj.twinRegistry(r);
                genStep = regEntry.genStep;
                misoriArray = regEntry.misoris;
                variantMap = regEntry.variantMap;
                sysPath = regEntry.sysPath;
                
                % Guard clause: Validate parent generation based on step size
                isValidParent = (obj.grainGen(parent_gid) == gen - genStep && ~isempty(obj.variantPath{parent_gid})) || ...
                                (gen == genStep && obj.grainGen(parent_gid) == 0);
                if ~isValidParent
                    continue;
                end
                
                O_ideal_parent = obj.idealOri(parent_gid);
                O_actual_parent = allOris(idToIndMap(parent_gid));
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
                % Append to path as [sysIdx, varIdx] pairs
                newPath = [obj.variantPath{parent_gid}; [sysPath(:), vars(:)]];
                success = true;
                return;
            end
        end
    end
    
    methods (Static)
        function [misoris, variantMap] = computeCompoundTwins(twinSysArray, path)
            tS = twinSysArray{path(1)};
            tS_sym = tS.symmetrise;
            misoris = tS_sym.parentTwinMisorientation;
            variantMap = (1:length(misoris))';
            
            for p = 2:length(path)
                tS_next = twinSysArray{path(p)};
                tS_next_sym = tS_next.symmetrise;
                mis_next = tS_next_sym.parentTwinMisorientation;
                
                newMisoris = [];
                newVariantMap = [];
                
                for i = 1:length(misoris)
                    for j = 1:length(mis_next)
                        newMisoris = [newMisoris; misoris(i) * mis_next(j)];
                        newVariantMap = [newVariantMap; variantMap(i, :), j];
                    end
                end
                misoris = newMisoris;
                variantMap = newVariantMap;
            end
        end

        function topo = generateTopologyStruct(gB, grains, settings, twinRegistry)
            % Identifies boundaries and builds the Equivalence and Quotient Graphs
            topo = struct();
            validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
            
            % Store the updated registry
            topo.twinRegistry = twinRegistry;
            topo.physicalEdges = [];
            topo.edgeToRegistryMap = containers.Map('KeyType', 'uint64', 'ValueType', 'any');
            
            activeGrains = [];
            
            % Find edges for each registry entry
            for r = 1:length(topo.twinRegistry)
                misoris = topo.twinRegistry(r).misoris;
                
                % Vectorized check against all misorientations for this sequence
                isBoundary = any(angle(gB.misorientation, misoris) < settings.twinThresh, 2);
                isBoundary = isBoundary & validEdges;
                
                edges = gB.grainId(isBoundary, :);
                uniqueEdges = unique(sort(edges, 2), 'rows');
                
                topo.twinRegistry(r).edges = uniqueEdges;
                activeGrains = [activeGrains; uniqueEdges(:)];
                topo.physicalEdges = [topo.physicalEdges; uniqueEdges];
                
                % Tag edges using bitshifted uint64 key
                if ~isempty(uniqueEdges)
                    allKeys = bitshift(uint64(max(uniqueEdges(:,1), uniqueEdges(:,2))), 32) + uint64(min(uniqueEdges(:,1), uniqueEdges(:,2)));
                    for e = 1:size(uniqueEdges, 1)
                        key = allKeys(e);
                        
                        if isKey(topo.edgeToRegistryMap, key)
                            topo.edgeToRegistryMap(key) = [topo.edgeToRegistryMap(key), r];
                        else
                            topo.edgeToRegistryMap(key) = r;
                        end
                    end
                end
            end
            
            topo.physicalEdges = unique(topo.physicalEdges, 'rows');
            topo.activeGrains = unique(activeGrains);
            
            if isempty(topo.activeGrains)
                topo.Q_G = graph();
                topo.virtualNodeBase = 0;
                topo.compBins = [];
                topo.uniqueComps = [];
                return;
            end
            maxGid = max(grains.id);
            
            % 4. Quotient Graph
            qPairs = [];
            topo.virtualNodeBase = maxGid;
            numVirtualNodes = 0;
            
            for r = 1:length(topo.twinRegistry)
                genStep = topo.twinRegistry(r).genStep;
                edges = topo.twinRegistry(r).edges;
                if isempty(edges)
                    continue; 
                end
                
                if genStep == 1
                    % Direct edges
                    pairs = edges;
                    isSelfLoop = pairs(:,1) == pairs(:,2);
                    pairs(isSelfLoop, :) = [];
                    qPairs = [qPairs; pairs];
                else
                    for i = 1:size(edges, 1)
                        numVirtualNodes = numVirtualNodes + 1;
                        vNode = topo.virtualNodeBase + numVirtualNodes;
                        u = edges(i, 1);
                        v = edges(i, 2);
                        qPairs = [qPairs; u, vNode; vNode, v];
                    end
                end
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
        end
    end
end
