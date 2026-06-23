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
        edgeToRegistryMap   % containers.Map: uint64 edge key -> array of registry indices that detected this edge
        edgeVariantMap      % containers.Map: uint64 edge key -> struct(registryIdx, variantIdx, fromNode) assigned once in buildTopology
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
            obj.edgeVariantMap = topo.edgeVariantMap;
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
            % Assigns twin generations, variant paths, and ideal orientations.
            %
            % ARCHITECTURE: Guess -> Compose -> Group -> Crown -> BFS Generations -> Assign Paths
            %
            % GUESS:   Largest grain in each component is the initial parent proxy.
            % COMPOSE: Traverse spanning tree using pre-assigned edgeVariantMap to compute
            %          exact composed misorientations from root (no orientation noise).
            %          For involution twins: mis_v6 * mis_v6 = identity -> T6,6 = parent.
            % GROUP:   Cluster grains with algebraically identical composed misorientations.
            % CROWN:   Entity with largest total grain area becomes the true parent set.
            % BFS:     distances() from true parents gives correct generation depths.
            % PATHS:   Build variant paths directly from edgeVariantMap (no refitting).
            
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
                activeGrainsInComp = nodesInComp(nodesInComp <= obj.virtualNodeBase);
                if isempty(activeGrainsInComp)
                    continue;
                end
                
                nComp = length(activeGrainsInComp);
                compGrainsInds = idToIndMap(activeGrainsInComp);
                
                % 1. GUESS: largest grain as initial parent proxy
                grainSizes_c = allSizes(compGrainsInds);
                [~, maxIdx_c] = max(grainSizes_c);
                guessParentGid = activeGrainsInComp(maxIdx_c);
                
                % Build local grain-ID -> component-position map
                localMap_c = zeros(max(activeGrainsInComp), 1);
                for li = 1:nComp
                    localMap_c(activeGrainsInComp(li)) = li;
                end
                
                % 2. COMPOSE: BFS through spanning tree, composing exact registry
                %    misorientations. Uses edgeVariantMap assigned in buildTopology.
                %    All compositions are exact math (no measurement noise).
                %    identity = rotation.id, so: T6,6 gets mis_v6 * mis_v6 = identity.
                composedMisoris_c = repmat(rotation.id, nComp, 1);
                
                T_c = shortestpathtree(obj.Q_G, guessParentGid);
                totalNodes = numnodes(obj.Q_G);
                visited_c = false(totalNodes, 1);
                visited_c(guessParentGid) = true;
                queue_c = guessParentGid;
                
                while ~isempty(queue_c)
                    u = queue_c(1); queue_c(1) = [];
                    childNodes_c = successors(T_c, u);
                    
                    for ni = 1:length(childNodes_c)
                        v = childNodes_c(ni);
                        if visited_c(v); continue; end
                        visited_c(v) = true;
                        queue_c(end+1) = v; %#ok<AGROW>
                        if v > obj.virtualNodeBase; continue; end
                        
                        u_li = 0;
                        if u <= max(activeGrainsInComp); u_li = localMap_c(u); end
                        v_li = localMap_c(v);
                        if u_li == 0 || v_li == 0; continue; end
                        
                        id1 = min(u, v); id2 = max(u, v);
                        key = bitshift(uint64(id2), 32) + uint64(id1);
                        
                        if isKey(obj.edgeVariantMap, key)
                            ev = obj.edgeVariantMap(key);
                            mis = obj.twinRegistry(ev.registryIdx).misoris(ev.variantIdx);
                            % If the variant was assigned in the opposite direction, use inverse
                            if ev.fromNode == v
                                mis = inv(mis);
                            end
                            composedMisoris_c(v_li) = composedMisoris_c(u_li) * mis;
                        end
                    end
                end
                
                % 3. GROUP: cluster by algebraic identity of composed misorientations
                entityBins = obj.groupKinematicEntities(activeGrainsInComp, composedMisoris_c);
                
                % 4. CROWN: entity with largest total grain area becomes true parents
                uniqueEntities = unique(entityBins);
                entityTotalArea = zeros(size(uniqueEntities));
                for ei = 1:length(uniqueEntities)
                    memberInds = compGrainsInds(entityBins == uniqueEntities(ei));
                    entityTotalArea(ei) = sum(allSizes(memberInds));
                end
                [~, crownIdx] = max(entityTotalArea);
                parentEntityBin = uniqueEntities(crownIdx);
                compParentGids = activeGrainsInComp(entityBins == parentEntityBin);
                
                % 5. BFS: compute generation depths from all crowned parents
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
                
                % Re-extract crowned parents (BFS may reassign based on distances)
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
                
                % 6. ASSIGN PATHS: BFS generation-by-generation from true parents.
                %    For each child, look up its edge to the parent in edgeVariantMap
                %    and build the variant path directly (no refitting needed).
                for gen = 1:maxGenComp
                    genGids = activeGrainsInComp(obj.grainGen(activeGrainsInComp) == gen);
                    
                    for i = 1:length(genGids)
                        child_gid = genGids(i);
                        child_ind = idToIndMap(child_gid);
                        
                        % Find a neighboring grain of lower generation (the parent)
                        childEdges = obj.physicalEdges(obj.physicalEdges(:,1) == child_gid | obj.physicalEdges(:,2) == child_gid, :);
                        
                        for e = 1:size(childEdges, 1)
                            id1 = childEdges(e,1); id2 = childEdges(e,2);
                            parent_gid = id1;
                            if parent_gid == child_gid; parent_gid = id2; end
                            
                            if ~isKey(obj.edgeVariantMap, key); continue; end
                            
                            ev = obj.edgeVariantMap(key);
                            r = ev.registryIdx;
                            m = ev.variantIdx;
                            mis = obj.twinRegistry(r).misoris(m);
                            % Determine topological step size
                            isCompound = obj.twinRegistry(r).genStep > 1;
                            expectedParentGen = gen - 1;
                            if isCompound
                                expectedParentGen = gen - 2;
                            end
                            
                            parentGen = obj.grainGen(parent_gid);
                            if ~isfinite(parentGen) || parentGen ~= expectedParentGen; continue; end
                            if expectedParentGen > 0 && isempty(obj.variantPath{parent_gid}); continue; end
                            
                            % If the stored direction is child->parent, use inverse variant
                            if ev.fromNode == child_gid
                                mis = inv(mis);
                                % Find closest matching variant index for the inverse
                                misoriArray = obj.twinRegistry(r).misoris;
                                fits = arrayfun(@(k) angle(misoriArray(k), mis), 1:length(misoriArray));
                                [~, m] = min(fits);
                            end
                            
                            vars = obj.twinRegistry(r).variantMap(m, :);
                            sysPath = obj.twinRegistry(r).sysPath;
                            
                            % Compute ideal child orientation via FZ-fixed parent
                            O_ideal_parent = obj.idealOri(parent_gid);
                            O_parent_fixed = allOris(idToIndMap(parent_gid)).project2FundamentalRegion(O_ideal_parent);
                            obj.idealOri(child_gid) = O_parent_fixed * mis;
                            obj.spawnVariant(child_ind) = vars(end);
                            obj.variantPath{child_gid} = [obj.variantPath{parent_gid}; [sysPath(:), vars(:)]];
                            break;
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
                    if g == 0
                        highlight(p, physNodes, 'MarkerSize', 10, 'NodeColor', [0.3 0.3 0.3]);
                    else
                        highlight(p, physNodes, 'MarkerSize', 6, 'NodeColor', cmapGen(g+1, :));
                    end
                end
                
                % Virtual nodes (faded by blending with white)
                virtNodes = find(nodeGen == g & (1:length(nodeGen))' > obj.virtualNodeBase);
                if ~isempty(virtNodes)
                    if g == 0
                        colorBlend = [0.3 0.3 0.3] * 0.5 + [1 1 1] * 0.5;
                    else
                        colorBlend = cmapGen(g+1, :) * 0.5 + [1 1 1] * 0.5;
                    end
                    highlight(p, virtNodes, 'NodeColor', colorBlend);
                end
            end
            
            title('Quotient Graph Topology (Node size = Parent, Colors = Generations, Faded = Virtual)');
        end
        
        function plotGenerations(obj)
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0.2); 
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            genProp = inf(max(obj.grains.id), 1);
            genProp(obj.activeGrains) = obj.grainGen(obj.activeGrains);
            maxGenPlot = max(genProp(genProp < inf));
            
            if isempty(maxGenPlot) || isnan(maxGenPlot)
                maxGenPlot = 0;
            end
            
            cmapGen = lines(maxGenPlot + 1);
            
            maxGid = max(obj.grains.id);
            idToIndMap = zeros(maxGid, 1);
            idToIndMap(obj.grains.id) = 1:length(obj.grains);
            
            for g = 0:maxGenPlot
                gids = find(genProp == g);
                if ~isempty(gids)
                    gInds = idToIndMap(gids);
                    if g == 0
                        hold on
                        plot(obj.grains(gInds), 'faceColor', [0.3 0.3 0.3], 'DisplayName', 'Parent (Gen 0)');
                    else
                        hold on
                        plot(obj.grains(gInds), 'faceColor', cmapGen(g+1, :), 'DisplayName', sprintf('Generation %d', g));
                    end
                end
            end
            
            title('Twin Generations');
            legend('show', 'Location', 'northeastoutside');
            hold off;
        end
        
        function plotTwinVariants(obj)
            % Plots twin variants with text labels
            figure;
            plot(obj.grains, obj.grains.meanOrientation, 'faceAlpha', 0);
            hold on;
            plot(obj.grains.boundary, 'lineColor', [0.5 0.5 0.5], 'lineWidth', 1);
            
            maxGid = max(obj.grains.id);
            idToIndMap = zeros(maxGid, 1);
            idToIndMap(obj.grains.id) = 1:length(obj.grains);
            
            labels = cell(length(obj.grains), 1);
            parentIndsPlot = [];
            twinIndsPlot = [];
            colorArray = zeros(length(obj.grains), 1);
            
            for i = 1:length(obj.activeGrains)
                gid = obj.activeGrains(i);
                gInd = idToIndMap(gid);
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
                hold on
                plot(obj.grains(parentIndsPlot), 'faceColor', [0.7 0.7 0.7], 'faceAlpha', 0.8, 'DisplayName', 'P');
            end
            
            if ~isempty(twinIndsPlot)
                % Unique final variants
                uVars = unique(colorArray(twinIndsPlot));
                cmap = lines(max(uVars));
                for v = uVars'
                    vInds = twinIndsPlot(colorArray(twinIndsPlot) == v);
                    hold on
                    plot(obj.grains(vInds), 'faceColor', cmap(v, :), 'faceAlpha', 0.8, 'DisplayName', sprintf('Variant %d', v));
                end
            end
            
            % Plot text labels (only non-empty ones are rendered)
            validLabels = ~cellfun(@isempty, labels);
            if any(validLabels)
                text(obj.grains(validLabels), labels(validLabels), 'FontSize', 8, 'BackgroundColor', 'w');
            end
            
            title('Twinning Variants (Colored by Final Spawned Variant)');
            legend('show', 'Location', 'northeastoutside');
            hold off;
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
        
        function entityBins = groupKinematicEntities(~, activeGrainsInComp, composedMisoris)
            % Clusters grains by algebraic identity of their composed rotations.
            % After composing exact registry misorientations from the root,
            % structurally equivalent grains (net-identity twin sequence, e.g. T6,6)
            % have mathematically identical composed rotations (angle = 0).
            % Epsilon is floating-point equality, NOT a twin detection threshold.
            
            nComp = length(activeGrainsInComp);
            entityBins = zeros(nComp, 1);
            nextBin = 1;
            epsilon = 1e-3 * degree; % floating-point equality for exact compositions
            
            for i = 1:nComp
                if entityBins(i) ~= 0; continue; end
                entityBins(i) = nextBin;
                for j = i+1:nComp
                    if entityBins(j) == 0 && angle(composedMisoris(i), composedMisoris(j)) < epsilon
                        entityBins(j) = nextBin;
                    end
                end
                nextBin = nextBin + 1;
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
            % Identifies boundaries and builds the Equivalence and Quotient Graphs.
            % Pre-computes exact structural assignments for edges and paths.
            
            topo = struct();
            validEdges = gB.grainId(:,1) > 0 & gB.grainId(:,2) > 0;
            
            % Store the updated registry
            topo.twinRegistry = twinRegistry;
            topo.physicalEdges = [];
            
            % Maps an unordered grain ID pair (key) to an array of valid registry indices
            topo.edgeToRegistryMap = containers.Map('KeyType', 'uint64', 'ValueType', 'any');
            
            activeGrains = [];
            
            % Step 1: Find valid boundary edges for each registry entry
            for r = 1:length(topo.twinRegistry)
                misoris = topo.twinRegistry(r).misoris;
                
                % Vectorized check against all ideal misorientations for this twin sequence
                isBoundary = any(angle(gB.misorientation, misoris) < settings.twinThresh, 2);
                isBoundary = isBoundary & validEdges;
                
                edges = gB.grainId(isBoundary, :);
                uniqueEdges = unique(sort(edges, 2), 'rows');
                
                topo.twinRegistry(r).edges = uniqueEdges;
                activeGrains = [activeGrains; uniqueEdges(:)];
                topo.physicalEdges = [topo.physicalEdges; uniqueEdges];
                
                % Tag unique physical boundaries using a bitshifted uint64 key
                % This makes lookup for any specific boundary instantly O(1).
                if ~isempty(uniqueEdges)
                    allKeys = bitshift(uint64(max(uniqueEdges(:,1), uniqueEdges(:,2))), 32) + uint64(min(uniqueEdges(:,1), uniqueEdges(:,2)));
                    for e = 1:size(uniqueEdges, 1)
                        key = allKeys(e);
                        % An edge might match multiple twin sequences (e.g. Primary and Double twin)
                        % so we store ALL valid registry indices for this boundary.
                        if isKey(topo.edgeToRegistryMap, key)
                            topo.edgeToRegistryMap(key) = [topo.edgeToRegistryMap(key), r];
                        else
                            topo.edgeToRegistryMap(key) = r;
                        end
                    end
                end
            end
            
            % Deduplicate global physical edges and extract all involved grains
            topo.physicalEdges = unique(topo.physicalEdges, 'rows');
            topo.activeGrains = unique(activeGrains);
            
            % Guard Clause: Exit cleanly if no twinning boundaries are found
            if isempty(topo.activeGrains)
                topo.Q_G = graph();
                topo.virtualNodeBase = 0;
                topo.compBins = [];
                topo.uniqueComps = [];
                return;
            end
            
            maxGid = max(grains.id);
            
            % Step 4: Build Quotient Graph (Topology)
            qPairs = [];
            
            % Virtual nodes start at maxGid + 1 so they don't collide with true MTEX grainIDs
            topo.virtualNodeBase = maxGid;
            numVirtualNodes = 0;
            
            for r = 1:length(topo.twinRegistry)
                genStep = topo.twinRegistry(r).genStep;
                edges = topo.twinRegistry(r).edges;
                if isempty(edges)
                    continue; 
                end
                
                if genStep == 1
                    % Primary twins (Generation step = 1): Add direct edge to graph
                    pairs = edges;
                    isSelfLoop = pairs(:,1) == pairs(:,2);
                    pairs(isSelfLoop, :) = [];
                    qPairs = [qPairs; pairs];
                else
                    % Double/Nested twins (Generation step > 1): Topologically separate them
                    % from direct boundaries by injecting a "virtual" ghost node in between.
                    % Example: Grain 1 ---[Virtual Node]--- Grain 5
                    % This ensures distance(1,5) = 2 in the graph, representing the 2-step variant sequence.
                    for i = 1:size(edges, 1)
                        numVirtualNodes = numVirtualNodes + 1;
                        vNode = topo.virtualNodeBase + numVirtualNodes;
                        u = edges(i, 1);
                        v = edges(i, 2);
                        qPairs = [qPairs; u, vNode; vNode, v];
                    end
                end
            end
            
            % Create MATLAB Graph
            if isempty(qPairs)
                topo.Q_G = graph();
            else
                topo.Q_G = graph(qPairs(:,1), qPairs(:,2));
            end
            
            % Ensure graph size is capable of fitting all assigned node IDs 
            % (handles case where high IDs are disjoint)
            maxNodeID = topo.virtualNodeBase + numVirtualNodes;
            if numnodes(topo.Q_G) < maxNodeID
                topo.Q_G = addnode(topo.Q_G, maxNodeID - numnodes(topo.Q_G));
            end
            
            % Step 5: Extract connected components from Graph
            topo.compBins = conncomp(topo.Q_G);
            topo.uniqueComps = unique(topo.compBins);
            
            % Step 6: Single-pass topological orientation fitting
            % For each connected component, find a spanning tree from the largest grain.
            % Traverse the tree and select the absolute mathematically best variant from the registry.
            % Stored in edgeVariantMap for perfect deterministic assignment in propagateReferences.
            
            % Pre-extract all data for fast vectorized lookups
            allOrisLocal = grains.meanOrientation;
            idToIndMapLocal = zeros(max(grains.id), 1);
            idToIndMapLocal(grains.id) = 1:length(grains);
            allSizesLocal = grains.grainSize;
            
            topo.edgeVariantMap = containers.Map('KeyType', 'uint64', 'ValueType', 'any');
            
            for c = 1:length(topo.uniqueComps)
                compID_c = topo.uniqueComps(c);
                nodesInComp_c = find(topo.compBins == compID_c);
                
                % Filter to true physical grains only
                physNodes_c = nodesInComp_c(nodesInComp_c <= topo.virtualNodeBase);
                if isempty(physNodes_c); continue; end
                
                % Guess initial parent: Largest physical grain in the component
                physGrainInds_c = idToIndMapLocal(physNodes_c);
                grainSizes_c = allSizesLocal(physGrainInds_c);
                [~, maxIdx_c] = max(grainSizes_c);
                guessParent_c = physNodes_c(maxIdx_c);
                
                % Local indexing map for this specific component
                maxPhysNode_c = max(physNodes_c);
                localMap_c = zeros(maxPhysNode_c, 1);
                for li_c = 1:length(physNodes_c)
                    localMap_c(physNodes_c(li_c)) = li_c;
                end
                
                % Setup for cascading ideal orientations
                idealOris_c = allOrisLocal(idToIndMapLocal(physNodes_c));
                
                % Build a directed Shortest Path Spanning Tree from the guess parent
                T_c = shortestpathtree(topo.Q_G, guessParent_c);
                visited_c = false(numnodes(topo.Q_G), 1);
                visited_c(guessParent_c) = true;
                queue_c = guessParent_c;
                
                % Map to remember which physical node spawned a given virtual node
                virtParent = zeros(numnodes(topo.Q_G), 1);
                
                % Breadth-First-Search traversal
                while ~isempty(queue_c)
                    u_c = queue_c(1); queue_c(1) = [];
                    childNodes_c = successors(T_c, u_c);
                    
                    for ni_c = 1:length(childNodes_c)
                        v_c = childNodes_c(ni_c);
                        if visited_c(v_c); continue; end
                        visited_c(v_c) = true;
                        queue_c(end+1) = v_c; %#ok<AGROW>
                        
                        % If reaching a virtual node, save where we came from and skip assignment.
                        % The edge assignment happens when we bridge from virtual -> physical.
                        if v_c > topo.virtualNodeBase
                            virtParent(v_c) = u_c;
                            continue;
                        end
                        
                        % Target v_c is a true physical grain. How did we get here?
                        if u_c > topo.virtualNodeBase
                            % Arrived via a virtual node (Compound double twin boundary)
                            phys_u = virtParent(u_c);
                            isDirect = false;
                        else
                            % Arrived directly from another physical grain (Primary twin boundary)
                            phys_u = u_c;
                            isDirect = true;
                        end
                        
                        % Lookup indices for orientation math
                        u_li_c = 0;
                        if phys_u <= maxPhysNode_c; u_li_c = localMap_c(phys_u); end
                        v_li_c = localMap_c(v_c);
                        if u_li_c == 0 || v_li_c == 0; continue; end
                        
                        % Retrieve all registry entries that hit this boundary
                        id1_c = min(phys_u, v_c); id2_c = max(phys_u, v_c);
                        key_c = bitshift(uint64(id2_c), 32) + uint64(id1_c);
                        if ~isKey(topo.edgeToRegistryMap, key_c); continue; end
                        regIndices_c = topo.edgeToRegistryMap(key_c);
                        
                        % TOPOLOGY ENFORCEMENT FILTER:
                        % A direct edge in the graph is strictly a 1-step twin.
                        % A virtual node bridge is strictly a >1-step nested twin.
                        validRegs = [];
                        for r_idx = 1:length(regIndices_c)
                            r = regIndices_c(r_idx);
                            if isDirect && topo.twinRegistry(r).genStep == 1
                                validRegs = [validRegs, r]; %#ok<AGROW>
                            elseif ~isDirect && topo.twinRegistry(r).genStep > 1
                                validRegs = [validRegs, r]; %#ok<AGROW>
                            end
                        end
                        
                        % If no valid topologies match, skip this edge
                        if isempty(validRegs); continue; end
                        
                        % Find the mathematically optimal variant out of the valid registries
                        O_u_ideal = idealOris_c(u_li_c);
                        O_v_actual = allOrisLocal(idToIndMapLocal(v_c));
                        
                        bestAngle_c = inf;
                        bestReg_c = validRegs(1);
                        bestVar_c = 1;
                        bestMis_c = topo.twinRegistry(bestReg_c).misoris(1);
                        
                        for k_c = 1:length(validRegs)
                            r_c = validRegs(k_c);
                            misArray_c = topo.twinRegistry(r_c).misoris;
                            
                            % Check every variant misorientation within this twin sequence
                            for m_c = 1:length(misArray_c)
                                % Compose ideal child orientation using parent's ideal frame
                                O_cand = O_u_ideal * misArray_c(m_c);
                                a_c = angle(O_cand, O_v_actual);
                                
                                % Keep the tightest fitting variant
                                if a_c < bestAngle_c
                                    bestAngle_c = a_c;
                                    bestReg_c = r_c;
                                    bestVar_c = m_c;
                                    bestMis_c = misArray_c(m_c);
                                end
                            end
                        end
                        
                        % Advance the perfect ideal orientation chain to the next node
                        idealOris_c(v_li_c) = O_u_ideal * bestMis_c;
                        
                        % Store the final assignment. 'fromNode' encodes graph direction
                        % so inverses can be computed if traversed backwards later.
                        topo.edgeVariantMap(key_c) = struct('registryIdx', bestReg_c, ...
                            'variantIdx', bestVar_c, 'fromNode', phys_u);
                    end
                end
            end
        end
    end
end
