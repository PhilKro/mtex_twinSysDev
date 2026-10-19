classdef twinNetwork < handle
    % TWINNETWORK Modernized class for identifying and analyzing twin networks in EBSD data.
    % This class implements a structured, pipeline-based approach to twin network reduction. 
    %
    % Pipeline Execution Example:
    % ---------------------------
    % The class utilizes a Stateful Handle paradigm allowing for method chaining.
    %
    %   network = twinNetwork(grains, gB, tS, ebsd) ...
    %       .buildInitialGraph() ...        % Builds raw topology graph (requires initialized object). State: 'GraphBuilt'
    %       .extractSpanningTree() ...      % Extracts MST from raw graph. State: 'TreeExtracted'
    %       .resolveNodeStates() ...        % Propagates FZ rules down the tree. State: 'StatesResolved'
    %       .pruneAndRewireOutliers() ...   % Optionally prune outliers and rewire. State: 'OutliersPruned'
    %       .reduceToQuotient();            % Reduces uniform subtrees to quotient nodes. State: 'QuotientReduced'
    %
    %   % Post-processing Example:
    %   [sf, iw] = network.analyzeMechanics(sigma); % Analyzes mechanics using quotient graph. State: 'MechanicsAnalyzed'
    %
    % Reverting Steps:
    % ----------------
    % Since intermediate graphs (Graph, mstGraph, QuotientGraph) are stored independently, 
    % you can revert pipeline steps by simply re-calling the method that overwrites the active state.
    % Example: To re-run outlier pruning with a different threshold, simply call `resolveNodeStates()` 
    % again (which regenerates a pristine `mstGraph` from the raw `Graph`), and then re-call pruning.
    
    properties
        grains              % MTEX grain object
        gB                  % MTEX grain boundary object
        tS                  % Twin symmetry (variant definitions)
        ebsd                % EBSD data
        twinThresh = 5*degree; % Threshold for considering a boundary as a twin
        maxGen = 2;         % Maximum generation of twins to consider
        
        Graph               % The primary raw topological digraph
        mstGraph            % The intermediate spanning tree graph
        QuotientGraph       % The final reduced quotient graph
        nodeGroups          % Mapping of original nodes to quotient groups
        pipelineState       % String indicating current pipeline status
        
        % Core properties for pipeline
        variants            % Array of twin variants
        numVariants         % Total number of variants
    end

    methods
        function obj = twinNetwork(grains, gB, tS, ebsd, varargin)
            % Constructor
            p = inputParser;
            addParameter(p, 'thresh', 5*degree);
            addParameter(p, 'twinThresh', []);
            addParameter(p, 'maxGen', 2);
            addParameter(p, 'maxGenerationLimit', []);
            parse(p, varargin{:});
            if ~isempty(p.Results.twinThresh)
                obj.twinThresh = p.Results.twinThresh;
            else
                obj.twinThresh = p.Results.thresh;
            end
            if ~isempty(p.Results.maxGenerationLimit)
                obj.maxGen = p.Results.maxGenerationLimit;
            else
                obj.maxGen = p.Results.maxGen;
            end

            obj.grains = grains;
            obj.gB = gB;
            obj.tS = tS;
            obj.ebsd = ebsd;

            % Precompute variants
            tS_sym = obj.tS.symmetrise;
            obj.variants = tS_sym.parentTwinMisorientation;
            if obj.variants.isImproper
                obj.variants = -obj.variants;
            end
            obj.numVariants = length(obj.variants);

            % Pipeline execution has been removed from constructor to allow step-by-step execution.
            obj.pipelineState = 'Uninitialized';
            disp('twinNetwork Initialized. Ready for step-by-step pipeline execution.');
        end
        
        function robust_ref = computeRobustCenter(obj, oris)
            if length(oris) <= 1
                robust_ref = oris;
                return;
            end
            
            scores = zeros(length(oris), 1);
            for i = 1:length(oris)
                close_idx = angle(oris(i), oris) < (5 * degree);
                scores(i) = sum(close_idx); % Majority vote!
            end
            
            [~, bestIdx] = max(scores);
            inlier_idx = angle(oris(bestIdx), oris) < (5 * degree);
            robust_ref = mean(oris(inlier_idx));
        end

        function obj = buildInitialGraph(obj, varargin)
            % Step A: Build pure topological graph with physical/virtual nodes
            % Pre-evaluate meanEbsdOri1/2 for edges. Variants are NaN.
            p = inputParser;
            addParameter(p, 'maxGen', obj.maxGen);
            parse(p, varargin{:});
            maxGen = p.Results.maxGen;
            
            G = digraph();
            
            % Add all physical grains
            numPhysNodes = max(obj.grains.id);
            if numPhysNodes == 0
                return;
            end
            G = addnode(G, numPhysNodes);
            
            % Store basic physical grain properties in Graph Nodes
            id2idx = zeros(numPhysNodes, 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);
            
            nodeTypes = repmat("Physical", numPhysNodes, 1);
            nodeGrains = (1:numPhysNodes)';
            nodeAreas = zeros(numPhysNodes, 1);
            nodeCentroids = zeros(numPhysNodes, 2);
            nodeMeanOris = orientation.nan(numPhysNodes, 1, obj.grains.CS);
            nodeComponentID = zeros(numPhysNodes, 1);
            nodeVariantPath = cell(numPhysNodes, 1);
            
            for i = 1:length(obj.grains)
                g_id = obj.grains(i).id;
                nodeAreas(g_id) = obj.grains(i).area;
                nodeCentroids(g_id, :) = [obj.grains(i).centroid.x, obj.grains(i).centroid.y];
                nodeMeanOris(g_id) = obj.grains(i).meanOrientation;
            end
            
            G.Nodes.Type = nodeTypes;
            G.Nodes.GrainId = nodeGrains;
            G.Nodes.Area = nodeAreas;
            G.Nodes.Centroid = nodeCentroids;
            G.Nodes.MeanOrientation = nodeMeanOris;
            G.Nodes.ComponentID = nodeComponentID;
            G.Nodes.VariantPath = nodeVariantPath;
            
            % Find valid boundary interfaces
            validMask = all(obj.gB.grainId ~= 0, 2) & all(obj.gB.ebsdId ~= 0, 2);
            gB_grainId = obj.gB.grainId(validMask, :);
            gB_ebsdId = obj.gB.ebsdId(validMask, :);
            
            % Extract raw EBSD data safely
            ind1 = obj.ebsd.id2ind(gB_ebsdId(:,1));
            ind2 = obj.ebsd.id2ind(gB_ebsdId(:,2));
            
            validInds = (ind1 > 0) & (ind2 > 0) & ~isnan(ind1) & ~isnan(ind2);
            if isprop(obj.ebsd, 'isIndexed')
                validInds = validInds & obj.ebsd.isIndexed(ind1) & obj.ebsd.isIndexed(ind2);
            end
            
            gB_grainId = gB_grainId(validInds, :);
            ind1 = ind1(validInds);
            ind2 = ind2(validInds);
            
            if isempty(gB_grainId)
                obj.Graph = G;
                obj.pipelineState = 'GraphBuilt';
                return;
            end
            
            ebsdOri1 = obj.ebsd(ind1).orientations;
            ebsdOri2 = obj.ebsd(ind2).orientations;
            
            % Group all boundary pixels by their unique physical grain pair
            [pairs, ~, dirTags] = unique(gB_grainId, 'rows');
            
            if ~all(pairs(:,1) < pairs(:,2))
                error('There is some flippage of grain ids in the gb segments going on')
            end
            
            if isempty(pairs)
                return;
            end
            
            % Average all EBSD pixels for each unique grain pair
            meanE1_cell = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(ebsdOri1(idx))});
            meanE2_cell = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(ebsdOri2(idx))});
            meanE1_array = [meanE1_cell{:}];
            meanE2_array = [meanE2_cell{:}];
            
            % Project each pair's mean orientation into the physical grain's fundamental region
            grainOri1 = obj.grains(obj.grains.id2ind(pairs(:,1))).meanOrientation;
            grainOri2 = obj.grains(obj.grains.id2ind(pairs(:,2))).meanOrientation;
            
            mE1_list = meanE1_array(:).project2FundamentalRegion(grainOri1(:));
            mE2_list = meanE2_array(:).project2FundamentalRegion(grainOri2(:));
            
            % Generate variant sequences up to max generation
            allSeqs = num2cell(1:obj.numVariants)';
            allMori = obj.variants(:);
            genMori = ones(obj.numVariants, 1);
            
            currentSeqs = allSeqs;
            currentMori = allMori;
            for g = 2:maxGen
                M = length(currentSeqs);
                idx_rep = repmat((1:M)', obj.numVariants, 1);
                idx_repel = repelem((1:obj.numVariants)', M, 1);
                
                nextSeqs = cellfun(@(x, y) [x, y], currentSeqs(idx_rep), num2cell(idx_repel), 'UniformOutput', false);
                repMori = currentMori(idx_rep);
                newMori = obj.variants(idx_repel);
                nextMori = repMori(:) .* newMori(:);
                
                allSeqs = [allSeqs; nextSeqs];
                allMori = [allMori; nextMori];
                genMori = [genMori; repmat(g, length(idx_rep), 1)];
                
                currentSeqs = nextSeqs;
                currentMori = nextMori;
            end
            
            % Calculate physical orientation angles
            % Evaluate: angle(physNode2, physNode1 * variant_sequence)
            angMatrix = zeros(size(pairs, 1), length(allMori));
            for p = 1:size(pairs, 1)
                % MTEX BUG FIX: `allMori(:)'` takes the INVERSE of the misorientations!
                % We must keep it as a column vector `allMori(:)` to avoid invoking inv().
                % The result of `angle` is a standard double column vector, which we then transpose with `.'`
                angMatrix(p, :) = angle(mE2_list(p), mE1_list(p) * allMori(:)).';
            end
            validMatches = angMatrix < obj.twinThresh;
            
            % Penalty to prefer shorter edges over longer ones (discourages virtual nodes)
            % Use repmat to ensure safe 2D addition for older MATLAB versions
            penaltyMatrix = repmat(genMori(:)' * 1000, size(angMatrix, 1), 1);
            scoreMatrix = angMatrix + penaltyMatrix;
            scoreMatrix(~validMatches) = inf;
            
            [bestScores, minIdxs] = min(scoreMatrix, [], 2);
            isTwinEdge = bestScores < inf;
            
            edgeSeqs = allSeqs(minIdxs(isTwinEdge));
            twinPairs = pairs(isTwinEdge, :);
            
            % Extract the precise angular deviations for the chosen variants
            linearIdxs = sub2ind(size(angMatrix), find(isTwinEdge), minIdxs(isTwinEdge));
            twinDevs = angMatrix(linearIdxs);
            twinMoriIdx = minIdxs(isTwinEdge);
            
            % Now we have valid physical twin pairs and their optimal sequences.
            % We will add Virtual Nodes to the graph where length(seq) > 1.
            s_nodes = [];
            t_nodes = [];
            e_seqs = {};
            e_devs = [];
            e_variant = [];
            e_meanE1 = orientation.empty(0, 1);
            e_meanE2 = orientation.empty(0, 1);
            
            vNodeIdx = numPhysNodes;
            debugging = false; % Set to true to enable verbose virtual node printouts
            
            for p_idx = 1:size(twinPairs, 1)
                global_pair_idx = find(isTwinEdge);
                global_pair_idx = global_pair_idx(p_idx); % Index into `pairs`
                
                physNode1 = twinPairs(p_idx, 1);
                physNode2 = twinPairs(p_idx, 2);
                
                % Simply grab the fully averaged EBSD orientation for this specific physical pair
                physE1 = mE1_list(global_pair_idx);
                physE2 = mE2_list(global_pair_idx);
                
                seq = edgeSeqs{p_idx};
                dev = twinDevs(p_idx);
                mori_idx = twinMoriIdx(p_idx);
                
                prevNode = physNode1;
                coord1 = nodeCentroids(physNode1, :);
                coord2 = nodeCentroids(physNode2, :);
                
                currTheoOri = nodeMeanOris(physNode1);
                endRealOri = nodeMeanOris(physNode2);
                
                if debugging && length(seq) > 1
                    fprintf('\n--- VIRTUAL NODE CHAIN CREATION ---\n');
                    fprintf('Connecting physNode1 (Grain %d) to physNode2 (Grain %d)\n', physNode1, physNode2);
                    fprintf('Variant sequence: [%s] (Index in allMori: %d)\n', num2str(seq), mori_idx);
                    
                    e1 = currTheoOri.Euler;
                    fprintf('physNode1 bulk MeanOri (rad): [%.4f, %.4f, %.4f]\n', e1(1), e1(2), e1(3));
                    e2 = endRealOri.Euler;
                    fprintf('physNode2 bulk MeanOri (rad): [%.4f, %.4f, %.4f]\n', e2(1), e2(2), e2(3));
                    
                    p1 = physE1.Euler;
                    fprintf('Edge EBSD physE1 (rad): [%.4f, %.4f, %.4f]\n', p1(1), p1(2), p1(3));
                    p2 = physE2.Euler;
                    fprintf('Edge EBSD physE2 (rad): [%.4f, %.4f, %.4f]\n', p2(1), p2(2), p2(3));
                    
                    fprintf('  -> physE1 deviation from bulk Grain %d: %.2f deg\n', physNode1, angle(physE1, currTheoOri)/degree);
                    fprintf('  -> physE2 deviation from bulk Grain %d: %.2f deg\n', physNode2, angle(physE2, endRealOri)/degree);
                    
                    % Print the actual verified angle from angMatrix
                    verified_ang = angle(physE2, physE1 * allMori(mori_idx))/degree;
                    fprintf('  -> Verified EBSD match angle (physE2 vs physE1*variants): %.2f deg\n', verified_ang);
                    
                    % Debug the MTEX vectorization mismatch
                    ang_matrix_val = twinDevs(p_idx) / degree;
                    fprintf('  -> CRITICAL DEBUG: angMatrix stored this angle as: %.2f deg\n', ang_matrix_val);
                    fprintf('  -> CRITICAL DEBUG: obj.twinThresh is set to: %.2f deg\n', obj.twinThresh / degree);
                end
                
                for step = 1:length(seq)
                    is_last_step = (step == length(seq));
                    
                    v_variant = seq(step);
                    
                    % Source theoretical orientation (before variant multiplication)
                    srcTheoOri = currTheoOri;
                    
                    currTheoOri = currTheoOri * obj.variants(v_variant);
                                    
                    if is_last_step
                        nxtNode = physNode2;
                    else
                        vNodeIdx = vNodeIdx + 1;
                        nxtNode = vNodeIdx;
                        
                        % Add virtual node to Graph
                        alpha = step / length(seq);
                        vCoord = coord1 + alpha * (coord2 - coord1);
                        
                        vProps = table("Virtual", NaN, 0, vCoord, currTheoOri, 0, {[]}, ...
                            'VariableNames', {'Type', 'GrainId', 'Area', 'Centroid', 'MeanOrientation', 'ComponentID', 'VariantPath'});
                        G = addnode(G, vProps);
                        
                        if debugging && length(seq) > 1
                            e_curr = currTheoOri.Euler;
                            fprintf('  Step %d (Variant %d): Virtual Node created with Euler (rad): [%.4f, %.4f, %.4f]\n', step, v_variant, e_curr(1), e_curr(2), e_curr(3));
                        end
                    end
                    
                    s_nodes(end+1, 1) = prevNode;
                    t_nodes(end+1, 1) = nxtNode;
                    e_seqs{end+1, 1} = num2str(seq(step));
                    
                    if is_last_step
                        e_devs(end+1, 1) = dev;
                    else
                        e_devs(end+1, 1) = 0; % Intermediate theoretical steps have 0 deviation
                    end
                    
                    e_variant(end+1, 1) = NaN; % Explicitly NaN as requested!
                    
                    % Only the FIRST sub-edge gets physE1, otherwise it gets the virtual node's orientation
                    if step == 1
                        e_meanE1 = [e_meanE1; physE1];
                    else
                        e_meanE1 = [e_meanE1; srcTheoOri];
                    end
                    
                    % Only the LAST sub-edge gets physE2, otherwise it gets the virtual node's orientation
                    if is_last_step
                        e_meanE2 = [e_meanE2; physE2];
                        
                        if debugging && length(seq) > 1
                            e_curr = currTheoOri.Euler;
                            fprintf('  Step %d (Variant %d): Final TheoOri Euler (rad): [%.4f, %.4f, %.4f]\n', step, v_variant, e_curr(1), e_curr(2), e_curr(3));
                            fprintf('  -> Final deviation to physNode2: %.2f deg\n', angle(endRealOri, currTheoOri)/degree);
                            angle(nodeMeanOris(physNode1), nodeMeanOris(physNode2)*allMori(mori_idx))/degree;
                            angle(nodeMeanOris(physNode2), nodeMeanOris(physNode1)*allMori(mori_idx))/degree;
                        end
                    else
                        e_meanE2 = [e_meanE2; currTheoOri];
                    end
                    
                    prevNode = nxtNode;
                end
            end
            
            if ~isempty(s_nodes)
                EdgeTable = table([s_nodes, t_nodes], e_seqs, e_devs, e_variant, e_meanE1, e_meanE2, ...
                    'VariableNames', {'EndNodes', 'Sequence', 'AngularDeviation', 'Variant', 'MeanEbsdOri1', 'MeanEbsdOri2'});
                G = addedge(G, EdgeTable);
            end
            
            obj.Graph = G;
            obj.pipelineState = 'GraphBuilt';
        end

        function obj = extractSpanningTree(obj)
            % Step B: Extract MST for each component
            if isempty(obj.Graph) || ~isa(obj.Graph, 'digraph')
                obj.mstGraph = graph();
                obj.pipelineState = 'TreeExtracted';
                return;
            end
            
            mstG = graph();
            mstG = addnode(mstG, obj.Graph.Nodes);
            
            EdgeTable = obj.Graph.Edges;
            if isempty(EdgeTable) || height(EdgeTable) == 0
                fprintf('No identified twin relationships!\n');
                obj.mstGraph = mstG;
                obj.pipelineState = 'TreeExtracted';
                return;
            end
            EdgeTable.Weight = EdgeTable.AngularDeviation;
            
            % Add penalty to virtual edges so spanning tree prefers physical chains
            isVirtualEdge = isnan(EdgeTable.Variant); % Variant is NaN for all edges right now!
            % Wait, we can identify virtual edges if they connect to a Virtual Node
            isVirtualS = obj.Graph.Nodes.Type(EdgeTable.EndNodes(:,1)) == "Virtual";
            isVirtualT = obj.Graph.Nodes.Type(EdgeTable.EndNodes(:,2)) == "Virtual";
            isVirtualEdge = isVirtualS | isVirtualT;
            EdgeTable.Weight(isVirtualEdge) = EdgeTable.Weight(isVirtualEdge) + 1000;
            
            G_undir = graph(EdgeTable, obj.Graph.Nodes);
            bins = conncomp(G_undir);
            numComps = max(bins);
            
            mstS = [];
            mstT = [];
            mstSeq = {};
            mstDev = [];
            mstVar = [];
            mstE1 = orientation.empty(0, 1);
            mstE2 = orientation.empty(0, 1);
            
            for c = 1:numComps
                compNodeIdx = find(bins == c);
                compNodeIdx = compNodeIdx(:); % Force column vector
                if length(compNodeIdx) == 1
                    continue;
                end
                
                [localS, localT, localSeq, localDev, localVar, localE1, localE2] = obj.extractComponentMST(G_undir, compNodeIdx);
                
                if ~isempty(localS)
                    mstS = [mstS; compNodeIdx(localS)];
                    mstT = [mstT; compNodeIdx(localT)];
                    mstSeq = [mstSeq; localSeq];
                    mstDev = [mstDev; localDev];
                    mstVar = [mstVar; localVar];
                    mstE1 = [mstE1; localE1];
                    mstE2 = [mstE2; localE2];
                end
            end
            
            if ~isempty(mstS)
                mstEdgeT = table([mstS, mstT], mstSeq, mstDev, mstVar, mstE1, mstE2, ...
                    'VariableNames', {'EndNodes', 'Sequence', 'AngularDeviation', 'Variant', 'MeanEbsdOri1', 'MeanEbsdOri2'});
                mstG = addedge(mstG, mstEdgeT);
            end
            
            obj.mstGraph = mstG;
            obj.pipelineState = 'TreeExtracted';
        end
        
        function [localS, localT, localSeq, localDev, localVar, localE1, localE2, localWeight] = extractComponentMST(obj, G_undir, compNodeIdx)
            subG = subgraph(G_undir, compNodeIdx);
            [subMST, ~] = minspantree(subG);
            
            if numedges(subMST) > 0
                localS = subMST.Edges.EndNodes(:, 1);
                localT = subMST.Edges.EndNodes(:, 2);
                localSeq = subMST.Edges.Sequence;
                localDev = subMST.Edges.AngularDeviation;
                localVar = subMST.Edges.Variant;
                localE1 = subMST.Edges.MeanEbsdOri1;
                localE2 = subMST.Edges.MeanEbsdOri2;
                localWeight = subMST.Edges.Weight;
            else
                localS = []; localT = []; localSeq = {}; localDev = [];
                localVar = []; localE1 = orientation.empty(0, 1); localE2 = orientation.empty(0, 1);
                localWeight = [];
            end
        end

        function obj = resolveNodeStates(obj)
            % Step C: Downward FZ propagation and variant locking
            if isempty(obj.mstGraph) || ~isa(obj.mstGraph, 'graph')
                obj.mstGraph = graph();
                obj.pipelineState = 'StatesResolved';
                return;
            end
            
            mstG = obj.mstGraph;
            numNodes = numnodes(mstG);
            if numNodes == 0
                obj.pipelineState = 'StatesResolved';
                return;
            end
            
            mstG.Nodes.Depth = zeros(numNodes, 1);
            mstG.Nodes.ParentNode = zeros(numNodes, 1);
            
            % Find components to locate roots
            bins = conncomp(mstG);
            numComps = max(bins);
            
            for c = 1:numComps
                compNodeIdx = find(bins == c);
                if length(compNodeIdx) == 1
                    mstG.Nodes.ComponentID(compNodeIdx(1)) = c;
                    mstG.Nodes.VariantPath{compNodeIdx(1)} = [];
                    continue;
                end
                
                % Find root (largest physical grain)
                physNodes = compNodeIdx(mstG.Nodes.Type(compNodeIdx) == "Physical");
                if ~isempty(physNodes)
                    [maxArea, maxAreaIdx] = max(mstG.Nodes.Area(physNodes));
                    rootNode = physNodes(maxAreaIdx);
                    fprintf('  [DEBUG-ROOT] Component %d: Chose root Grain %d with area %.2f.\n', c, rootNode, maxArea);
                    % fprintf('               Physical nodes in this component: %s\n', mat2str(physNodes'));
                else
                    rootNode = compNodeIdx(1);
                    fprintf('  [DEBUG-ROOT] Component %d: No physical nodes. Chose root node %d.\n', c, rootNode);
                end
                
                mstG.Nodes.ComponentID(rootNode) = c;
                mstG.Nodes.VariantPath{rootNode} = [];
                
                % BFS Traversal
                visited = false(numNodes, 1);
                visited(rootNode) = true;
                queue = rootNode;
                
                % State array to hold numerical sequence for each node
                nodeStateArr = cell(numNodes, 1);
                nodeStateArr{rootNode} = [];
                
                while ~isempty(queue)
                    currNode = queue(1);
                    queue(1) = [];
                    currState = nodeStateArr{currNode};
                    currOri = mstG.Nodes.MeanOrientation(currNode);
                    
                    % Find neighbors (undirected edges in mstG)
                    allEdges = outedges(mstG, currNode);
                    
                    for e = 1:length(allEdges)
                        edgeIdx = allEdges(e);
                        s = mstG.Edges.EndNodes(edgeIdx, 1);
                        t = mstG.Edges.EndNodes(edgeIdx, 2);
                        
                        if s == currNode
                            nxtNode = t;
                            isForward = true;
                        else
                            nxtNode = s;
                            isForward = false;
                        end
                        
                        if ~visited(nxtNode)
                            visited(nxtNode) = true;
                            
                            mstG.Nodes.Depth(nxtNode) = mstG.Nodes.Depth(currNode) + 1;
                            mstG.Nodes.ParentNode(nxtNode) = currNode;
                            
                            % Get base sequence ID from Step A
                            base_v = str2double(mstG.Edges.Sequence{edgeIdx});
                            
                            % Get the actual orientation of the child (measured for Physical, theoretical for Virtual)
                            nxtOriOld = mstG.Nodes.MeanOrientation(nxtNode);
                            
                            % Determine the true variant that minimizes deviation
                            % Note: obj.variants ALWAYS maps from EndNodes(:,1) to EndNodes(:,2)
                            if isForward
                                % currNode is s, nxtNode is t: ori_t = ori_s * v
                                devs = angle(nxtOriOld, currOri * obj.variants(:));
                            else
                                % currNode is t, nxtNode is s: ori_s = ori_t * inv(v)
                                devs = angle(nxtOriOld, currOri * inv(obj.variants(:)));
                            end
                            
                            [min_dev, actual_v] = min(devs);
                            
                            if isForward
                                nxtCenter = currOri * obj.variants(actual_v);
                            else
                                nxtCenter = currOri * inv(obj.variants(actual_v));
                            end
                            
                            % Project orientation into the FZ closest to nxtCenter!
                            nxtOriNew = nxtOriOld.project2FundamentalRegion(nxtCenter);
                            mstG.Nodes.MeanOrientation(nxtNode) = nxtOriNew;
                            
                            mstG.Edges.Variant(edgeIdx) = actual_v;
                            mstG.Edges.AngularDeviation(edgeIdx) = min_dev;
                            % FZ Propagate ALL EBSD boundary data associated with nxtNode!
                            allNxtEdges = outedges(mstG, nxtNode);
                            for e_nxt = 1:length(allNxtEdges)
                                edgeIdx_nxt = allNxtEdges(e_nxt);
                                if mstG.Edges.EndNodes(edgeIdx_nxt, 1) == nxtNode
                                    if ~isnan(mstG.Edges.MeanEbsdOri1(edgeIdx_nxt))
                                        mstG.Edges.MeanEbsdOri1(edgeIdx_nxt) = mstG.Edges.MeanEbsdOri1(edgeIdx_nxt).project2FundamentalRegion(nxtCenter);
                                    end
                                end
                                if mstG.Edges.EndNodes(edgeIdx_nxt, 2) == nxtNode
                                    if ~isnan(mstG.Edges.MeanEbsdOri2(edgeIdx_nxt))
                                        mstG.Edges.MeanEbsdOri2(edgeIdx_nxt) = mstG.Edges.MeanEbsdOri2(edgeIdx_nxt).project2FundamentalRegion(nxtCenter);
                                    end
                                end
                            end
                            
                            % Update algebraic string state using the new variant
                            nxtState = currState;
                            if ~isempty(nxtState) && nxtState(end) == actual_v
                                nxtState(end) = []; % Self-annihilation (1_2_2 -> 1)
                            else
                                nxtState(end+1) = actual_v;
                            end
                            
                            nodeStateArr{nxtNode} = nxtState;
                            mstG.Nodes.ComponentID(nxtNode) = c;
                            mstG.Nodes.VariantPath{nxtNode} = nxtState;
                            
                            queue(end+1) = nxtNode;
                        end
                    end
                end
            end
            
            obj.mstGraph = mstG;
            obj.pipelineState = 'StatesResolved';
        end

        obj = pruneAndRewireOutliers(obj, varargin)
        function obj = reduceToQuotient(obj)
            % REDUCETOQUOTIENT Contracts physical and virtual nodes of identical twin generations into a quotient graph.
            %
            % Description:
            %   Groups nodes in the Spanning Tree (mstGraph) based on their ComponentID and VariantPath.
            %   Nodes with identical states are collapsed into a single quotient node representing that twin generation.
            %   Calculates state-level properties (robust mean orientation, centroid, area) from the constituent grains.
            %   Constructs undirected edges linking parent and child twin generations.
            %   
            % Output:
            %   Updates obj.QuotientGraph with the contracted undirected graph.
            %   Updates obj.nodeGroups with the mapping from physical/virtual nodes to quotient nodes.
            if isempty(obj.mstGraph) || ~isa(obj.mstGraph, 'graph') || numnodes(obj.mstGraph) == 0
                obj.QuotientGraph = graph();
                obj.nodeGroups = [];
                obj.pipelineState = 'QuotientReduced';
                return;
            end
            mstG = obj.mstGraph;
            Q = graph();
            
            % Convert VariantPath to string just for grouping purposes
            % First, normalize all empty paths to 0x0 double so mat2str('[]') is consistent
            for i = 1:numnodes(mstG)
                if isempty(mstG.Nodes.VariantPath{i})
                    mstG.Nodes.VariantPath{i} = [];
                end
            end
            
            pathStrs = cellfun(@mat2str, mstG.Nodes.VariantPath, 'UniformOutput', false);
            [unqStates, ~, groupIDs] = unique(table(mstG.Nodes.ComponentID, pathStrs, 'VariableNames', {'ComponentID', 'pathStrs'}));
            numGroups = height(unqStates);
            
            Q = addnode(Q, numGroups);
            Q.Nodes.ComponentID = unqStates.ComponentID;
            Q.Nodes.VariantPath = cellfun(@eval, unqStates.pathStrs, 'UniformOutput', false);
            
            % Compute mean orientation for each quotient node
            nodeMeanOris = orientation.nan(numGroups, 1, obj.grains.CS);
            nodeCentroids = zeros(numGroups, 2);
            nodeAreas = zeros(numGroups, 1);
            
            for g = 1:numGroups
                g_nodes = find(groupIDs == g);
                
                oris = mstG.Nodes.MeanOrientation(g_nodes);
                areas = mstG.Nodes.Area(g_nodes);
                
                % Compute Robust Center (both physical and virtual orientations)
                if length(oris) <= 1
                    nodeMeanOris(g) = oris;
                else
                    nodeMeanOris(g) = obj.computeRobustCenter(oris);
                end
                
                total_area = sum(areas);
                
                % Compute Centroid (Virtual nodes have Area=0, so they inherently contribute 0 to the weighted sum)
                if total_area > 0
                    weighted_x = sum(mstG.Nodes.Centroid(g_nodes, 1) .* areas) / total_area;
                    weighted_y = sum(mstG.Nodes.Centroid(g_nodes, 2) .* areas) / total_area;
                    nodeCentroids(g, :) = [weighted_x, weighted_y];
                else
                    % Fallback for purely virtual groups (Area = 0)
                    nodeCentroids(g, :) = mean(mstG.Nodes.Centroid(g_nodes, :), 1);
                end
                
                nodeAreas(g) = total_area;
            end
            
            Q.Nodes.robustMeanOrientation = nodeMeanOris;
            Q.Nodes.Centroid = nodeCentroids;
            Q.Nodes.Area = nodeAreas;
            
            % Create Edges
            s_nodes = [];
            t_nodes = [];
            e_variants = [];
            
            for e = 1:numedges(mstG)
                u = mstG.Edges.EndNodes(e, 1);
                v = mstG.Edges.EndNodes(e, 2);
                
                group_u = groupIDs(u);
                group_v = groupIDs(v);
                
                if group_u ~= group_v
                    s_nodes(end+1, 1) = group_u;
                    t_nodes(end+1, 1) = group_v;
                    e_variants(end+1, 1) = mstG.Edges.Variant(e);
                end
            end
            
            if ~isempty(s_nodes)
                % Unique quotient edges
                edgeMat = [s_nodes, t_nodes, e_variants];
                [unqEdges, ~, ~] = unique(edgeMat, 'rows');
                
                EdgeT = table(unqEdges(:, 1:2), unqEdges(:, 3), 'VariableNames', {'EndNodes', 'Variant'});
                Q = addedge(Q, EdgeT);
            end
            
            obj.QuotientGraph = Q;
            obj.pipelineState = 'QuotientReduced';
            obj.nodeGroups = groupIDs;
        end
        

        
        function str = formatState(~, comp, path)
            if isempty(path)
                str = sprintf('C%d_ROOT', comp);
            else
                path_str = sprintf('%d_', path);
                str = sprintf('C%d_%s', comp, path_str(1:end-1));
            end
        end
        
        plot(obj, varargin)
        plotSigmaSections(obj, varargin)
        plotGraph(obj, G, varargin)     
        
        function [grain_SF, grain_IW] = analyzeMechanics(obj, sigma, varargin)
            % ANALYZEMECHANICS Evaluates Schmid Factor and Interaction Work for twin generations
            %
            % Syntax:
            %   [grain_SF, grain_IW] = obj.analyzeMechanics(sigma)
            %   [grain_SF, grain_IW] = obj.analyzeMechanics(sigma, 'relative')
            %
            % Input:
            %   sigma - stressTensor (in specimen coordinates)
            %   'relative' - flag to compute SF relative to CRSS
            %
            % Output:
            %   grain_SF - (n x 1) array mapping the Schmid factor directly to obj.grains
            %   grain_IW - (n x 1) array mapping the Interaction Work directly to obj.grains
            
            Q = obj.QuotientGraph;
            if isempty(Q) || numnodes(Q) == 0
                error('twinNetwork:QuotientGraphMissing', 'Quotient Graph is empty. Please run step 5 (reduceToQuotient) first.');
            end
            
            % Initialize properties on the Quotient nodes if not already present
            numQ = numnodes(Q);
            qSF = nan(numQ, 1);
            qIW = nan(numQ, 1);
            
            % Symmetrize the twin system to access individual variants
            tS_sym = obj.tS.symmetrise;
            
            % Perform a BFS tree traversal starting from the roots
            % In an undirected quotient graph, roots have an empty VariantPath
            isRoot = cellfun(@isempty, Q.Nodes.VariantPath);
            queue = find(isRoot);
            queue = queue(:)'; % Row vector for queue
            
            visited = false(numQ, 1);
            visited(queue) = true;
            
            while ~isempty(queue)
                curr_node = queue(1);
                queue(1) = [];
                
                currOri = Q.Nodes.robustMeanOrientation(curr_node);
                
                % Get incident edges for the current node (outedges returns all incident edges for undirected graphs)
                edge_idxs = outedges(Q, curr_node);
                
                for i = 1:length(edge_idxs)
                    edge_idx = edge_idxs(i);
                    s = Q.Edges.EndNodes(edge_idx, 1);
                    t = Q.Edges.EndNodes(edge_idx, 2);
                    
                    % Identify which node is the child
                    if s == curr_node
                        child_node = t;
                    else
                        child_node = s;
                    end
                    
                    if ~visited(child_node)
                        visited(child_node) = true;
                        
                        variant_idx = Q.Edges.Variant(edge_idx);
                        
                        % We only calculate if the parent has a valid orientation
                        if ~isnan(currOri) && variant_idx > 0 && variant_idx <= obj.numVariants
                            % Calculate Schmid Factor
                            sf_val = SchmidFactor(tS_sym(variant_idx), sigma, currOri, varargin{:});
                            % Calculate Interaction Work
                            iw_val = interactionWork(tS_sym(variant_idx), sigma, currOri, varargin{:});
                            
                            % Store in child node (the twin generation)
                            qSF(child_node) = max(sf_val);
                            qIW(child_node) = max(iw_val);
                        end
                        
                        queue(end+1) = child_node;
                    end
                end
            end
            
            % Attach the calculated values permanently to the QuotientGraph Nodes table
            obj.QuotientGraph.Nodes.SchmidFactor = qSF;
            obj.QuotientGraph.Nodes.InteractionWork = qIW;
            
            % Now map these values back to the individual physical grains
            numGrains = length(obj.grains);
            grain_SF = nan(numGrains, 1);
            grain_IW = nan(numGrains, 1);
            
            % Create an inverse mapping from physical ID to grain linear index (1 to numGrains)
            id2idx = zeros(max(obj.grains.id), 1);
            id2idx(obj.grains.id) = 1:numGrains;
            
            % For each mapped physical node, look up the quotient node's SF/IW
            max_phys_id = max(obj.grains.id);
            for phys_id = 1:min(length(obj.nodeGroups), max_phys_id)
                idx = id2idx(phys_id);
                if idx > 0
                    q_group = obj.nodeGroups(phys_id);
                    if q_group > 0
                        grain_SF(idx) = qSF(q_group);
                        grain_IW(idx) = qIW(q_group);
                    end
                end
            end
            
            % Also store on the raw Graph for consistency (aligning with its nodes)
            gSF = nan(numnodes(obj.Graph), 1);
            gIW = nan(numnodes(obj.Graph), 1);
            % obj.Graph nodes correspond directly to physical grains by GrainId
            if ismember('GrainId', obj.Graph.Nodes.Properties.VariableNames)
                for i = 1:numnodes(obj.Graph)
                    gid = obj.Graph.Nodes.GrainId(i);
                    if ~isnan(gid) && gid > 0 && gid <= length(obj.nodeGroups)
                        q_grp = obj.nodeGroups(gid);
                        if q_grp > 0
                            gSF(i) = qSF(q_grp);
                            gIW(i) = qIW(q_grp);
                        end
                    end
                end
            end
            obj.Graph.Nodes.SchmidFactor = gSF;
            obj.Graph.Nodes.InteractionWork = gIW;
            
            obj.pipelineState = 'MechanicsAnalyzed';
        end
    end
end
