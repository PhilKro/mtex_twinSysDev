classdef TwinNetwork < handle
    % TwinNetwork Modernized class for identifying and analyzing twin networks in EBSD data.
    % This class implements a structured, pipeline-based approach to twin network reduction.

    properties
        grains              % MTEX grain object
        gB                  % MTEX grain boundary object
        tS                  % Twin symmetry (variant definitions)
        ebsd                % EBSD data
        twinThresh = 5*degree; % Threshold for considering a boundary as a twin

        Graph               % The primary raw topological digraph
        QuotientGraph       % The final reduced quotient graph
        nodeGroups          % Mapping of original nodes to quotient groups

        % Core properties for pipeline
        variants            % Array of twin variants
        numVariants         % Total number of variants
    end

    methods
        function obj = TwinNetwork(grains, gB, tS, ebsd, varargin)
            % TwinNetwork Constructor for the TwinNetwork class.
            %
            % Parameters:
            %   grains   - MTEX grain object containing the physical grains
            %   gB       - MTEX grain boundary object
            %   tS       - MTEX twin symmetry (variant definitions)
            %   ebsd     - EBSD data underlying the grains
            %   varargin - Optional parameters (e.g., 'thresh', default 5*degree)
            %
            % Returns:
            %   obj - Initialized TwinNetwork object
            if nargin == 0, return; end

            p = inputParser;
            addParameter(p, 'thresh', 5*degree);
            parse(p, varargin{:});
            obj.twinThresh = p.Results.thresh;

            obj.grains = grains;
            obj.gB = gB;
            obj.tS = tS;
            obj.ebsd = ebsd;

            % Precompute variants
            tS_sym = symmetrise(obj.tS);
            obj.variants = parentTwinMisorientation(tS_sym);
            improper = isImproper(obj.variants);
            if any(improper(:))
                obj.variants(improper) = -obj.variants(improper);
            end
            obj.numVariants = length(obj.variants);

            % Pipeline execution has been removed from constructor to allow step-by-step execution.
            disp('TwinNetwork Initialized. Ready for step-by-step pipeline execution.');
        end

        function robust_ref = computeRobustCenter(obj, oris)
            % COMPUTEROBUSTCENTER Computes a robust central orientation from a list using majority voting.
            % Identifies the central orientation by scoring each orientation based on how many
            % other orientations are within a 5-degree tolerance, explicitly rejecting area-weighted means.
            %
            % Parameters:
            %   oris - Array of MTEX orientations
            %
            % Returns:
            %   robust_ref - The robust mean orientation
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

        function G = buildInitialGraph(obj, varargin)

            % Workflow Step 1: Build pure topological graph with physical/virtual nodes
            % Pre-evaluate meanEbsdOri1/2 for edges. Variants are NaN.
            %
            % Parameters:
            %   varargin - Optional parameters:
            %              'maxGenStep' - Maximum twin generation to consider for virtual nodes (default: 2)

            % BUILDINITIALGRAPH Constructs the initial topological graph of physical and virtual nodes.
            % "Physical" nodes represent grains and "virtual" nodes represent possibly consumed parent grains if maxGenStep>1
            % (i.e. secondary/tertiary/.. twins are recognized even if they consumed the parent grain).
            %
            % Edges represent grain boundaries with weights corresponding to
            % angular deviations from theoretical twin relationships.
            %
            % Returns:
            %   G - The initialized topological digraph

            p = inputParser;
            addParameter(p, 'maxGenStep', 1, @isnumeric);
            parse(p, varargin{:});
            maxGenStep = p.Results.maxGenStep;

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
            nodeStates = repmat("[]", numPhysNodes, 1);

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
            G.Nodes.State = nodeStates;

            % Find valid boundary interfaces
            validMask = all(obj.gB.grainId ~= 0, 2);
            gB_grainId = obj.gB.grainId(validMask, :);
            gB_ebsdId = obj.gB.ebsdId(validMask, :);

            % Extract raw EBSD data
            ebsdOri1 = obj.ebsd(obj.ebsd.id2ind(gB_ebsdId(:,1))).orientations;
            ebsdOri2 = obj.ebsd(obj.ebsd.id2ind(gB_ebsdId(:,2))).orientations;

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
            if isscalar(meanE1_cell)
                meanE1_array = meanE1_cell{1};
                meanE2_array = meanE2_cell{1};
            else
                meanE1_array = [meanE1_cell{:}];
                meanE2_array = [meanE2_cell{:}];
            end

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
            for g = 2:maxGenStep
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
            e_meanE1 = {};
            e_meanE2 = {};

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

                        vProps = table("Virtual", NaN, 0, vCoord, currTheoOri, "[]", ...
                            'VariableNames', {'Type', 'GrainId', 'Area', 'Centroid', 'MeanOrientation', 'State'});
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
                        e_meanE1{end+1, 1} = physE1;
                    else
                        e_meanE1{end+1, 1} = srcTheoOri;
                    end

                    % Only the LAST sub-edge gets physE2, otherwise it gets the virtual node's orientation
                    if is_last_step
                        e_meanE2{end+1, 1} = physE2;

                        if debugging && length(seq) > 1
                            e_curr = currTheoOri.Euler;
                            fprintf('  Step %d (Variant %d): Final TheoOri Euler (rad): [%.4f, %.4f, %.4f]\n', step, v_variant, e_curr(1), e_curr(2), e_curr(3));
                            fprintf('  -> Final deviation to physNode2: %.2f deg\n', angle(endRealOri, currTheoOri)/degree);
                            angle(nodeMeanOris(physNode1), nodeMeanOris(physNode2)*allMori(mori_idx))/degree;
                            angle(nodeMeanOris(physNode2), nodeMeanOris(physNode1)*allMori(mori_idx))/degree;
                        end
                    else
                        e_meanE2{end+1, 1} = currTheoOri;
                    end

                    prevNode = nxtNode;
                end
            end

            if ~isempty(s_nodes)
                e_meanE1 = catOrientationCells(e_meanE1);
                e_meanE2 = catOrientationCells(e_meanE2);
                EdgeTable = table([s_nodes, t_nodes], e_seqs, e_devs, e_variant, e_meanE1, e_meanE2, ...
                    'VariableNames', {'EndNodes', 'Sequence', 'AngularDeviation', 'Variant', 'MeanEbsdOri1', 'MeanEbsdOri2'});
                G = digraph(EdgeTable,G.Nodes);
            end

            obj.Graph = G;
        end

        function mstG = extractSpanningTree(obj, G)
            % Workflow Step 2: Extract MST for each component
            %
            % EXTRACTSPANNINGTREE Extracts the Minimum Spanning Tree (MST) from the component graph.
            % Decomposes the initial graph into connected components (nodes that are connected by twin boundary edges) and finds the MST for each,
            % minimizing the overall angular deviation from twin relationships (edge weights).
            %
            % Parameters:
            %   G - The initial topological graph
            %
            % Returns:
            %   mstG - An undirected graph containing only the MST edges
            EdgeTable = G.Edges;
            if isempty(EdgeTable)
                emptyEdges = table(zeros(0,2),'VariableNames',{'EndNodes'});
                mstG = graph(emptyEdges,G.Nodes);
                fprintf('No identified twin relationships!');
                return
            end
            EdgeTable.Weight = EdgeTable.AngularDeviation;

            % Add penalty to virtual edges so spanning tree prefers physical chains
            isVirtualEdge = isnan(EdgeTable.Variant); % Variant is NaN for all edges right now!
            % Wait, we can identify virtual edges if they connect to a Virtual Node
            isVirtualS = G.Nodes.Type(EdgeTable.EndNodes(:,1)) == "Virtual";
            isVirtualT = G.Nodes.Type(EdgeTable.EndNodes(:,2)) == "Virtual";
            isVirtualEdge = isVirtualS | isVirtualT;
            EdgeTable.Weight(isVirtualEdge) = EdgeTable.Weight(isVirtualEdge) + 1000;

            G_undir = graph(EdgeTable, G.Nodes);
            bins = conncomp(G_undir);
            numComps = max(bins);

            mstS = [];
            mstT = [];
            mstSeq = {};
            mstDev = [];
            mstVar = [];
            mstE1 = {};
            mstE2 = {};

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
                    mstE1{end+1, 1} = localE1;
                    mstE2{end+1, 1} = localE2;
                end
            end

            if ~isempty(mstS)
                mstE1 = catOrientationCells(mstE1);
                mstE2 = catOrientationCells(mstE2);
                mstEdgeT = table([mstS, mstT], mstSeq, mstDev, mstVar, mstE1, mstE2, ...
                    'VariableNames', {'EndNodes', 'Sequence', 'AngularDeviation', 'Variant', 'MeanEbsdOri1', 'MeanEbsdOri2'});
                mstG = graph(mstEdgeT,G.Nodes);
            else
                mstG = graph();
                mstG = addnode(mstG, G.Nodes);
            end
        end

        function [localS, localT, localSeq, localDev, localVar, localE1, localE2, localWeight] = extractComponentMST(obj, G_undir, compNodeIdx)
            % EXTRACTCOMPONENTMST Helper function to compute the MST for a specific graph component.
            %
            % Parameters:
            %   G_undir - The full undirected graph
            %   compNodeIdx - Array of node indices belonging to the component
            %
            % Returns:
            %   Arrays describing the edges (source, target, deviation, variants, etc.) of the component MST.
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

        function mstG = resolveNodeStates(obj, mstG)
            % Workflow Step 3: Downward FZ propagation and variant locking

            % RESOLVENODESTATES Propagates variant states from component roots down the MSTs.
            % Uses Breadth-First Search (BFS) to traverse each component from its root, assigning
            % variant paths (e.g., 'C1_ROOT', 'C1_3_5'). For the consistency of the twin variants the
            % orientation of the next depth's nodes is projected into the fundamental zone of the determined twin variants at the present depth
            %
            % Parameters:
            %   mstG - Minimum Spanning Tree graph
            %
            % Returns:
            %   mstG - Updated graph with 'State', 'Depth', and 'ParentNode' properties


            numNodes = numnodes(mstG);
            if numNodes == 0
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
                    mstG.Nodes.State(compNodeIdx(1)) = sprintf('C%d_ROOT', c);
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

                mstG.Nodes.State(rootNode) = sprintf('C%d_ROOT', c);

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

                            if isempty(nxtState)
                                mstG.Nodes.State(nxtNode) = sprintf('C%d_ROOT', c);
                            else
                                str_s = sprintf('%d_', nxtState);
                                mstG.Nodes.State(nxtNode) = sprintf('C%d_%s', c, str_s(1:end-1));
                            end

                            queue(end+1) = nxtNode;
                        end
                    end
                end
            end
        end

        function mstG = pruneAndRewireOutliers(obj, mstG, varargin)

            % Workflow Step 4 (optional): Detect and rewire outliers and isolated nodes based on empirical FZ centers
            %
            % PRUNEANDREWIREOUTLIERS Detects and rewires nodes that deviate from theoretical orientations.
            % Uses either 'empirical' (robust center of nodes with a common variant path) or
            % 'theoretical' (theoretical twin orientations calculated from roots of each component) modes to identify outliers.

            % Attempts to rewire outliers by substituting alternative edges from the graph into the component MST. The MST is recalculated using resolveNodeStates.
            %
            % Parameters:
            %   mstG             - Minimum Spanning Tree graph
            %   varargin         - 'outlierMode'        - 'empirical' or 'theoretical' (default: 'empirical')
            %                      'maxErrThreshold'    - maximum angular deviation threshold to be recoginzed as outlier (default: 15 * degree)
            %                      'debug'              - enable debug mode (default: false)
            %                      'rejectGenIncrease'  - reject rewiring attempts that would increase the generation number of the outlier node (default: false)
            %
            % Returns:
            %   mstG - Updated MST graph with outliers rewired or removed
            p = inputParser;
            addParameter(p, 'targetNode', 0);
            addParameter(p, 'debug', false);
            addParameter(p, 'outlierMode', 'empirical', @(x) any(validatestring(x, {'empirical', 'theoretical'})));
            addParameter(p, 'maxErrThreshold', 15 * degree);
            addParameter(p, 'rejectGenIncrease', true, @islogical);
            parse(p, varargin{:});
            target_node_override = p.Results.targetNode;
            debug_mode = p.Results.debug;
            outlierMode = p.Results.outlierMode;
            max_err_threshold = p.Results.maxErrThreshold;
            reject_gen_increase = p.Results.rejectGenIncrease;

            numNodes = numnodes(mstG);
            if numNodes == 0
                return;
            end

            % The physical topology is inherently undirected. We must search an undirected graph!
            search_G_dir = obj.Graph;
            search_G = graph(search_G_dir.Edges,search_G_dir.Nodes);

            if ismember('AngularDeviation', search_G.Edges.Properties.VariableNames)
                search_G.Edges.Weight = search_G.Edges.AngularDeviation;
            else
                search_G.Edges.Weight = ones(numedges(search_G), 1);
            end

            loop_changed = true;
            max_iters = 100;
            iter = 0;
            unfixable_nodes = [];

            while loop_changed && iter < max_iters
                iter = iter + 1;
                loop_changed = false;

                % 1. Extract Topological Data from mstG
                bins = conncomp(mstG);
                numComps = max(bins);

                depths = mstG.Nodes.Depth;
                parent_nodes = mstG.Nodes.ParentNode;
                comp_roots = zeros(numComps, 1);

                for c = 1:numComps
                    compNodeIdx = find(bins == c);
                    if length(compNodeIdx) == 1
                        comp_roots(c) = compNodeIdx(1);
                        continue;
                    end

                    % The root is intrinsically defined by the BFS in resolveNodeStates
                    rootNodes = compNodeIdx(depths(compNodeIdx) == 0);
                    if isempty(rootNodes)
                        comp_roots(c) = compNodeIdx(1);
                    else
                        comp_roots(c) = rootNodes(1);
                    end
                end

                % 2. Compute group reference centers based on selected mode
                [unqStates, ~, groupIDs] = unique(mstG.Nodes.State);
                numStates = length(unqStates);
                groupRefOris = cell(numStates, 1);

                if strcmp(outlierMode, 'empirical')
                    for s = 1:numStates
                        s_nodes = find(groupIDs == s);
                        phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");

                        if isempty(phys_nodes)
                            continue;
                        end

                        oris = mstG.Nodes.MeanOrientation(phys_nodes);

                        if length(oris) <= 1
                            groupRefOris{s} = oris;
                        else
                            groupRefOris{s} = obj.computeRobustCenter(oris);
                        end
                    end
                else
                    % Theoretical mode: Roots first
                    root_centers = containers.Map();
                    safe_roots = containers.Map();

                    % Pass 1: Extract and store root centers and safe nodes
                    for s = 1:numStates
                        state_str = char(unqStates{s});
                        if endsWith(state_str, 'ROOT')
                            s_nodes = find(groupIDs == s);
                            phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");
                            if ~isempty(phys_nodes)
                                oris = mstG.Nodes.MeanOrientation(phys_nodes);
                                if length(oris) <= 1
                                    r_ori = oris;
                                else
                                    r_ori = obj.computeRobustCenter(oris);
                                end

                                parts = strsplit(state_str, '_');
                                comp_idx_str = parts{1};
                                root_centers(comp_idx_str) = r_ori;

                                devs = angle(oris, r_ori);
                                s_safe = phys_nodes(devs < max_err_threshold);
                                if isempty(s_safe)
                                    s_safe = phys_nodes;
                                end
                                safe_roots(comp_idx_str) = s_safe;
                            end
                        end
                    end
                end

                % 3. Identify Outliers and Isolated Nodes
                all_problems = []; % [node_id, dev_val, is_outlier, depth, closest_sr]

                if strcmp(outlierMode, 'empirical')
                    for s = 1:numStates
                        s_nodes = find(groupIDs == s);
                        phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");

                        % Empirical mode: Single nodes cannot form a consensus center, flag as isolated
                        if length(phys_nodes) == 1
                            v_node = phys_nodes(1);
                            if parent_nodes(v_node) ~= 0 && ~ismember(v_node, unfixable_nodes)
                                all_problems = [all_problems; v_node, 0, 0, depths(v_node), 0];
                            end
                        end

                        if length(phys_nodes) <= 1 || isempty(groupRefOris{s})
                            continue;
                        end

                        ref_ori = groupRefOris{s};
                        oris = mstG.Nodes.MeanOrientation(phys_nodes);
                        devs = angle(oris, ref_ori);

                        outlier_idx = find(devs > max_err_threshold);
                        for k = 1:length(outlier_idx)
                            v_node = phys_nodes(outlier_idx(k));
                            if parent_nodes(v_node) ~= 0 && ~ismember(v_node, unfixable_nodes)
                                all_problems = [all_problems; v_node, devs(outlier_idx(k)), 1, depths(v_node), 0];
                            end
                        end
                    end
                else
                    % Theoretical mode: No such thing as isolated! Every state has a theoretical center.
                    for s = 1:numStates
                        s_nodes = find(groupIDs == s);
                        phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");

                        if isempty(phys_nodes)
                            continue;
                        end

                        state_str = char(unqStates{s});
                        parts = strsplit(state_str, '_');
                        comp_idx_str = parts{1};

                        if ~safe_roots.isKey(comp_idx_str)
                            continue;
                        end

                        s_safe = safe_roots(comp_idx_str);

                        for p_idx = 1:length(phys_nodes)
                            v_node = phys_nodes(p_idx);
                            if ismember(v_node, unfixable_nodes) || parent_nodes(v_node) == 0
                                continue;
                            end

                            if length(s_safe) == 1
                                closest_sr = s_safe(1);
                            else
                                v_cent = mstG.Nodes.Centroid(v_node, :);
                                s_cents = mstG.Nodes.Centroid(s_safe, :);
                                D_sq = (s_cents(:,1) - v_cent(1)).^2 + (s_cents(:,2) - v_cent(2)).^2;
                                [~, min_idx] = min(D_sq);
                                closest_sr = s_safe(min_idx);
                            end

                            theoOri = mstG.Nodes.MeanOrientation(closest_sr);
                            for v_idx = 2:length(parts)
                                if strcmp(parts{v_idx}, 'ROOT')
                                    continue;
                                end
                                variant_idx = str2double(parts{v_idx});
                                theoOri = theoOri * obj.variants(variant_idx);
                            end

                            if v_node == 15
                                disp('Is 15 a problem? Its root node is: ')
                                disp(closest_sr)
                                disp('. All safe roots are: ')
                                disp(s_safe)
                                % keyboard; % DEBUG HALT: Inspect theoOri and parts before it gets flagged!
                            end

                            dev_val = angle(mstG.Nodes.MeanOrientation(v_node), theoOri);
                            if dev_val > max_err_threshold
                                all_problems = [all_problems; v_node, dev_val, 1, depths(v_node), closest_sr];
                            end
                        end
                    end
                end

                if isempty(all_problems)
                    break; % Graph is stable!
                end

                if target_node_override > 0
                    all_problems = all_problems(all_problems(:, 1) == target_node_override, :);
                    if isempty(all_problems)
                        break;
                    end
                end

                % 4. Sort problems by depth (shallowest first)
                [~, sort_idx] = sort(all_problems(:, 4));
                all_problems = all_problems(sort_idx, :);

                if debug_mode
                    fprintf('\n  [DEBUG] Iteration %d - Identified %d problem nodes:\n', iter, size(all_problems, 1));
                    for p = 1:size(all_problems, 1)
                        p_id = all_problems(p, 1);
                        p_dev = all_problems(p, 2);
                        p_is_outlier = all_problems(p, 3) == 1;
                        p_depth = all_problems(p, 4);
                        if p_is_outlier
                            fprintf('    - Grain %d: OUTLIER (Dev: %.2f deg, Depth: %d)\n', p_id, p_dev / degree, p_depth);
                        else
                            fprintf('    - Grain %d: ISOLATED (Depth: %d)\n', p_id, p_depth);
                        end
                    end
                end

                for p_idx = 1:size(all_problems, 1)
                    v_node = all_problems(p_idx, 1);
                    is_outlier = all_problems(p_idx, 3) == 1;

                if is_outlier
                    if strcmp(outlierMode, 'theoretical')
                        fprintf('\nAttempting to rewire OUTLIER Grain %d (State: %s) at depth %d with deviation %.2f deg to its theoretical orientation as calculated from the closest safe root Grain %d...\n', ...
                            v_node, unqStates{groupIDs(v_node)}, all_problems(p_idx, 4), all_problems(p_idx, 2)/degree, all_problems(p_idx, 5));
                    else
                        fprintf('\nAttempting to rewire OUTLIER Grain %d (State: %s) at depth %d...\n', v_node, unqStates{groupIDs(v_node)}, all_problems(p_idx, 4));
                    end
                else
                    fprintf('\nAttempting to rewire ISOLATED Grain %d (State: %s) at depth %d...\n', v_node, unqStates{groupIDs(v_node)}, all_problems(p_idx, 4));
                end

                if v_node == 16
                    disp('Halt stop')
                end

                % 5. Subtree Collection & Component Identification
                p_node = parent_nodes(v_node);
                c_idx = bins(v_node);
                comp_root = comp_roots(c_idx);
                compNodeIdx = find(bins == c_idx);
                compNodeIdx = compNodeIdx(:); % Force column vector

                subtree_nodes = v_node;
                sq = v_node;
                while ~isempty(sq)
                    curr = sq(1); sq(1) = [];
                    chs = find(parent_nodes == curr);
                    if ~isempty(chs)
                        subtree_nodes = [subtree_nodes; chs];
                        sq = [sq; chs];
                    end
                end

                % Temporarily mask the old edge in search_G
                edgeIdx1 = findedge(search_G, p_node, v_node);
                edgeIdx2 = findedge(search_G, v_node, p_node);
                old_w1 = Inf; old_w2 = Inf;
                if edgeIdx1 > 0, old_w1 = search_G.Edges.Weight(edgeIdx1); search_G.Edges.Weight(edgeIdx1) = inf; end
                if edgeIdx2 > 0, old_w2 = search_G.Edges.Weight(edgeIdx2); search_G.Edges.Weight(edgeIdx2) = inf; end

                fallback_success = false;

                masked_fallback_edges = [];
                masked_fallback_weights = [];

                for fallback_iter = 1:10
                    % Bias search_G to strongly prefer the CURRENT topology (mstG)
                    % This prevents zero-weight edges from arbitrarily shuffling and destroying geometric safe roots!
                    bias_val = 1e-6;
                    search_G.Edges.Weight = search_G.Edges.Weight + bias_val;

                    for e_idx = 1:numedges(search_G)
                        s_g = search_G.Edges.EndNodes(e_idx, 1);
                        t_g = search_G.Edges.EndNodes(e_idx, 2);
                        if findedge(mstG, s_g, t_g) > 0 && search_G.Edges.Weight(e_idx) < inf
                            search_G.Edges.Weight(e_idx) = search_G.Edges.Weight(e_idx) - bias_val;
                        end
                    end

                    % Extract the MST for this component using the stabilized search_G
                    [localS, localT, localSeq, localDev, localVar, localE1, localE2, localWeight] = obj.extractComponentMST(search_G, compNodeIdx);

                    % Remove the bias to restore search_G mathematically
                    for e_idx = 1:numedges(search_G)
                        s_g = search_G.Edges.EndNodes(e_idx, 1);
                        t_g = search_G.Edges.EndNodes(e_idx, 2);
                        if findedge(mstG, s_g, t_g) > 0 && search_G.Edges.Weight(e_idx) < inf
                            search_G.Edges.Weight(e_idx) = search_G.Edges.Weight(e_idx) + bias_val;
                        end
                    end
                    search_G.Edges.Weight = search_G.Edges.Weight - bias_val;

                    if isempty(localS)
                        if debug_mode
                            fprintf('    [DEBUG] Iter %d: Subtree is fundamentally disconnected (no finite paths left). Aborting fallback.\n', fallback_iter);
                        end
                        break;
                    end

                    % Find the edge bridging the subtree and main tree in the NEW MST
                    bridge_s = 0; bridge_t = 0; bridge_w = 0;
                    for ei = 1:length(localS)
                        g_s = compNodeIdx(localS(ei));
                        g_t = compNodeIdx(localT(ei));
                        is_s_in_sub = ismember(g_s, subtree_nodes);
                        is_t_in_sub = ismember(g_t, subtree_nodes);
                        if is_s_in_sub ~= is_t_in_sub
                            bridge_s = g_s;
                            bridge_t = g_t;
                            bridge_w = localWeight(ei);
                            break;
                        end
                    end

                    if isinf(bridge_w)
                        if debug_mode
                            fprintf('    [DEBUG] Iter %d: Subtree is fundamentally disconnected (no finite paths left). Aborting fallback.\n', fallback_iter);
                        end
                        break;
                    end

                    % Create a temporary graph to test the new MST
                    temp_mstG = mstG;

                    % Remove all existing edges for this component from temp_mstG
                    subG_old = subgraph(temp_mstG, compNodeIdx);
                    if numedges(subG_old) > 0
                        sub_s = compNodeIdx(subG_old.Edges.EndNodes(:,1));
                        sub_t = compNodeIdx(subG_old.Edges.EndNodes(:,2));
                        temp_mstG = rmedge(temp_mstG, sub_s, sub_t);
                    end

                    % Add the new MST edges
                    newEdgesTable = table([compNodeIdx(localS), compNodeIdx(localT)], localSeq, localDev, localVar, localE1, localE2, ...
                        'VariableNames', {'EndNodes', 'Sequence', 'AngularDeviation', 'Variant', 'MeanEbsdOri1', 'MeanEbsdOri2'});
                    temp_mstG = addedge(temp_mstG, newEdgesTable);

                    % Re-propagate FZ states perfectly down the new topology!
                    temp_mstG.Nodes.State(:) = "[]";
                    temp_mstG = obj.resolveNodeStates(temp_mstG);

                    % Validate the new orientation of v_node
                    new_state = temp_mstG.Nodes.State(v_node);

                    % Compute robust center for the new state to check deviation
                    % [temp_unqStates, ~, temp_groupIDs] = unique(temp_mstG.Nodes.State);
                    % s_idx = find(temp_unqStates == new_state);
                    phys_nodes = find(temp_mstG.Nodes.State == new_state & temp_mstG.Nodes.Type == "Physical");

                    dev_val = inf;
                    if strcmp(outlierMode, 'empirical')
                        if length(phys_nodes) > 1
                            oris = temp_mstG.Nodes.MeanOrientation(phys_nodes);
                            robust_ref = obj.computeRobustCenter(oris);

                            v_ori = temp_mstG.Nodes.MeanOrientation(v_node);
                            dev_val = angle(v_ori, robust_ref);
                        end
                    else
                        % Theoretical mode: Calculate deviation against closest safe root
                        new_state_char = char(new_state);
                        parts = strsplit(new_state_char, '_');
                        comp_idx_str = parts{1};
                        if safe_roots.isKey(comp_idx_str)
                            s_safe = safe_roots(comp_idx_str);
                            if length(s_safe) == 1
                                closest_sr = s_safe(1);
                            else
                                v_cent = temp_mstG.Nodes.Centroid(v_node, :);
                                s_cents = temp_mstG.Nodes.Centroid(s_safe, :);
                                D_sq = (s_cents(:,1) - v_cent(1)).^2 + (s_cents(:,2) - v_cent(2)).^2;
                                [~, min_idx] = min(D_sq);
                                closest_sr = s_safe(min_idx);
                            end

                            theoOri = temp_mstG.Nodes.MeanOrientation(closest_sr);
                            for v_idx = 2:length(parts)
                                if strcmp(parts{v_idx}, 'ROOT')
                                    continue;
                                end
                                variant_idx = str2double(parts{v_idx});
                                theoOri = theoOri * obj.variants(variant_idx);
                            end
                            v_ori = temp_mstG.Nodes.MeanOrientation(v_node);
                            dev_val = angle(v_ori, theoOri);
                        end
                    end

                    % Determine generation for the old and new states
                    old_state = mstG.Nodes.State(v_node);
                    old_parts = strsplit(char(old_state), '_');
                    old_gen = length(old_parts) - 1;
                    if strcmp(old_parts{end}, 'ROOT')
                        old_gen = 0;
                    end

                    new_parts = strsplit(char(new_state), '_');
                    new_gen = length(new_parts) - 1;
                    if strcmp(new_parts{end}, 'ROOT')
                        new_gen = 0;
                    end

                    gen_increased = new_gen > old_gen;

                    if dev_val < max_err_threshold && ~(reject_gen_increase && gen_increased)
                        fallback_success = true;

                        if debug_mode
                            fprintf('    [DEBUG] Iter %d SUCCESS for Grain %d:\n', fallback_iter, v_node);
                            fprintf('      - Severed Edge:  Grain %d <-> Grain %d\n', v_node, p_node);
                            fprintf('      - New State:     %s\n', new_state);
                            fprintf('      - New Deviation: %.2f deg\n', dev_val / degree);
                        else
                            fprintf('  -> SUCCESS: Rewired %d optimally via Component MST.\n', v_node);
                        end

                        mstG = temp_mstG;
                        break;
                    else
                        if debug_mode
                            if reject_gen_increase && gen_increased && dev_val < max_err_threshold
                                fprintf('    [DEBUG] Iter %d FAILED validation (Generation increase rejected!):\n', fallback_iter);
                                fprintf('      - Original State: %s (Gen %d)\n', old_state, old_gen);
                                fprintf('      - Resulting State: %s (Gen %d)\n', new_state, new_gen);
                                fprintf('      - Deviation: %.2f deg (Passed %.2f deg threshold)\n', dev_val/degree, max_err_threshold/degree);
                            else
                                fprintf('    [DEBUG] Iter %d FAILED validation:\n', fallback_iter);
                                fprintf('      - Resulting State: %s\n', new_state);
                                fprintf('      - Deviation: %.2f deg (Exceeds %.2f deg threshold!)\n', dev_val/degree, max_err_threshold/degree);
                            end
                        end

                        if bridge_s > 0
                            fail_edge = findedge(search_G, bridge_s, bridge_t);
                            if fail_edge == 0, fail_edge = findedge(search_G, bridge_t, bridge_s); end
                            if fail_edge > 0
                                masked_fallback_edges(end+1) = fail_edge;
                                masked_fallback_weights(end+1) = search_G.Edges.Weight(fail_edge);
                                search_G.Edges.Weight(fail_edge) = inf;
                            end
                            if debug_mode
                                fprintf('      - Masking fallback bridge edge: %d <-> %d\n', bridge_s, bridge_t);
                            end
                        else
                            break; % Cannot find a bridge edge, abort fallback
                        end
                    end
                end

                % Restore all temporarily masked edges in search_G
                for k = 1:length(masked_fallback_edges)
                    search_G.Edges.Weight(masked_fallback_edges(k)) = masked_fallback_weights(k);
                end

                if edgeIdx1 > 0, search_G.Edges.Weight(edgeIdx1) = old_w1; end
                if edgeIdx2 > 0, search_G.Edges.Weight(edgeIdx2) = old_w2; end

                if ~fallback_success
                    fprintf('  -> FAILED: Exhausted all alternative component MSTs for Grain %d.\n', v_node);
                    unfixable_nodes(end+1) = v_node;
                else
                    loop_changed = true;
                    break; % Break the p_idx loop to restart the graph evaluation
                end
                end
            end

            % --- POST-PRUNING THEORETICAL VALIDATION ---
            if debug_mode
                fprintf('\n--- Post-Pruning Theoretical Validation ---\n');
                max_err_threshold = 15 * degree;

                % Extract root robust centers for each component
                root_centers = containers.Map();
                for s = 1:numStates
                    state_str = char(unqStates{s});
                    if endsWith(state_str, 'ROOT') && ~isempty(groupRefOris{s})
                        parts = strsplit(state_str, '_');
                        comp_idx_str = parts{1}; % e.g., 'C1'
                        root_centers(comp_idx_str) = groupRefOris{s};
                    end
                end

                % Validate every non-root group
                for s = 1:numStates
                    state_str = char(unqStates{s});
                    if endsWith(state_str, 'ROOT') || isempty(groupRefOris{s})
                        continue;
                    end

                    parts = strsplit(state_str, '_');
                    comp_idx_str = parts{1};

                    if ~root_centers.isKey(comp_idx_str)
                        continue; % Shouldn't happen, but safe fallback
                    end

                    % Reconstruct theoretical orientation purely from the root and the twin chain
                    theoOri = root_centers(comp_idx_str);
                    for v_idx = 2:length(parts)
                        variant_idx = str2double(parts{v_idx});
                        theoOri = theoOri * obj.variants(variant_idx);
                    end

                    empiricalOri = groupRefOris{s};
                    dev_val = angle(empiricalOri, theoOri);

                    if dev_val > max_err_threshold
                        fprintf('  [WARNING] Group %s deviates from theoretical chain by %.2f deg!\n', state_str, dev_val / degree);
                    else
                        fprintf('  [OK] Group %s aligns theoretically (Dev: %.2f deg)\n', state_str, dev_val / degree);
                    end
                end
                fprintf('-------------------------------------------\n\n');
            end
        end
        function Q = reduceToQuotient(obj, mstG)
            % Workflow Step 5: Contract nodes with same state into quotient graph
            %
            % REDUCETOQUOTIENT Reduces the MST to a quotient graph.
            % Groups nodes with identical state strings into single
            % quotient nodes. Quotient nodes represent different twin variants
            % and twin generations (primary/secondary/tertiary/..).
            %
            % Parameters:
            %   mstG - Minimum Spanning Tree graph
            %
            % Returns:
            %   Q - The reduced quotient digraph

            Q = digraph();
            if numnodes(mstG) == 0
                obj.QuotientGraph = Q;
                return;
            end

            [unqStates, ~, groupIDs] = unique(mstG.Nodes.State);
            numGroups = length(unqStates);

            Q = addnode(Q, numGroups);
            Q.Nodes.State = unqStates;

            % Compute mean orientation for each quotient node
            nodeMeanOris = orientation.nan(numGroups, 1, obj.grains.CS);
            nodeCentroids = zeros(numGroups, 2);
            nodeAreas = zeros(numGroups, 1);

            for g = 1:numGroups
                g_nodes = find(groupIDs == g);
                phys_nodes = g_nodes(mstG.Nodes.Type(g_nodes) == "Physical");

                if ~isempty(phys_nodes)
                    oris = mstG.Nodes.MeanOrientation(phys_nodes);
                    areas = mstG.Nodes.Area(phys_nodes);

                    if length(oris) <= 1
                        nodeMeanOris(g) = oris;
                    else
                        nodeMeanOris(g) = obj.computeRobustCenter(oris);
                    end

                    total_area = sum(areas);
                    if total_area > 0
                        weighted_x = sum(mstG.Nodes.Centroid(phys_nodes, 1) .* areas) / total_area;
                        weighted_y = sum(mstG.Nodes.Centroid(phys_nodes, 2) .* areas) / total_area;
                        nodeCentroids(g, :) = [weighted_x, weighted_y];
                    end
                    nodeAreas(g) = total_area;
                else
                    % Purely virtual group! (Represents a twin generation with no physical grains)
                    virt_nodes = g_nodes(mstG.Nodes.Type(g_nodes) == "Virtual");
                    if ~isempty(virt_nodes)
                        oris = mstG.Nodes.MeanOrientation(virt_nodes);
                        nodeMeanOris(g) = mean(oris); % Fallback to virtual orientations
                        nodeCentroids(g, :) = mean(mstG.Nodes.Centroid(virt_nodes, :), 1);
                        nodeAreas(g) = 0; % Virtual nodes have 0 physical area
                    end
                end
            end

            Q.Nodes.MeanOrientation = nodeMeanOris;
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
            obj.nodeGroups = groupIDs;
        end

        function plot(obj, varargin)
            % PLOT High-level plotting function to visualize the grains and the quotient graph.
            % Displays the MTEX physical grains with the reduced quotient network overlaid.
            %
            % Parameters:
            %   varargin - Additional plotting options passed to MTEX plot functions
            p = inputParser;
            addParameter(p, 'twinsOnly', true);
            addParameter(p, 'noText', false);
            addParameter(p, 'lineWidth', 1);
            parse(p, varargin{:});

            twinsOnly = p.Results.twinsOnly;
            noText = p.Results.noText;
            lineWidth = p.Results.lineWidth;

            if isempty(obj.QuotientGraph)
                warning('TwinNetwork:plot:noQuotient', 'Quotient graph not found. Please call reduceToQuotient() first.');
                return;
            end

            figure;
            % Initialize white base for all grains
            numGrains = length(obj.grains);
            faceColors = ones(numGrains, 3);

            labelCentersX = [];
            labelCentersY = [];
            labelTexts = {};

            cmap = lines();
            numCmap = size(cmap, 1);

            maxPhysId = max(obj.grains.id);
            id2idx = zeros(maxPhysId, 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);

            toPlot = true(numGrains, 1);
            Q_deg = indegree(obj.QuotientGraph) + outdegree(obj.QuotientGraph);

            for i = 1:maxPhysId
                idx = id2idx(i);
                if idx > 0
                    grp = obj.nodeGroups(i);
                    isIsolated = Q_deg(grp) == 0;

                    if twinsOnly && isIsolated
                        toPlot(idx) = false;
                        continue;
                    end

                    if isIsolated
                        continue; % Remains white, no text overlay
                    end

                    % Parse state string, e.g., "C1_ROOT", "C1_1", "C1_1_2"
                    stateStr = char(obj.QuotientGraph.Nodes.State(grp));
                    parts = strsplit(stateStr, '_');

                    if length(parts) >= 2 && strcmp(parts{2}, 'ROOT')
                        % Root Generation
                        gen = 0;
                        path = [];
                    else
                        % Has variants
                        gen = length(parts) - 1;
                        path = str2double(parts(2:end));
                    end

                    % Determine color based on generation and first variant
                    if gen == 0
                        faceColors(idx, :) = [0.7 0.7 0.7]; % Grey
                        lbl = 'P';
                    elseif gen == 1
                        variant = path(1);
                        cIdx = mod(variant - 1, numCmap) + 1;
                        faceColors(idx, :) = cmap(cIdx, :);
                        lbl = sprintf('V%d', variant);
                    elseif gen == 2
                        variant = path(1);
                        cIdx = mod(variant - 1, numCmap) + 1;
                        baseColor = cmap(cIdx, :);
                        faceColors(idx, :) = baseColor * 0.6; % Darker

                        pathStr = sprintf('%d,', path);
                        lbl = sprintf('V%s', pathStr(1:end-1));
                    else
                        faceColors(idx, :) = [0.50 1.00 0.00]; % Green

                        pathStr = sprintf('%d,', path);
                        lbl = sprintf('V%s', pathStr(1:end-1));
                    end

                    c = obj.grains(idx).centroid;
                    labelCentersX(end+1) = c.x;
                    labelCentersY(end+1) = c.y;
                    labelTexts{end+1} = lbl;
                end
            end

            % Plot the grains
            if twinsOnly
                plot(obj.grains(toPlot), faceColors(toPlot, :), 'micronbar', 'off','lineWidth', lineWidth);
            else
                plot(obj.grains, faceColors, 'micronbar', 'off', 'lineWidth', lineWidth);
            end
            hold on;

            % Plot text
            if ~isempty(labelTexts) && ~noText
                text(labelCentersX, labelCentersY, labelTexts, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 10, ...
                    'Color', 'k');
            end
            hold off;
        end

        function plotSigmaSections(obj, mstG, varargin)
            % PLOTSIGMASECTIONS Plots the sigma sections of components of the twin network individually.
            % Useful for isolating specific twin clusters in a complex microstructure.
            %
            % Parameters:
            %   mstG     - The Minimum Spanning Tree graph containing component bins
            %   varargin - Additional plotting options
            % Plots orientations of groups in sigma sections to verify clustering
            % % Plot ALL groups with their robust centers (as stars)
            % network.plotSigmaSections(mstG);
            % % Plot only a specific group without the legend
            % network.plotSigmaSections(mstG, 'group', 'C1_1_2', 'showLegend', false);
            % % Plot all groups but turn off the robust center stars
            % network.plotSigmaSections(mstG, 'robustCenter', false);
            p = inputParser;
            addParameter(p, 'group', 'all');
            addParameter(p, 'showLegend', true);
            addParameter(p, 'robustCenter', true);
            parse(p, varargin{:});

            targetGroup = p.Results.group;
            showLegend = p.Results.showLegend;
            plotRobust = p.Results.robustCenter;

            if isempty(obj.QuotientGraph)
                warning('TwinNetwork:plotSigmaSections:noQuotient', 'Quotient graph not found. Please call reduceToQuotient() first.');
                return;
            end

            states = obj.QuotientGraph.Nodes.State;
            numStates = length(states);

            % Determine which states to plot
            if strcmpi(targetGroup, 'all')
                plotIdxs = 1:numStates;
            else
                matchIdx = find(states == string(targetGroup));
                if isempty(matchIdx)
                    warning('Group %s not found!', targetGroup);
                    return;
                end
                plotIdxs = matchIdx;
            end

            figure;
            cmap = lines();
            numCmap = size(cmap, 1);

            h_leg = [];
            l_leg = {};
            isFirstPlot = true;

            for k = 1:length(plotIdxs)
                s = plotIdxs(k);
                state_str = states(s);
                color = cmap(mod(s-1, numCmap) + 1, :);

                % Find physical nodes in this group
                s_nodes = find(obj.nodeGroups == s);
                phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");

                if isempty(phys_nodes)
                    continue;
                end

                % Get FZ-propagated orientations from mstG
                oris = mstG.Nodes.MeanOrientation(phys_nodes);

                try
                    plotSection(oris, 'sigma', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'Marker', 'o', 'MarkerSize', 5);
                catch
                    plotSection(oris, 'sigma');
                end

                if isFirstPlot
                    hold on;
                    isFirstPlot = false;
                end

                % Dummy handle for legend
                h_dummy = plot(NaN, NaN, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
                h_leg(end+1) = h_dummy;
                l_leg{end+1} = char(state_str);

                if plotRobust
                    % Robust center is stored in QuotientGraph
                    ref_ori = obj.QuotientGraph.Nodes.MeanOrientation(s);
                    try
                        plotSection(ref_ori, 'sigma', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'Marker', 'p', 'MarkerSize', 14, 'LineWidth', 1.5);
                        plotSection(ref_ori, 'sigma', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'Marker', 'p', 'MarkerSize', 10);
                    catch
                    end
                end
            end

            if showLegend && ~isempty(h_leg)
                legend(h_leg, l_leg, 'Interpreter', 'none', 'Location', 'best');
            end
            hold off;
        end

        function plotGraph(obj, G, varargin)
            % PLOTGRAPH Low-level plotting utility for displaying TwinNetwork graphs on top of the grain plot.
            % Overlays nodes and edges on the current axis based on spatial centroids.
            % Edges are color-coded by sequence length or deviation.
            %
            % Parameters:
            %   G        - Any digraph or graph from the TwinNetwork pipeline
            %   varargin - Additional plotting options
            %              'noGrainsBG' - logical, if true, suppresses background grains and uses force layout

            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'noGrainsBG', false, @islogical);
            parse(p, varargin{:});
            noGrainsBG = p.Results.noGrainsBG;

            % Basic plotting utility matching original MTEX overlay style
            if nargin < 2 || isempty(G) || numnodes(G) == 0
                disp('Graph is empty or not provided.');
                return;
            end

            figure;
            tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
            nexttile;

            if ~noGrainsBG
                plot(obj.grains, 'micronbar', 'off', 'FaceColor', [0.85 0.85 0.85]);
                hold on;
                plot(obj.grains.boundary)

                if ismember('Centroid', G.Nodes.Properties.VariableNames)
                    X = G.Nodes.Centroid(:, 1);
                    Y = G.Nodes.Centroid(:, 2);
                    pPlot = plot(G, 'XData', X, 'YData', Y);
                else
                    pPlot = plot(G);
                end
            else
                % When background is off, use the force layout to cleanly display the graph structure
                pPlot = plot(G, 'Layout', 'force');
                hold on;
            end

            pPlot.MarkerSize = 5;
            pPlot.NodeColor = 'k';
            pPlot.LineWidth = 1.5;

            % Labeling logic
            if ismember('State', G.Nodes.Properties.VariableNames) && any(G.Nodes.State ~= "[]")
                pPlot.NodeLabel = G.Nodes.State;
            elseif ismember('GrainId', G.Nodes.Properties.VariableNames)
                nodeLabels = cell(numnodes(G), 1);
                for i = 1:numnodes(G)
                    if isnan(G.Nodes.GrainId(i))
                        nodeLabels{i} = ['v', num2str(i)];
                    else
                        nodeLabels{i} = num2str(G.Nodes.GrainId(i));
                    end
                end
                pPlot.NodeLabel = nodeLabels;
            end

            % Edge coloring
            if ismember('AngularDeviation', G.Edges.Properties.VariableNames)
                pPlot.EdgeCData = G.Edges.AngularDeviation ./ degree;
                colormap(gca, 'jet');
                cb = colorbar;
                cb.Label.String = 'Angular Deviation (degrees)';
            else
                pPlot.EdgeColor = 'k';
            end

            % Edge labeling
            if ismember('Variant', G.Edges.Properties.VariableNames)
                % Show variant. If NaN (like step 1), it will just show NaN.
                % Or we can show sequence if variant is NaN.
                labels = cell(numedges(G), 1);
                for i = 1:numedges(G)
                    v = G.Edges.Variant(i);
                    if isnan(v) && ismember('Sequence', G.Edges.Properties.VariableNames)
                        labels{i} = G.Edges.Sequence{i};
                    else
                        labels{i} = num2str(v);
                    end
                end
                labeledge(pPlot, 1:numedges(G), labels);
            end

        end
        
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
                
                % Use robustMeanOrientation if present, fallback to MeanOrientation
                if ismember('robustMeanOrientation', Q.Nodes.Properties.VariableNames)
                    currOri = Q.Nodes.robustMeanOrientation(curr_node);
                else
                    currOri = Q.Nodes.MeanOrientation(curr_node);
                end
                
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
        end
    end
end

function ori = catOrientationCells(values)

if isscalar(values)
    ori = values{1};
else
    ori = vertcat(values{:});
end

end
