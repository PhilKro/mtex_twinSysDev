classdef twinNetwork < handle
    properties
        grains
        gB
        tS
        maxGenerationLimit
        twinThresh
        Graph
        GraphUndirected
        virtualNodes
    end

    methods
        function obj = twinNetwork(grains, gB, tS, varargin)
            obj.grains = grains;
            obj.gB = gB;
            obj.tS = tS;

            p = inputParser;
            addParameter(p, 'maxGenerationLimit', 2);
            addParameter(p, 'twinThresh', 5*degree);
            parse(p, varargin{:});

            obj.maxGenerationLimit = p.Results.maxGenerationLimit;
            obj.twinThresh = p.Results.twinThresh;

            % Symmetrize twin system to get variants
            tS_sym = obj.tS.symmetrise;
            variants = tS_sym.parentTwinMisorientation;

            % Precompute possible sequences
            allSeqs = {};
            allMori = orientation.empty;

            N = length(variants);

            % Generation 1
            allSeqs = num2cell(1:N)';
            allMori = variants(:);
            genMori = ones(N,1);

            currentSeqs = allSeqs;
            currentMori = allMori;

            % Generations 2 to maxGenerationLimit
            for g = 2:obj.maxGenerationLimit
                M = length(currentSeqs);
                idx_rep = repmat((1:M)', N, 1);
                idx_repel = repelem((1:N)', M, 1);

                nextSeqs = cellfun(@(x, y) [x, y], currentSeqs(idx_rep), num2cell(idx_repel), 'UniformOutput', false);

                repMori = currentMori(idx_rep);
                newMori = variants(idx_repel);
                nextMori = repMori(:) .* newMori(:);

                allSeqs = [allSeqs; nextSeqs];
                allMori = [allMori; nextMori];
                genMori = [genMori; repmat(g, length(idx_rep), 1)];

                currentSeqs = nextSeqs;
                currentMori = nextMori;
            end

            % 1. Extract and sort grain IDs (ensures [1, 2] and [2, 1] are treated as the same pair)
            grainIds = sort(obj.gB.grainId, 2);

            % 2. Filter out specimen boundaries (where grainId is 0) first
            % Filtering early keeps your tags perfectly aligned with your valid data
            validMask = all(grainIds ~= 0, 2);
            validGrainIds = grainIds(validMask, :);
            validMisorientations = obj.gB.misorientation(validMask);

            % 3. Find unique pairs AND get the grouping tag for every segment
            % pairTags is a vector assigning an integer (1 to numPairs) to every valid boundary segment
            [pairs, ~, pairTags] = unique(validGrainIds, 'rows');

            if isempty(pairs)
                obj.Graph = digraph();
                obj.virtualNodes = table([], [], {}, 'VariableNames', {'x', 'y', 'Sequence'});
                return;
            end

            % Map IDs to indices (Needed for BFS later)
            id2idx = zeros(max(obj.grains.id), 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);

            % 4. Calculate the mean misorientation per pair without a loop
            % accumarray groups the row indices by their pairTags and applies the mean function
            rowIndices = (1:length(pairTags))';
            meanMoriCell = accumarray(pairTags, rowIndices, [], @(idx) {mean(validMisorientations(idx))});

            % 5. Vectorized angle calculation across all pairs and variants
            % Concatenate cell array into a single MTEX misorientation array
            mori_pairs = [meanMoriCell{:}];

            % Ensure mori_pairs is a column and allMori is a row for broadcasting
            angMatrix = angle(mori_pairs(:), allMori(:)'); % Returns P x N matrix

            % Favor lower generation twins by adding a penalty to higher generation angles
            validMatches = angMatrix < obj.twinThresh;

            penaltyMatrix = repmat(genMori(:)' * 1000, size(angMatrix, 1), 1);
            scoreMatrix = angMatrix + penaltyMatrix;
            scoreMatrix(~validMatches) = inf;

            [bestScores, minIdxs] = min(scoreMatrix, [], 2);
            isTwinEdge = bestScores < inf;

            edgeSeqs = allSeqs(minIdxs);
            twinPairs = pairs(isTwinEdge, :);



            % Create undirected graph for connected components
            numPhysNodes = max(obj.grains.id);
            G_undir = graph();
            G_undir = addnode(G_undir, numPhysNodes);
            if ~isempty(twinPairs)
                G_undir = addedge(G_undir, twinPairs(:,1), twinPairs(:,2));
            end

            bins = conncomp(G_undir);

            visited = false(numPhysNodes, 1);
            edgeProcessed = false(size(pairs, 1), 1);

            s_nodes = [];
            t_nodes = [];
            variants_list = {};

            vNodeIdx = numPhysNodes;
            vNodeCoords = [];
            vNodeGens = {};

            for c = 1:max(bins)
                compNodes = find(bins == c);
                if length(compNodes) == 1
                    continue;
                end

                validNodes = compNodes(id2idx(compNodes) > 0);
                if isempty(validNodes)
                    continue;
                end
                areas = obj.grains(id2idx(validNodes)).area;
                [~, maxAreaIdx] = max(areas);
                rootNode = validNodes(maxAreaIdx);

                queue = rootNode;
                visited(rootNode) = true;

                while ~isempty(queue)
                    currNode = queue(1);
                    queue(1) = [];

                    nbrs = neighbors(G_undir, currNode);
                    for n = 1:length(nbrs)
                        nxtNode = nbrs(n);
                        if ~visited(nxtNode)
                            visited(nxtNode) = true;
                            queue(end+1) = nxtNode;
                        end

                        % Identify the specific undirected edge
                        pair_id = [min(currNode, nxtNode), max(currNode, nxtNode)];
                        [~, loc] = ismember(pair_id, pairs, 'rows');

                        % Ensure every edge is traced exactly once in the directed graph
                        if edgeProcessed(loc)
                            continue;
                        end
                        edgeProcessed(loc) = true;

                        % Determine directed variant sequence using boundary misorientation
                        is_forward = currNode < nxtNode;
                        m_boundary = mori_pairs(loc);

                        % MTEX boundary misorientation typically goes from grainId(2) to grainId(1)
                        if is_forward
                            mori_directed = inv(m_boundary);
                        else
                            mori_directed = m_boundary;
                        end

                        angList = angle(mori_directed, allMori(:)');
                        validMatches = angList < obj.twinThresh;

                        scoreList = angList + genMori(:)' * 1000;
                        scoreList(~validMatches) = inf;

                        [bestScore, minMoriIdx] = min(scoreList);

                        if bestScore < inf
                            seq = allSeqs{minMoriIdx};
                        else
                            % Should not happen since we already checked undirected, but just in case
                            seq = 1; % fallback
                        end

                        prevNode = currNode;
                        for step = 1:length(seq)
                            v_step = seq(step);
                            if step == length(seq)
                                s_nodes(end+1, 1) = prevNode;
                                t_nodes(end+1, 1) = nxtNode;
                                variants_list{end+1, 1} = num2str(v_step);
                            else
                                vNodeIdx = vNodeIdx + 1;

                                coord1 = obj.grains(id2idx(currNode)).centroid;
                                coord2 = obj.grains(id2idx(nxtNode)).centroid;
                                alpha = step / length(seq);
                                vCoord = coord1 + alpha * (coord2 - coord1);

                                vNodeCoords(end+1, :) = vCoord;
                                vNodeGens{end+1, 1} = seq(1:step);

                                s_nodes(end+1, 1) = prevNode;
                                t_nodes(end+1, 1) = vNodeIdx;
                                variants_list{end+1, 1} = num2str(v_step);

                                prevNode = vNodeIdx;
                            end
                        end
                    end
                end
            end
            
            numTotalNodes = vNodeIdx;

        finalG = digraph();
        if numTotalNodes > 0
            finalG = addnode(finalG, numTotalNodes);
        end

        if ~isempty(s_nodes)
            EdgeTable = table([s_nodes, t_nodes], variants_list, 'VariableNames', {'EndNodes', 'Variant'});
            finalG = addedge(finalG, EdgeTable);
        end

        obj.Graph = finalG;
        obj.GraphUndirected = G_undir;

        if ~isempty(vNodeCoords)
            obj.virtualNodes = table(vNodeCoords(:,1), vNodeCoords(:,2), vNodeGens, 'VariableNames', {'x', 'y', 'Sequence'});
        else
            obj.virtualNodes = table([], [], {}, 'VariableNames', {'x', 'y', 'Sequence'});
        end
    end

    function plot(obj, varargin)
        p = inputParser;
        addOptional(p, 'type', 'directed');
        parse(p, varargin{:});
        plotType = p.Results.type;

        figure;
        plot(obj.grains, 'micronbar', 'off');
        hold on;

        if strcmpi(plotType, 'undirected')
            G = obj.GraphUndirected;
        else
            G = obj.Graph;
        end

        if isempty(G) || numnodes(G) == 0
            hold off;
            return;
        end

        numNodes = numnodes(G);
        X = zeros(numNodes, 1);
        Y = zeros(numNodes, 1);

        id2idx = zeros(max(obj.grains.id), 1);
        id2idx(obj.grains.id) = 1:length(obj.grains);

        for i = 1:max(obj.grains.id)
            if id2idx(i) > 0
                c = obj.grains(id2idx(i)).centroid;
                X(i) = c.x;
                Y(i) = c.y;
            else
                X(i) = NaN;
                Y(i) = NaN;
            end
        end

        if strcmpi(plotType, 'directed') && ~isempty(obj.virtualNodes) && size(obj.virtualNodes, 1) > 0
            vStart = max(obj.grains.id) + 1;
            X(vStart:end) = obj.virtualNodes.x;
            Y(vStart:end) = obj.virtualNodes.y;
        end

        p = plot(G, 'XData', X, 'YData', Y);
        p.MarkerSize = 5;
        p.NodeColor = 'k';
        p.EdgeColor = 'k';
        if strcmpi(plotType, 'directed')
            p.ArrowSize = 10;
        end
        p.LineWidth = 1.5;

        if ismember('Variant', G.Edges.Properties.VariableNames)
            labeledge(p, 1:numedges(G), G.Edges.Variant);
        end

        if strcmpi(plotType, 'directed') && ~isempty(obj.virtualNodes) && size(obj.virtualNodes, 1) > 0
            vIdx = (max(obj.grains.id) + 1):numNodes;
            highlight(p, vIdx, 'NodeColor', 'r', 'Marker', 'o', 'MarkerSize', 8);
        end

        nodeLabels = cell(numNodes, 1);
        for i = 1:numNodes
            if i <= max(obj.grains.id) && id2idx(i) > 0
                nodeLabels{i} = num2str(i);
            else
                nodeLabels{i} = '';
            end
        end
        p.NodeLabel = nodeLabels;

        hold off;
    end
end
end
