classdef twinNetwork < handle
    properties
        grains
        gB
        tS
        maxGenerationLimit
        twinThresh
        Graph
        untaggedGraph
        virtualNodes
        ebsd
        edgeData
        QuotientGraph
        nodeGroups
        mstGraph
        groupGenerations
        groupVariantPaths
    end
    methods
        function obj = twinNetwork(grains, gB, tS, ebsd, varargin)
            obj.grains = grains;
            obj.gB = gB;
            obj.tS = tS;
            obj.ebsd = ebsd;

            p = inputParser;
            addParameter(p, 'maxGenerationLimit', 2);
            addParameter(p, 'twinThresh', 3*degree);
            parse(p, varargin{:});

            obj.maxGenerationLimit = p.Results.maxGenerationLimit;
            obj.twinThresh = p.Results.twinThresh;

            % Symmetrize twin system to get variants
            tS_sym = obj.tS.symmetrise;
            variants = tS_sym.parentTwinMisorientation;
            if variants.isImproper
                variants = -variants;
            end

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
                nextMori = newMori(:) .* repMori(:);
                warning('twinNetwork:twinNetwork:moriConfusion','This feels like the wrong way round, but like this the variant order along the virtual edges seems correct. I expected nextMori = repMori(:) .* newMori(:); however then its flipped, but should be correct since its ori*mori1*mori2 and not ori*mori2*mori1. I dont understand.')


                allSeqs = [allSeqs; nextSeqs];
                allMori = [allMori; nextMori];
                genMori = [genMori; repmat(g, length(idx_rep), 1)];

                currentSeqs = nextSeqs;
                currentMori = nextMori;
            end

            % 1. Create edgeData struct to cache valid boundary orientations
            validMask = all(obj.gB.grainId ~= 0, 2);

            obj.edgeData = struct();
            obj.edgeData.gB_index = find(validMask);

            % Using original grain IDs for alignment with ebsdId
            gB_grainId = obj.gB.grainId(validMask, :);
            gB_ebsdId = obj.gB.ebsdId(validMask, :);

            obj.edgeData.grainId = gB_grainId;
            obj.edgeData.ebsdId = gB_ebsdId;

            % Orientations of EBSD pixels on both sides
            obj.edgeData.ebsdOri1 = obj.ebsd(ebsd.id2ind(gB_ebsdId(:,1))).orientations;
            obj.edgeData.ebsdOri2 = obj.ebsd(ebsd.id2ind(gB_ebsdId(:,2))).orientations;

            % Compute meanEbsdOri grouped by boundaries
            [~, ~, dirTags] = unique(gB_grainId, 'rows');
            meanE1 = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(obj.edgeData.ebsdOri1(idx))});
            meanE2 = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(obj.edgeData.ebsdOri2(idx))});

            meanE1_array = [meanE1{:}];
            meanE2_array = [meanE2{:}];

            % Extract local grain orientations to align EBSD data to the fundamental region
            grainOri1 = obj.grains(grains.id2ind(gB_grainId(:,1))).meanOrientation;
            grainOri2 = obj.grains(grains.id2ind(gB_grainId(:,2))).meanOrientation;

            obj.edgeData.meanEbsdOri1 = meanE1_array(dirTags);
            obj.edgeData.meanEbsdOri1 = obj.edgeData.meanEbsdOri1(:).project2FundamentalRegion(grainOri1);

            obj.edgeData.meanEbsdOri2 = meanE2_array(dirTags);
            obj.edgeData.meanEbsdOri2 = obj.edgeData.meanEbsdOri2(:).project2FundamentalRegion(grainOri2);

            % Explicit boundary misorientation (ebsdOri2 -> ebsdOri1)
            obj.edgeData.misorientation = inv(obj.edgeData.ebsdOri2) .* obj.edgeData.ebsdOri1;

            % 2. Extract and sort grain IDs for unique boundary pair grouping
            validGrainIds = sort(gB_grainId, 2);

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

            % 4. Calculate the mean misorientation per pair using edgeData
            % accumarray groups the row indices by their pairTags and applies the mean function
            rowIndices = (1:length(pairTags))';
            meanMoriCell = accumarray(pairTags, rowIndices, [], @(idx) {mean(obj.edgeData.misorientation(idx))});

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
            angular_deviation_list = [];

            vNodeIdx = numPhysNodes;
            vNodeCoords = [];
            vNodeGens = {};
            vNodeOris = orientation.empty(0, 1);
            warning('twinNetwork:twinNetwork:initialGrainErrorReduction','Maybe start the graph from the grain with the lowest GROD to size ratio? or sth similar')
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
                warning('a better code has been commented out here in order to troubleshoot the outlier logic')
                % --- REVERTED SWEEP LOGIC (COMMENTED OUT) ---
                % for current_max_gen = 1:max(genMori)
                %     queue = compNodes(visited(compNodes))';
                
                [~, maxAreaIdx] = max(areas);
                rootNode = validNodes(maxAreaIdx);
                
                % old solution (b)
                queue = rootNode;
                % old solution (e)
                visited(rootNode) = true;
                while ~isempty(queue)
                    currNode = queue(1);
                    queue(1) = [];
                    nbrs = neighbors(G_undir, currNode);
                    for n = 1:length(nbrs)
                        nxtNode = nbrs(n);
                        % old solution (b)
                        if ~visited(nxtNode)
                            visited(nxtNode) = true;
                            queue(end+1) = nxtNode;
                            is_newly_visited = true;
                        else
                            is_newly_visited = false;
                        end
                        % old solution (e)
                        % Identify the specific undirected edge
                        pair_id = [min(currNode, nxtNode), max(currNode, nxtNode)];
                        [~, loc] = ismember(pair_id, pairs, 'rows');
                        % Ensure every edge is traced exactly once in the directed graph
                        if edgeProcessed(loc)
                            continue;
                        end
                        
                        % --- REVERTED SWEEP LOGIC ---
                        % initial_gen = length(allSeqs{minIdxs(loc)});
                        % if initial_gen > current_max_gen
                        %     continue;
                        % end
                        
                        % old solution (b)
                        edgeProcessed(loc) = true;
                        % old solution (e)
                        % Determine directed variant sequence using dynamic boundary misorientation
                        pairIdxs = find(pairTags == loc);

                        firstIdx = pairIdxs(1);

                        ebsd1ori = mean(obj.edgeData.meanEbsdOri1(pairIdxs));
                        ebsd2ori = mean(obj.edgeData.meanEbsdOri2(pairIdxs));

                        % Rigorously determine which orientation belongs to currNode (parent)
                        if obj.edgeData.grainId(firstIdx, 1) == currNode
                            currOri = ebsd1ori;
                            nxtOri  = ebsd2ori;
                        else
                            currOri = ebsd2ori;
                            nxtOri  = ebsd1ori;
                        end

                        % Direct orientation angle evaluation (Angle between nxtOri and currOri * variant)
                        angList = angle(nxtOri, currOri * allMori(:)');

                        angList = angList(:);
                        validMatches = angList < obj.twinThresh;

                        scoreList = angList + genMori(:) * 1000;
                        scoreList(~validMatches) = inf;

                        [bestScore, minMoriIdx] = min(scoreList);

                        if bestScore < inf
                            seq = allSeqs{minMoriIdx};
                            dev = angList(minMoriIdx);
                        else
                            % Should not happen since we already checked undirected, but just in case
                            seq = 1; % fallback
                            dev = inf;
                        end

                        if is_newly_visited
                            % Calculate target fundamental region center based purely on the newly matched variant
                            variant_mori = allMori(minMoriIdx);
                            center_ori = currOri * variant_mori;

                            % Find all boundary segments containing nxtNode
                            idx1 = obj.edgeData.grainId(:,1) == nxtNode;
                            idx2 = obj.edgeData.grainId(:,2) == nxtNode;

                            % Project EBSD and meanEBSD orientations (Removed grainOri1/2 projection)
                            if any(idx1)
                                obj.edgeData.ebsdOri1(idx1)  = obj.edgeData.ebsdOri1(idx1).project2FundamentalRegion(center_ori);
                                obj.edgeData.meanEbsdOri1(idx1) = obj.edgeData.meanEbsdOri1(idx1).project2FundamentalRegion(center_ori);
                            end
                            if any(idx2)
                                obj.edgeData.ebsdOri2(idx2)  = obj.edgeData.ebsdOri2(idx2).project2FundamentalRegion(center_ori);
                                obj.edgeData.meanEbsdOri2(idx2) = obj.edgeData.meanEbsdOri2(idx2).project2FundamentalRegion(center_ori);
                            end

                            % Recompute misorientations for affected segments
                            idx_update = idx1 | idx2;
                            obj.edgeData.misorientation(idx_update) = inv(obj.edgeData.ebsdOri2(idx_update)) .* obj.edgeData.ebsdOri1(idx_update);
                        end

                        v_step_ori = currOri;
                        prevNode = currNode;
                        for step = 1:length(seq)
                            v_step = seq(step);
                            v_step_ori = v_step_ori * variants(v_step);
                            
                            if step == length(seq)
                                s_nodes(end+1, 1) = prevNode;
                                t_nodes(end+1, 1) = nxtNode;
                                variants_list{end+1, 1} = num2str(v_step);
                                angular_deviation_list(end+1, 1) = dev;
                            else
                                vNodeIdx = vNodeIdx + 1;
                                vNodeOris = [vNodeOris; v_step_ori];

                                coord1 = obj.grains(id2idx(currNode)).centroid;
                                coord2 = obj.grains(id2idx(nxtNode)).centroid;
                                alpha = step / length(seq);
                                vCoord = coord1 + alpha * (coord2 - coord1);

                                vNodeCoords(end+1, :) = vCoord;
                                vNodeGens{end+1, 1} = seq(1:step);

                                s_nodes(end+1, 1) = prevNode;
                                t_nodes(end+1, 1) = vNodeIdx;
                                variants_list{end+1, 1} = num2str(v_step);
                                angular_deviation_list(end+1, 1) = dev;

                                prevNode = vNodeIdx;
                            end
                        end
                    end
                end
                % end % (for current_max_gen)
            end

            numTotalNodes = vNodeIdx;

            finalG = graph();
            if numTotalNodes > 0
                finalG = addnode(finalG, numTotalNodes);
            end

            if ~isempty(s_nodes)
                EdgeTable = table([s_nodes, t_nodes], variants_list, angular_deviation_list, 'VariableNames', {'EndNodes', 'Variant', 'AngularDeviation'});
                finalG = addedge(finalG, EdgeTable);
            end

            obj.Graph = finalG;
            obj.untaggedGraph = G_undir;

            if ~isempty(vNodeCoords)
                obj.virtualNodes = table(vNodeCoords(:,1), vNodeCoords(:,2), vNodeGens, vNodeOris, 'VariableNames', {'x', 'y', 'Sequence', 'meanOrientation'});
            else
                obj.virtualNodes = table([], [], {}, orientation.empty(0, 1), 'VariableNames', {'x', 'y', 'Sequence', 'meanOrientation'});
            end
            
            % Give ALL nodes their mathematically exact meanOrientation directly on the Graph
            % We append in a loop to be 100% MTEX-safe.
                all_node_oris = orientation.empty(0, 1);
                for n = 1:numTotalNodes
                    if n <= numPhysNodes && id2idx(n) > 0
                        all_node_oris = [all_node_oris; obj.grains(id2idx(n)).meanOrientation(:)];
                    elseif n > numPhysNodes
                        vIdx = n - numPhysNodes;
                        all_node_oris = [all_node_oris; obj.virtualNodes.meanOrientation(vIdx)];
                    else
                        % If node doesn't exist, append a NaN orientation or a dummy one
                        all_node_oris = [all_node_oris; orientation.nan(obj.grains.CS)];
                    end
                end
                obj.Graph.Nodes.meanOrientation = all_node_oris;
        end

        function plotGraph(obj, varargin)
            p = inputParser;
            addOptional(p, 'type', 'tagged');
            parse(p, varargin{:});
            plotType = p.Results.type;

            if strcmpi(plotType, 'untagged')
                G = obj.untaggedGraph;
            elseif strcmpi(plotType, 'mst')
                if isempty(obj.mstGraph)
                    warning('twinNetwork:plotGraph:emptyMST', 'MST is empty. Compute the quotient graph first.');
                    return
                end
                G = obj.mstGraph;
            elseif strcmpi(plotType, 'quotient')
                if isempty(obj.QuotientGraph)
                    warning('twinNetwork:plotGraph:emptyQuotient', 'Quotient is empty. Compute the quotient graph first.');
                    return
                end
                G = obj.Graph;
            else
                G = obj.Graph;
            end

            if isempty(G) || numnodes(G) == 0
                warning('TwinNetwork:plot:emptyGraph', 'Graph is empty or non-existent.');
                return;
            end

            figure;
            plot(obj.grains, 'micronbar', 'off');
            hold on;

            numNodesG = numnodes(obj.Graph);
            X_G = zeros(numNodesG, 1);
            Y_G = zeros(numNodesG, 1);

            id2idx = zeros(max(obj.grains.id), 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);

            for i = 1:max(obj.grains.id)
                if id2idx(i) > 0
                    c = obj.grains(id2idx(i)).centroid;
                    X_G(i) = c.x;
                    Y_G(i) = c.y;
                else
                    X_G(i) = NaN;
                    Y_G(i) = NaN;
                end
            end

            if ~isempty(obj.virtualNodes) && size(obj.virtualNodes, 1) > 0
                vStart = max(obj.grains.id) + 1;
                X_G(vStart:end) = obj.virtualNodes.x;
                Y_G(vStart:end) = obj.virtualNodes.y;
            end

            numPlotNodes = numnodes(G);
            X_plot = X_G(1:numPlotNodes);
            Y_plot = Y_G(1:numPlotNodes);

            pPlot = plot(G, 'XData', X_plot, 'YData', Y_plot);
            pPlot.MarkerSize = 5;
            pPlot.NodeColor = 'k';

            if (strcmpi(plotType, 'tagged') || strcmpi(plotType, 'mst')) && ismember('AngularDeviation', G.Edges.Properties.VariableNames)
                pPlot.EdgeCData = G.Edges.AngularDeviation ./ degree;
                colormap(gca, 'jet');
                cb = colorbar;
                cb.Label.String = 'Angular Deviation (degrees)';
            else
                pPlot.EdgeColor = 'k';
            end
            pPlot.LineWidth = 1.5;

            if ismember('Variant', G.Edges.Properties.VariableNames)
                labeledge(pPlot, 1:numedges(G), G.Edges.Variant);
            end

            nodeLabels = cell(numPlotNodes, 1);
            for i = 1:numPlotNodes
                if (i <= max(obj.grains.id) && id2idx(i) > 0) || i > max(obj.grains.id)
                    nodeLabels{i} = num2str(i);
                else
                    nodeLabels{i} = '';
                end
            end
            pPlot.NodeLabel = nodeLabels;

            if strcmpi(plotType, 'quotient')
                cmap = lines();
                numCmap = size(cmap, 1);

                bins = conncomp(G);
                for c = 1:max(bins)
                    nodesInComp = find(bins == c);
                    groupsInComp = unique(obj.nodeGroups(nodesInComp));

                    for k = 1:length(groupsInComp)
                        grp = groupsInComp(k);
                        colorIdx = mod(k - 1, numCmap) + 1;
                        nodesInThisGroup = nodesInComp(obj.nodeGroups(nodesInComp) == grp);

                        highlight(pPlot, nodesInThisGroup, 'NodeColor', cmap(colorIdx, :), 'MarkerSize', 8, 'Marker', 's');
                    end
                end
            elseif strcmpi(plotType, 'tagged') && ~isempty(obj.virtualNodes) && size(obj.virtualNodes, 1) > 0
                vIdx = (max(obj.grains.id) + 1):numPlotNodes;
                highlight(pPlot, vIdx, 'NodeColor', 'r', 'Marker', 'o', 'MarkerSize', 8);
            end

            hold off;
        end
        
        function plot(obj, varargin)
            p = inputParser;
            addParameter(p, 'twinsOnly', true);
            parse(p, varargin{:});
            twinsOnly = p.Results.twinsOnly;
            
            if isempty(obj.groupGenerations)
                warning('TwinNetwork:plot:noParents', 'Parent identification not run. Please call identifyParents() first.');
                return;
            end
            
            figure;
            % Initialize white base for all grains
            numGrains = length(obj.grains);
            faceColors = ones(numGrains, 3); 
            
            % We will also collect labels to plot text
            labelCentersX = [];
            labelCentersY = [];
            labelTexts = {};
            
            cmap = lines();
            numCmap = size(cmap, 1);
            
            maxPhysId = max(obj.grains.id);
            id2idx = zeros(maxPhysId, 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);
            
            numNodesInG = length(obj.nodeGroups);
            maxPhysIdSafe = min(maxPhysId, numNodesInG);
            physNodeGroups = obj.nodeGroups(1:maxPhysIdSafe);
            
            toPlot = true(numGrains, 1);
            bins = conncomp(obj.Graph);
            compSizes = histcounts(bins, 1:max(bins)+1);
            
            for i = 1:maxPhysIdSafe
                idx = id2idx(i);
                if idx > 0
                    c = bins(i);
                    isIsolated = compSizes(c) <= 1;
                    
                    if twinsOnly && isIsolated
                        toPlot(idx) = false;
                        continue;
                    end
                    
                    if isIsolated
                        % Isolated grains stay default white, no text overlay
                        continue;
                    end
                    
                    g = physNodeGroups(i);
                    if g > 0 && g <= length(obj.groupGenerations) && ~isnan(obj.groupGenerations(g))
                        gen = obj.groupGenerations(g);
                        path = obj.groupVariantPaths{g};
                        
                        % Determine color
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
                            
                            if isempty(path)
                                lbl = 'V';
                            else
                                pathStr = sprintf('%d,', path);
                                lbl = sprintf('V%s', pathStr(1:end-1));
                            end
                        end
                        
                        c = obj.grains(idx).centroid;
                        labelCentersX(end+1) = c.x;
                        labelCentersY(end+1) = c.y;
                        labelTexts{end+1} = lbl;
                    end
                end
            end
            
            % Plot the grains
            if twinsOnly
                plot(obj.grains(toPlot), faceColors(toPlot, :), 'micronbar', 'off');
            else
                plot(obj.grains, faceColors, 'micronbar', 'off');
            end
            hold on;
            
            % Plot text
            if ~isempty(labelTexts)
                text(labelCentersX, labelCentersY, labelTexts, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 10, ...
                    'Color', 'k');
            end
            hold off;
        end
        
        function identifyParents(obj, varargin)
            p = inputParser;
            addParameter(p, 'criterion', 'grainArea');
            parse(p, varargin{:});
            criterion = p.Results.criterion;
            
            if isempty(obj.QuotientGraph)
                warning('TwinNetwork:identifyParents:emptyQuotient', 'QuotientGraph is empty. Compute quotient graph first.');
                return
            end
            
            Q = obj.QuotientGraph;
            if isempty(Q) || numnodes(Q) == 0
                return;
            end
            
            numGroups = numnodes(Q);
            obj.groupGenerations = nan(numGroups, 1);
            obj.groupVariantPaths = cell(numGroups, 1);
            
            bins = conncomp(Q);
            numComps = max(bins);
            
            maxPhysId = max(obj.grains.id);
            id2idx = zeros(maxPhysId, 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);
            
            numNodesInG = length(obj.nodeGroups);
            maxPhysIdSafe = min(maxPhysId, numNodesInG);
            physNodeGroups = obj.nodeGroups(1:maxPhysIdSafe);
            
            for c = 1:numComps
                groupsInComp = find(bins == c);
                
                if strcmpi(criterion, 'grainArea')
                    bestArea = -1;
                    parentGroup = groupsInComp(1);
                    
                    for i = 1:length(groupsInComp)
                        g = groupsInComp(i);
                        physIds = find(physNodeGroups == g);
                        
                        validPhysIds = physIds(physIds <= length(id2idx));
                        grainIndices = id2idx(validPhysIds);
                        validIndices = grainIndices(grainIndices > 0);
                        
                        if isempty(validIndices)
                            areaSum = 0;
                        else
                            areaSum = sum(obj.grains(validIndices).area);
                        end
                        if areaSum > bestArea
                            bestArea = areaSum;
                            parentGroup = g;
                        end
                    end
                else
                    warning('TwinNetwork:identifyParents:unimplementedCriterion', 'Criterion %s not implemented. Falling back to arbitrary parent.', criterion);
                    parentGroup = groupsInComp(1);
                end
                
                % BFS from parentGroup
                visited = false(numGroups, 1);
                queue = parentGroup;
                visited(parentGroup) = true;
                
                obj.groupGenerations(parentGroup) = 0;
                obj.groupVariantPaths{parentGroup} = [];
                
                while ~isempty(queue)
                    currGroup = queue(1);
                    queue(1) = [];
                    
                    currGen = obj.groupGenerations(currGroup);
                    currPath = obj.groupVariantPaths{currGroup};
                    
                    nbrs = neighbors(Q, currGroup);
                    for n = 1:length(nbrs)
                        nxtGroup = nbrs(n);
                        if ~visited(nxtGroup)
                            visited(nxtGroup) = true;
                            
                            edgeIdx = findedge(Q, currGroup, nxtGroup);
                            edgeIdx = edgeIdx(1); % In case of multiple edges with different variants
                            variantStr = Q.Edges.Variant{edgeIdx};
                            variantNum = str2double(variantStr);
                            
                            obj.groupGenerations(nxtGroup) = currGen + 1;
                            obj.groupVariantPaths{nxtGroup} = [currPath, variantNum];
                            
                            queue(end+1) = nxtGroup;
                        end
                    end
                end
            end
        end
        function SF = calcSF(obj, sigma)
            if ~isa(sigma, 'stressTensor')
                error('TwinNetwork:calcSF:invalidInput', 'Input must be a stressTensor object.');
            end
            
            warning('TwinNetwork:calcSF:primaryOnly', 'This method currently only supports primary twins.');
            
            if isempty(obj.groupGenerations)
                error('TwinNetwork:calcSF:noParents', 'Parent identification not run. Please call identifyParents() first.');
            end
            
            tS_sym = obj.tS.symmetrise;
            variants = tS_sym.parentTwinMisorientation;
            if variants.isImproper
                variants = -variants;
            end
            
            numGrains = length(obj.grains);
            SF = nan(numGrains, 1);
            
            maxPhysId = max(obj.grains.id);
            id2idx = zeros(maxPhysId, 1);
            id2idx(obj.grains.id) = 1:numGrains;
            
            numNodesInG = length(obj.nodeGroups);
            maxPhysIdSafe = min(maxPhysId, numNodesInG);
            physNodeGroups = obj.nodeGroups(1:maxPhysIdSafe);
            
            for i = 1:maxPhysIdSafe
                idx = id2idx(i);
                if idx > 0
                    g = physNodeGroups(i);
                    if g > 0 && g <= length(obj.groupGenerations) && ~isnan(obj.groupGenerations(g))
                        gen = obj.groupGenerations(g);
                        if gen == 1
                            % Primary twin
                            path = obj.groupVariantPaths{g};
                            V = path(1);
                            
                            % Find corrected grain orientation natively from Graph
                            correctedTwinOri = obj.Graph.Nodes.meanOrientation(idx);
                            
                            variant_mori = variants(V);
                            
                            % Virtual parent
                            virtualParentOri = correctedTwinOri * inv(variant_mori);
                            
                            % Transform stress tensor to virtual parent crystal coordinates
                            sigma_crystal = inv(virtualParentOri) * sigma;
                            
                            % Calculate SF for all variants
                            sf_all = SchmidFactor(tS_sym, sigma_crystal);
                            
                            % Extract the SF for this specific variant
                            SF(idx) = sf_all(V);
                        end
                    end
                end
            end
        end
        function [Q, groupIDs] = computeQuotientGraph(obj, varargin)
            p = inputParser;
            addParameter(p, 'stable', false);
            parse(p, varargin{:});
            stable_val = p.Results.stable;
            if ischar(stable_val) || isstring(stable_val)
                stable = strcmpi(stable_val, 'all') || strcmpi(stable_val, 'true');
            else
                stable = logical(stable_val);
            end

            G = obj.Graph;
            search_G = G;
            if isempty(G) || numnodes(G) == 0
                Q = graph();
                groupIDs = [];
                obj.QuotientGraph = Q;
                obj.nodeGroups = groupIDs;
                return;
            end

            numNodes = numnodes(G);
            maxPhysId = max(obj.grains.id);
            id2idx = zeros(maxPhysId, 1);
            id2idx(obj.grains.id) = 1:length(obj.grains);

            tS_sym = obj.tS.symmetrise;
            variants = tS_sym.parentTwinMisorientation;
            if variants.isImproper, variants = -variants; end
            outlierThresh = 15 * degree;

            bins = conncomp(G);
            numComps = max(bins);

            % Step 1: Initial MST construction
            mstS = [];
            mstT = [];
            mstVar = {};
            mstDev = [];
            hasDeviation = ismember('AngularDeviation', G.Edges.Properties.VariableNames);

            for c = 1:numComps
                compNodeIdx = find(bins == c);
                compNodeIdx = compNodeIdx(:);
                if length(compNodeIdx) == 1
                    continue;
                end

                subG = subgraph(G, compNodeIdx);

                if hasDeviation
                    subG.Edges.Weight = subG.Edges.AngularDeviation;
                else
                    warning('TwinNetwork:computeQuotientGraph:noWeight', ...
                        'No ''AngularDeviation'' variable found in graph edges. Using default connectivity.');
                    subG.Edges.Weight = ones(numedges(subG), 1);
                end

                virtualLocalIdx = find(compNodeIdx > maxPhysId); 
                if ~isempty(virtualLocalIdx)
                    endNodes = subG.Edges.EndNodes;
                    virtualEdges = find(ismember(endNodes(:, 1), virtualLocalIdx) | ...
                                        ismember(endNodes(:, 2), virtualLocalIdx));
                    virtualPenalty = 1000; 
                    subG.Edges.Weight(virtualEdges) = subG.Edges.Weight(virtualEdges) + virtualPenalty;
                end
                
                [subMST, ~] = minspantree(subG);

                if numedges(subMST) > 0
                    localS = subMST.Edges.EndNodes(:, 1);
                    localT = subMST.Edges.EndNodes(:, 2);
                    mstS = [mstS; compNodeIdx(localS)];
                    mstT = [mstT; compNodeIdx(localT)];
                    mstVar = [mstVar; subMST.Edges.Variant];
                    if hasDeviation
                        mstDev = [mstDev; subMST.Edges.AngularDeviation];
                    end
                end
            end

            mstG = graph();
            mstG = addnode(mstG, numNodes);
            if ~isempty(mstS)
                if hasDeviation
                    EdgeT = table([mstS, mstT], mstVar, mstDev, 'VariableNames', {'EndNodes', 'Variant', 'AngularDeviation'});
                else
                    EdgeT = table([mstS, mstT], mstVar, 'VariableNames', {'EndNodes', 'Variant'});
                end
                mstG = addedge(mstG, EdgeT);
            end

            % Step 2: Global Refinement Loop
            maxIters = stable * 5 + 1;
            nodeStates = cell(numNodes, 1);
            groupIDs = zeros(numNodes, 1);
            
            for iter = 1:maxIters
                parentNode = zeros(numNodes, 1);
                
                for c = 1:numComps
                    compNodeIdx = find(bins == c);
                    compNodeIdx = compNodeIdx(:);
                    if length(compNodeIdx) == 1
                        nodeStates{compNodeIdx(1)} = sprintf('C%d_ROOT', c);
                        continue;
                    end
                    
                    physNodes = compNodeIdx(compNodeIdx <= maxPhysId);
                    validNodes = physNodes(id2idx(physNodes) > 0);
                    if ~isempty(validNodes)
                        areas = obj.grains(id2idx(validNodes)).area;
                        [~, maxAreaIdx] = max(areas);
                        rootNode = validNodes(maxAreaIdx);
                    else
                        rootNode = compNodeIdx(1);
                    end
                    
                    nodeStates{rootNode} = sprintf('C%d_ROOT', c);
                    
                    subStates = cell(numNodes, 1);
                    subStates{rootNode} = [];
                    
                    visited = false(numNodes, 1);
                    visited(rootNode) = true;
                    
                    queue = rootNode;
                    
                    while ~isempty(queue)
                        currNode = queue(1);
                        queue(1) = [];
                        
                        currState = subStates{currNode};
                        
                        nbrs = neighbors(mstG, currNode);
                        for n = 1:length(nbrs)
                            nxtNode = nbrs(n);
                            if bins(nxtNode) == c && ~visited(nxtNode)
                                visited(nxtNode) = true;
                                parentNode(nxtNode) = currNode;
                                
                                edgeIdx = findedge(mstG, currNode, nxtNode);
                                v_num = str2double(mstG.Edges.Variant{edgeIdx});
                                
                                nxtState = currState;
                                if ~isempty(nxtState) && nxtState(end) == v_num
                                    nxtState(end) = [];
                                else
                                    nxtState(end+1) = v_num;
                                end
                                
                                subStates{nxtNode} = nxtState;
                                
                                if isempty(nxtState)
                                    nodeStates{nxtNode} = sprintf('C%d_ROOT', c);
                                else
                                    str_s = sprintf('%d_', nxtState);
                                    nodeStates{nxtNode} = sprintf('C%d_%s', c, str_s(1:end-1));
                                end
                                
                                queue(end+1) = nxtNode;
                            end
                        end
                    end
                end
                
                [unqStates, ~, groupIDs] = unique(string(nodeStates));
                
                changed = false;
                if stable
                    % Calculate depths for hierarchy
                    depths = zeros(maxPhysId, 1);
                    for v = 1:maxPhysId
                        if id2idx(v) > 0
                            d = 0;
                            curr = v;
                            while curr ~= 0 && parentNode(curr) ~= 0
                                curr = parentNode(curr);
                                d = d + 1;
                            end
                            depths(v) = d;
                        end
                    end
                    
                    outliers_loop_changed = true;
                    while outliers_loop_changed
                        outliers_loop_changed = false;
                        
                        % Recompute groups
                        [unqStates, ~, groupIDs] = unique(string(nodeStates));
                        numGroups = length(unqStates);
                        groupRefOris = cell(numGroups, 1);
                        for g = 1:numGroups
                            g_nodes = find(groupIDs == g);
                            phys_g_nodes = g_nodes(g_nodes <= maxPhysId);
                            valid_phys = phys_g_nodes(id2idx(phys_g_nodes) > 0);
                            if ~isempty(valid_phys)
                                oris = search_G.Nodes.meanOrientation(valid_phys);
                                grp_areas = obj.grains(id2idx(valid_phys)).area;
                                
                                num_nodes = length(oris);
                                if num_nodes == 1
                                    bestIdx = 1;
                                    robust_ref_ori = oris;
                                else
                                    scores = zeros(num_nodes, 1);
                                    for i = 1:num_nodes
                                        % Find how many grains are within a tight cluster (5 degrees)
                                        close_idx = angle(oris(i), oris) < (5 * degree);
                                        % Score by cluster density, use normalized area to break ties
                                        scores(i) = sum(close_idx) + (grp_areas(i) / max(grp_areas)) * 0.1;
                                    end
                                    [~, bestIdx] = max(scores);
                                    
                                    % Calculate the true robust mean of this dense inlier cluster
                                    inlier_idx = angle(oris(bestIdx), oris) < (5 * degree);
                                    robust_ref_ori = mean(oris(inlier_idx));
                                end
                                
                                groupRefOris{g} = robust_ref_ori;
                            end
                        end
                        
                        % Collect all outliers
                        all_outliers = [];
                        for g = 1:numGroups
                            g_nodes = find(groupIDs == g);
                            phys_g_nodes = g_nodes(g_nodes <= maxPhysId);
                            valid_phys = phys_g_nodes(id2idx(phys_g_nodes) > 0);
                            
                            if length(valid_phys) <= 1
                                continue;
                            end
                            
                            ref_ori = groupRefOris{g};
                            oris = search_G.Nodes.meanOrientation(valid_phys);
                            devs = angle(oris, ref_ori);
                            outlier_local_idx = find(devs > outlierThresh);
                            
                            for k = 1:length(outlier_local_idx)
                                v_global = valid_phys(outlier_local_idx(k));
                                if parentNode(v_global) == 0
                                    continue;
                                end
                                all_outliers = [all_outliers; v_global, devs(outlier_local_idx(k))];
                            end
                        end
                        
                        if ~isempty(all_outliers)
                            % Sort by depth (ascending) to fix ancestors first
                            outlier_depths = depths(all_outliers(:, 1));
                            [~, sort_idx] = sort(outlier_depths);
                            all_outliers = all_outliers(sort_idx, :);
                            
                            for k = 1:size(all_outliers, 1)
                                v_global = all_outliers(k, 1);
                                current_dev = all_outliers(k, 2);
                                v_ori = search_G.Nodes.meanOrientation(v_global);
                                
                                v_group_id = groupIDs(v_global);
                                v_children_count = sum(parentNode == v_global);
                                g_nodes = find(groupIDs == v_group_id);
                                phys_g_nodes = g_nodes(g_nodes <= maxPhysId);
                                valid_phys = phys_g_nodes(id2idx(phys_g_nodes) > 0);
                                if ~isempty(valid_phys)
                                    grp_areas = obj.grains(id2idx(valid_phys)).area;
                                    [~, sortIdx] = sort(grp_areas, 'descend');
                                    sorted_nodes = valid_phys(sortIdx);
                                    nodes_str = sprintf('%d,', sorted_nodes);
                                    nodes_str = nodes_str(1:end-1);
                                else
                                    nodes_str = 'none';
                                end
                                
                                warning('TwinNetwork:computeQuotientGraph:outlierFound', ...
                                    '\n=== Outlier grain %d found ===\nHas %d children node(s).\nCurrent nodes in the same group: %s (listed largest to smallest area).\nAng deviation from robust group center: %.1f deg', ...
                                    v_global, v_children_count, nodes_str, current_dev/degree);
                                
                                % FUSED OUTLIER SEARCH LOOP
                                % Collect subtree nodes
                                subtree_nodes = v_global;
                                sq = v_global;
                                while ~isempty(sq)
                                    curr = sq(1); sq(1) = [];
                                    chs = find(parentNode == curr);
                                    if isempty(subtree_nodes)
                                        subtree_nodes = chs;
                                    else
                                        subtree_nodes = [subtree_nodes; chs];
                                    end
                                    
                                    if isempty(sq)
                                        sq = chs;
                                    else
                                        sq = [sq; chs];
                                    end
                                end
                                
                                v_comp = bins(v_global);
                                compNodes = find(bins == v_comp);
                                comp_root = 0;
                                for ci = 1:length(compNodes)
                                    if parentNode(compNodes(ci)) == 0
                                        comp_root = compNodes(ci);
                                        break;
                                    end
                                end
                                
                                fallback_success = false;
                                if comp_root ~= 0
                                    % Do NOT reset search_G = G because search_G contains the progressively updated Orientations!
                                    if ismember('AngularDeviation', search_G.Edges.Properties.VariableNames)
                                        search_G.Edges.Weight = search_G.Edges.AngularDeviation;
                                    else
                                        search_G.Edges.Weight = ones(numedges(search_G), 1);
                                    end
                                    
                                    p_global = parentNode(v_global);
                                    edgeIdx1 = findedge(search_G, p_global, v_global);
                                    edgeIdx2 = findedge(search_G, v_global, p_global);
                                    if edgeIdx1 > 0, search_G.Edges.Weight(edgeIdx1) = inf; end
                                    if edgeIdx2 > 0, search_G.Edges.Weight(edgeIdx2) = inf; end
                                    
                                    warning('TwinNetwork:computeQuotientGraph:fusedSearch', '  -> Searching for optimal reconnection for node %d (and its subtree).', v_global);
                                    
                                    for fallback_iter = 1:10
                                        [P, D] = shortestpathtree(search_G, comp_root, 'Method', 'positive');
                                        subtree_dists = D(subtree_nodes);
                                        [min_dist, min_idx] = min(subtree_dists);
                                        
                                        if isinf(min_dist)
                                            warning('TwinNetwork:computeQuotientGraph:fusedSearch', '     * Attempt %d: No path found to the subtree.', fallback_iter);
                                            break;
                                        end
                                        
                                        target_node = subtree_nodes(min_idx);
                                        % Use shortestpath to get the actual path
                                        path_nodes = shortestpath(search_G, comp_root, target_node, 'Method', 'positive');
                                        
                                        if length(path_nodes) < 2
                                            break;
                                        end
                                        new_anchor = path_nodes(end-1);
                                        edgeIdx = findedge(search_G, new_anchor, target_node);
                                        if edgeIdx == 0, edgeIdx = findedge(search_G, target_node, new_anchor); end
                                        
                                        % 1. Create a temporary MST
                                        temp_mstG = mstG;
                                        temp_parentNode = parentNode;
                                        temp_mstG = rmedge(temp_mstG, p_global, v_global);
                                        
                                        % Invert subtree connections from target_node up to v_global
                                        if target_node ~= v_global
                                            curr_inv = target_node;
                                            prev_inv = 0;
                                            while curr_inv ~= 0
                                                next_inv = parentNode(curr_inv);
                                                temp_parentNode(curr_inv) = prev_inv;
                                                if curr_inv == v_global
                                                    break;
                                                end
                                                prev_inv = curr_inv;
                                                curr_inv = next_inv;
                                            end
                                        end
                                        
                                        % Finally, connect target_node to new_anchor
                                        temp_parentNode(target_node) = new_anchor;
                                        
                                        % 2. Get effective orientation of new_anchor natively from the Graph Nodes
                                        eff_center = search_G.Nodes.meanOrientation(new_anchor);
                                        
                                        % 3. Evaluate the connection to target_node using actual variant
                                        target_ori = search_G.Nodes.meanOrientation(target_node);
                                        angList = angle(target_ori, eff_center * variants(:)');
                                        [conn_dev, best_v_num] = min(angList);
                                        target_center = eff_center * variants(best_v_num);
                                        
                                        % Add the new edge to temp_mstG using best_v_num
                                        hasDev = ismember('AngularDeviation', temp_mstG.Edges.Properties.VariableNames);
                                        if hasDev
                                            newEdgeT = table([new_anchor, target_node], {num2str(best_v_num)}, conn_dev, 'VariableNames', {'EndNodes', 'Variant', 'AngularDeviation'});
                                        else
                                            newEdgeT = table([new_anchor, target_node], {num2str(best_v_num)}, 'VariableNames', {'EndNodes', 'Variant'});
                                        end
                                        temp_mstG = addedge(temp_mstG, newEdgeT);
                                        
                                        % 4. Propagate down the subtree and find max deviation
                                        sq = target_node;
                                        sq_centers = target_center;
                                        max_subtree_dev = 0;
                                        worst_node = -1;
                                        
                                        while ~isempty(sq)
                                            curr = sq(1); sq(1) = [];
                                            curr_c = sq_centers(1); sq_centers(1) = [];
                                            
                                            curr_ori = search_G.Nodes.meanOrientation(curr);
                                            dev = angle(curr_ori, curr_c);
                                            if dev > max_subtree_dev
                                                max_subtree_dev = dev;
                                                worst_node = curr;
                                            end
                                            
                                            chs = find(temp_parentNode == curr);
                                            for c_idx = 1:length(chs)
                                                child = chs(c_idx);
                                                eIdx = findedge(temp_mstG, curr, child);
                                                if eIdx == 0, eIdx = findedge(temp_mstG, child, curr); end
                                                v_num = str2double(temp_mstG.Edges.Variant{eIdx});
                                                
                                                if isempty(sq)
                                                    sq = child;
                                                else
                                                    sq = [sq; child];
                                                end
                                                
                                                new_c = curr_c * variants(v_num);
                                                if isempty(sq_centers)
                                                    sq_centers = new_c;
                                                else
                                                    sq_centers = [sq_centers; new_c];
                                                end
                                            end
                                        end
                                        % Trace full node path from outlier node (v_global) to root in temp_parentNode
                                        curr_trace = v_global;
                                        full_node_path = curr_trace;
                                        while curr_trace ~= 0
                                            p_node = temp_parentNode(curr_trace);
                                            if p_node == 0, break; end
                                            full_node_path = [p_node, full_node_path];
                                            curr_trace = p_node;
                                        end
                                        node_path_str = sprintf('%d-', full_node_path);
                                        node_path_str = node_path_str(1:end-1);
                                        
                                        % Propagate down the path to correctly evaluate variants in the Fundamental Zone
                                        full_var_path = [];
                                        root_node = full_node_path(1);
                                        eff_center = search_G.Nodes.meanOrientation(root_node);
                                        
                                        for path_idx = 1:length(full_node_path)-1
                                            n1 = full_node_path(path_idx);
                                            n2 = full_node_path(path_idx+1);
                                            
                                            target_ori = search_G.Nodes.meanOrientation(n2);
                                            
                                            eIdx_trace = findedge(temp_mstG, n1, n2);
                                            if eIdx_trace == 0, eIdx_trace = findedge(temp_mstG, n2, n1); end
                                            if eIdx_trace > 0
                                                orig_v_num = str2double(temp_mstG.Edges.Variant{eIdx_trace});
                                            else
                                                eIdx_search = findedge(search_G, n1, n2);
                                                if eIdx_search == 0, eIdx_search = findedge(search_G, n2, n1); end
                                                orig_v_num = str2double(search_G.Edges.Variant{eIdx_search});
                                            end
                                            
                                            angList = angle(target_ori, eff_center * variants(:)');
                                            [~, best_v_num] = min(angList);
                                            
                                            % Project the target into the fundamental region exactly as in line 283
                                                center_ori = eff_center * variants(best_v_num);
                                                target_ori_proj = target_ori.project2FundamentalRegion(center_ori);
                                                
                                                disp(['--- Tracing edge ' num2str(n1) ' -> ' num2str(n2) ' ---']);
                                                disp(['Original Variant in graph: ' num2str(orig_v_num)]);
                                                disp(['New Recalculated Variant: ' num2str(best_v_num)]);
                                                disp('Euler angles BEFORE recasting to FZ (deg):');
                                                try
                                                    disp(target_ori.Euler ./ degree);
                                                catch
                                                    disp(target_ori);
                                                end
                                                disp('Euler angles AFTER recasting to FZ (deg):');
                                                try
                                                    disp(target_ori_proj.Euler ./ degree);
                                                catch
                                                    disp(target_ori_proj);
                                                end
                                                
                                                eff_center = target_ori_proj;
                                            
                                            full_var_path(end+1) = best_v_num;
                                        end
                                        
                                        if isempty(full_var_path)
                                            full_var_str = '[]';
                                        else
                                            full_var_str = ['[', num2str(full_var_path), ']'];
                                        end
                                        
                                        % Evaluate reduced group state for v_global
                                        v_global_state = [];
                                        for v_idx_reduce = 1:length(full_var_path)
                                            if ~isempty(v_global_state) && v_global_state(end) == full_var_path(v_idx_reduce)
                                                v_global_state(end) = [];
                                            else
                                                v_global_state(end+1) = full_var_path(v_idx_reduce);
                                            end
                                        end
                                        
                                        c_idx = bins(v_global);
                                        if isempty(v_global_state)
                                            v_global_state_str = sprintf('C%d_ROOT', c_idx);
                                            var_str = '[]';
                                        else
                                            str_s = sprintf('%d_', v_global_state);
                                            v_global_state_str = sprintf('C%d_%s', c_idx, str_s(1:end-1));
                                            var_str = ['[', num2str(v_global_state), ']'];
                                        end
                                        
                                        % Determine the theoretical global expected orientation
                                        expected_global_ori = search_G.Nodes.meanOrientation(full_node_path(1));
                                        for v = v_global_state
                                            expected_global_ori = expected_global_ori * variants(v);
                                        end
                                        
                                        % Find if this group state already exists and has an empirical center
                                        state_idx = find(unqStates == string(v_global_state_str));
                                        if ~isempty(state_idx) && ~isempty(groupRefOris{state_idx})
                                            ref_ori = groupRefOris{state_idx};
                                            is_empirical = true;
                                        else
                                            ref_ori = expected_global_ori;
                                            is_empirical = false;
                                        end
                                        
                                        global_dev = angle(search_G.Nodes.meanOrientation(v_global), ref_ori);
                                        
                                        warning('TwinNetwork:computeQuotientGraph:fusedSearch', '     * Attempt %d: Found path via %d->%d. Subtree max dev %.1f deg (node %d). Global dev %.1f deg (Empirical Ref: %d). Temp group: %s (reduced group state: %s, full path: %s).', fallback_iter, new_anchor, target_node, max_subtree_dev/degree, worst_node, global_dev/degree, is_empirical, v_global_state_str, var_str, full_var_str);
                                        
                                        if max_subtree_dev < outlierThresh && global_dev < outlierThresh
                                            fallback_success = true;
                                            best_u_global = new_anchor;
                                            best_v_num_final = best_v_num;
                                            best_u_eff_ori = eff_center;
                                            final_target_node = target_node;
                                            final_max_dev = max_subtree_dev;
                                            final_state_str = v_global_state_str;
                                            final_var_str = var_str;
                                            final_full_var_str = full_var_str;
                                            final_root_node = full_node_path(1);
                                            final_node_path_str = node_path_str;
                                            mstG = temp_mstG;
                                            parentNode = temp_parentNode;
                                            break;
                                        else
                                            search_G.Edges.Weight(edgeIdx) = inf;
                                            warning('TwinNetwork:computeQuotientGraph:fusedSearch', '     * Path invalid. Masking edge %d-%d and retrying.', new_anchor, target_node);
                                        end
                                    end
                                end
                                
                                    if ~fallback_success
                                        warning('TwinNetwork:computeQuotientGraph:unfixableOutlier', ...
                                            'Outlier grain %d (dev %.1f deg) could not be reassigned. No valid alternative path found.', ...
                                            v_global, current_dev/degree);
                                    else
                                        warning('TwinNetwork:computeQuotientGraph:outlierCheck', ...
                                            'Outlier grain %d (old max dev %.1f deg) reassigned to neighbor %d (new max dev %.1f deg) via variant %d. Now in temp group: %s (reduced group state: %s, full variant path: %s, root: %d, nodes: %s).', ...
                                            v_global, current_dev/degree, best_u_global, final_max_dev/degree, best_v_num_final, final_state_str, final_var_str, final_full_var_str, final_root_node, final_node_path_str);
                                            
                                        v_center_ori = best_u_eff_ori * variants(best_v_num_final);
                                        new_v_state = subStates{best_u_global};
                                        if ~isempty(new_v_state) && new_v_state(end) == best_v_num_final
                                            new_v_state(end) = [];
                                        else
                                            new_v_state(end+1) = best_v_num_final;
                                        end
                                        
                                        queue_nodes = final_target_node;
                                        queue_centers = v_center_ori;
                                        queue_states = {new_v_state};
                                        
                                        % Final downward FZ propagation
                                        while ~isempty(queue_nodes)
                                            curr = queue_nodes(1); queue_nodes(1) = [];
                                            curr_center = queue_centers(1); queue_centers(1) = [];
                                            curr_state = queue_states{1}; queue_states(1) = [];
                                            
                                            subStates{curr} = curr_state;
                                            if isempty(curr_state)
                                                nodeStates{curr} = sprintf('C%d_ROOT', bins(curr));
                                            else
                                                str_s = sprintf('%d_', curr_state);
                                                nodeStates{curr} = sprintf('C%d_%s', bins(curr), str_s(1:end-1));
                                            end
                                            
                                            % Propagate FZ to the graph node natively!
                                            search_G.Nodes.meanOrientation(curr) = search_G.Nodes.meanOrientation(curr).project2FundamentalRegion(curr_center);
                                            curr_ori = search_G.Nodes.meanOrientation(curr);
                                            
                                            % FZ Propagate EBSD data ONLY IF it's a physical node
                                            if curr <= maxPhysId && id2idx(curr) > 0
                                                idx1 = obj.edgeData.grainId(:,1) == curr;
                                                idx2 = obj.edgeData.grainId(:,2) == curr;
                                                if any(idx1)
                                                    obj.edgeData.ebsdOri1(idx1)  = obj.edgeData.ebsdOri1(idx1).project2FundamentalRegion(curr_center);
                                                    obj.edgeData.meanEbsdOri1(idx1) = obj.edgeData.meanEbsdOri1(idx1).project2FundamentalRegion(curr_center);
                                                end
                                                if any(idx2)
                                                    obj.edgeData.ebsdOri2(idx2)  = obj.edgeData.ebsdOri2(idx2).project2FundamentalRegion(curr_center);
                                                    obj.edgeData.meanEbsdOri2(idx2) = obj.edgeData.meanEbsdOri2(idx2).project2FundamentalRegion(curr_center);
                                                end
                                                idx_update = idx1 | idx2;
                                                obj.edgeData.misorientation(idx_update) = inv(obj.edgeData.ebsdOri2(idx_update)) .* obj.edgeData.ebsdOri1(idx_update);
                                            end
                                            
                                            children = find(parentNode == curr);
                                            for c_idx_child = 1:length(children)
                                                child = children(c_idx_child);
                                                
                                                edgeIdx = findedge(mstG, curr, child);
                                                if edgeIdx == 0
                                                    edgeIdx = findedge(mstG, child, curr);
                                                end
                                                
                                                % Uniform variant assignment for physical and virtual nodes!
                                                child_ori_old = search_G.Nodes.meanOrientation(child);
                                                angList = angle(child_ori_old, curr_ori * variants(:)');
                                                [~, new_v_num_child] = min(angList);
                                                
                                                mstG.Edges.Variant{edgeIdx} = num2str(new_v_num_child);
                                                child_center = curr_ori * variants(new_v_num_child);
                                                
                                                next_child_state = curr_state;
                                                if ~isempty(next_child_state) && next_child_state(end) == new_v_num_child
                                                    next_child_state(end) = [];
                                                else
                                                    next_child_state(end+1) = new_v_num_child;
                                                end
                                                
                                                if isempty(queue_nodes)
                                                    queue_nodes = child;
                                                else
                                                    queue_nodes = [queue_nodes; child];
                                                end
                                                
                                                if isempty(queue_centers)
                                                    queue_centers = child_center;
                                                else
                                                    queue_centers = [queue_centers; child_center];
                                                end
                                                
                                                if isempty(queue_states)
                                                    queue_states = {next_child_state};
                                                else
                                                    queue_states = [queue_states; {next_child_state}];
                                                end
                                            end
                                        end
                                        
                                        outliers_loop_changed = true;
                                        break; % Restart the loop to recompute groupRefOris and sort remaining outliers!
                                    end
                            end
                        end
                        changed = outliers_loop_changed; % ensure parent scope knows if anything changed
                    end
                end
                
                if ~changed
                    break;
                end
            end

            % Step 3: Build the final Quotient Graph
            numGroups = length(unqStates);
            Q = graph();
            Q = addnode(Q, numGroups);

            if numedges(G) > 0
                u = groupIDs(G.Edges.EndNodes(:, 1));
                v = groupIDs(G.Edges.EndNodes(:, 2));
                vars = G.Edges.Variant;

                qEdges = [u, v];
                qEdges = sort(qEdges, 2);

                strKeys = cell(size(qEdges, 1), 1);
                for i = 1:size(qEdges, 1)
                    strKeys{i} = sprintf('%d_%d_%s', qEdges(i,1), qEdges(i,2), vars{i});
                end

                [~, unqIdx] = unique(strKeys);

                validEdges = qEdges(unqIdx, 1) ~= qEdges(unqIdx, 2);

                EdgeT = table(qEdges(unqIdx(validEdges), :), vars(unqIdx(validEdges)), ...
                    'VariableNames', {'EndNodes', 'Variant'});

                if ~isempty(EdgeT)
                    Q = addedge(Q, EdgeT);
                end
            end

            obj.Graph = search_G;
            obj.QuotientGraph = Q;
            obj.nodeGroups = groupIDs;
            obj.mstGraph = mstG;
        end
        function exportGraphData(obj, filename)
            % Default filename if none is provided
            if nargin < 2
                filename = 'twinGraphData.json';
            end

            G = obj.Graph;
            if isempty(G) || numnodes(G) == 0
                warning('Graph is empty. Nothing to export.');
                return;
            end

            % 1. Extract Node Data
            numNodes = numnodes(G);
            maxPhysId = max(obj.grains.id);

            % Initialize a struct array for nodes
            nodeData = struct('ID', cell(numNodes, 1), 'IsVirtual', cell(numNodes, 1));

            for i = 1:numNodes
                nodeData(i).ID = i;
                % Any node ID greater than the max physical grain ID is a virtual sequence node
                nodeData(i).IsVirtual = (i > maxPhysId);
            end

            % 2. Extract Edge Data
            numEdgesG = numedges(G);
            edgeData = struct('Source', cell(numEdgesG, 1), 'Target', cell(numEdgesG, 1));

            hasVariant = ismember('Variant', G.Edges.Properties.VariableNames);
            hasDeviation = ismember('AngularDeviation', G.Edges.Properties.VariableNames);

            for i = 1:numEdgesG
                edgeData(i).Source = G.Edges.EndNodes(i, 1);
                edgeData(i).Target = G.Edges.EndNodes(i, 2);

                if hasVariant
                    edgeData(i).Variant = G.Edges.Variant{i};
                end
                if hasDeviation
                    edgeData(i).AngularDeviation = G.Edges.AngularDeviation(i);
                end
            end

            % 3. Package and Export
            exportStruct = struct('Nodes', nodeData, 'Edges', edgeData);

            % Convert to JSON string (PrettyPrint makes it highly readable)
            jsonStr = jsonencode(exportStruct, 'PrettyPrint', true);

            % Write to file
            fid = fopen(filename, 'w');
            if fid == -1
                error('Could not open file for writing.');
            end
            fprintf(fid, '%s', jsonStr);
            fclose(fid);

            fprintf('Graph topology successfully exported to %s\n', filename);
        end
    end

end
