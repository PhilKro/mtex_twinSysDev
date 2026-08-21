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

            % Mean Orientations of grains on both sides
            obj.edgeData.grainOri1 = obj.grains(grains.id2ind(gB_grainId(:,1))).meanOrientation;
            obj.edgeData.grainOri2 = obj.grains(grains.id2ind(gB_grainId(:,2))).meanOrientation;

            % Orientations of EBSD pixels on both sides
            obj.edgeData.ebsdOri1 = obj.ebsd(ebsd.id2ind(gB_ebsdId(:,1))).orientations;
            obj.edgeData.ebsdOri2 = obj.ebsd(ebsd.id2ind(gB_ebsdId(:,2))).orientations;

            % Compute meanEbsdOri grouped by boundaries
            [~, ~, dirTags] = unique(gB_grainId, 'rows');
            meanE1 = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(obj.edgeData.ebsdOri1(idx))});
            meanE2 = accumarray(dirTags, (1:length(dirTags))', [], @(idx) {mean(obj.edgeData.ebsdOri2(idx))});

            meanE1_array = [meanE1{:}];
            meanE2_array = [meanE2{:}];

            obj.edgeData.meanEbsdOri1 = meanE1_array(dirTags);
            obj.edgeData.meanEbsdOri1 = obj.edgeData.meanEbsdOri1(:).project2FundamentalRegion(obj.edgeData.grainOri1);

            obj.edgeData.meanEbsdOri2 = meanE2_array(dirTags);
            obj.edgeData.meanEbsdOri2 = obj.edgeData.meanEbsdOri2(:).project2FundamentalRegion(obj.edgeData.grainOri2);

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
            angMatrix = angle_outer(mori_pairs(:), allMori(:)'); % Returns P x N matrix

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
                            is_newly_visited = true;
                        else
                            is_newly_visited = false;
                        end

                        % Identify the specific undirected edge
                        pair_id = [min(currNode, nxtNode), max(currNode, nxtNode)];
                        [~, loc] = ismember(pair_id, pairs, 'rows');

                        % Ensure every edge is traced exactly once in the directed graph
                        if edgeProcessed(loc)
                            continue;
                        end
                        edgeProcessed(loc) = true;

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

                            % Project grain, EBSD, and meanEBSD orientations
                            if any(idx1)
                                obj.edgeData.grainOri1(idx1) = obj.edgeData.grainOri1(idx1).project2FundamentalRegion(center_ori);
                                obj.edgeData.ebsdOri1(idx1)  = obj.edgeData.ebsdOri1(idx1).project2FundamentalRegion(center_ori);
                                obj.edgeData.meanEbsdOri1(idx1) = obj.edgeData.meanEbsdOri1(idx1).project2FundamentalRegion(center_ori);
                            end
                            if any(idx2)
                                obj.edgeData.grainOri2(idx2) = obj.edgeData.grainOri2(idx2).project2FundamentalRegion(center_ori);
                                obj.edgeData.ebsdOri2(idx2)  = obj.edgeData.ebsdOri2(idx2).project2FundamentalRegion(center_ori);
                                obj.edgeData.meanEbsdOri2(idx2) = obj.edgeData.meanEbsdOri2(idx2).project2FundamentalRegion(center_ori);
                            end

                            % Recompute misorientations for affected segments
                            idx_update = idx1 | idx2;
                            obj.edgeData.misorientation(idx_update) = inv(obj.edgeData.ebsdOri2(idx_update)) .* obj.edgeData.ebsdOri1(idx_update);
                        end

                        prevNode = currNode;
                        for step = 1:length(seq)
                            v_step = seq(step);
                            if step == length(seq)
                                s_nodes(end+1, 1) = prevNode;
                                t_nodes(end+1, 1) = nxtNode;
                                variants_list{end+1, 1} = num2str(v_step);
                                angular_deviation_list(end+1, 1) = dev;
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
                                angular_deviation_list(end+1, 1) = dev;

                                prevNode = vNodeIdx;
                            end
                        end
                    end
                end
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
                obj.virtualNodes = table(vNodeCoords(:,1), vNodeCoords(:,2), vNodeGens, 'VariableNames', {'x', 'y', 'Sequence'});
            else
                obj.virtualNodes = table([], [], {}, 'VariableNames', {'x', 'y', 'Sequence'});
            end
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
                if i <= max(obj.grains.id) && id2idx(i) > 0
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
                            
                            % Find corrected grain orientation
                            edgeIdx1 = find(obj.edgeData.grainId(:, 1) == i, 1);
                            edgeIdx2 = find(obj.edgeData.grainId(:, 2) == i, 1);
                            
                            if ~isempty(edgeIdx1)
                                correctedTwinOri = obj.edgeData.grainOri1(edgeIdx1);
                            elseif ~isempty(edgeIdx2)
                                correctedTwinOri = obj.edgeData.grainOri2(edgeIdx2);
                            else
                                correctedTwinOri = obj.grains(idx).meanOrientation; % fallback
                            end
                            
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
        function [Q, groupIDs] = computeQuotientGraph(obj)
            G = obj.Graph;
            if isempty(G) || numnodes(G) == 0
                Q = graph();
                groupIDs = [];
                obj.QuotientGraph = Q;
                obj.nodeGroups = groupIDs;
                return;
            end

            numNodes = numnodes(G);
            groupIDs = zeros(numNodes, 1);

            % Store the "reduced word" state for every node
            nodeStates = cell(numNodes, 1);

            % Isolate connected components
            bins = conncomp(G);
            numComps = max(bins);

            mstS = [];
            mstT = [];
            mstVar = {};
            mstDev = [];
            hasDeviation = ismember('AngularDeviation', G.Edges.Properties.VariableNames);

            for c = 1:numComps
                compNodeIdx = find(bins == c);
                compNodeIdx = compNodeIdx(:);

                % Base case for isolated nodes
                if length(compNodeIdx) == 1
                    nodeStates{compNodeIdx(1)} = sprintf('C%d_ROOT', c);
                    continue;
                end

                subG = subgraph(G, compNodeIdx);

                % 1. Create Minimum Spanning Tree to avoid spurious cycles
                if ismember('AngularDeviation', subG.Edges.Properties.VariableNames)
                    subG.Edges.Weight = subG.Edges.AngularDeviation;
                else
                    warning('TwinNetwork:computeQuotientGraph:noWeight', ...
                        'No ''AngularDeviation'' variable found in graph edges. Using default connectivity.');
                    subG.Edges.Weight = ones(numedges(subG), 1);
                end

            % VIRTUAL NODE PENALTY
            % Virtual nodes are defined as having an ID greater than the max physical grain ID
            maxPhysId = max(obj.grains.id);
            
            % compNodeIdx maps the local subG nodes to their global G IDs
            virtualLocalIdx = find(compNodeIdx > maxPhysId); 
            
            if ~isempty(virtualLocalIdx)
                % Extract all edge endpoints in the current subgraph
                endNodes = subG.Edges.EndNodes;
                
                % Find any edge where either the source OR the target is a virtual node
                virtualEdges = find(ismember(endNodes(:, 1), virtualLocalIdx) | ...
                                    ismember(endNodes(:, 2), virtualLocalIdx));
                
                % Add a massive penalty to these edges
                virtualPenalty = 1000; 
                subG.Edges.Weight(virtualEdges) = subG.Edges.Weight(virtualEdges) + virtualPenalty;
            end
            
                % Compute the MST (it will use subG.Edges.Weight if we just created it)
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

                % 2. Breadth-First Search to assign algebraic states
                rootNode = 1;

                subStates = cell(numnodes(subMST), 1);
                subStates{rootNode} = []; % The identity state

                visited = false(numnodes(subMST), 1);
                visited(rootNode) = true;

                queue = rootNode;

                while ~isempty(queue)
                    currNode = queue(1);
                    queue(1) = [];

                    currState = subStates{currNode};

                    nbrs = neighbors(subMST, currNode);
                    for n = 1:length(nbrs)
                        nxtNode = nbrs(n);
                        if ~visited(nxtNode)
                            visited(nxtNode) = true;

                            % Find the edge variant in the MST
                            edgeIdx = findedge(subMST, currNode, nxtNode);
                            v_str = subMST.Edges.Variant{edgeIdx};
                            v_num = str2double(v_str);

                            % 3. The Algebraic Contraction Rule
                            nxtState = currState;
                            if ~isempty(nxtState) && nxtState(end) == v_num
                                % Identical adjacent variants contract (cancel out)
                                nxtState(end) = [];
                            else
                                % Otherwise, append to the sequence
                                nxtState(end+1) = v_num;
                            end

                            subStates{nxtNode} = nxtState;
                            queue(end+1) = nxtNode;
                        end
                    end
                end

                % 4. Convert numeric arrays to strings for unique grouping
                for i = 1:length(compNodeIdx)
                    globalID = compNodeIdx(i);
                    s = subStates{i};

                    if isempty(s)
                        nodeStates{globalID} = sprintf('C%d_ROOT', c);
                    else
                        str_s = sprintf('%d_', s);
                        % Example state: "C1_2_1"
                        nodeStates{globalID} = sprintf('C%d_%s', c, str_s(1:end-1));
                    end
                end
            end

            % 5. Group nodes that resolved to the exact same reduced word
            [unqStates, ~, groupIDs] = unique(string(nodeStates));

            % 6. Build the final Quotient Graph
            numGroups = length(unqStates);
            Q = graph();
            Q = addnode(Q, numGroups);

            if numedges(G) > 0
                u = groupIDs(G.Edges.EndNodes(:, 1));
                v = groupIDs(G.Edges.EndNodes(:, 2));
                vars = G.Edges.Variant;

                qEdges = [u, v];
                qEdges = sort(qEdges, 2); % Standardize undirected edges

                strKeys = cell(size(qEdges, 1), 1);
                for i = 1:size(qEdges, 1)
                    strKeys{i} = sprintf('%d_%d_%s', qEdges(i,1), qEdges(i,2), vars{i});
                end

                [~, unqIdx] = unique(strKeys);

                % Filter out self-loops (edges where both nodes contract to the same group)
                validEdges = qEdges(unqIdx, 1) ~= qEdges(unqIdx, 2);

                EdgeT = table(qEdges(unqIdx(validEdges), :), vars(unqIdx(validEdges)), ...
                    'VariableNames', {'EndNodes', 'Variant'});

                if ~isempty(EdgeT)
                    Q = addedge(Q, EdgeT);
                end
            end

            obj.QuotientGraph = Q;
            obj.nodeGroups = groupIDs;

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
            obj.mstGraph = mstG;
        warning('TwinNetwork:computeQuotientGraph:sanityCheck','Here we can introduce a sanity check, via grain orientation, we check for every group in the quotient graph if there are significant outliers (>15°) and then see if any other edges attached to these outliers allows for a different grouping, for which the orientation gets checked again. If non of the edges allow for a more fitting grouping, the original grouping in the quotient graph is kept and a warning is displayed. If a more fitting grouping is found (all groupings have to be possible in the tagged graph, the orinentation is only the backcheck not the source of the grouping), then this node gets moved to another group in the quotient and the necessary edge-edits will be made in the MST. Maybe this behavior can be toggled with a "stable" true/false flag.')
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
