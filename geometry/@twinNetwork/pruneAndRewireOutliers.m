function obj = pruneAndRewireOutliers(obj, varargin)
    % Step D: Detect and rewire outliers and isolated nodes based on empirical FZ centers
    p = inputParser;
    addParameter(p, 'targetNode', 0);
    addParameter(p, 'debug', false);
    addParameter(p, 'outlierMode', 'empirical', @(x) any(validatestring(x, {'empirical', 'theoretical'})));
    addParameter(p, 'maxErrThreshold', 15 * degree);
    addParameter(p, 'rejectGenIncrease', false, @islogical);
    parse(p, varargin{:});
    target_node_override = p.Results.targetNode;
    debug_mode = p.Results.debug;
    outlierMode = p.Results.outlierMode;
    max_err_threshold = p.Results.maxErrThreshold;
    reject_gen_increase = p.Results.rejectGenIncrease;
    
    if isempty(obj.mstGraph) || ~isa(obj.mstGraph, 'graph')
        return;
    end
    mstG = obj.mstGraph;
    numNodes = numnodes(mstG);
    if numNodes == 0 || numedges(mstG) == 0
        return;
    end
    
    % The physical topology is inherently undirected. We must search an undirected graph!
    search_G_dir = obj.Graph;
    search_G = graph();
    search_G = addnode(search_G, search_G_dir.Nodes);
    search_G = addedge(search_G, search_G_dir.Edges);
    
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
        % Convert VariantPath to string just for grouping purposes
        pathStrs = cellfun(@mat2str, mstG.Nodes.VariantPath, 'UniformOutput', false);
        [unqStates, ~, groupIDs] = unique(table(mstG.Nodes.ComponentID, pathStrs, 'VariableNames', {'ComponentID', 'pathStrs'}));
        numStates = height(unqStates);
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
                
                if isempty(eval(char(unqStates.pathStrs(s))))
                    s_nodes = find(groupIDs == s);
                    phys_nodes = s_nodes(mstG.Nodes.Type(s_nodes) == "Physical");
                    if ~isempty(phys_nodes)
                        oris = mstG.Nodes.MeanOrientation(phys_nodes);
                        if length(oris) <= 1
                            r_ori = oris;
                        else
                            r_ori = obj.computeRobustCenter(oris);
                        end
                        
                        comp_idx = unqStates.ComponentID(s);
                        root_centers(num2str(comp_idx)) = r_ori;
                        
                        devs = angle(oris, r_ori);
                        s_safe = phys_nodes(devs < max_err_threshold);
                        if isempty(s_safe)
                            s_safe = phys_nodes;
                        end
                        safe_roots(num2str(comp_idx)) = s_safe;
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
                
                
                comp_idx_str = num2str(unqStates.ComponentID(s));
                
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
                    var_path = eval(char(unqStates.pathStrs(s)));
                    for v_idx = 1:length(var_path)
                        theoOri = theoOri * obj.variants(var_path(v_idx));
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
                fprintf('\nAttempting to rewire OUTLIER Grain %d at depth %d with deviation %.2f deg to its theoretical orientation as calculated from the closest safe root Grain %d...\n', ...
                    v_node, all_problems(p_idx, 4), all_problems(p_idx, 2)/degree, all_problems(p_idx, 5));
            else
                fprintf('\nAttempting to rewire OUTLIER Grain %d at depth %d...\n', v_node, all_problems(p_idx, 4));
            end
        else
            fprintf('\nAttempting to rewire ISOLATED Grain %d at depth %d...\n', v_node, all_problems(p_idx, 4));
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
            temp_mstG.Nodes.ComponentID(:) = 0;
            temp_mstG.Nodes.VariantPath(:) = {[]};
            temp_mstG = obj.resolveNodeStates(temp_mstG);
            
            % Validate the new orientation of v_node
            new_state_comp = temp_mstG.Nodes.ComponentID(v_node);
            new_state_path = temp_mstG.Nodes.VariantPath{v_node};
            
            % Compute robust center for the new state to check deviation
            % [temp_unqStates, ~, temp_groupIDs] = unique(temp_mstG.Nodes.State);
            % s_idx = find(temp_unqStates == new_state);
            phys_nodes = find(temp_mstG.Nodes.ComponentID == new_state_comp & cellfun(@(x) isequal(x, new_state_path), temp_mstG.Nodes.VariantPath) & temp_mstG.Nodes.Type == "Physical");
            
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
                comp_idx_str = num2str(new_state_comp);
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
                    var_path = eval(char(unqStates.pathStrs(s)));
                    for v_idx = 1:length(var_path)
                        theoOri = theoOri * obj.variants(var_path(v_idx));
                    end
                    v_ori = temp_mstG.Nodes.MeanOrientation(v_node);
                    dev_val = angle(v_ori, theoOri);
                end
            end
            
            % Determine generation for the old and new states
            old_gen = length(mstG.Nodes.VariantPath{v_node});
            new_gen = length(new_state_path);
            
            gen_increased = new_gen > old_gen;
            
            if dev_val < max_err_threshold && ~(reject_gen_increase && gen_increased)
                fallback_success = true;
                
                if debug_mode
                    fprintf('    [DEBUG] Iter %d SUCCESS for Grain %d:\n', fallback_iter, v_node);
                    fprintf('      - Severed Edge:  Grain %d <-> Grain %d\n', v_node, p_node);
                    fprintf('      - New State Gen: %d\n', new_gen);
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
                        fprintf('      - Original State Gen: %d\n', old_gen);
                        fprintf('      - Resulting State Gen: %d\n', new_gen);
                        fprintf('      - Deviation: %.2f deg (Passed %.2f deg threshold)\n', dev_val/degree, max_err_threshold/degree);
                    else
                        fprintf('    [DEBUG] Iter %d FAILED validation:\n', fallback_iter);
                        fprintf('      - Resulting State Gen: %d\n', new_gen);
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
            
            if isempty(eval(char(unqStates.pathStrs(s)))) && ~isempty(groupRefOris{s})
                comp_idx_str = num2str(unqStates.ComponentID(s));
                root_centers(comp_idx_str) = groupRefOris{s};
            end
        end
        
        % Validate every non-root group
        for s = 1:numStates
            
            if isempty(eval(char(unqStates.pathStrs(s)))) || isempty(groupRefOris{s})
                continue;
            end
            
            comp_idx_str = num2str(unqStates.ComponentID(s));
            var_path = eval(char(unqStates.pathStrs(s)));
            
            if ~root_centers.isKey(comp_idx_str)
                continue; % Shouldn't happen, but safe fallback
            end
            
            % Reconstruct theoretical orientation purely from the root and the twin chain
            theoOri = root_centers(comp_idx_str);
            for v_idx = 1:length(var_path)
                variant_idx = var_path(v_idx);
                theoOri = theoOri * obj.variants(variant_idx);
            end
            
            empiricalOri = groupRefOris{s};
            dev_val = angle(empiricalOri, theoOri);
            
            if dev_val > max_err_threshold
                fprintf('  [WARNING] Group %d deviates from theoretical chain by %.2f deg!\n', s, dev_val / degree);
            else
                fprintf('  [OK] Group %d aligns theoretically (Dev: %.2f deg)\n', s, dev_val / degree);
            end
        end
        fprintf('-------------------------------------------\n\n');
    end
    
    obj.mstGraph = mstG;
    obj.pipelineState = 'OutliersPruned';
end
