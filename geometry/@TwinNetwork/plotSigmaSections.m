function plotSigmaSections(obj, varargin)
    % Plots orientations of groups in sigma sections to verify clustering
    % % Plot ALL groups with their robust centers (as stars)
    % network.plotSigmaSections();
    % % Plot only a specific group without the legend
    % network.plotSigmaSections('group', 'C1_1_2', 'showLegend', false);
    % % Plot all groups but turn off the robust center stars
    % network.plotSigmaSections('robustCenter', false);
    
    if isempty(obj.mstGraph) || ~isa(obj.mstGraph, 'graph') || numnodes(obj.mstGraph) == 0
        warning('twinNetwork:plotSigmaSections:noMstGraph', 'Spanning tree not found. Please run the pipeline steps first.');
        return;
    end
    mstG = obj.mstGraph;
    
    p = inputParser;
    addParameter(p, 'group', 'all');
    addParameter(p, 'showLegend', true);
    addParameter(p, 'robustCenter', true);
    parse(p, varargin{:});
    
    targetGroup = p.Results.group;
    showLegend = p.Results.showLegend;
    plotRobust = p.Results.robustCenter;
    
    if isempty(obj.QuotientGraph) || ~isa(obj.QuotientGraph, 'graph') || numnodes(obj.QuotientGraph) == 0
        warning('twinNetwork:plotSigmaSections:noQuotient', 'Quotient graph not found. Please call reduceToQuotient() first.');
        return;
    end
    
    numStates = numnodes(obj.QuotientGraph);
    states = strings(numStates, 1);
    for i = 1:numStates
        states(i) = string(obj.formatState(obj.QuotientGraph.Nodes.ComponentID(i), obj.QuotientGraph.Nodes.VariantPath{i}));
    end
    
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
            ref_ori = obj.QuotientGraph.Nodes.robustMeanOrientation(s);
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
