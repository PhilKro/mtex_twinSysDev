function plotGraph(obj, G, varargin)
    % PLOTGRAPH Visualizes any of the topological graphs over the EBSD grain map.
    %
    % Syntax:
    %   network.plotGraph(network.Graph)
    %   network.plotGraph(network.mstGraph)
    %   network.plotGraph(network.QuotientGraph)
    %
    % Description:
    %   This is a robust plotting utility that automatically adapts its visualization
    %   based on the properties present in the supplied graph `G`:
    %
    %   Node Labeling:
    %     - If G has `VariantPath`, nodes are labeled with their FZ generation path (e.g., C1_1_2)
    %     - If G has `GrainId`, nodes are labeled with their physical MTEX grain IDs.
    %
    %   Edge Coloring & Labeling:
    %     - If G has `AngularDeviation`, edges are colored using a 'jet' colormap representing the 
    %       misorientation deviation from the theoretical twin boundary (in degrees).
    %     - If G has `Variant`, edges are labeled with the twin variant index (1..numVariants).
    %     - If G has `Sequence` (and Variant is NaN), edges are labeled with the theoretical generation sequence.
    %
    %   Options:
    %     'showBackground' - Logical flag to toggle the grain map overlay. 
    %                        Set to false to plot the graph cleanly on a white background. (Default: true)
    %     'layout'         - Graph layout method. 'spatial' uses real grain centroids. 
    %                        'layered' is excellent for hierarchical generation trees. (Default: 'spatial')
    %     'scaleNodesByArea' - Logical flag to scale node size relative to physical area. (Default: false)
    if nargin < 2 || isempty(G) || (~isa(G, 'graph') && ~isa(G, 'digraph')) || numnodes(G) == 0
        disp('Graph is empty or not provided.');
        return;
    end
    
    p = inputParser;
    addParameter(p, 'showBackground', true, @islogical);
    addParameter(p, 'noGrainsBG', false, @islogical);
    addParameter(p, 'layout', 'spatial', @ischar);
    addParameter(p, 'scaleNodesByArea', false, @islogical);
    parse(p, varargin{:});
    show_bg = p.Results.showBackground && ~p.Results.noGrainsBG;
    graph_layout = p.Results.layout;
    scale_nodes = p.Results.scaleNodesByArea;
    
    figure;
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    nexttile;
    
    if show_bg
        plot(obj.grains, 'micronbar', 'off', 'FaceColor', [0.85 0.85 0.85]);
        hold on;
        plot(obj.grains.boundary);
    else
        hold on;
        axis equal;
        axis off;
    end
    
    if strcmpi(graph_layout, 'spatial') && ismember('Centroid', G.Nodes.Properties.VariableNames)
        X = G.Nodes.Centroid(:, 1);
        Y = G.Nodes.Centroid(:, 2);
        pPlot = plot(G, 'XData', X, 'YData', Y);
    else
        pPlot = plot(G, 'Layout', graph_layout);
    end
    
    if scale_nodes && ismember('Area', G.Nodes.Properties.VariableNames)
        areas = G.Nodes.Area;
        minA = min(areas); maxA = max(areas);
        if maxA > minA
            pPlot.MarkerSize = 5 + 20 * (areas - minA) / (maxA - minA);
        else
            pPlot.MarkerSize = 10;
        end
    else
        pPlot.MarkerSize = 5;
    end
    
    pPlot.NodeColor = 'k';
    pPlot.LineWidth = 1.5;
    
    % Labeling logic
    if ismember('VariantPath', G.Nodes.Properties.VariableNames)
        nodeLabels = cell(numnodes(G), 1);
        for i = 1:numnodes(G)
            nodeLabels{i} = obj.formatState(G.Nodes.ComponentID(i), G.Nodes.VariantPath{i});
        end
        pPlot.NodeLabel = nodeLabels;
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
    
    hold off;
end
