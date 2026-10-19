function plot(obj, varargin)
    p = inputParser;
    addParameter(p, 'twinsOnly', true);
    addParameter(p, 'noText', false, @islogical);
    addParameter(p, 'showText', true, @islogical);
    addParameter(p, 'lineWidth', 1);
    addParameter(p, 'showBoundary', true, @islogical);
    parse(p, varargin{:});
    twinsOnly = p.Results.twinsOnly;
    noText = p.Results.noText || ~p.Results.showText;
    lineWidth = p.Results.lineWidth;
    showBoundary = p.Results.showBoundary;
    
    if isempty(obj.QuotientGraph) || ~isa(obj.QuotientGraph, 'graph')
        warning('twinNetwork:plot:noQuotient', 'Quotient graph not found. Please call reduceToQuotient() first.');
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
    
    if numnodes(obj.QuotientGraph) > 0 && ~isempty(obj.nodeGroups)
        Q_deg = degree(obj.QuotientGraph);
        
        for i = 1:maxPhysId
            idx = id2idx(i);
            if idx > 0 && i <= length(obj.nodeGroups)
                grp = obj.nodeGroups(i);
                if grp <= 0 || grp > numnodes(obj.QuotientGraph)
                    continue;
                end
                isIsolated = Q_deg(grp) == 0;
                
                if twinsOnly && isIsolated
                    toPlot(idx) = false;
                    continue;
                end
                
                if isIsolated
                    continue; % Remains white, no text overlay
                end
                
                var_path = obj.QuotientGraph.Nodes.VariantPath{grp};
                
                if isempty(var_path)
                    % Root Generation
                    gen = 0;
                    path = [];
                else
                    % Has variants
                    gen = length(var_path);
                    path = var_path;
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
    else
        if twinsOnly
            toPlot(:) = false;
        end
    end
    
    % Plot the grains
    if twinsOnly
        if any(toPlot)
            plot(obj.grains(toPlot), faceColors(toPlot, :), 'micronbar', 'off');
        else
            plot(obj.grains, 'FaceColor', [0.85 0.85 0.85], 'micronbar', 'off');
        end
    else
        plot(obj.grains, faceColors, 'micronbar', 'off');
    end
    hold on;
    
    if showBoundary
        plot(obj.grains.boundary, 'linewidth', lineWidth, 'linecolor', 'k');
    end
    
    % Plot text
    if ~noText && ~isempty(labelTexts)
        text(labelCentersX, labelCentersY, labelTexts, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 10, ...
            'Color', 'k');
    end
    hold off;
end
