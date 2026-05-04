function [final_mask, p] = interactive_specimen_mask(ci_matrix, init_p)
% INTERACTIVE_SPECIMEN_MASK Interactive GUI for tracing and snapping specimen boundaries
%
% Inputs:
%   ci_matrix - 2D matrix of EBSD Confidence Index (or similar property)
%   init_p    - Struct containing initialization parameters
%
% Outputs:
%   final_mask - 2D boolean mask indicating the specimen area
%   p          - Struct containing the updated parameters

    fig = figure('Name', 'Trace Specimen Boundaries', 'Position', [100 100 1100 650], ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'UserData', 'waiting');
    
    ax = axes('Parent', fig, 'Position', [0.35 0.1 0.6 0.8]);
    pnl = uipanel('Parent', fig, 'Position', [0.02 0.1 0.3 0.8], 'Title', 'Edge Tracing Parameters');
    
    % Make init_p optional by providing robust defaults based on the image size
    if nargin < 2 || isempty(init_p)
        p.sigma = 2.0;       % Smoothing for gradient calculation
        p.disp_min = 0.0;    % Display contrast min
        p.disp_max = 0.1;    % Display contrast max
        p.snap_window = 20;  % Snapping window radius in pixels
        
        p.left_edge_x = [size(ci_matrix,2)*0.2; size(ci_matrix,2)*0.2];
        p.left_edge_y = [1; size(ci_matrix,1)];
        
        p.right_edge_x = [size(ci_matrix,2)*0.8; size(ci_matrix,2)*0.8];
        p.right_edge_y = [1; size(ci_matrix,1)];
    else
        p = init_p;
    end
    
    final_mask = false(size(ci_matrix));
    left_roi = [];
    right_roi = [];
    img_handle = [];
    bound_handles = [];
    
    % Initial image setup
    axes(ax);
    img_handle = imagesc(ax, zeros(size(ci_matrix)), [p.disp_min, p.disp_max]);
    colormap(ax, hot); hold(ax, 'on');
    axis(ax, 'image'); axis(ax, 'off');
    set(ax, 'Color', 'k');
    set(fig, 'Color', [0.15 0.15 0.15]);
    set(pnl, 'BackgroundColor', [0.9 0.9 0.9]);
    
    function [sld, lbl] = make_slider(name, bounds, ynorm, val)
        lbl = uicontrol('Parent', pnl, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.05, ynorm+0.04, 0.9, 0.04], 'String', sprintf('%s: %.3f', name, val), ...
            'HorizontalAlignment', 'left', 'FontSize', 10);
        sld = uicontrol('Parent', pnl, 'Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.05, ynorm, 0.9, 0.04], 'Min', bounds(1), 'Max', bounds(2), 'Value', val);
        addlistener(sld, 'ContinuousValueChange', @(s,e) set(lbl, 'String', sprintf('%s: %.3f', name, s.Value)));
    end
    
    [sld_sigma, lbl_sigma] = make_slider('Smoothing Sigma', [0.1, 10], 0.82, p.sigma);
    [sld_min, lbl_min] = make_slider('Edge Display Min', [0, 0.2], 0.70, p.disp_min);
    [sld_max, lbl_max] = make_slider('Edge Display Max', [0.01, 0.5], 0.58, p.disp_max);
    [sld_snap, lbl_snap] = make_slider('Snap Window (px)', [0, 100], 0.46, p.snap_window);
    
    uicontrol('Parent', pnl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.05, 0.32, 0.9, 0.08], 'String', 'Trace LEFT Edge', 'Callback', @draw_left_cb, 'FontSize', 10, 'ForegroundColor', [0.8 0.8 0]);
    uicontrol('Parent', pnl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.05, 0.22, 0.9, 0.08], 'String', 'Trace RIGHT Edge', 'Callback', @draw_right_cb, 'FontSize', 10, 'ForegroundColor', [0 0.8 0]);
    
    uicontrol('Parent', pnl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.05, 0.03, 0.4, 0.08], 'String', 'Apply / Update View', 'Callback', @apply_cb, 'FontSize', 10, 'FontWeight', 'bold');
    uicontrol('Parent', pnl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.5, 0.03, 0.45, 0.08], 'String', 'Save & Close', 'Callback', @save_cb, 'FontSize', 10, 'FontWeight', 'bold');
        
    function draw_left_cb(~, ~)
        title(ax, 'DRAWING LEFT EDGE: Click waypoints along the edge. Double-click to finish.', 'Color', 'y', 'FontSize', 11, 'FontWeight', 'bold');
        if ~isempty(left_roi) && isvalid(left_roi), delete(left_roi); end
        try
            roi = drawpolyline(ax, 'Color', 'y', 'LineWidth', 2);
            if isvalid(roi) 
                p.left_edge_x = roi.Position(:,1);
                p.left_edge_y = roi.Position(:,2);
                left_roi = roi;
                addlistener(left_roi, 'ROIMoved', @apply_cb);
            end
        catch ME
            error('Drawing failed: %s', ME.message);
        end
        apply_cb();
    end

    function draw_right_cb(~, ~)
        title(ax, 'DRAWING RIGHT EDGE: Click waypoints along the edge. Double-click to finish.', 'Color', 'g', 'FontSize', 11, 'FontWeight', 'bold');
        if ~isempty(right_roi) && isvalid(right_roi), delete(right_roi); end
        try
            roi = drawpolyline(ax, 'Color', 'g', 'LineWidth', 2);
            if isvalid(roi) 
                p.right_edge_x = roi.Position(:,1);
                p.right_edge_y = roi.Position(:,2);
                right_roi = roi;
                addlistener(right_roi, 'ROIMoved', @apply_cb);
            end
        catch ME
            error('Drawing failed: %s', ME.message);
        end
        apply_cb();
    end

    function apply_cb(~, ~)
        p.sigma = get(sld_sigma, 'Value');
        p.disp_min = get(sld_min, 'Value');
        p.disp_max = get(sld_max, 'Value');
        p.snap_window = round(get(sld_snap, 'Value'));
        
        if ~isempty(left_roi) && isvalid(left_roi)
            p.left_edge_x = left_roi.Position(:,1);
            p.left_edge_y = left_roi.Position(:,2);
        end
        if ~isempty(right_roi) && isvalid(right_roi)
            p.right_edge_x = right_roi.Position(:,1);
            p.right_edge_y = right_roi.Position(:,2);
        end
        
        % 1. Generate Mask from Traces and get calculated Edge Map
        [final_mask, Gmag] = compute_mask(ci_matrix, p);
        
        % 2. Update Display
        if isvalid(img_handle)
            set(img_handle, 'CData', Gmag);
            set(ax, 'CLim', [p.disp_min, p.disp_max]);
        end
        
        if ~isempty(bound_handles)
            delete(bound_handles(isvalid(bound_handles)));
        end
        if any(final_mask(:))
            bound_handles = visboundaries(ax, final_mask, 'Color', 'r', 'LineWidth', 1.5);
        end
        
        title(ax, 'Edge Map Preview (Red = Snapped Boundary)', 'Color', 'w', 'FontSize', 12);
    end
    
    function save_cb(~, ~)
        set(fig, 'UserData', 'done');
    end
    
    % Instantiate ROIs from initialization parameters before first render
    if isfield(p, 'left_edge_x') && ~isempty(p.left_edge_x)
        left_roi = drawpolyline(ax, 'Position', [p.left_edge_x, p.left_edge_y], 'Color', 'y', 'LineWidth', 2);
        addlistener(left_roi, 'ROIMoved', @apply_cb);
    end
    if isfield(p, 'right_edge_x') && ~isempty(p.right_edge_x)
        right_roi = drawpolyline(ax, 'Position', [p.right_edge_x, p.right_edge_y], 'Color', 'g', 'LineWidth', 2);
        addlistener(right_roi, 'ROIMoved', @apply_cb);
    end
    
    apply_cb();
    
    % Block script execution using waitfor instead of uiwait
    waitfor(fig, 'UserData', 'done');
    
    if ishandle(fig)
        close(fig);
    end
end

function [mask_2d, Gmag] = compute_mask(ci_matrix, p)
    mask_2d = false(size(ci_matrix));
    
    % Always calculate the edge map for both snapping and display
    ci_smooth = imgaussfilt(ci_matrix, p.sigma);
    [Gmag, ~] = imgradient(ci_smooth, 'Sobel');
    
    if ~isfield(p, 'left_edge_x') || isempty(p.left_edge_x) || ...
       ~isfield(p, 'right_edge_x') || isempty(p.right_edge_x)
        return;
    end
    
    % Image dimensions
    num_rows = size(ci_matrix, 1);
    num_cols = size(ci_matrix, 2);
    y_query = (1:num_rows)';
    
    % Prepare Left Edge data (must have unique Y values for interpolation)
    [uy_left, idx_left] = unique(p.left_edge_y);
    ux_left = p.left_edge_x(idx_left);
    min_y_left = min(uy_left);
    max_y_left = max(uy_left);
    
    % Prepare Right Edge data
    [uy_right, idx_right] = unique(p.right_edge_y);
    ux_right = p.right_edge_x(idx_right);
    min_y_right = min(uy_right);
    max_y_right = max(uy_right);
    
    % Interpolate the X coordinates strictly within the drawn bounds
    if length(uy_left) >= 2
        x_left_full = interp1(uy_left, ux_left, y_query, 'linear');
    else
        x_left_full = NaN(num_rows, 1);
        if length(uy_left) == 1
            x_left_full(uy_left) = ux_left;
        end
    end
    
    if length(uy_right) >= 2
        x_right_full = interp1(uy_right, ux_right, y_query, 'linear');
    else
        x_right_full = NaN(num_rows, 1);
        if length(uy_right) == 1
            x_right_full(uy_right) = ux_right;
        end
    end
    
    % FINISH LINES TO BOUNDARIES
    % If the drawing stops before the top/bottom, bridge it to the left/right canvas borders
    x_left_full(y_query < min_y_left) = 1;
    x_left_full(y_query > max_y_left) = 1;
    
    x_right_full(y_query < min_y_right) = num_cols;
    x_right_full(y_query > max_y_right) = num_cols;
    
    % Constrain interpolated coordinates to image bounds before snapping
    x_left_full = round(max(1, min(num_cols, x_left_full)));
    x_right_full = round(max(1, min(num_cols, x_right_full)));
    
    % MAGNETIC SNAPPING LOGIC
    snap_w = round(p.snap_window);
    if snap_w > 0
        for r = 1:num_rows
            % Snap Left Edge ONLY within the user-drawn vertical bounds
            if r >= min_y_left && r <= max_y_left
                c_left = x_left_full(r);
                c_min = max(1, c_left - snap_w);
                c_max = min(num_cols, c_left + snap_w);
                [~, max_idx] = max(Gmag(r, c_min:c_max));
                x_left_full(r) = c_min + max_idx - 1;
            end
            
            % Snap Right Edge ONLY within the user-drawn vertical bounds
            if r >= min_y_right && r <= max_y_right
                c_right = x_right_full(r);
                c_min = max(1, c_right - snap_w);
                c_max = min(num_cols, c_right + snap_w);
                [~, max_idx] = max(Gmag(r, c_min:c_max));
                x_right_full(r) = c_min + max_idx - 1;
            end
        end
    end
    
    % Fill the mask horizontally between the two traces
    for r = 1:num_rows
        left_col = min(x_left_full(r), x_right_full(r));
        right_col = max(x_left_full(r), x_right_full(r));
        mask_2d(r, left_col:right_col) = true;
    end
end