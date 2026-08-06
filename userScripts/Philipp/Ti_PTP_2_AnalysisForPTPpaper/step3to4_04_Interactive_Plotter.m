function step3to4_04_Interactive_Plotter()
% =========================================================================
% step3to4_04_Interactive_Plotter.m
% =========================================================================
% INTERACTIVE REAL-TIME EBSD & DIC PLOTTING STUDIO (ZERO-LAG ENGINE)
%
% This tool provides a real-time interactive studio to blend greyscale EBSD
% and DIC diagnostic background maps with Inverse Pole Figure (IPF) crystal 
% orientation overlays.
%
% Features:
%   1. Real-Time Zero-Lag Rendering: Utilizes pre-gridded spatial buffers 
%      and dynamic RGB CData injection to ensure silky-smooth ~60 FPS update
%      rates during live slider dragging (via PostSet listeners).
%   2. IPF Overlay Controls: Toggle between MTEX coloring schemes 
%      (ipfColorKey, ipfHKLKey, ipfHSVKey, ipfSpotKey, ipfTSLKey) and vary
%      linear alpha blending from 0.0 (pure background) to 1.0 (pure IPF).
%   3. Advanced Greyscale Background Mastery: Select from any scalar field
%      in ebsd.prop (band contrast, confidence index, strain/displacement 
%      gradients) with precise Black Point, White Point, and Gamma sliders.
%   4. Publication-Ready Export Suite: On saving, generates a clean, non-GUI
%      high-resolution export via exportScaledFigure, along with an audit
%      MAT parameter archive and an interactive Markdown summary report.
%
% Author: Philipp / Antigravity
% Date:   2026-08-05
% =========================================================================

    clear variables; close all; clc;
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');

    % Ensure standard MTEX spatial convention
    plotx2east;
    plotzIntoPlane;

    %% 1. Locate & Load Default Dataset
    fprintf('Initializing Interactive EBSD Plotting Studio...\n');
    try
        [~, defaults] = a00_Dataset_Registry();
        inDir  = fullfile(defaults.dir, 'Processed_Output_Step3to4');
        matFile = fullfile(inDir, 'step3to4_EBSD_withDIC.mat');
        if ~exist(matFile, 'file')
            matFile = fullfile(defaults.dir, 'Processed_Output_Undeformed', 'step3to4_EBSD_withDIC.mat');
        end
    catch
        inDir = pwd;
        matFile = '';
    end

    if ~exist(matFile, 'file')
        [fname, fpath] = uigetfile('*.mat', 'Select EBSD Data MAT File (*.mat)');
        if isequal(fname, 0), error('No EBSD dataset selected. Exiting.'); end
        matFile = fullfile(fpath, fname);
        inDir   = fpath;
    end

    fprintf('Loading EBSD Dataset:\n  %s\n', matFile);
    dataLoaded = load(matFile);
    
    if isfield(dataLoaded, 'newEBSD')
        ebsd = dataLoaded.newEBSD;
    elseif isfield(dataLoaded, 'ebsd')
        ebsd = dataLoaded.ebsd;
    else
        error('Dataset must contain an EBSD object named "newEBSD" or "ebsd".');
    end
    
    if isfield(dataLoaded, 'grains')
        grains = dataLoaded.grains;
    else
        fprintf('Calculating grains dynamically (>300 px)...\n');
        [grains, ~] = calcGrains(ebsd, 'minPixel', 300, 'angle', 5*degree);
    end

    % Setup Default Export Routing
    exportDir = fullfile(inDir, 'Interactive_Exports');
    if ~exist(exportDir, 'dir'), mkdir(exportDir); end

    %% 2. Precompute Spatial Grid & Extract Property Fields
    fprintf('Gridifying dataset for ultra-fast GPU matrix blending...\n');
    ebsd_grid = ebsd.gridify;
    
    % Identify valid numeric fields in ebsd.prop for background selection
    prop_names = fields(ebsd_grid.prop);
    valid_fields = {};
    for idx = 1:length(prop_names)
        fn = prop_names{idx};
        val = ebsd_grid.prop.(fn);
        if isnumeric(val) || isa(val, 'double') || isa(val, 'vector3d') || isstruct(val)
            valid_fields{end+1} = fn; %#ok<AGROW>
        end
    end
    
    % Ensure base diagnostic qualities are listed first
    base_props = {'ci', 'bc', 'imagequality', 'iq', 'band_contrast'};
    for idx = 1:length(base_props)
        if isfield(ebsd_grid.prop, base_props{idx}) && ~any(strcmp(valid_fields, base_props{idx}))
            valid_fields = [base_props{idx}, valid_fields];
        elseif isprop(ebsd_grid, base_props{idx}) && ~any(strcmp(valid_fields, base_props{idx}))
            valid_fields = [base_props{idx}, valid_fields];
        end
    end
    if isempty(valid_fields), valid_fields = {'ci'}; end

    % Supported IPF Color Keys in MTEX
    ipf_keys = {'ipfColorKey', 'ipfHKLKey', 'ipfHSVKey', 'ipfSpotKey', 'ipfTSLKey'};

    %% 3. Build Interactive Studio GUI Layout
    fig_w = 1200; fig_h = 820;
    scr = get(0, 'ScreenSize');
    pos_x = max(50, round((scr(3) - fig_w)/2));
    pos_y = max(50, round((scr(4) - fig_h)/2));
    
    fig = figure('Name', 'Interactive EBSD & DIC Plotting Studio', ...
                 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure', ...
                 'Position', [pos_x, pos_y, fig_w, fig_h], 'Color', [0.93 0.93 0.93], ...
                 'Resize', 'on');

    % Left Area: Map Canvas
    ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05, 0.06, 0.65, 0.90]);
    axis(ax, 'image'); hold(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'x [\mum]', 'FontWeight', 'bold'); ylabel(ax, 'y [\mum]', 'FontWeight', 'bold');
    title(ax, 'Interactive Composite Map', 'FontSize', 14, 'FontWeight', 'bold');

    % Right Area: Control Sidebar Panel
    pan = uipanel('Parent', fig, 'Title', 'STUDIO CONTROLS & EXPORT', ...
                  'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [0.95 0.95 0.95], ...
                  'Units', 'normalized', 'Position', [0.72, 0.02, 0.27, 0.96]);

    % --- SECTION 1: IPF OVERLAY CONTROLS ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 1. IPF COLOUR OVERLAY ===', ...
              'Units', 'normalized', 'Position', [0.05, 0.94, 0.90, 0.03], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.85 0.85 0.90]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'IPF Color Key Scheme:', ...
              'Units', 'normalized', 'Position', [0.05, 0.89, 0.90, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    pu_ipf = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', ipf_keys, ...
                       'Units', 'normalized', 'Position', [0.05, 0.85, 0.90, 0.04], 'FontSize', 10);

    lbl_alpha = uicontrol('Parent', pan, 'Style', 'text', 'String', 'IPF Overlay Alpha: 0.50', ...
                          'Units', 'normalized', 'Position', [0.05, 0.81, 0.90, 0.03], ...
                          'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_alpha  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.0, 'Max', 1.0, 'Value', 0.50, ...
                          'Units', 'normalized', 'Position', [0.05, 0.78, 0.90, 0.03]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Grain Boundary Style:', ...
              'Units', 'normalized', 'Position', [0.05, 0.73, 0.90, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    gb_styles = {'White (1.2 px)', 'Black (1.2 px)', 'White (2.0 px)', 'Black (2.0 px)', 'Hidden / None'};
    pu_gb     = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', gb_styles, ...
                          'Units', 'normalized', 'Position', [0.05, 0.69, 0.90, 0.04], 'FontSize', 10);

    % --- SECTION 2: GREYSCALE BACKGROUND CONTROLS ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 2. GREYSCALE BACKGROUND ===', ...
              'Units', 'normalized', 'Position', [0.05, 0.63, 0.90, 0.03], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.85 0.90 0.85]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Select Background Field:', ...
              'Units', 'normalized', 'Position', [0.05, 0.58, 0.90, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    pu_bg = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', valid_fields, ...
                      'Units', 'normalized', 'Position', [0.05, 0.54, 0.90, 0.04], 'FontSize', 10);

    lbl_bp = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Black Point: 0.00', ...
                       'Units', 'normalized', 'Position', [0.05, 0.49, 0.90, 0.03], ...
                       'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_bp  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.0, 'Max', 0.98, 'Value', 0.00, ...
                       'Units', 'normalized', 'Position', [0.05, 0.46, 0.90, 0.03]);

    lbl_wp = uicontrol('Parent', pan, 'Style', 'text', 'String', 'White Point: 1.00', ...
                       'Units', 'normalized', 'Position', [0.05, 0.41, 0.90, 0.03], ...
                       'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_wp  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.02, 'Max', 1.00, 'Value', 1.00, ...
                       'Units', 'normalized', 'Position', [0.05, 0.38, 0.90, 0.03]);

    lbl_gam = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Gamma (\gamma): 1.00', ...
                        'Units', 'normalized', 'Position', [0.05, 0.33, 0.90, 0.03], ...
                        'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_gam  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.10, 'Max', 3.50, 'Value', 1.00, ...
                        'Units', 'normalized', 'Position', [0.05, 0.30, 0.90, 0.03]);

    % --- SECTION 3: EXPORT SUITE ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 3. EXPORT CONFIGURATION ===', ...
              'Units', 'normalized', 'Position', [0.05, 0.24, 0.90, 0.03], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.90 0.85 0.85]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Base Filename:', ...
              'Units', 'normalized', 'Position', [0.05, 0.19, 0.40, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    ed_name = uicontrol('Parent', pan, 'Style', 'edit', 'String', 'PrettyMap_Export', ...
                        'Units', 'normalized', 'Position', [0.45, 0.19, 0.50, 0.035], 'BackgroundColor', 'w');

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Shorter Edge (px):', ...
              'Units', 'normalized', 'Position', [0.05, 0.14, 0.40, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    ed_px = uicontrol('Parent', pan, 'Style', 'edit', 'String', '1500', ...
                      'Units', 'normalized', 'Position', [0.45, 0.14, 0.22, 0.035], 'BackgroundColor', 'w');

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'DPI:', ...
              'Units', 'normalized', 'Position', [0.70, 0.14, 0.10, 0.03], ...
              'HorizontalAlignment', 'center', 'BackgroundColor', pan.BackgroundColor);
    ed_dpi = uicontrol('Parent', pan, 'Style', 'edit', 'String', '150', ...
                       'Units', 'normalized', 'Position', [0.80, 0.14, 0.15, 0.035], 'BackgroundColor', 'w');

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'File Format:', ...
              'Units', 'normalized', 'Position', [0.05, 0.09, 0.40, 0.03], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    pu_fmt = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', {'.jpg', '.png', '.tif'}, ...
                       'Units', 'normalized', 'Position', [0.45, 0.08, 0.50, 0.04]);

    btn_save = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'SAVE PRETTY PLOT + MAT & MD', ...
                         'Units', 'normalized', 'Position', [0.05, 0.02, 0.90, 0.05], ...
                         'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [0.30 0.75 0.35], 'ForegroundColor', 'w');

    %% 4. Initialize State Archive & Callbacks
    app = struct();
    app.fig = fig; app.ax = ax; app.pan = pan;
    app.ebsd_grid = ebsd_grid; app.grains = grains;
    app.exportDir = exportDir;
    app.x_range   = [min(ebsd_grid.x(:)), max(ebsd_grid.x(:))];
    app.y_range   = [min(ebsd_grid.y(:)), max(ebsd_grid.y(:))];
    app.hImage    = []; app.hBoundaries = [];

    % UI Widget handles
    ui = struct('pu_ipf', pu_ipf, 'sl_alpha', sl_alpha, 'lbl_alpha', lbl_alpha, ...
                'pu_gb', pu_gb, 'pu_bg', pu_bg, 'sl_bp', sl_bp, 'lbl_bp', lbl_bp, ...
                'sl_wp', sl_wp, 'lbl_wp', lbl_wp, 'sl_gam', sl_gam, 'lbl_gam', lbl_gam, ...
                'ed_name', ed_name, 'ed_px', ed_px, 'ed_dpi', ed_dpi, 'pu_fmt', pu_fmt, 'btn_save', btn_save);
    app.ui = ui;
    
    fig.UserData = app;

    % Hook up Standard Callbacks for dropdowns and actions
    set(pu_ipf, 'Callback', @(s, e) onIpfChange(fig));
    set(pu_bg,  'Callback', @(s, e) onBgChange(fig));
    set(pu_gb,  'Callback', @(s, e) onGbChange(fig));
    set(btn_save, 'Callback', @(s, e) onSaveExport(fig));

    % Attach Real-Time Continuous PostSet Listeners to all sliders (Zero-Lag scrubbing)
    addlistener(sl_alpha, 'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_bp,    'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_wp,    'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_gam,   'Value', 'PostSet', @(s, e) onSliderChange(fig));

    % Perform First Initialization
    onIpfChange(fig, false);
    onBgChange(fig, true);
    onGbChange(fig);
    fprintf('Studio Ready! Adjust sliders and menus for instantaneous rendering.\n');
end

%% ========================================================================
%  GUI CALLBACKS & ZERO-LAG RENDER ENGINE
%  ========================================================================

function onIpfChange(fig, doUpdate)
    if nargin < 2, doUpdate = true; end
    app = fig.UserData;
    idx = get(app.ui.pu_ipf, 'Value');
    keys = get(app.ui.pu_ipf, 'String');
    sel_key = keys{idx};
    
    try
        oM = feval(sel_key, app.ebsd_grid);
        try, oM.inversePoleFigureDirection = zvector; catch, end
        ipf_rgb = oM.orientation2color(app.ebsd_grid.orientations);
    catch ME
        fprintf(2, 'Could not instantiate %s (%s). Using ipfColorKey.\n', sel_key, ME.message);
        oM = ipfColorKey(app.ebsd_grid);
        ipf_rgb = oM.orientation2color(app.ebsd_grid.orientations);
    end
    
    % Guarantee 3D RGB array [M x N x 3]
    if size(ipf_rgb, 3) ~= 3 && size(ipf_rgb, 2) == 3
        ipf_rgb = reshape(ipf_rgb, [size(app.ebsd_grid, 1), size(app.ebsd_grid, 2), 3]);
    end
    
    app.ipf_rgb = double(ipf_rgb);
    fig.UserData = app;
    if doUpdate, updatePlot(fig); end
end

function onBgChange(fig, doUpdate)
    if nargin < 2, doUpdate = true; end
    app = fig.UserData;
    idx = get(app.ui.pu_bg, 'Value');
    fields_list = get(app.ui.pu_bg, 'String');
    sel_field = fields_list{idx};
    
    if isfield(app.ebsd_grid.prop, sel_field)
        val = app.ebsd_grid.prop.(sel_field);
    elseif isprop(app.ebsd_grid, sel_field)
        val = app.ebsd_grid.(sel_field);
    else
        val = app.ebsd_grid.ci;
    end
    
    % Handle structural / vector data gracefully
    if isa(val, 'vector3d')
        val = val.norm;
    elseif isstruct(val) && isfield(val, 'x') && isfield(val, 'y')
        val = sqrt(val.x.^2 + val.y.^2);
    end
    
    val = double(val);
    if ismatrix(val) && ~isequal(size(val), [size(app.ebsd_grid,1), size(app.ebsd_grid,2)])
        try, val = reshape(val, [size(app.ebsd_grid,1), size(app.ebsd_grid,2)]); catch, end
    end
    
    app.bg_raw = val;
    
    % Compute baseline 1st and 99th percentiles for robust contrast bounds
    valid_vals = val(~isnan(val) & ~isinf(val));
    if isempty(valid_vals)
        app.bg_min = 0; app.bg_max = 1;
    else
        app.bg_min = prctile(valid_vals, 1);
        app.bg_max = prctile(valid_vals, 99);
        if app.bg_min == app.bg_max, app.bg_max = app.bg_min + 1e-5; end
    end
    
    fig.UserData = app;
    if doUpdate, updatePlot(fig); end
end

function onGbChange(fig)
    app = fig.UserData;
    idx = get(app.ui.pu_gb, 'Value');
    styles = get(app.ui.pu_gb, 'String');
    sel_style = styles{idx};
    
    % Delete existing boundaries if present
    if ~isempty(app.hBoundaries) && isvalid(app.hBoundaries)
        delete(app.hBoundaries);
        app.hBoundaries = [];
    end
    
    if contains(sel_style, 'Hidden') || contains(sel_style, 'None')
        fig.UserData = app;
        return;
    end
    
    if contains(sel_style, 'White'), col = 'white'; else, col = 'black'; end
    if contains(sel_style, '1.2'), lw = 1.2; else, lw = 2.0; end
    
    hold(app.ax, 'on');
    app.hBoundaries = plot(app.grains.boundary, 'Parent', app.ax, 'lineColor', col, 'lineWidth', lw, 'micronbar', 'off');
    hold(app.ax, 'off');
    
    fig.UserData = app;
end

function onSliderChange(fig)
    updatePlot(fig);
end

function updatePlot(fig)
    app = fig.UserData;
    
    % Retrieve active slider parameters
    alpha_val = get(app.ui.sl_alpha, 'Value');
    bp_norm   = get(app.ui.sl_bp,    'Value');
    wp_norm   = get(app.ui.sl_wp,    'Value');
    gamma_val = get(app.ui.sl_gam,   'Value');
    
    % Ensure valid threshold order
    if wp_norm <= bp_norm
        wp_norm = min(1.0, bp_norm + 0.02);
        set(app.ui.sl_wp, 'Value', wp_norm);
    end
    
    % Update Text Labels in real-time
    set(app.ui.lbl_alpha, 'String', sprintf('IPF Overlay Alpha: %.2f', alpha_val));
    set(app.ui.lbl_bp,    'String', sprintf('Black Point: %.2f', bp_norm));
    set(app.ui.lbl_wp,    'String', sprintf('White Point: %.2f', wp_norm));
    set(app.ui.lbl_gam,   'String', sprintf('Gamma (\\gamma): %.2f', gamma_val));
    
    % Map normalized slider bounds to physical data values
    rng = app.bg_max - app.bg_min;
    phys_bp = app.bg_min + bp_norm * rng;
    phys_wp = app.bg_min + wp_norm * rng;
    if phys_wp <= phys_bp, phys_wp = phys_bp + 1e-5; end
    
    % Normalize and apply gamma correction to greyscale background
    norm_bg = (app.bg_raw - phys_bp) / (phys_wp - phys_bp);
    norm_bg = min(max(norm_bg, 0), 1);
    gamma_bg = norm_bg .^ (1.0 / max(0.01, gamma_val));
    
    bg_rgb = repmat(gamma_bg, [1, 1, 3]);
    
    % Handle NaN / Unindexed regions cleanly
    nan_mask = isnan(app.bg_raw) | isinf(app.bg_raw);
    nan_mask_3d = repmat(nan_mask, [1, 1, 3]);
    bg_rgb(nan_mask_3d) = 0.85; % neutral light grey background
    
    % Blend background with IPF overlay
    final_rgb = (1 - alpha_val) * bg_rgb + alpha_val * app.ipf_rgb;
    final_rgb(nan_mask_3d) = 0.85; % retain clean background in unindexed gaps
    
    % ZERO-LAG INJECTION: Update image array buffer directly without object rebuilding
    if isempty(app.hImage) || ~isvalid(app.hImage)
        app.hImage = image(app.ax, app.x_range, app.y_range, final_rgb);
        set(app.ax, 'YDir', 'normal');
        axis(app.ax, 'image');
        % Ensure image sits behind grain boundaries
        uistack(app.hImage, 'bottom');
    else
        app.hImage.CData = final_rgb;
    end
    
    % Update Figure Title with currently active layers
    idx_bg = get(app.ui.pu_bg, 'Value'); bg_names = get(app.ui.pu_bg, 'String');
    idx_ipf = get(app.ui.pu_ipf, 'Value'); ipf_names = get(app.ui.pu_ipf, 'String');
    title(app.ax, sprintf('Composite: %s + %s (\alpha=%.2f, \gamma=%.2f)', ...
          strrep(bg_names{idx_bg}, '_', '\_'), ipf_names{idx_ipf}, alpha_val, gamma_val), ...
          'FontSize', 13, 'FontWeight', 'bold');
          
    app.curr_rgb = final_rgb;
    fig.UserData = app;
end

%% ========================================================================
%  EXPORT SUITE (SAVING HIGH-RES IMAGE, MAT & MD DOCUMENTATION)
%  ========================================================================

function onSaveExport(fig)
    app = fig.UserData;
    ui  = app.ui;
    
    % Read UI configurations
    base_name = strtrim(get(ui.ed_name, 'String'));
    if isempty(base_name), base_name = 'PrettyMap'; end
    
    target_px = str2double(get(ui.ed_px, 'String'));
    if isnan(target_px) || target_px <= 0, target_px = 1500; end
    
    dpi_val = str2double(get(ui.ed_dpi, 'String'));
    if isnan(dpi_val) || dpi_val <= 0, dpi_val = 150; end
    
    idx_fmt = get(ui.pu_fmt, 'Value'); fmts = get(ui.pu_fmt, 'String');
    ext_fmt = fmts{idx_fmt};
    
    idx_bg  = get(ui.pu_bg, 'Value');  bg_list  = get(ui.pu_bg, 'String');  sel_bg  = bg_list{idx_bg};
    idx_ipf = get(ui.pu_ipf, 'Value'); ipf_list = get(ui.pu_ipf, 'String'); sel_ipf = ipf_list{idx_ipf};
    idx_gb  = get(ui.pu_gb, 'Value');  gb_list  = get(ui.pu_gb, 'String');  sel_gb  = gb_list{idx_gb};
    
    alpha_val = get(ui.sl_alpha, 'Value');
    bp_val    = get(ui.sl_bp,    'Value');
    wp_val    = get(ui.sl_wp,    'Value');
    gamma_val = get(ui.sl_gam,   'Value');
    
    % Build filename tag
    timestamp_str = datestr(now, 'yyyy-mm-dd_HHMMSS');
    file_tag  = sprintf('%s_%s_%s_a%03d_g%03d_%s', base_name, sel_bg, sel_ipf, ...
                        round(alpha_val*100), round(gamma_val*100), timestamp_str);
    clean_tag = strrep(strrep(file_tag, '.', '_'), ' ', '');
    export_path = fullfile(app.exportDir, clean_tag);
    
    fprintf('\n========================================================================\n');
    fprintf('SAVING PUBLICATION-GRADE COMPOSITE EXPORT...\n');
    fprintf('  Target path: %s%s\n', export_path, ext_fmt);
    
    % 1. Create clean, dedicated export window (without interactive sidebar clutter)
    figExp = figure('Name', 'Export Renderer', 'Visible', 'off', 'Position', [100, 100, 1000, 850], 'Color', 'w');
    axExp  = axes('Parent', figExp, 'Units', 'normalized', 'Position', [0.08, 0.08, 0.88, 0.86]);
    
    image(axExp, app.x_range, app.y_range, app.curr_rgb);
    set(axExp, 'YDir', 'normal');
    axis(axExp, 'image'); box(axExp, 'on');
    xlabel(axExp, 'x [\mum]', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel(axExp, 'y [\mum]', 'FontWeight', 'bold', 'FontSize', 12);
    title(axExp, sprintf('EBSD Composite: %s + %s (\\alpha=%.2f, \\gamma=%.2f)', ...
          strrep(sel_bg, '_', '\_'), sel_ipf, alpha_val, gamma_val), ...
          'FontSize', 14, 'FontWeight', 'bold');
          
    if ~contains(sel_gb, 'Hidden') && ~contains(sel_gb, 'None')
        if contains(sel_gb, 'White'), c_col = 'white'; else, c_col = 'black'; end
        if contains(sel_gb, '1.2'), c_lw = 1.2; else, c_lw = 2.0; end
        hold(axExp, 'on');
        plot(app.grains.boundary, 'Parent', axExp, 'lineColor', c_col, 'lineWidth', c_lw, 'micronbar', 'off');
        hold(axExp, 'off');
    end
    
    % Export figure via exportScaledFigure
    img_filepath = [export_path, ext_fmt];
    exportScaledFigure(figExp, img_filepath, 'target_px', target_px, 'dpi', dpi_val, 'nodatestring');
    close(figExp);
    fprintf('  -> Successfully exported visualization image: %s\n', img_filepath);

    % 2. Save MAT Parameter Archive
    params = struct();
    params.base_name     = base_name;
    params.background_field = sel_bg;
    params.ipf_color_key = sel_ipf;
    params.alpha_overlay = alpha_val;
    params.black_point   = bp_val;
    params.white_point   = wp_val;
    params.gamma_val     = gamma_val;
    params.grain_boundary_style = sel_gb;
    params.target_shorter_edge_px = target_px;
    params.dpi_val       = dpi_val;
    params.export_time   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    params.final_rgb_map = app.curr_rgb;
    
    mat_filepath = [export_path, '_Params.mat'];
    save(mat_filepath, 'params', '-v7.3');
    fprintf('  -> Successfully archived parameters and RGB buffer: %s\n', mat_filepath);

    % 3. Generate Interactive Markdown Documentation Report
    md_filepath = [export_path, '_Report.md'];
    fid = fopen(md_filepath, 'w');
    fprintf(fid, '# Interactive EBSD & DIC Visualization Export Report\n\n');
    fprintf(fid, 'This report details the exact visualization parameters and colormap blending algorithms used to generate the composite publication figure.\n\n');
    fprintf(fid, '## Exported Visualization Preview\n');
    fprintf(fid, '![Composite Map](%s%s)\n\n---\n\n', clean_tag, ext_fmt);
    
    fprintf(fid, '## Layer & Blending Configuration\n\n');
    fprintf(fid, '| Parameter | Selected Value | Functional Description |\n');
    fprintf(fid, '| :--- | :--- | :--- |\n');
    fprintf(fid, '| **Background Field** | `%s` | Extracted from EBSD dataset properties |\n', sel_bg);
    fprintf(fid, '| **Black Point** | `%.3f` | Normalized scalar floor threshold ($0.00 = 1\\%%$ percentile) |\n', bp_val);
    fprintf(fid, '| **White Point** | `%.3f` | Normalized scalar ceiling threshold ($1.00 = 99\\%%$ percentile) |\n', wp_val);
    fprintf(fid, '| **Gamma ($\\\\gamma$)** | `%.2f` | Non-linear luminance midpoint scaling |\n', gamma_val);
    fprintf(fid, '| **IPF Color Key** | `%s` | MTEX Inverse Pole Figure coloring scheme |\n', sel_ipf);
    fprintf(fid, '| **IPF Alpha Overlay** | `%.2f` | Blending weight ($0.00 = \\\\text{pure grey}, 1.00 = \\\\text{pure IPF}$) |\n', alpha_val);
    fprintf(fid, '| **Grain Boundaries** | `%s` | Boundary plotting style |\n\n---\n\n', sel_gb);
    
    fprintf(fid, '## Image Resolution & Export Specifications\n');
    fprintf(fid, '- **Target Shorter Edge:** `%d` px\n', target_px);
    fprintf(fid, '- **Raster DPI:** `%d`\n', dpi_val);
    fprintf(fid, '- **File Format:** `%s`\n', ext_fmt);
    fprintf(fid, '- **MAT Parameter File:** [`%s`](%s_Params.mat)\n', [clean_tag, '_Params.mat'], clean_tag);
    fprintf(fid, '- **Export Timestamp:** `%s`\n', params.export_time);
    fclose(fid);
    fprintf('  -> Successfully formatted Markdown report: %s\n', md_filepath);
    fprintf('========================================================================\n');
    
    try
        msgbox(sprintf('Successfully saved:\n1) %s\n2) %s\n3) %s\n\nin:\n%s', ...
               [clean_tag, ext_fmt], [clean_tag, '_Params.mat'], [clean_tag, '_Report.md'], app.exportDir), ...
               'Export Complete');
    catch, end
end
