function varargout = Interactive_PrettyPlotter(ebsd, grains, varargin)
% =========================================================================
% Interactive_PrettyPlotter.m
% =========================================================================
% INTERACTIVE REAL-TIME EBSD & TKD PRETTY PLOTTER FUNCTIONAL STUDIO
%
% Modern, modular visualization component following clean separation of
% concerns. Accepts a domain EBSD object and renders an interactive, zero-lag
% plotting studio with customizable greyscale background curves, IPF color
% keys, independent property thresholding, non-destructive specimen masking,
% transparency blending, and publication-ready batch export engines.
%
% SYNTAX:
%   studio = Interactive_PrettyPlotter(ebsd, grains, options) % Programmatic handle
%   [status, opts] = Interactive_PrettyPlotter(ebsd, grains, options) % Batch loop mode
%
% INPUTS:
%   ebsd    : Required MTEX EBSD object (gridified or ungridified). All 
%             imported scalar layers (BS, BC, FSE, DIC gradients) should 
%             exist in ebsd.prop.
%   grains  : Optional MTEX grain2d object for boundary overlays. If empty
%             or omitted, grains (>300 px, 5°) are computed lazily.
%   options : Structural struct or Name-Value pairs for defaults & persistence:
%             - datasetName          : Tag used in title and export filenames.
%             - exportDir            : Folder destination for visualizations.
%             - defaultBgField       : Initial greyscale property (e.g. 'BS').
%             - defaultIpfKey        : Initial color scheme (default 'ipfColorKey').
%             - defaultIpfDir        : Initial reference projection ('Y', 'Z', or 'X').
%             - defaultIpfThreshProp : Initial IPF masking property (default 'ci').
%             - mask_spec / crop_spec: Enable specimenlocation mask and crop toggles.
%             - batchMode            : Boolean flag to pause via uiwait for batch loops.
%
% OUTPUTS:
%   status  : Workflow signaling string ('next', 'skip', or 'exit').
%   opts    : Updated parameter structure preserving user visual settings across batches.
%
% Author: Philipp / Antigravity
% Date:   2026-08-06
% =========================================================================

    if nargin < 1 || isempty(ebsd)
        error('[Interactive_PrettyPlotter] Error: A valid MTEX EBSD object is required as the first argument.');
    end
    if nargin < 2
        grains = [];
    end

    % --- Parse Configuration Options ---
    opts = struct();
    opts.datasetName = 'EBSD_Composite_Map';
    opts.exportDir   = fullfile(pwd, 'PrettyPlots_Interactive');
    opts.defaultBgField       = 'ci';
    opts.defaultIpfKey        = 'ipfColorKey';
    opts.defaultIpfDir        = 'Y';
    opts.defaultIpfThreshProp = 'ci';
    opts.mask_spec = true;
    opts.crop_spec = true;
    opts.batchMode = false;
    
    if nargin > 2
        if isstruct(varargin{1})
            flds = fields(varargin{1});
            for k = 1:length(flds), opts.(flds{k}) = varargin{1}.(flds{k}); end
        else
            for k = 1:2:length(varargin)
                if k+1 <= length(varargin), opts.(varargin{k}) = varargin{k+1}; end
            end
        end
    end

    if ~exist(opts.exportDir, 'dir')
        mkdir(opts.exportDir);
    end

    % Ensure standard MTEX spatial plotting conventions
    plotx2east;
    plotzIntoPlane;

    % Probe active MTEX plotting convention to determine correct MATLAB axes directions
    ydir_str = 'reverse'; % Default for standard microscopy convention (y↓→x)
    xdir_str = 'normal';
    try
        if plottingConvention.default.north.y < 0
            ydir_str = 'reverse';
        else
            ydir_str = 'normal';
        end
        if plottingConvention.default.east.x < 0
            xdir_str = 'reverse';
        else
            xdir_str = 'normal';
        end
    catch
    end
    fprintf('[Interactive_PrettyPlotter] Aligned spatial axes with active MTEX plotting convention (XDir: %s, YDir: %s).\n', xdir_str, ydir_str);

    % --- Ensure Gridified Spatial Data & Grain Topology ---
    if ~isa(ebsd, 'EBSDgrid')
        try, ebsd_grid = ebsd.gridify; catch, ebsd_grid = ebsd; end
    else
        ebsd_grid = ebsd;
    end

    if isempty(grains)
        fprintf('[Interactive_PrettyPlotter] Lazy grain boundary computation (>300 px, 5°)...\n');
        try, [grains, ~] = calcGrains(ebsd_grid, 'minPixel', 300, 'angle', 5*degree); catch, grains = []; end
    end

    % Check for non-destructive specimen location mask property
    if isfield(ebsd_grid.prop, 'specimen_location')
        spec_loc = logical(ebsd_grid.prop.specimen_location);
    elseif isfield(ebsd_grid.prop, 'specimenlocation')
        spec_loc = logical(ebsd_grid.prop.specimenlocation);
    else
        spec_loc = [];
    end

    % Identify available numerical property layers in ebsd_grid.prop
    prop_names = fields(ebsd_grid.prop);
    valid_fields = {};
    for k = 1:length(prop_names)
        fn = prop_names{k}; val = ebsd_grid.prop.(fn);
        if isnumeric(val) || isa(val, 'double') || isa(val, 'vector3d') || isstruct(val)
            valid_fields{end+1} = fn; %#ok<AGROW>
        end
    end

    % Organize top diagnostic layer names near top of dropdown list (pure string operations)
    priority_order = {'BS', 'BC', 'ci', 'bc', 'FSE', 'BSE', 'MSE', 'GradMag', 'LapVar', 'imagequality', 'iq'};
    for k = length(priority_order):-1:1
        pn = priority_order{k};
        if any(strcmpi(valid_fields, pn))
            valid_fields = [pn, valid_fields(~strcmpi(valid_fields, pn))];
        end
    end
    if isempty(valid_fields), valid_fields = {'ci'}; end

    % Set initial background default if available
    default_bg_idx = find(strcmpi(valid_fields, opts.defaultBgField), 1);
    if isempty(default_bg_idx), default_bg_idx = 1; end

    %% 1. Build Studio GUI Layout
    fig_w = 1280; fig_h = 880;
    scr = get(0, 'ScreenSize');
    pos_x = max(20, round((scr(3) - fig_w)/2));
    pos_y = max(20, round((scr(4) - fig_h)/2));

    title_str = sprintf('Interactive Pretty Plotter Studio | %s', opts.datasetName);
    if isfield(opts, 'currentStep') && isfield(opts, 'totalSteps')
        title_str = sprintf('[Batch %d/%d] Studio | %s', opts.currentStep, opts.totalSteps, opts.datasetName);
    end

    fig = figure('Name', title_str, ...
                 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure', ...
                 'Position', [pos_x, pos_y, fig_w, fig_h], 'Color', [0.93 0.93 0.93], ...
                 'Resize', 'on');

    % Left Map Canvas
    ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.04, 0.06, 0.65, 0.90]);
    axis(ax, 'image'); hold(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'x [\mum]', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel(ax, 'y [\mum]', 'FontWeight', 'bold', 'FontSize', 12);

    % Right Sidebar Studio Controls
    pan = uipanel('Parent', fig, 'Title', 'VISUALIZATION & BATCH STUDIO', ...
                  'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [0.96 0.96 0.97], ...
                  'Units', 'normalized', 'Position', [0.71, 0.02, 0.28, 0.96]);

    % Dataset metadata panel display
    try, min_name = ebsd_grid.CS.mineral; catch, min_name = 'Unknown'; end
    idx_cnt = sum(ebsd_grid.isIndexed(:));
    info_str = sprintf('Name: %s\nMineral: %s | Indexed Px: %s\nProperties: %d fields ready', ...
               opts.datasetName, min_name, num2str(idx_cnt), length(valid_fields));
    uicontrol('Parent', pan, 'Style', 'text', 'String', info_str, ...
              'Units', 'normalized', 'Position', [0.04, 0.915, 0.92, 0.075], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'FontWeight', 'bold', 'FontSize', 9.5);

    btn_import = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'IMPORT OPTIONS (.MAT / .MD)...', ...
                           'Units', 'normalized', 'Position', [0.04, 0.875, 0.92, 0.034], ...
                           'FontSize', 9.5, 'FontWeight', 'bold', 'BackgroundColor', [0.35 0.45 0.65], 'ForegroundColor', 'w', ...
                           'TooltipString', 'Load layer parameters and threshold configurations from a previously exported _Params.mat archive or _Report.md report.');

    % --- SECTION 1: COLOR OVERLAY CONTROLS ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 1. COLOR OVERLAY & MASKS ===', ...
              'Units', 'normalized', 'Position', [0.04, 0.845, 0.92, 0.025], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.85 0.85 0.90]);

    % Overlay Mode Selection & Colormap Limits
    ov_modes = {'1. IPF Orientation Color', '2. Scalar Property / Strain Map'};
    pu_ov_mode = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', ov_modes, ...
                           'Units', 'normalized', 'Position', [0.04, 0.815, 0.42, 0.028], 'FontSize', 9.5, ...
                           'TooltipString', 'Select between IPF orientation maps or scalar property colormaps (e.g. DIC strain fields).');
    
    lbl_cmin = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Min:', ...
                         'Units', 'normalized', 'Position', [0.45, 0.815, 0.07, 0.024], ...
                         'HorizontalAlignment', 'right', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'Visible', 'off');
    ed_cmin  = uicontrol('Parent', pan, 'Style', 'edit', 'String', '0.00', ...
                         'Units', 'normalized', 'Position', [0.53, 0.815, 0.10, 0.026], 'BackgroundColor', 'w', 'Visible', 'off');
    lbl_cmax = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Max:', ...
                         'Units', 'normalized', 'Position', [0.64, 0.815, 0.07, 0.024], ...
                         'HorizontalAlignment', 'right', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'Visible', 'off');
    ed_cmax  = uicontrol('Parent', pan, 'Style', 'edit', 'String', '1.00', ...
                         'Units', 'normalized', 'Position', [0.72, 0.815, 0.10, 0.026], 'BackgroundColor', 'w', 'Visible', 'off');
    btn_prop_auto = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'AUTO', ...
                              'Units', 'normalized', 'Position', [0.84, 0.813, 0.13, 0.029], ...
                              'FontSize', 9.0, 'FontWeight', 'bold', 'BackgroundColor', [0.25 0.55 0.75], 'ForegroundColor', 'w', 'Visible', 'off', ...
                              'TooltipString', 'Automatically crop colormap limits to the 2nd and 98th percentiles.');

    % Row 2: IPF controls (visible when Mode == 1) vs Property controls (visible when Mode == 2)
    ipf_keys = {'ipfColorKey', 'ipfHKLKey', 'ipfHSVKey', 'ipfSpotKey', 'ipfTSLKey'};
    pu_ipf   = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', ipf_keys, ...
                         'Units', 'normalized', 'Position', [0.04, 0.785, 0.55, 0.030], 'FontSize', 10);
    
    ipf_dirs = {'Y (Loading)', 'Z (Normal)', 'X (East)'};
    pu_dir   = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', ipf_dirs, ...
                         'Units', 'normalized', 'Position', [0.61, 0.785, 0.35, 0.030], 'FontSize', 10);
    
    if strcmpi(opts.defaultIpfDir, 'Z'), set(pu_dir, 'Value', 2);
    elseif strcmpi(opts.defaultIpfDir, 'X'), set(pu_dir, 'Value', 3);
    else, set(pu_dir, 'Value', 1); end

    pu_prop_ov = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', valid_fields, ...
                           'Units', 'normalized', 'Position', [0.04, 0.785, 0.52, 0.030], 'FontSize', 10, 'Visible', 'off', ...
                           'TooltipString', 'Select property or strain field to color as overlay map.');
    cmap_list  = {'viridis', 'jet', 'parula', 'turbo', 'inferno', 'magma', 'plasma', 'hot', 'cool', 'hsv', 'gray'};
    pu_cmap    = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', cmap_list, ...
                           'Units', 'normalized', 'Position', [0.58, 0.785, 0.38, 0.030], 'FontSize', 10, 'Visible', 'off', ...
                           'TooltipString', 'Select colormap for scalar property overlay.');
    
    if isfield(opts, 'defaultOverlayMode') && (contains(lower(opts.defaultOverlayMode), 'scalar') || contains(lower(opts.defaultOverlayMode), 'strain') || contains(lower(opts.defaultOverlayMode), 'property') || (isnumeric(opts.defaultOverlayMode) && opts.defaultOverlayMode == 2))
        set(pu_ov_mode, 'Value', 2);
        set(pu_ipf, 'Visible', 'off'); set(pu_dir, 'Visible', 'off');
        set(pu_prop_ov, 'Visible', 'on'); set(pu_cmap, 'Visible', 'on');
        set(lbl_cmin, 'Visible', 'on'); set(ed_cmin, 'Visible', 'on');
        set(lbl_cmax, 'Visible', 'on'); set(ed_cmax, 'Visible', 'on');
        set(btn_prop_auto, 'Visible', 'on');
    end
    if isfield(opts, 'defaultPropOverlay')
        idx_p = find(strcmpi(valid_fields, opts.defaultPropOverlay), 1);
        if ~isempty(idx_p), set(pu_prop_ov, 'Value', idx_p); end
    elseif any(strcmpi(valid_fields, 'geff'))
        set(pu_prop_ov, 'Value', find(strcmpi(valid_fields, 'geff'), 1));
    end
    if isfield(opts, 'defaultColormap')
        idx_c = find(strcmpi(cmap_list, opts.defaultColormap), 1);
        if ~isempty(idx_c), set(pu_cmap, 'Value', idx_c); end
    end
    if isfield(opts, 'cmin_val') && ~isnan(opts.cmin_val), try, set(ed_cmin, 'String', num2str(opts.cmin_val)); catch, end; end
    if isfield(opts, 'cmax_val') && ~isnan(opts.cmax_val), try, set(ed_cmax, 'String', num2str(opts.cmax_val)); catch, end; end

    lbl_alpha = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Overlay Alpha: 0.55', ...
                          'Units', 'normalized', 'Position', [0.04, 0.760, 0.54, 0.026], ...
                          'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    chk_above_bp = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Show > BlackPt', ...
                             'Units', 'normalized', 'Position', [0.60, 0.760, 0.38, 0.026], ...
                             'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9, ...
                             'TooltipString', 'Only display overlay colors where background is above the Black Point threshold.');

    sl_alpha  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.0, 'Max', 1.0, 'Value', 0.55, ...
                          'Units', 'normalized', 'Position', [0.04, 0.735, 0.54, 0.022]);
    chk_below_wp = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Show < WhitePt', ...
                             'Units', 'normalized', 'Position', [0.60, 0.735, 0.38, 0.026], ...
                             'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9, ...
                             'TooltipString', 'Only display overlay colors where background is below the White Point threshold.');

    % Independent Overlay Property Thresholding (Does not affect background)
    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Overlay Property Filter:', ...
              'Units', 'normalized', 'Position', [0.04, 0.695, 0.42, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    pu_ipf_thresh = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', valid_fields, ...
                              'Units', 'normalized', 'Position', [0.47, 0.695, 0.49, 0.032], 'FontSize', 10, ...
                              'TooltipString', 'Select property to threshold overlay colors independently of greyscale background.');
    idx_th = find(strcmpi(valid_fields, opts.defaultIpfThreshProp), 1);
    if isempty(idx_th), idx_th = find(strcmpi(valid_fields, 'ci'), 1); end
    if ~isempty(idx_th), set(pu_ipf_thresh, 'Value', idx_th); end

    lbl_ipf_thresh = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Min IPF Thresh: 0.000', ...
                               'Units', 'normalized', 'Position', [0.04, 0.665, 0.92, 0.025], ...
                               'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_ipf_thresh  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.0, 'Max', 1.0, 'Value', 0.0, ...
                               'Units', 'normalized', 'Position', [0.04, 0.64, 0.92, 0.022]);

    % Non-Destructive Specimen Location Mask & Map Cropping
    chk_mask_spec = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Mask IPF to Specimen', ...
                              'Units', 'normalized', 'Position', [0.04, 0.60, 0.47, 0.028], ...
                              'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9, ...
                              'Value', (opts.mask_spec && ~isempty(spec_loc)), ...
                              'TooltipString', 'Mask out IPF colors outside specimen location while preserving greyscale background.');
    chk_crop_spec = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Crop Map to Specimen', ...
                              'Units', 'normalized', 'Position', [0.52, 0.60, 0.44, 0.028], ...
                              'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9, ...
                              'Value', (opts.crop_spec && ~isempty(spec_loc)), ...
                              'TooltipString', 'Dynamically crop spatial view box to specimen location.');
    if isempty(spec_loc)
        set(chk_mask_spec, 'Enable', 'off', 'String', 'Mask (No Specimen Mask)');
        set(chk_crop_spec, 'Enable', 'off', 'String', 'Crop (No Specimen Mask)');
    end

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'GB Style:', ...
              'Units', 'normalized', 'Position', [0.04, 0.56, 0.20, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    chk_micronbar = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Micron Bar', ...
                              'Units', 'normalized', 'Position', [0.24, 0.56, 0.25, 0.028], ...
                              'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9, ...
                              'Value', 1, 'TooltipString', 'Display a dynamic scale bar.');
    gb_styles = {'White (1.2 px)', 'Black (1.2 px)', 'White (2.0 px)', 'Black (2.0 px)', 'Hidden / None'};
    pu_gb     = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', gb_styles, ...
                          'Units', 'normalized', 'Position', [0.50, 0.555, 0.46, 0.032], 'FontSize', 10);
    if isfield(opts, 'defaultGbStyle')
        gb_idx = find(strcmpi(gb_styles, opts.defaultGbStyle), 1);
        if ~isempty(gb_idx), set(pu_gb, 'Value', gb_idx); end
    end

    % --- SECTION 2: GREYSCALE BACKGROUND CONTROLS ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 2. GREYSCALE BACKGROUND ===', ...
              'Units', 'normalized', 'Position', [0.04, 0.515, 0.92, 0.025], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.85 0.90 0.85]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Background Property:', ...
              'Units', 'normalized', 'Position', [0.04, 0.48, 0.55, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    pu_bg   = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', valid_fields, 'Value', default_bg_idx, ...
                        'Units', 'normalized', 'Position', [0.04, 0.45, 0.58, 0.032], 'FontSize', 10);
    chk_inv = uicontrol('Parent', pan, 'Style', 'checkbox', 'String', 'Invert Colors', ...
                        'Units', 'normalized', 'Position', [0.65, 0.45, 0.32, 0.032], ...
                        'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold', 'FontSize', 9);

    lbl_bp = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Black Point: 0.05', ...
                       'Units', 'normalized', 'Position', [0.04, 0.415, 0.73, 0.025], ...
                       'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_bp  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.0, 'Max', 0.95, 'Value', 0.05, ...
                       'Units', 'normalized', 'Position', [0.04, 0.39, 0.73, 0.022]);

    lbl_wp = uicontrol('Parent', pan, 'Style', 'text', 'String', 'White Point: 0.98', ...
                       'Units', 'normalized', 'Position', [0.04, 0.355, 0.73, 0.025], ...
                       'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_wp  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.05, 'Max', 1.00, 'Value', 0.98, ...
                       'Units', 'normalized', 'Position', [0.04, 0.33, 0.73, 0.022]);

    btn_auto = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'AUTO (98%)', ...
                         'Units', 'normalized', 'Position', [0.79, 0.33, 0.17, 0.082], ...
                         'FontSize', 9.5, 'FontWeight', 'bold', 'BackgroundColor', [0.25 0.55 0.75], 'ForegroundColor', 'w', ...
                         'TooltipString', 'Automatically stretch greyscale background contrast between the 2nd and 98th percentiles.');

    lbl_gam = uicontrol('Parent', pan, 'Style', 'text', 'String', 'Gamma (\gamma): 1.50', ...
                        'Units', 'normalized', 'Position', [0.04, 0.295, 0.92, 0.025], ...
                        'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor, 'FontWeight', 'bold');
    sl_gam  = uicontrol('Parent', pan, 'Style', 'slider', 'Min', 0.10, 'Max', 3.50, 'Value', 1.50, ...
                        'Units', 'normalized', 'Position', [0.04, 0.27, 0.92, 0.022]);

    % --- SECTION 3: EXPORT & BATCH SUITE ---
    uicontrol('Parent', pan, 'Style', 'text', 'String', '=== 3. EXPORT & BATCH SUITE ===', ...
              'Units', 'normalized', 'Position', [0.04, 0.23, 0.92, 0.025], ...
              'FontWeight', 'bold', 'BackgroundColor', [0.90 0.85 0.85]);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Base Filename:', ...
              'Units', 'normalized', 'Position', [0.04, 0.19, 0.38, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    ed_name = uicontrol('Parent', pan, 'Style', 'edit', 'String', opts.datasetName, ...
                        'Units', 'normalized', 'Position', [0.43, 0.19, 0.53, 0.03], 'BackgroundColor', 'w', 'FontSize', 10);

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'Shorter Edge (px):', ...
              'Units', 'normalized', 'Position', [0.04, 0.15, 0.38, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    ed_px = uicontrol('Parent', pan, 'Style', 'edit', 'String', '1500', ...
                      'Units', 'normalized', 'Position', [0.43, 0.15, 0.20, 0.03], 'BackgroundColor', 'w');

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'DPI:', ...
              'Units', 'normalized', 'Position', [0.65, 0.15, 0.12, 0.028], ...
              'HorizontalAlignment', 'center', 'BackgroundColor', pan.BackgroundColor);
    ed_dpi = uicontrol('Parent', pan, 'Style', 'edit', 'String', '150', ...
                       'Units', 'normalized', 'Position', [0.78, 0.15, 0.18, 0.03], 'BackgroundColor', 'w');

    uicontrol('Parent', pan, 'Style', 'text', 'String', 'File Format:', ...
              'Units', 'normalized', 'Position', [0.04, 0.11, 0.38, 0.028], ...
              'HorizontalAlignment', 'left', 'BackgroundColor', pan.BackgroundColor);
    pu_fmt = uicontrol('Parent', pan, 'Style', 'popupmenu', 'String', {'.jpg', '.png', '.tif', '.pdf'}, ...
                       'Units', 'normalized', 'Position', [0.43, 0.105, 0.53, 0.032]);

    % Batch Navigation & Save Buttons
    btn_save_curr = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'SAVE CURRENT PRETTY PLOT', ...
                              'Units', 'normalized', 'Position', [0.04, 0.065, 0.92, 0.038], ...
                              'FontSize', 10.5, 'FontWeight', 'bold', 'BackgroundColor', [0.28 0.75 0.35], 'ForegroundColor', 'w', ...
                              'TooltipString', 'Export current plot while staying on this specimen.');

    btn_save_next = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'SAVE & NEXT ->', ...
                              'Units', 'normalized', 'Position', [0.04, 0.015, 0.44, 0.042], ...
                              'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.60 0.85], 'ForegroundColor', 'w', ...
                              'TooltipString', 'Export current plot and advance to next dataset in batch loop.');

    btn_skip      = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'SKIP ->', ...
                              'Units', 'normalized', 'Position', [0.50, 0.015, 0.23, 0.042], ...
                              'FontSize', 9.5, 'FontWeight', 'bold', 'BackgroundColor', [0.85 0.65 0.20], 'ForegroundColor', 'w', ...
                              'TooltipString', 'Skip exporting this dataset and move to next.');

    btn_exit      = uicontrol('Parent', pan, 'Style', 'pushbutton', 'String', 'EXIT', ...
                              'Units', 'normalized', 'Position', [0.75, 0.015, 0.21, 0.042], ...
                              'FontSize', 9.5, 'FontWeight', 'bold', 'BackgroundColor', [0.85 0.30 0.30], 'ForegroundColor', 'w', ...
                              'TooltipString', 'Terminate semi-automated batch processing loop.');

    %% 2. Initialize App State Architecture
    app = struct();
    app.fig = fig; app.ax = ax; app.pan = pan;
    app.ebsd_grid = ebsd_grid; app.grains = grains;
    app.opts = opts; app.datasetName = opts.datasetName;
    app.exportDir = opts.exportDir;
    app.spec_loc  = spec_loc;
    app.ipf_thresh_data = [];
    app.ipf_thresh_min  = 0;
    app.ipf_thresh_max  = 1;
    app.batch_status    = 'next';
    app.xdir            = xdir_str;
    app.ydir            = ydir_str;
    
    valid_x = ebsd_grid.x(~isnan(ebsd_grid.x) & ~isinf(ebsd_grid.x));
    valid_y = ebsd_grid.y(~isnan(ebsd_grid.y) & ~isinf(ebsd_grid.y));
    if isempty(valid_x), valid_x = [0; 100]; end
    if isempty(valid_y), valid_y = [0; 100]; end
    min_x = double(min(valid_x(:))); max_x = double(max(valid_x(:)));
    min_y = double(min(valid_y(:))); max_y = double(max(valid_y(:)));
    if max_x <= min_x, max_x = min_x + 1.0; end
    if max_y <= min_y, max_y = min_y + 1.0; end
    app.x_range = [min_x, max_x];
    app.y_range = [min_y, max_y];
    app.hImage = []; app.hBoundaries = []; app.latestExport = struct();

    ui = struct('pu_ov_mode', pu_ov_mode, 'lbl_cmin', lbl_cmin, 'ed_cmin', ed_cmin, 'lbl_cmax', lbl_cmax, 'ed_cmax', ed_cmax, 'btn_prop_auto', btn_prop_auto, ...
                'pu_ipf', pu_ipf, 'pu_dir', pu_dir, 'pu_prop_ov', pu_prop_ov, 'pu_cmap', pu_cmap, ...
                'sl_alpha', sl_alpha, 'lbl_alpha', lbl_alpha, 'chk_above_bp', chk_above_bp, 'chk_below_wp', chk_below_wp, ...
                'pu_ipf_thresh', pu_ipf_thresh, 'lbl_ipf_thresh', lbl_ipf_thresh, 'sl_ipf_thresh', sl_ipf_thresh, ...
                'chk_mask_spec', chk_mask_spec, 'chk_crop_spec', chk_crop_spec, 'chk_micronbar', chk_micronbar, ...
                'pu_gb', pu_gb, 'pu_bg', pu_bg, 'chk_inv', chk_inv, 'sl_bp', sl_bp, 'lbl_bp', lbl_bp, ...
                'sl_wp', sl_wp, 'lbl_wp', lbl_wp, 'sl_gam', sl_gam, 'lbl_gam', lbl_gam, 'btn_auto', btn_auto, ...
                'ed_name', ed_name, 'ed_px', ed_px, 'ed_dpi', ed_dpi, 'pu_fmt', pu_fmt, ...
                'btn_import', btn_import, 'btn_save_curr', btn_save_curr, 'btn_save_next', btn_save_next, 'btn_skip', btn_skip, 'btn_exit', btn_exit);
    app.ui = ui;
    fig.UserData = app;

    % Connect Callbacks & Continuous Zero-Lag PostSet Listeners
    set(btn_import,    'Callback', @(s, e) onImportOptions(fig));
    set(pu_ov_mode,    'Callback', @(s, e) onOverlayModeChange(fig));
    set(pu_ipf,        'Callback', @(s, e) onIpfChange(fig));
    set(pu_dir,        'Callback', @(s, e) onIpfChange(fig));
    set(pu_prop_ov,    'Callback', @(s, e) onIpfChange(fig));
    set(pu_cmap,       'Callback', @(s, e) onIpfChange(fig));
    set(ed_cmin,       'Callback', @(s, e) onIpfChange(fig));
    set(ed_cmax,       'Callback', @(s, e) onIpfChange(fig));
    set(btn_prop_auto, 'Callback', @(s, e) onPropAuto(fig));
    set(chk_above_bp,  'Callback', @(s, e) onSliderChange(fig));
    set(chk_below_wp,  'Callback', @(s, e) onSliderChange(fig));
    set(pu_ipf_thresh, 'Callback', @(s, e) onIpfThreshPropChange(fig));
    set(chk_mask_spec, 'Callback', @(s, e) onSliderChange(fig));
    set(chk_crop_spec, 'Callback', @(s, e) onSliderChange(fig));
    set(chk_micronbar, 'Callback', @(s, e) onSliderChange(fig));
    set(pu_bg,         'Callback', @(s, e) onBgChange(fig));
    set(chk_inv,       'Callback', @(s, e) onSliderChange(fig));
    set(pu_gb,         'Callback', @(s, e) onGbChange(fig));
    set(btn_auto,      'Callback', @(s, e) onAutoContrast(fig));
    
    set(btn_save_curr, 'Callback', @(s, e) onSaveExport(fig, false));
    set(btn_save_next, 'Callback', @(s, e) onSaveExport(fig, true));
    set(btn_skip,      'Callback', @(s, e) onBatchNavigate(fig, 'skip'));
    set(btn_exit,      'Callback', @(s, e) onBatchNavigate(fig, 'exit'));

    addlistener(sl_alpha,      'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_ipf_thresh, 'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_bp,         'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_wp,         'Value', 'PostSet', @(s, e) onSliderChange(fig));
    addlistener(sl_gam,        'Value', 'PostSet', @(s, e) onSliderChange(fig));

    % Initial rendering invocation & Persistent parameter restore
    onIpfChange(fig, false);
    onBgChange(fig, false);
    onIpfThreshPropChange(fig, false);
    
    % Restore persistent UI parameter settings across sequential batch executions
    if isfield(opts, 'alpha_val'), set(sl_alpha, 'Value', opts.alpha_val); end
    if isfield(opts, 'bp_val'),    set(sl_bp,    'Value', opts.bp_val);    end
    if isfield(opts, 'wp_val'),    set(sl_wp,    'Value', opts.wp_val);    end
    if isfield(opts, 'gamma_val'), set(sl_gam,   'Value', opts.gamma_val); end
    if isfield(opts, 'invert_bg'), set(chk_inv,  'Value', opts.invert_bg); end
    if isfield(opts, 'restrict_above_bp'), set(chk_above_bp, 'Value', opts.restrict_above_bp); end
    if isfield(opts, 'restrict_below_wp'), set(chk_below_wp, 'Value', opts.restrict_below_wp); end
    if isfield(opts, 'ipf_thresh_val')
        try
            min_v = get(sl_ipf_thresh, 'Min'); max_v = get(sl_ipf_thresh, 'Max');
            if opts.ipf_thresh_val >= min_v && opts.ipf_thresh_val <= max_v
                set(sl_ipf_thresh, 'Value', opts.ipf_thresh_val);
            end
        catch
        end
    end
    
    updatePlot(fig);
    onGbChange(fig);
    
    % --- Handle Workflow Execution Mode ---
    if nargout >= 2 || (isfield(opts, 'batchMode') && opts.batchMode)
        uiwait(fig);
        if ~isvalid(fig)
            status = 'exit';
            ret_opts = opts;
        else
            app = fig.UserData;
            status = app.batch_status;
            ret_opts = saveCurrentSettingsToOpts(fig);
            close(fig);
        end
        varargout{1} = status;
        if nargout > 1, varargout{2} = ret_opts; end
    elseif nargout == 1
        out = struct();
        out.fig = fig;
        out.ax  = ax;
        out.appData = fig.UserData;
        out.updatePlot = @() updatePlot(fig);
        out.saveExport = @() onSaveExport(fig, false);
        varargout{1} = out;
    end
end

%% ========================================================================
%  INTERNAL ZERO-LAG GPU RENDER ENGINE & UI CALLBACKS
%  ========================================================================

function updated_opts = saveCurrentSettingsToOpts(fig)
    app = fig.UserData;
    ui  = app.ui;
    updated_opts = app.opts;
    
    ov_modes = get(ui.pu_ov_mode, 'String'); updated_opts.defaultOverlayMode = ov_modes{get(ui.pu_ov_mode, 'Value')};
    prop_ov_list = get(ui.pu_prop_ov, 'String'); updated_opts.defaultPropOverlay = prop_ov_list{get(ui.pu_prop_ov, 'Value')};
    cmap_list = get(ui.pu_cmap, 'String'); updated_opts.defaultColormap = cmap_list{get(ui.pu_cmap, 'Value')};
    updated_opts.cmin_val = str2double(get(ui.ed_cmin, 'String'));
    updated_opts.cmax_val = str2double(get(ui.ed_cmax, 'String'));
    
    idx_bg  = get(ui.pu_bg,  'Value'); bg_list = get(ui.pu_bg, 'String'); updated_opts.defaultBgField = bg_list{idx_bg};
    idx_ipf = get(ui.pu_ipf, 'Value'); ipf_list = get(ui.pu_ipf, 'String'); updated_opts.defaultIpfKey = ipf_list{idx_ipf};
    idx_dir = get(ui.pu_dir, 'Value'); dir_list = get(ui.pu_dir, 'String'); updated_opts.defaultIpfDir = dir_list{idx_dir};
    idx_gb  = get(ui.pu_gb,  'Value'); gb_list = get(ui.pu_gb, 'String'); updated_opts.defaultGbStyle = gb_list{idx_gb};
    
    updated_opts.alpha_val = get(ui.sl_alpha, 'Value');
    updated_opts.bp_val    = get(ui.sl_bp,    'Value');
    updated_opts.wp_val    = get(ui.sl_wp,    'Value');
    updated_opts.gamma_val = get(ui.sl_gam,   'Value');
    updated_opts.invert_bg = (get(ui.chk_inv, 'Value') == 1);
    
    updated_opts.restrict_above_bp = (get(ui.chk_above_bp, 'Value') == 1);
    updated_opts.restrict_below_wp = (get(ui.chk_below_wp, 'Value') == 1);
    
    idx_th = get(ui.pu_ipf_thresh, 'Value'); th_list = get(ui.pu_ipf_thresh, 'String');
    updated_opts.defaultIpfThreshProp = th_list{idx_th};
    updated_opts.ipf_thresh_val       = get(ui.sl_ipf_thresh, 'Value');
    
    updated_opts.mask_spec = (get(ui.chk_mask_spec, 'Value') == 1);
    updated_opts.crop_spec = (get(ui.chk_crop_spec, 'Value') == 1);
    updated_opts.show_micronbar = (get(ui.chk_micronbar, 'Value') == 1);
end

function onBatchNavigate(fig, action)
    app = fig.UserData;
    app.batch_status = action;
    app.opts = saveCurrentSettingsToOpts(fig);
    fig.UserData = app;
    uiresume(fig);
end

function onImportOptions(fig)
    app = fig.UserData;
    [filename, pathname] = uigetfile({'*_Params.mat;*_Report.md', 'Export Archives & Reports (*_Params.mat, *_Report.md)'; ...
                                      '*.mat', 'MAT Parameter Archives (*.mat)'; ...
                                      '*.md', 'Markdown Reports (*.md)'; ...
                                      '*.*', 'All Files (*.*)'}, ...
                                     'Select Previous Visualization Parameters to Import', app.exportDir);
    if isequal(filename, 0) || isequal(pathname, 0)
        return;
    end
    filepath = fullfile(pathname, filename);
    p = struct();
    
    try
        if endsWith(lower(filename), '.mat')
            data = load(filepath);
            if isfield(data, 'params')
                p = data.params;
            elseif isfield(data, 'opts')
                p = data.opts;
            else
                p = data;
            end
        elseif endsWith(lower(filename), '.md')
            txt = fileread(filepath);
            tokens = regexp(txt, '\| \*\*Background Property\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.background_property = tokens{1}{1}; end
            tokens = regexp(txt, '\| \*\*Color Inverted\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.is_background_inverted = strcmpi(tokens{1}{1}, 'true'); end
            tokens = regexp(txt, '\| \*\*Black Point\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.black_point = str2double(tokens{1}{1}); end
            tokens = regexp(txt, '\| \*\*White Point\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.white_point = str2double(tokens{1}{1}); end
            tokens = regexp(txt, '\| \*\*Gamma \(\\\gamma\)\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.gamma_val = str2double(tokens{1}{1}); end
            tokens = regexp(txt, '\| \*\*IPF Color Key\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.ipf_color_key = tokens{1}{1}; end
            tokens = regexp(txt, '\| \*\*IPF Reference Direction\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.ipf_direction = tokens{1}{1}; end
            tokens = regexp(txt, '\| \*\*IPF Alpha Overlay\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.alpha_overlay = str2double(tokens{1}{1}); end
            tokens = regexp(txt, '\| \*\*IPF Restrict > BlackPt\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.restrict_above_bp = strcmpi(tokens{1}{1}, 'true'); end
            tokens = regexp(txt, '\| \*\*IPF Restrict < WhitePt\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.restrict_below_wp = strcmpi(tokens{1}{1}, 'true'); end
            tokens = regexp(txt, '\| \*\*IPF Property Filter\*\* \| `([^ ]+) >= ([^`]+)`', 'tokens');
            if ~isempty(tokens)
                p.ipf_thresh_property = tokens{1}{1};
                p.ipf_thresh_value = str2double(tokens{1}{2});
            end
            tokens = regexp(txt, '\| \*\*Specimen Mask Applied\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.mask_to_specimen = strcmpi(tokens{1}{1}, 'true'); end
            tokens = regexp(txt, '\| \*\*Map Cropped to Specimen\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.crop_to_specimen = strcmpi(tokens{1}{1}, 'true'); end
            tokens = regexp(txt, '\| \*\*Grain Boundaries\*\* \| `([^`]+)`', 'tokens');
            if ~isempty(tokens), p.grain_boundary_style = tokens{1}{1}; end
        end
    catch ME
        fprintf(2, '[Import Error] Failed to read parameter file %s: %s\n', filename, ME.message);
        errordlg(sprintf('Could not parse file:\n%s', ME.message), 'Import Error');
        return;
    end
    
    % --- Security Validation & Application against current EBSD dataset ---
    ui = app.ui;
    
    % 1. Background Property (with security check)
    bg_fld = '';
    if isfield(p, 'background_property'), bg_fld = p.background_property;
    elseif isfield(p, 'defaultBgField'), bg_fld = p.defaultBgField; end
    if ~isempty(bg_fld)
        bg_list = get(ui.pu_bg, 'String');
        idx = find(strcmpi(bg_list, bg_fld), 1);
        if ~isempty(idx)
            set(ui.pu_bg, 'Value', idx);
            onBgChange(fig, false);
        else
            fprintf('[Import Warning] Property "%s" not found in current dataset properties; keeping current background.\n', bg_fld);
        end
    end
    
    % 2. Invert Toggle
    if isfield(p, 'is_background_inverted'), set(ui.chk_inv, 'Value', p.is_background_inverted);
    elseif isfield(p, 'invert_bg'), set(ui.chk_inv, 'Value', p.invert_bg); end
    
    % 3. Sliders: BP, WP, Gamma
    if isfield(p, 'black_point'), try, set(ui.sl_bp, 'Value', min(0.95, max(0.0, p.black_point))); catch, end;
    elseif isfield(p, 'bp_val'), try, set(ui.sl_bp, 'Value', min(0.95, max(0.0, p.bp_val))); catch, end; end
    
    if isfield(p, 'white_point'), try, set(ui.sl_wp, 'Value', min(1.0, max(0.05, p.white_point))); catch, end;
    elseif isfield(p, 'wp_val'), try, set(ui.sl_wp, 'Value', min(1.0, max(0.05, p.wp_val))); catch, end; end
    
    if isfield(p, 'gamma_val'), try, set(ui.sl_gam, 'Value', min(3.5, max(0.1, p.gamma_val))); catch, end; end
    
    % 4. Overlay Mode & Settings
    ov_mode_str = '';
    if isfield(p, 'overlay_mode'), ov_mode_str = p.overlay_mode;
    elseif isfield(p, 'defaultOverlayMode'), ov_mode_str = p.defaultOverlayMode; end
    
    if ~isempty(ov_mode_str) && (contains(lower(ov_mode_str), 'scalar') || contains(lower(ov_mode_str), 'strain') || contains(lower(ov_mode_str), 'property') || (isnumeric(ov_mode_str) && ov_mode_str == 2))
        set(ui.pu_ov_mode, 'Value', 2);
        set(ui.pu_ipf, 'Visible', 'off'); set(ui.pu_dir, 'Visible', 'off');
        set(ui.pu_prop_ov, 'Visible', 'on'); set(ui.pu_cmap, 'Visible', 'on');
        set(ui.lbl_cmin, 'Visible', 'on'); set(ui.ed_cmin, 'Visible', 'on');
        set(ui.lbl_cmax, 'Visible', 'on'); set(ui.ed_cmax, 'Visible', 'on');
        
        prop_ov = '';
        if isfield(p, 'property_overlay'), prop_ov = p.property_overlay;
        elseif isfield(p, 'defaultPropOverlay'), prop_ov = p.defaultPropOverlay; end
        if ~isempty(prop_ov)
            plist = get(ui.pu_prop_ov, 'String'); idx_po = find(strcmpi(plist, prop_ov), 1);
            if ~isempty(idx_po), set(ui.pu_prop_ov, 'Value', idx_po); end
        end
        
        cmap_str = '';
        if isfield(p, 'colormap'), cmap_str = p.colormap;
        elseif isfield(p, 'defaultColormap'), cmap_str = p.defaultColormap; end
        if ~isempty(cmap_str)
            clist = get(ui.pu_cmap, 'String'); idx_cm = find(strcmpi(clist, cmap_str), 1);
            if ~isempty(idx_cm), set(ui.pu_cmap, 'Value', idx_cm); end
        end
        
        if isfield(p, 'cmin_val') && ~isnan(p.cmin_val), try, set(ui.ed_cmin, 'String', num2str(p.cmin_val)); catch, end;
        elseif isfield(p, 'cmin') && ~isnan(p.cmin), try, set(ui.ed_cmin, 'String', num2str(p.cmin)); catch, end; end
        if isfield(p, 'cmax_val') && ~isnan(p.cmax_val), try, set(ui.ed_cmax, 'String', num2str(p.cmax_val)); catch, end;
        elseif isfield(p, 'cmax') && ~isnan(p.cmax), try, set(ui.ed_cmax, 'String', num2str(p.cmax)); catch, end; end
    else
        set(ui.pu_ov_mode, 'Value', 1);
        set(ui.pu_ipf, 'Visible', 'on'); set(ui.pu_dir, 'Visible', 'on');
        set(ui.pu_prop_ov, 'Visible', 'off'); set(ui.pu_cmap, 'Visible', 'off');
        set(ui.lbl_cmin, 'Visible', 'off'); set(ui.ed_cmin, 'Visible', 'off');
        set(ui.lbl_cmax, 'Visible', 'off'); set(ui.ed_cmax, 'Visible', 'off');
        
        ipf_key = '';
        if isfield(p, 'ipf_color_key'), ipf_key = p.ipf_color_key;
        elseif isfield(p, 'defaultIpfKey'), ipf_key = p.defaultIpfKey; end
        if ~isempty(ipf_key)
            key_list = get(ui.pu_ipf, 'String');
            idx = find(strcmpi(key_list, ipf_key), 1);
            if ~isempty(idx), set(ui.pu_ipf, 'Value', idx); end
        end
        
        ipf_dir = '';
        if isfield(p, 'ipf_direction'), ipf_dir = p.ipf_direction;
        elseif isfield(p, 'defaultIpfDir'), ipf_dir = p.defaultIpfDir; end
        if ~isempty(ipf_dir)
            dir_list = get(ui.pu_dir, 'String');
            idx = find(strcmpi(dir_list, ipf_dir), 1);
            if isempty(idx)
                if strcmpi(ipf_dir, 'Z'), idx = 2; elseif strcmpi(ipf_dir, 'X'), idx = 3; else, idx = 1; end
            end
            if ~isempty(idx), set(ui.pu_dir, 'Value', idx); end
        end
    end
    onIpfChange(fig, false);
    
    % 5. Alpha & Restrictions
    if isfield(p, 'alpha_overlay'), try, set(ui.sl_alpha, 'Value', min(1.0, max(0.0, p.alpha_overlay))); catch, end;
    elseif isfield(p, 'alpha_val'), try, set(ui.sl_alpha, 'Value', min(1.0, max(0.0, p.alpha_val))); catch, end; end
    
    if isfield(p, 'restrict_above_bp'), set(ui.chk_above_bp, 'Value', p.restrict_above_bp); end
    if isfield(p, 'restrict_below_wp'), set(ui.chk_below_wp, 'Value', p.restrict_below_wp); end
    
    % 6. IPF Threshold Property & Value (with security check)
    th_fld = '';
    if isfield(p, 'ipf_thresh_property'), th_fld = p.ipf_thresh_property;
    elseif isfield(p, 'defaultIpfThreshProp'), th_fld = p.defaultIpfThreshProp; end
    if ~isempty(th_fld)
        th_list = get(ui.pu_ipf_thresh, 'String');
        idx = find(strcmpi(th_list, th_fld), 1);
        if ~isempty(idx)
            set(ui.pu_ipf_thresh, 'Value', idx);
            onIpfThreshPropChange(fig, false);
        else
            fprintf('[Import Warning] IPF Threshold property "%s" not found in dataset; using existing property.\n', th_fld);
        end
    end
    if isfield(p, 'ipf_thresh_value')
        try
            min_v = get(ui.sl_ipf_thresh, 'Min'); max_v = get(ui.sl_ipf_thresh, 'Max');
            val = min(max_v, max(min_v, p.ipf_thresh_value));
            set(ui.sl_ipf_thresh, 'Value', val);
        catch
        end
    elseif isfield(p, 'ipf_thresh_val')
        try
            min_v = get(ui.sl_ipf_thresh, 'Min'); max_v = get(ui.sl_ipf_thresh, 'Max');
            val = min(max_v, max(min_v, p.ipf_thresh_val));
            set(ui.sl_ipf_thresh, 'Value', val);
        catch
        end
    end
    
    % 7. Specimen Masking & Cropping Toggles
    if strcmpi(get(ui.chk_mask_spec, 'Enable'), 'on')
        if isfield(p, 'mask_to_specimen'), set(ui.chk_mask_spec, 'Value', p.mask_to_specimen);
        elseif isfield(p, 'mask_spec'), set(ui.chk_mask_spec, 'Value', p.mask_spec); end
        
        if isfield(p, 'crop_to_specimen'), set(ui.chk_crop_spec, 'Value', p.crop_to_specimen);
        elseif isfield(p, 'crop_spec'), set(ui.chk_crop_spec, 'Value', p.crop_spec); end
    end
    
    if isfield(p, 'show_micronbar'), set(ui.chk_micronbar, 'Value', p.show_micronbar); end
    
    % 8. Grain Boundary Style
    gb_style = '';
    if isfield(p, 'grain_boundary_style'), gb_style = p.grain_boundary_style;
    elseif isfield(p, 'defaultGbStyle'), gb_style = p.defaultGbStyle; end
    if ~isempty(gb_style)
        gb_list = get(ui.pu_gb, 'String');
        idx = find(strcmpi(gb_list, gb_style), 1);
        if ~isempty(idx)
            set(ui.pu_gb, 'Value', idx);
        end
    end
    
    % 9. Resolution settings
    if isfield(p, 'target_shorter_edge_px'), try, set(ui.ed_px, 'String', num2str(p.target_shorter_edge_px)); catch, end;
    elseif isfield(p, 'target_px'), try, set(ui.ed_px, 'String', num2str(p.target_px)); catch, end; end
    if isfield(p, 'dpi'), try, set(ui.ed_dpi, 'String', num2str(p.dpi)); catch, end;
    elseif isfield(p, 'dpi_val'), try, set(ui.ed_dpi, 'String', num2str(p.dpi_val)); catch, end; end
    
    updatePlot(fig);
    onGbChange(fig);
    fprintf('[Interactive_PrettyPlotter] Successfully imported and applied parameters from: %s\n', filename);
    try, msgbox(sprintf('Successfully applied settings from:\n%s', filename), 'Import Complete'); catch, end
end

function onOverlayModeChange(fig)
    app = fig.UserData;
    ui  = app.ui;
    mode_val = get(ui.pu_ov_mode, 'Value');
    if mode_val == 1
        set(ui.pu_ipf, 'Visible', 'on');  set(ui.pu_dir, 'Visible', 'on');
        set(ui.pu_prop_ov, 'Visible', 'off'); set(ui.pu_cmap, 'Visible', 'off');
        set(ui.lbl_cmin, 'Visible', 'off'); set(ui.ed_cmin, 'Visible', 'off');
        set(ui.lbl_cmax, 'Visible', 'off'); set(ui.ed_cmax, 'Visible', 'off');
        set(ui.btn_prop_auto, 'Visible', 'off');
    else
        set(ui.pu_ipf, 'Visible', 'off');  set(ui.pu_dir, 'Visible', 'off');
        set(ui.pu_prop_ov, 'Visible', 'on'); set(ui.pu_cmap, 'Visible', 'on');
        set(ui.lbl_cmin, 'Visible', 'on'); set(ui.ed_cmin, 'Visible', 'on');
        set(ui.lbl_cmax, 'Visible', 'on'); set(ui.ed_cmax, 'Visible', 'on');
        set(ui.btn_prop_auto, 'Visible', 'on');
    end
    onIpfChange(fig, true);
end

function onPropAuto(fig)
    app = fig.UserData;
    idx_p = get(app.ui.pu_prop_ov, 'Value'); props = get(app.ui.pu_prop_ov, 'String'); sel_p = props{idx_p};
    
    if isfield(app.ebsd_grid.prop, sel_p)
        p_data = app.ebsd_grid.prop.(sel_p);
    elseif isprop(app.ebsd_grid, sel_p)
        p_data = app.ebsd_grid.(sel_p);
    else
        p_data = double(app.ebsd_grid.ci);
    end
    p_data = double(p_data);
    
    valid_vals = p_data(~isnan(p_data) & ~isinf(p_data));
    if all(p_data(~isnan(p_data)) >= 0) && any(p_data(:) > 0)
        valid_vals = p_data(~isnan(p_data) & ~isinf(p_data) & p_data > 0);
    end
    if isempty(valid_vals), return; end
    
    % Crop colormap limits to the 2nd and 98th percentiles
    p_low  = prctile(valid_vals, 2);
    p_high = prctile(valid_vals, 98);
    if p_low >= p_high
        p_high = p_low + 1e-4;
    end
    
    set(app.ui.ed_cmin, 'String', num2str(p_low, '%.4f'));
    set(app.ui.ed_cmax, 'String', num2str(p_high, '%.4f'));
    
    onIpfChange(fig, true);
end

function onIpfChange(fig, doUpdate)
    if nargin < 2, doUpdate = true; end
    app = fig.UserData;
    mode_val = get(app.ui.pu_ov_mode, 'Value');
    
    if mode_val == 1
        % --- 1. IPF ORIENTATION OVERLAY MODE ---
        idx_k = get(app.ui.pu_ipf, 'Value'); keys = get(app.ui.pu_ipf, 'String'); sel_k = keys{idx_k};
        idx_d = get(app.ui.pu_dir, 'Value');
        if idx_d == 1, dir_vec = yvector; elseif idx_d == 2, dir_vec = zvector; else, dir_vec = xvector; end
        
        try
            oM = feval(sel_k, app.ebsd_grid);
            try, oM.inversePoleFigureDirection = dir_vec; catch, end
            ipf_rgb = oM.orientation2color(app.ebsd_grid.orientations);
        catch ME
            fprintf(2, '[Warning] Could not construct %s (%s). Falling back to ipfColorKey.\n', sel_k, ME.message);
            oM = ipfColorKey(app.ebsd_grid);
            oM.inversePoleFigureDirection = dir_vec;
            ipf_rgb = oM.orientation2color(app.ebsd_grid.orientations);
        end
        
        if size(ipf_rgb, 3) ~= 3 && size(ipf_rgb, 2) == 3
            ipf_rgb = reshape(ipf_rgb, [size(app.ebsd_grid, 1), size(app.ebsd_grid, 2), 3]);
        end
        app.ipf_rgb = double(ipf_rgb);
        try, colorbar(app.ax, 'off'); catch, end
    else
        % --- 2. SCALAR PROPERTY / STRAIN COLORMAP MODE ---
        idx_p = get(app.ui.pu_prop_ov, 'Value'); props = get(app.ui.pu_prop_ov, 'String'); sel_p = props{idx_p};
        idx_c = get(app.ui.pu_cmap, 'Value'); cmaps = get(app.ui.pu_cmap, 'String'); sel_c = cmaps{idx_c};
        
        if isfield(app.ebsd_grid.prop, sel_p)
            p_data = app.ebsd_grid.prop.(sel_p);
        elseif isprop(app.ebsd_grid, sel_p)
            p_data = app.ebsd_grid.(sel_p);
        else
            p_data = double(app.ebsd_grid.ci);
        end
        p_data = double(p_data);
        
        cmin = str2double(get(app.ui.ed_cmin, 'String'));
        cmax = str2double(get(app.ui.ed_cmax, 'String'));
        if isnan(cmin) || isnan(cmax) || cmin >= cmax
            valid_vals = p_data(~isnan(p_data) & ~isinf(p_data) & p_data ~= 0);
            if ~isempty(valid_vals)
                cmin = min(valid_vals(:));
                cmax = max(valid_vals(:));
                if cmin == cmax, cmax = cmin + 0.1; end
            else
                cmin = 0.0; cmax = 1.0;
            end
            set(app.ui.ed_cmin, 'String', num2str(cmin, '%.4f'));
            set(app.ui.ed_cmax, 'String', num2str(cmax, '%.4f'));
        end
        
        % Build colormap matrix
        try
            if strcmpi(sel_c, 'viridis') || strcmpi(sel_c, 'magma') || strcmpi(sel_c, 'inferno') || strcmpi(sel_c, 'plasma')
                fig_tmp = figure('Visible', 'off');
                try, mtexColorMap(fig_tmp, sel_c); cmap_mat = colormap(fig_tmp);
                catch, colormap(fig_tmp, sel_c); cmap_mat = colormap(fig_tmp); end
                close(fig_tmp);
            else
                cmap_mat = feval(sel_c, 256);
            end
        catch
            cmap_mat = parula(256);
        end
        
        % Normalize property values into colormap indices
        norm_p = (p_data - cmin) ./ (cmax - cmin);
        norm_p(norm_p < 0) = 0; norm_p(norm_p > 1) = 1;
        num_c = size(cmap_mat, 1);
        col_idx = round(norm_p .* (num_c - 1)) + 1;
        col_idx(isnan(col_idx) | isinf(col_idx)) = 1;
        
        rgb_flat = cmap_mat(col_idx(:), :);
        app.ipf_rgb = reshape(rgb_flat, [size(app.ebsd_grid, 1), size(app.ebsd_grid, 2), 3]);
        
        % Enable live publication-grade colorbar on GUI axes without altering raster RGB colors
        try
            colormap(app.ax, cmap_mat);
            try, clim(app.ax, [cmin, cmax]); catch, caxis(app.ax, [cmin, cmax]); end
            cb = colorbar(app.ax, 'Location', 'eastoutside');
            ylabel(cb, strrep(sel_p, '_', '\_'), 'FontWeight', 'bold', 'FontSize', 11);
        catch
        end
    end
    
    fig.UserData = app;
    if doUpdate, updatePlot(fig); end
end

function onIpfThreshPropChange(fig, doUpdate)
    if nargin < 2, doUpdate = true; end
    app = fig.UserData;
    
    idx = get(app.ui.pu_ipf_thresh, 'Value'); flds = get(app.ui.pu_ipf_thresh, 'String'); sel_fld = flds{idx};
    
    try
        val = app.ebsd_grid.prop.(sel_fld);
    catch
        try, val = app.ebsd_grid.(sel_fld); catch, try, val = app.ebsd_grid.ci; catch, val = zeros(size(app.ebsd_grid)); end; end
    end
    if isa(val, 'vector3d'), val = val.norm; elseif isstruct(val) && isfield(val, 'x') && isfield(val, 'y'), val = sqrt(val.x.^2 + val.y.^2); end
    
    val = double(val(:));
    if ismatrix(val) && ~isequal(size(val), [size(app.ebsd_grid,1), size(app.ebsd_grid,2)])
        try, val = reshape(val, [size(app.ebsd_grid,1), size(app.ebsd_grid,2)]); catch, end
    end
    
    app.ipf_thresh_data = val;
    valid_vals = val(~isnan(val) & ~isinf(val));
    if isempty(valid_vals)
        app.ipf_thresh_min = 0; app.ipf_thresh_max = 1;
    else
        app.ipf_thresh_min = min(valid_vals(:));
        app.ipf_thresh_max = max(valid_vals(:));
        if app.ipf_thresh_min == app.ipf_thresh_max, app.ipf_thresh_max = app.ipf_thresh_min + 1e-5; end
    end
    
    set(app.ui.sl_ipf_thresh, 'Min', app.ipf_thresh_min, 'Max', app.ipf_thresh_max);
    curr_val = get(app.ui.sl_ipf_thresh, 'Value');
    if curr_val < app.ipf_thresh_min || curr_val > app.ipf_thresh_max
        if strcmpi(sel_fld, 'ci') && app.ipf_thresh_min <= 0.10 && app.ipf_thresh_max >= 0.10
            set(app.ui.sl_ipf_thresh, 'Value', 0.10);
        else
            set(app.ui.sl_ipf_thresh, 'Value', app.ipf_thresh_min);
        end
    end
    
    fig.UserData = app;
    if doUpdate, updatePlot(fig); end
end

function onBgChange(fig, doUpdate)
    if nargin < 2, doUpdate = true; end
    app = fig.UserData;
    
    idx = get(app.ui.pu_bg, 'Value'); flds = get(app.ui.pu_bg, 'String'); sel_fld = flds{idx};
    
    try
        val = app.ebsd_grid.prop.(sel_fld);
    catch
        try
            val = app.ebsd_grid.(sel_fld);
        catch
            try, val = app.ebsd_grid.ci; catch, val = ones(size(app.ebsd_grid)); end
        end
    end
    
    if isa(val, 'vector3d'), val = val.norm;
    elseif isstruct(val) && isfield(val, 'x') && isfield(val, 'y'), val = sqrt(val.x.^2 + val.y.^2); end
    
    val = double(val);
    if ismatrix(val) && ~isequal(size(val), [size(app.ebsd_grid,1), size(app.ebsd_grid,2)])
        try, val = reshape(val, [size(app.ebsd_grid,1), size(app.ebsd_grid,2)]); catch, end
    end
    
    app.bg_raw = val;
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

function onAutoContrast(fig)
    app = fig.UserData;
    if ~isfield(app, 'bg_raw') || isempty(app.bg_raw), return; end
    
    valid_vals = app.bg_raw(~isnan(app.bg_raw) & ~isinf(app.bg_raw));
    if isempty(valid_vals), return; end
    
    % Crop contrast to robust 2nd and 98th percentiles
    p_low  = prctile(valid_vals, 2);
    p_high = prctile(valid_vals, 98);
    
    rng = app.bg_max - app.bg_min;
    if rng <= 0 || isnan(rng), rng = 1e-5; end
    
    new_bp_norm = (p_low - app.bg_min) / rng;
    new_wp_norm = (p_high - app.bg_min) / rng;
    
    new_bp_norm = min(0.95, max(0.0, new_bp_norm));
    new_wp_norm = min(1.0, max(0.05, new_wp_norm));
    
    if new_wp_norm <= new_bp_norm
        new_bp_norm = 0.02;
        new_wp_norm = 0.98;
    end
    
    set(app.ui.sl_bp, 'Value', new_bp_norm);
    set(app.ui.sl_wp, 'Value', new_wp_norm);
    
    updatePlot(fig);
end

function onGbChange(fig)
    app = fig.UserData;
    if isempty(app.grains), return; end
    
    idx = get(app.ui.pu_gb, 'Value'); styles = get(app.ui.pu_gb, 'String'); sel_style = styles{idx};
    
    % Robustly purge ALL grain boundary line, patch, and graphic overlay objects from axes
    all_children = get(app.ax, 'Children');
    if ~isempty(all_children)
        types = get(all_children, 'Type');
        if ischar(types), types = {types}; end
        non_image_idx = ~strcmpi(types, 'image') & ~strcmpi(types, 'text');
        delete(all_children(non_image_idx));
    end
    app.hBoundaries = [];
    
    if contains(sel_style, 'Hidden') || contains(sel_style, 'None')
        fig.UserData = app; return;
    end
    
    if contains(sel_style, 'White'), col = 'white'; else, col = 'black'; end
    if contains(sel_style, '1.2'), lw = 1.2; else, lw = 2.0; end
    
    hold(app.ax, 'on');
    app.hBoundaries = plot(app.grains.boundary, 'Parent', app.ax, 'lineColor', col, 'lineWidth', lw, 'micronbar', 'off');
    hold(app.ax, 'off');
    
    % Ensure raster map image remains strictly beneath newly plotted boundaries
    if ~isempty(app.hImage) && isvalid(app.hImage)
        uistack(app.hImage, 'bottom');
    end
    fig.UserData = app;
end

function onSliderChange(fig)
    updatePlot(fig);
end

function updatePlot(fig)
    app = fig.UserData;
    if ~isfield(app, 'ipf_rgb'), return; end
    
    alpha_val = get(app.ui.sl_alpha, 'Value');
    bp_norm   = get(app.ui.sl_bp,    'Value');
    wp_norm   = get(app.ui.sl_wp,    'Value');
    gamma_val = get(app.ui.sl_gam,   'Value');
    t_val     = get(app.ui.sl_ipf_thresh, 'Value');
    
    if wp_norm <= bp_norm
        wp_norm = min(1.0, bp_norm + 0.02);
        set(app.ui.sl_wp, 'Value', wp_norm);
    end
    
    set(app.ui.lbl_alpha, 'String', sprintf('IPF Overlay Alpha: %.2f', alpha_val));
    set(app.ui.lbl_bp,    'String', sprintf('Black Point: %.2f', bp_norm));
    set(app.ui.lbl_wp,    'String', sprintf('White Point: %.2f', wp_norm));
    set(app.ui.lbl_gam,   'String', sprintf('Gamma (\\gamma): %.2f', gamma_val));
    
    idx_th = get(app.ui.pu_ipf_thresh, 'Value'); flds = get(app.ui.pu_ipf_thresh, 'String'); sel_th_fld = flds{idx_th};
    set(app.ui.lbl_ipf_thresh, 'String', sprintf('Min IPF Thresh (%s): %.3f', sel_th_fld, t_val));
    
    rng = app.bg_max - app.bg_min;
    phys_bp = app.bg_min + bp_norm * rng;
    phys_wp = app.bg_min + wp_norm * rng;
    if phys_wp <= phys_bp, phys_wp = phys_bp + 1e-5; end
    
    norm_bg = (app.bg_raw - phys_bp) / (phys_wp - phys_bp);
    norm_bg = min(max(norm_bg, 0), 1);
    gamma_bg = norm_bg .^ (1.0 / max(0.01, gamma_val));
    if get(app.ui.chk_inv, 'Value') == 1
        gamma_bg = 1.0 - gamma_bg;
    end
    
    bg_rgb = repmat(gamma_bg, [1, 1, 3]);
    nan_mask = isnan(app.bg_raw) | isinf(app.bg_raw);
    nan_mask_3d = repmat(nan_mask, [1, 1, 3]);
    bg_rgb(nan_mask_3d) = 0.85; 
    
    % --- Assemble Composite IPF Visibility Mask ---
    ipf_mask = true(size(app.bg_raw));
    
    % 1. Black/White point band-pass restrictions
    if get(app.ui.chk_above_bp, 'Value') == 1
        ipf_mask = ipf_mask & (app.bg_raw >= phys_bp);
    end
    if get(app.ui.chk_below_wp, 'Value') == 1
        ipf_mask = ipf_mask & (app.bg_raw <= phys_wp);
    end
    
    % 2. Independent property thresholding (e.g., CI < 0.10 -> notIndexed without affecting background)
    if ~isempty(app.ipf_thresh_data) && t_val > app.ipf_thresh_min
        ipf_mask = ipf_mask & (app.ipf_thresh_data >= t_val);
    end
    
    % 3. Specimen Location Masking
    if get(app.ui.chk_mask_spec, 'Value') == 1 && ~isempty(app.spec_loc)
        ipf_mask = ipf_mask & app.spec_loc;
    end
    
    % Exclude NaN IPF colors from the mask
    ipf_nan_mask = isnan(app.ipf_rgb(:,:,1));
    ipf_mask = ipf_mask & ~ipf_nan_mask;
    
    alpha_map = repmat(alpha_val * double(ipf_mask), [1, 1, 3]);
    
    % Sanitize IPF RGB to prevent 0.0 * NaN = NaN poisoning in the background
    safe_ipf_rgb = app.ipf_rgb;
    safe_ipf_rgb(isnan(safe_ipf_rgb)) = 0;
    
    final_rgb = (1.0 - alpha_map) .* bg_rgb + alpha_map .* safe_ipf_rgb;
    final_rgb(nan_mask_3d) = 0.85; 
    
    if isempty(app.hImage) || ~isvalid(app.hImage)
        app.hImage = image(app.ax, app.x_range, app.y_range, final_rgb);
        set(app.ax, 'YDir', app.ydir, 'XDir', app.xdir); axis(app.ax, 'image');
        uistack(app.hImage, 'bottom');
    else
        app.hImage.CData = final_rgb;
        set(app.ax, 'YDir', app.ydir, 'XDir', app.xdir);
    end
    
    % 4. Dynamic Spatial Bounding Box Cropping to Specimen
    if get(app.ui.chk_crop_spec, 'Value') == 1 && ~isempty(app.spec_loc)
        valid_px = app.spec_loc;
        if any(valid_px(:))
            try
                xs = app.ebsd_grid.x(valid_px); ys = app.ebsd_grid.y(valid_px);
                xs = xs(~isnan(xs) & ~isinf(xs)); ys = ys(~isnan(ys) & ~isinf(ys));
                min_xs = double(min(xs(:))); max_xs = double(max(xs(:)));
                min_ys = double(min(ys(:))); max_ys = double(max(ys(:)));
                if max_xs <= min_xs, max_xs = min_xs + 1.0; end
                if max_ys <= min_ys, max_ys = min_ys + 1.0; end
                pad_x = 0.02 * max(1.0, max_xs - min_xs);
                pad_y = 0.02 * max(1.0, max_ys - min_ys);
                set(app.ax, 'XLim', [min_xs - pad_x, max_xs + pad_x], 'YLim', [min_ys - pad_y, max_ys + pad_y]);
            catch
                set(app.ax, 'XLim', [min(app.x_range(:)), max(app.x_range(:))], 'YLim', [min(app.y_range(:)), max(app.y_range(:))]);
            end
        else
            set(app.ax, 'XLim', [min(app.x_range(:)), max(app.x_range(:))], 'YLim', [min(app.y_range(:)), max(app.y_range(:))]);
        end
    else
        set(app.ax, 'XLim', [min(app.x_range(:)), max(app.x_range(:))], 'YLim', [min(app.y_range(:)), max(app.y_range(:))]);
    end
    
    % 5. Render Custom Dynamic Micron Bar
    hBar = findobj(app.ax, 'Tag', 'CustomScaleBar');
    delete(hBar);
    
    if get(app.ui.chk_micronbar, 'Value') == 1
        x_lim = get(app.ax, 'XLim');
        y_lim = get(app.ax, 'YLim');
        width_um = diff(x_lim);
        
        % Target about 20% of the view width
        order = 10^floor(log10(width_um * 0.2));
        candidates = [1, 2, 5, 10] * order;
        [~, idx_c] = min(abs(candidates - width_um * 0.2));
        bar_width = candidates(idx_c);
        
        x_end = x_lim(2) - 0.05 * width_um;
        x_start = x_end - bar_width;
        
        % Y position depends on YDir
        if strcmpi(get(app.ax, 'YDir'), 'reverse')
            y_pos = y_lim(2) - 0.05 * diff(y_lim); % bottom of image
            y_text = y_pos - 0.01 * diff(y_lim); % slightly above line (smaller y)
            valign = 'bottom';
        else
            y_pos = y_lim(1) + 0.05 * diff(y_lim); % bottom of image
            y_text = y_pos + 0.01 * diff(y_lim); % slightly above line (larger y)
            valign = 'bottom';
        end
        
        % Draw thick black border line, then white inner line
        line(app.ax, [x_start, x_end], [y_pos, y_pos], 'Color', 'k', 'LineWidth', 6, 'Tag', 'CustomScaleBar');
        line(app.ax, [x_start, x_end], [y_pos, y_pos], 'Color', 'w', 'LineWidth', 3, 'Tag', 'CustomScaleBar');
        
        % Draw text with black background/edge so it's always readable
        txt = sprintf('%g \\mum', bar_width);
        text(app.ax, mean([x_start, x_end]), y_text, txt, ...
            'Color', 'k', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', valign, 'Tag', 'CustomScaleBar', ...
            'BackgroundColor', 'w', 'Margin', 1);
    end
    
    idx_bg = get(app.ui.pu_bg, 'Value'); bg_names = get(app.ui.pu_bg, 'String');
    sel_bg_title = bg_names{idx_bg};
    if get(app.ui.chk_inv, 'Value') == 1, sel_bg_title = [sel_bg_title, ' (Inv)']; end
    if get(app.ui.pu_ov_mode, 'Value') == 1
        idx_ipf = get(app.ui.pu_ipf, 'Value'); ipf_names = get(app.ui.pu_ipf, 'String');
        idx_dir = get(app.ui.pu_dir, 'Value'); dir_names = get(app.ui.pu_dir, 'String');
        sel_ov_title = sprintf('%s (%s)', ipf_names{idx_ipf}, dir_names{idx_dir});
    else
        idx_p = get(app.ui.pu_prop_ov, 'Value'); prop_names = get(app.ui.pu_prop_ov, 'String');
        idx_cm = get(app.ui.pu_cmap, 'Value'); cmap_names = get(app.ui.pu_cmap, 'String');
        sel_ov_title = sprintf('%s (%s)', prop_names{idx_p}, cmap_names{idx_cm});
    end
    title(app.ax, sprintf('%s: %s + %s (\\alpha=%.2f, \\gamma=%.2f)', ...
          strrep(app.datasetName, '_', '\_'), strrep(sel_bg_title, '_', '\_'), strrep(sel_ov_title, '_', '\_'), alpha_val, gamma_val), ...
          'FontSize', 13, 'FontWeight', 'bold');
          
    app.curr_rgb = final_rgb;
    fig.UserData = app;
end

%% ========================================================================
%  EXPORT SUITE (PUBLICATION FIGURE, MAT DICTIONARY & MARKDOWN REPORT)
%  ========================================================================

function onSaveExport(fig, isNext)
    if nargin < 2, isNext = false; end
    app = fig.UserData;
    ui  = app.ui;
    
    base_name = strtrim(get(ui.ed_name, 'String'));
    if isempty(base_name), base_name = app.datasetName; end
    
    target_px = str2double(get(ui.ed_px, 'String'));
    if isnan(target_px) || target_px < 100, target_px = 1500; end
    dpi_val = str2double(get(ui.ed_dpi, 'String'));
    if isnan(dpi_val) || dpi_val <= 0, dpi_val = 150; end
    
    idx_fmt = get(ui.pu_fmt, 'Value'); fmts = get(ui.pu_fmt, 'String'); ext_fmt = fmts{idx_fmt};
    idx_bg  = get(ui.pu_bg,  'Value'); bg_list = get(ui.pu_bg, 'String'); sel_bg = bg_list{idx_bg};
    is_inverted = (get(ui.chk_inv, 'Value') == 1);
    idx_ipf = get(ui.pu_ipf, 'Value'); ipf_list = get(ui.pu_ipf, 'String'); sel_ipf = ipf_list{idx_ipf};
    idx_dir = get(ui.pu_dir, 'Value'); dir_list = get(ui.pu_dir, 'String'); sel_dir = dir_list{idx_dir};
    idx_gb  = get(ui.pu_gb,  'Value'); gb_list = get(ui.pu_gb, 'String'); sel_gb = gb_list{idx_gb};
    idx_th  = get(ui.pu_ipf_thresh, 'Value'); th_list = get(ui.pu_ipf_thresh, 'String'); sel_th_prop = th_list{idx_th};
    
    mode_val = get(ui.pu_ov_mode, 'Value'); ov_mode_str = 'IPF Orientation Map';
    idx_p = get(ui.pu_prop_ov, 'Value'); prop_list = get(ui.pu_prop_ov, 'String'); sel_prop = prop_list{idx_p};
    idx_c = get(ui.pu_cmap, 'Value'); cmap_list = get(ui.pu_cmap, 'String'); sel_cmap = cmap_list{idx_c};
    cmin_val = str2double(get(ui.ed_cmin, 'String')); cmax_val = str2double(get(ui.ed_cmax, 'String'));
    
    alpha_val = get(ui.sl_alpha, 'Value');
    thresh_val = get(ui.sl_ipf_thresh, 'Value');
    restrict_above_bp = (get(ui.chk_above_bp, 'Value') == 1);
    restrict_below_wp = (get(ui.chk_below_wp, 'Value') == 1);
    mask_to_specimen  = (get(ui.chk_mask_spec, 'Value') == 1);
    crop_to_specimen  = (get(ui.chk_crop_spec, 'Value') == 1);
    
    bp_val    = get(ui.sl_bp,    'Value');
    wp_val    = get(ui.sl_wp,    'Value');
    gamma_val = get(ui.sl_gam,   'Value');
    
    timestamp_str = datestr(now, 'yyyy-mm-dd_HHMMSS');
    inv_tag = ''; if is_inverted, inv_tag = '_Inv'; end
    if mode_val == 1
        ov_mode_str = 'IPF Orientation Map';
        ov_name = sel_ipf; ov_sub = sel_dir;
        file_tag = sprintf('%s_%s%s_%s_a%03d_g%03d_%s', base_name, sel_bg, inv_tag, sel_ipf, round(alpha_val*100), round(gamma_val*100), timestamp_str);
    else
        ov_mode_str = 'Scalar Property / Strain Map';
        ov_name = sel_prop; ov_sub = sprintf('%s [%.3f, %.3f]', sel_cmap, cmin_val, cmax_val);
        file_tag = sprintf('%s_%s%s_Prop-%s-%s_a%03d_g%03d_%s', base_name, sel_bg, inv_tag, sel_prop, sel_cmap, round(alpha_val*100), round(gamma_val*100), timestamp_str);
    end
    clean_tag = strrep(strrep(file_tag, '.', '_'), ' ', '');
    export_path = fullfile(app.exportDir, clean_tag);
    
    fprintf('\n========================================================================\n');
    fprintf('SAVING PUBLICATION-GRADE COMPOSITE EXPORT...\n');
    
    % Use copyobj to directly clone the perfect, visible GUI axes without triggering MTEX axes reinitialization
    figExp = figure('Name', 'Export Renderer', 'Visible', 'on', 'Renderer', 'opengl', 'Position', [100, 100, 1050, 900], 'Color', 'w');
    axExp  = copyobj(app.ax, figExp);
    
    if mode_val == 2
        set(axExp, 'Units', 'normalized', 'Position', [0.08, 0.08, 0.76, 0.86]);
        try
            % Recreate colorbar for export (copyobj doesn't copy the colorbar object automatically)
            cmap_mat = colormap(app.ax);
            colormap(axExp, cmap_mat);
            try, clim(axExp, [cmin_val, cmax_val]); catch, caxis(axExp, [cmin_val, cmax_val]); end
            cbExp = colorbar(axExp, 'Location', 'eastoutside');
            ylabel(cbExp, strrep(sel_prop, '_', '\_'), 'FontWeight', 'bold', 'FontSize', 12);
        catch
        end
    else
        set(axExp, 'Units', 'normalized', 'Position', [0.08, 0.08, 0.84, 0.86]);
    end
    
    % Ensure base composite image layer remains strictly beneath all boundary vector overlays
    img_cloned = findobj(axExp, 'Type', 'image');
    if ~isempty(img_cloned), uistack(img_cloned, 'bottom'); end
    
    sel_bg_exp = sel_bg; if is_inverted, sel_bg_exp = [sel_bg, ' (Inverted)']; end
    title(axExp, sprintf('EBSD Composite: %s + %s (%s) | \\alpha=%.2f, \\gamma=%.2f', ...
          strrep(sel_bg_exp, '_', '\_'), strrep(ov_name, '_', '\_'), strrep(ov_sub, '_', '\_'), alpha_val, gamma_val), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Remove labeled X and Y axes on sides of the exported visualization while preserving title
    axis(axExp, 'off');
    set(axExp, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    xlabel(axExp, ''); ylabel(axExp, '');
          
    drawnow; % Ensure OpenGL buffer realization before rasterization
    img_filepath = [export_path, ext_fmt];
    exportScaledFigure(figExp, img_filepath, 'target_px', target_px, 'dpi', dpi_val, 'nodatestring');
    close(figExp);
    fprintf('  -> Successfully exported visualization image: %s\n', img_filepath);

    params = struct('dataset', app.datasetName, 'base_filename', base_name, 'background_property', sel_bg, ...
                    'is_background_inverted', is_inverted, 'overlay_mode', ov_mode_str, 'ipf_color_key', sel_ipf, 'ipf_direction', sel_dir, ...
                    'property_overlay', sel_prop, 'colormap', sel_cmap, 'cmin_val', cmin_val, 'cmax_val', cmax_val, ...
                    'alpha_overlay', alpha_val, 'restrict_above_bp', restrict_above_bp, 'restrict_below_wp', restrict_below_wp, ...
                    'ipf_thresh_property', sel_th_prop, 'ipf_thresh_value', thresh_val, ...
                    'mask_to_specimen', mask_to_specimen, 'crop_to_specimen', crop_to_specimen, 'show_micronbar', (get(ui.chk_micronbar, 'Value') == 1), ...
                    'black_point', bp_val, 'white_point', wp_val, 'gamma_val', gamma_val, ...
                    'grain_boundary_style', sel_gb, 'target_shorter_edge_px', target_px, 'dpi', dpi_val, ...
                    'export_timestamp', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'final_rgb_buffer', app.curr_rgb, ...
                    'image_filepath', img_filepath);
    mat_filepath = [export_path, '_Params.mat'];
    save(mat_filepath, 'params', '-v7.3');
    fprintf('  -> Successfully archived parameters and RGB buffer: %s\n', mat_filepath);

    md_filepath = [export_path, '_Report.md'];
    fid = fopen(md_filepath, 'w');
    fprintf(fid, '# Interactive EBSD & TKD Visualization Export Report\n\n');
    fprintf(fid, '- **Domain Dataset:** `%s`\n', app.datasetName);
    fprintf(fid, '## Exported Visualization Preview\n');
    fprintf(fid, '![Composite Map](%s%s)\n\n---\n\n', clean_tag, ext_fmt);
    
    fprintf(fid, '## Layer & Blending Configuration\n\n');
    fprintf(fid, '| Parameter | Selected Value | Functional Description |\n');
    fprintf(fid, '| :--- | :--- | :--- |\n');
    fprintf(fid, '| **Overlay Mode** | `%s` | Active visualization layer algorithm |\n', ov_mode_str);
    fprintf(fid, '| **Background Property** | `%s` | Extracted from domain EBSD properties |\n', sel_bg);
    fprintf(fid, '| **Color Inverted** | `%s` | Greyscale inversion ($1.0 - \\\\text{val}$) |\n', mat2str(is_inverted));
    fprintf(fid, '| **Black Point** | `%.3f` | Normalized scalar floor threshold ($0.00 = 1\\%%$ percentile) |\n', bp_val);
    fprintf(fid, '| **White Point** | `%.3f` | Normalized scalar ceiling threshold ($1.00 = 99\\%%$ percentile) |\n', wp_val);
    fprintf(fid, '| **Gamma ($\\\\gamma$)** | `%.2f` | Non-linear luminance midpoint scaling |\n', gamma_val);
    if mode_val == 1
        fprintf(fid, '| **IPF Color Key** | `%s` | MTEX Inverse Pole Figure coloring algorithm |\n', sel_ipf);
        fprintf(fid, '| **IPF Reference Direction** | `%s` | Projection vector axis |\n', sel_dir);
    else
        fprintf(fid, '| **Scalar Property / Strain** | `%s` | Property field projected as color layer |\n', sel_prop);
        fprintf(fid, '| **Colormap & Limits** | `%s [%.4f to %.4f]` | Colormap table and clipping boundaries |\n', sel_cmap, cmin_val, cmax_val);
    end
    fprintf(fid, '| **Overlay Alpha** | `%.2f` | Blending weight ($0.00 = \\\\text{pure grey}, 1.00 = \\\\text{pure overlay}$) |\n', alpha_val);
    fprintf(fid, '| **Restrict > BlackPt** | `%s` | Mask overlay exclusively to pixels with background $\\\\geq$ Black Point |\n', mat2str(restrict_above_bp));
    fprintf(fid, '| **Restrict < WhitePt** | `%s` | Mask overlay exclusively to pixels with background $\\\\leq$ White Point |\n', mat2str(restrict_below_wp));
    fprintf(fid, '| **Overlay Property Filter** | `%s >= %.4f` | Independent domain thresholding without affecting background |\n', sel_th_prop, thresh_val);
    fprintf(fid, '| **Specimen Mask Applied** | `%s` | Mask out overlay coloring outside specimen boundary |\n', mat2str(mask_to_specimen));
    fprintf(fid, '| **Map Cropped to Specimen** | `%s` | Dynamically fit axes bounding box around specimen boundary |\n', mat2str(crop_to_specimen));
    fprintf(fid, '| **Show Micron Bar** | `%s` | Render dynamic scale bar in bottom right corner |\n', mat2str(get(ui.chk_micronbar, 'Value') == 1));
    fprintf(fid, '| **Grain Boundaries** | `%s` | Boundary plotting style |\n\n---\n\n', sel_gb);
    
    fprintf(fid, '## Image Resolution & Hardware Specifications\n');
    fprintf(fid, '- **Target Shorter Edge:** `%d` px\n', target_px);
    fprintf(fid, '- **Raster DPI:** `%d`\n', dpi_val);
    fprintf(fid, '- **File Format:** `%s`\n', ext_fmt);
    fprintf(fid, '- **MAT Parameter Archive:** [`%s_Params.mat`](%s_Params.mat)\n', clean_tag, clean_tag);
    fprintf(fid, '- **Export Timestamp:** `%s`\n', params.export_timestamp);
    fclose(fid);
    fprintf('  -> Successfully formatted Markdown report: %s\n', md_filepath);
    fprintf('========================================================================\n');
    
    app.latestExport = params;
    app.opts = saveCurrentSettingsToOpts(fig);
    fig.UserData = app;
    
    if isNext
        app.batch_status = 'next';
        fig.UserData = app;
        uiresume(fig);
    else
        try
            msgbox(sprintf('Successfully exported:\n1) %s\n2) %s\n3) %s\n\nto:\n%s', ...
                   [clean_tag, ext_fmt], [clean_tag, '_Params.mat'], [clean_tag, '_Report.md'], app.exportDir), ...
                   'Export Complete');
        catch, end
    end
end
