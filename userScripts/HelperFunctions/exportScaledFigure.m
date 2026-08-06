function exportScaledFigure(fig, fname, varargin)
% exportScaledFigure Export a figure with scaled shorter side
%
%   exportScaledFigure(fig, fname)
%       exports the figure FIG to file FNAME with default settings:
%       (shorter side = 1300 px at 300 dpi), and appends a datetime string
%       before the file extension.
%
%   exportScaledFigure(fig, fname, 'target_px', N, 'dpi', D, 'noDateString')
%       allows you to override target size and dpi, and disable the date string.
%
%   Example:
%       f = figure;
%       plot(rand(10,1));
%       exportScaledFigure(f, 'myplot.png') 
%       % → myplot_09_04_15_30_12.png
%
%       exportScaledFigure(f, 'myplot.png', 'noDateString')
%       % → myplot.png

    % --- Modern Argument Validation & Default Setup ---
    p = inputParser;
    p.CaseSensitive = false;
    addRequired(p, 'fig', @(x) isgraphics(x, 'figure') || isgraphics(x, 'axes'));
    addRequired(p, 'fname', @(x) ischar(x) || isstring(x));
    addParameter(p, 'target_px', 1500, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'dpi', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'nodatestring', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'extraheight', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'extrawidth', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    
    % Support keyword-only flags like 'nodatestring', 'extraheight' passed as single arguments
    customArgs = varargin;
    k = 1;
    while k <= numel(customArgs)
        if ischar(customArgs{k}) || isstring(customArgs{k})
            key = lower(char(customArgs{k}));
            if ismember(key, {'nodatestring', 'extraheight', 'extrawidth'})
                if k == numel(customArgs) || ischar(customArgs{k+1}) || isstring(customArgs{k+1})
                    customArgs = [customArgs(1:k), {true}, customArgs(k+1:end)];
                end
            end
        end
        k = k + 1;
    end
    parse(p, fig, fname, customArgs{:});
    opts = p.Results;

    % --- Resolve Filepath & Ensure Directory Structure Exists ---
    fname = char(opts.fname);
    [filepath, name, ext] = fileparts(fname);
    if isempty(ext), ext = '.png'; end
    if ~opts.nodatestring
        dtstr = string(datetime('now', 'Format', 'MM_dd_HH_mm_ss'));
        name = char(name + "_" + dtstr);
    end
    if ~isempty(filepath) && ~exist(filepath, 'dir')
        mkdir(filepath);
    end
    fname = fullfile(filepath, [name, ext]);

    % --- Exception-Safe Figure Unit Preservation ---
    oldUnits = fig.Units;
    cleanupObj = onCleanup(@() set(fig, 'Units', oldUnits));
    fig.Units = 'inches';
    pos = fig.Position; % [x, y, width, height] in inches

    % --- MATHEMATICAL SOUNDNESS FIX: Apply Aspect Ratio Modifiers FIRST ---
    % Any modification to height or width must precede scale calculation,
    % otherwise changing dimensions can flip which side is shortest, invalidating target_in.
    if opts.extraheight
        pos(4) = pos(4) * 1.3;
    elseif opts.extrawidth
        pos(3) = pos(3) * 1.3;
    end
    
    % --- Calculate Precise Scale Factor for Shorter Side ---
    target_in = opts.target_px / opts.dpi;
    shorter_side = min(pos(3:4));
    scale = target_in / shorter_side;
    
    fig.Position(3:4) = pos(3:4) * scale;

    % --- Compute actual export dimensions in pixels ---
    width_px  = round(fig.Position(3) * opts.dpi);
    height_px = round(fig.Position(4) * opts.dpi);

    % --- Export via modern high-resolution engine ---
    drawnow; % Realize window transformations before capturing raster
    exportgraphics(fig, fname, 'Resolution', opts.dpi);

    % --- Report export telemetry (cleanupObj automatically restores fig.Units upon return) ---
    fprintf('Exported %s at %dx%d px (%d dpi)\n', fname, width_px, height_px, opts.dpi);
end
