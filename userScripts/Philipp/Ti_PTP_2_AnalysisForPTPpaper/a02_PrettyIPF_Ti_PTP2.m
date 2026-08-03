%% 02_PrettyIPF_Ti_PTP2.m
% Script to loop over Ti_PTP2 .ang files and generate pretty IPF Y maps
% overlaid on a tunable greyscale Confidence Index (CI) background,
% cropping to the interactive specimen mask.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Tuning Parameters
baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2';
outDir = fullfile(baseDir, 'PrettyPlots');
maskDir = fullfile(baseDir, 'Masks'); % Where masks from step 01 are saved

if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output directory: %s\n', outDir);
end

% --- "PHOTOSHOP CURVES" TUNING PARAMETERS ---
% 1. Black Point (0.0 to 1.0): CI values below this will be pure black.
blackPoint = 0.05; 

% 2. White Point (0.0 to 1.0): CI values above this will be pure white.
whitePoint = 0.98; 

% 3. Gamma (Midtones): 
%    < 1.0 makes the mid-range lighter (good if the map is too dark)
%    > 1.0 makes the mid-range darker (good if the map looks washed out)
%    = 1.0 applies no gamma curve.
gammaVal = 1.5; 

% 4. IPF Overlay Opacity (0.0 to 1.0)
ipf_alpha = 0.55;
% --------------------------------------------

%% 3. Load Data & Loop over .ang files
% Define Crystal Symmetry for Titanium
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

tkd_tilt = -7 * degree; % TKD tilt angle for rotation correction
thresh_ci = 0.1;

% Load Dataset Registry to inspect exclusion flags
[registry, ~] = a00_Dataset_Registry();
excluded_tkd_fnames = {registry([registry.exclude] & ~cellfun(@isempty, {registry.tkd_fname})).tkd_fname};
clean_excluded = cellfun(@(x) strrep(x, '.ang', ''), excluded_tkd_fnames, 'UniformOutput', false);

% Find all .ang files matching the pattern
angFiles = dir(fullfile(baseDir, '**', '*.ang'));

if isempty(angFiles)
    error('No .ang files found in %s', baseDir);
end

for i = 1:length(angFiles)
    currentFile = angFiles(i).name;
    ang_path = fullfile(angFiles(i).folder, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    if any(strcmpi(clean_excluded, baseName))
        fprintf('Skipping (%d/%d): [%s] -> Marked as EXCLUDED in a00_Dataset_Registry().\n', i, length(angFiles), baseName);
        continue;
    end
    
    fprintf('Processing %d/%d: %s\n', i, length(angFiles), currentFile);
    
    % Load EBSD Data
    ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % Load optional extra properties from .h5oina if extracted
    prop_file = fullfile(maskDir, sprintf('%s_h5oina_props.mat', baseName));
    if exist(prop_file, 'file')
        props = load(prop_file);
        if isfield(props, 'BS') && length(props.BS) == length(ebsd), ebsd.prop.BS = props.BS; end
        if isfield(props, 'BC') && length(props.BC) == length(ebsd), ebsd.prop.BC = props.BC; end
        if isfield(props, 'FSE') && length(props.FSE) == length(ebsd), ebsd.prop.FSE = props.FSE; end
        if isfield(props, 'BSE') && length(props.BSE) == length(ebsd), ebsd.prop.BSE = props.BSE; end
        if isfield(props, 'MSE') && length(props.MSE) == length(ebsd), ebsd.prop.MSE = props.MSE; end
        if isfield(props, 'GradMag') && length(props.GradMag) == length(ebsd), ebsd.prop.GradMag = props.GradMag; end
        if isfield(props, 'LapVar') && length(props.LapVar) == length(ebsd), ebsd.prop.LapVar = props.LapVar; end
    end
    
    % Apply standard rotations for TKD alignment
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Base filtering
    ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
    
    ebsd = ebsd.gridify;
    
    % Store a copy of the unmasked data for the FSE plots
    ebsd_full = ebsd;
    
    % --- Calculate Tuned CI Greyscale (Curves Adjustment) ---
    ci_raw = ebsd.ci;
    ci_min = min(ci_raw(:));
    ci_max = max(ci_raw(:));
    
    ci_norm = (ci_raw - ci_min) ./ (ci_max - ci_min);
    ci_tuned = (ci_norm - blackPoint) ./ (whitePoint - blackPoint);
    ci_tuned(ci_tuned < 0) = 0; % clip darks
    ci_tuned(ci_tuned > 1) = 1; % clip highlights
    ci_tuned = ci_tuned .^ gammaVal; % Apply gamma
    
    cibw = [ci_tuned(:), ci_tuned(:), ci_tuned(:)];
    cibw(~ebsd.isIndexed, :) = 0; % Mask out the background (black)

    % --- Apply Specimen Mask ---
    % Look for a specimen mask .mat file containing 'specimenlocation'
    mask_file = fullfile(maskDir, sprintf('%s_specimen_mask.mat', baseName));
    
    if exist(mask_file, 'file')
        fprintf('  -> Applying specimen mask...\n');
        load(mask_file, 'specimenlocation');
        ebsd(~specimenlocation) = 'notIndexed';
        
        % --- Crop to the bounding box of the specimen location ---
        valid_ebsd = ebsd(specimenlocation);
        if length(valid_ebsd) > 0
            xmin = min(valid_ebsd.x); xmax = max(valid_ebsd.x);
            ymin = min(valid_ebsd.y); ymax = max(valid_ebsd.y);
            ebsd = ebsd(inpolygon(ebsd, [xmin, ymin, xmax-xmin, ymax-ymin]));
        end
    else
        warning('  -> No specimen mask found at %s. Proceeding with unmasked data.', mask_file);
    end
    
    % --- Plotting IPF Y ---
    ipf = ipfColorKey(ebsd.CS);
    ipf.inversePoleFigureDirection = yvector;
    
    fig = figure('Name', currentFile, 'Visible', 'off');
    plot(ebsd, cibw, 'micronbar', 'on');
    hold on;

    colors = ipf.orientation2color(ebsd('indexed').orientations);
    plot(ebsd('indexed'), colors, 'FaceAlpha', ipf_alpha);
    
    title(['IPF Y - ', baseName], 'Interpreter', 'none');
    
    % --- Export Figure ---
    out_name = fullfile(outDir, sprintf('%s_Y_IPF_TunedCI.jpg', baseName)); 
    exportScaledFigure(fig, out_name, 'noDateString');
    close(fig);
    
    fprintf('  -> Saved plot to %s\n', out_name);
    
    % --- Optional: Plot IPF on FSE and Side-by-Side Plots ---
    if isfield(ebsd_full.prop, 'FSE') && isfield(ebsd_full.prop, 'BSE') && isfield(ebsd_full.prop, 'MSE')
        % 1. IPF over FSE (Unmasked)
        fse_raw = ebsd_full.prop.FSE;
        fse_min = min(fse_raw(:));
        fse_max = max(fse_raw(:));
        fse_norm = (fse_raw - fse_min) ./ (fse_max - fse_min);
        fsebw = [fse_norm(:), fse_norm(:), fse_norm(:)];
        fsebw(~ebsd_full.isIndexed, :) = 0;
        
        fig2 = figure('Name', currentFile, 'Visible', 'off');
        plot(ebsd_full, fsebw, 'micronbar', 'on');
        hold on;
        colors_full = ipf.orientation2color(ebsd_full('indexed').orientations);
        plot(ebsd_full('indexed'), colors_full, 'FaceAlpha', ipf_alpha);
        title(['IPF Y on FSE (Unmasked) - ', baseName], 'Interpreter', 'none');
        
        out_name2 = fullfile(outDir, sprintf('%s_Y_IPF_FSE.jpg', baseName)); 
        exportScaledFigure(fig2, out_name2, 'noDateString');
        close(fig2);
        fprintf('  -> Saved plot to %s\n', out_name2);
        
        % 2. Side-by-Side BSE, MSE, FSE, GradMag, LapVar (Unmasked)
        % Using mtexFigure layout to cleanly place them side-by-side, making it visible to the user
        fig3 = mtexFigure('layout', [1, 5], 'Visible', 'on');
        
        plot(ebsd_full, ebsd_full.prop.BSE, 'micronbar', 'off');
        title('BSE');
        
        nextAxis;
        plot(ebsd_full, ebsd_full.prop.MSE, 'micronbar', 'off');
        title('MSE');
        
        nextAxis;
        plot(ebsd_full, ebsd_full.prop.FSE, 'micronbar', 'on');
        % caxis([fse_min + 0.5*(fse_max - fse_min), fse_max]);
        title('FSE');
        
        nextAxis;
        if isfield(ebsd_full.prop, 'GradMag')
            plot(ebsd_full, ebsd_full.prop.GradMag, 'micronbar', 'off');
        end
        title('GradMag');
        
        nextAxis;
        if isfield(ebsd_full.prop, 'LapVar')
            plot(ebsd_full, ebsd_full.prop.LapVar, 'micronbar', 'off');
        end
        title('LapVar');
        
        % Apply greyscale colormap to the entire figure
        colormap(gcf, 'viridis');
        
        out_name3 = fullfile(outDir, sprintf('%s_BSE_MSE_FSE_Grad_Lap.jpg', baseName));
        % Standard saveas on current figure
        saveas(gcf, out_name3); 
        close(gcf);
        fprintf('  -> Saved plot to %s\n', out_name3);
        
        % 3. Side-by-Side CI, BS, BC, FSE (Unmasked) with Viridis colormap
        fig4 = mtexFigure('layout', [1, 4], 'Visible', 'on');
        
        plot(ebsd_full, ebsd_full.ci, 'micronbar', 'off');
        title('CI');
        
        nextAxis;
        if isfield(ebsd_full.prop, 'BS')
            plot(ebsd_full, ebsd_full.prop.BS, 'micronbar', 'off');
        end
        title('BS');
        
        nextAxis;
        if isfield(ebsd_full.prop, 'BC')
            plot(ebsd_full, ebsd_full.prop.BC, 'micronbar', 'off');
        end
        title('BC');
        
        nextAxis;
        plot(ebsd_full, ebsd_full.prop.FSE, 'micronbar', 'on');
        caxis([fse_min + 0.5*(fse_max - fse_min), fse_max]);
        title('FSE');
        
        % Apply viridis colormap (or fallback to parula if viridis is missing in older MATLAB)
        try
            colormap(gcf, 'viridis');
        catch
            colormap(gcf, 'parula');
        end
        
        out_name4 = fullfile(outDir, sprintf('%s_CI_BS_BC_FSE.jpg', baseName));
        saveas(gcf, out_name4);
        close(gcf);
        fprintf('  -> Saved plot to %s\n', out_name4);
    end
    
end

fprintf('\nFinished generating pretty IPF plots for Ti_PTP_2!\n');
