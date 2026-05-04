%% Re_PTP3_prettyIPF_loop.m
% Script to loop over Re_PTP_3 .ang files and generate pretty IPF Y maps
% overlaid on a tunable greyscale Confidence Index (CI) background.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Tuning Parameters
baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
outDir = fullfile(baseDir, 'PrettyPlots');

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
% Define Crystal Symmetry for Rhenium
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

tkd_tilt = -7 * degree; % Adjust this depending on your exact TKD geometry setup if needed
thresh_ci = 0.1;

% Find all .ang files matching the insitu Re sequence (searches recursively)
angFiles = dir(fullfile(baseDir, '**', 'Re_7deg_insitu Specimen*.ang'));

if isempty(angFiles)
    error('No .ang files found matching "Re_7deg_insitu Specimen*.ang" in %s', baseDir);
end

for i = 1:length(angFiles)
    currentFile = angFiles(i).name;
    ang_path = fullfile(angFiles(i).folder, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    fprintf('Processing %d/%d: %s\n', i, length(angFiles), currentFile);
    
    % Load EBSD Data
    ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % Apply standard rotations for TKD alignment
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Base filtering
    ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
    

    ebsd = ebsd.gridify;
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
    mask_file = fullfile(angFiles(i).folder, sprintf('%s_specimen_mask.mat', baseName));
    
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
    out_name = fullfile(outDir, sprintf('%s_Y_IPF_TunedCI.jpg', baseName)); %sprintf('%s_Y_IPF_TunedCI.jpg', baseName)
    exportScaledFigure(gcf, out_name, 'noDateString');
    close(fig);
    
    fprintf('  -> Saved plot to %s\n', out_name);
end

fprintf('\nFinished generating pretty IPF plots for Re_PTP_3!\n');