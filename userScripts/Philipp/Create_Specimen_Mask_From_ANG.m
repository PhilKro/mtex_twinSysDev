%% Create_Specimen_Mask_From_ANG.m
% A modular script to interactively create and save specimen masks for raw
% .ang files based on a folder and file pattern.

clear; close all; clc;

%% 1. Setup
% Add path to helper functions if not already in the MATLAB path
% Ensure 'interactive_specimen_mask.m' is in this path.
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

% MTEX preferences
plotx2east;
plotzIntoPlane;

%% 2. User Configuration
% --- Define the dataset to process ---
baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
filePattern = '**/*BW158.ang'; % e.g., '*.ang', '**/*_some_pattern.ang'

% --- Define Material and EBSD loading parameters ---
% Define Crystal Symmetry (change for your material)
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

tkd_tilt = -7 * degree; % TKD tilt angle for rotation correction
thresh_ci = 0.1;        % Confidence Index threshold for basic filtering

% --- Overwrite existing masks? ---
% 'ask': prompts for each file
% 'always': overwrites without asking
% 'never': skips files with existing masks
overwrite_mode = 'ask'; % 'ask', 'always', 'never'

%% 3. Processing Loop
fprintf('--- Starting Specimen Mask Creation ---\n');
fprintf('Searching for files matching "%s" in "%s"\n\n', filePattern, baseDir);

% Find all .ang files matching the pattern (searches recursively)
angFiles = dir(fullfile(baseDir, filePattern));

if isempty(angFiles)
    error('No .ang files found matching the pattern in the specified directory.');
end

for i = 1:length(angFiles)
    currentFile = angFiles(i).name;
    ang_path = fullfile(angFiles(i).folder, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    fprintf('Processing %d/%d: %s\n', i, length(angFiles), currentFile);
    
    % Define the output .mat filename
    output_mat_path = fullfile(angFiles(i).folder, sprintf('%s_specimen_mask.mat', baseName));
    
    % --- Handle existing files based on overwrite_mode ---
    if exist(output_mat_path, 'file')
        should_skip = false;
        switch lower(overwrite_mode)
            case 'ask'
                reply = input(sprintf('  -> Mask file already exists. Overwrite? (y/n): '), 's');
                if ~strcmpi(reply, 'y')
                    fprintf('  -> Skipping file.\n\n');
                    should_skip = true;
                end
            case 'never'
                fprintf('  -> Mask file exists. Skipping.\n\n');
                should_skip = true;
            case 'always'
                fprintf('  -> Mask file exists. Overwriting.\n');
        end
        if should_skip, continue; end
    end
    
    % --- Load and Process EBSD Data ---
    fprintf('  -> Loading EBSD data...\n');
    % Load the raw 1D EBSD data
    ebsd_1d = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % Apply standard rotations for TKD alignment
    ebsd_1d = rotate(ebsd_1d, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd_1d = rotate(ebsd_1d, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    % Base filtering
    ebsd_1d(ebsd_1d.ci < thresh_ci) = 'notIndexed';
    
    % Gridify the data to get a 2D matrix for the interactive tool
    ebsd_grid = ebsd_1d.gridify;
    
    % Get the CI matrix for the mask function
    ci_matrix = ebsd_grid.ci;
    ci_matrix(isnan(ci_matrix)) = 0; % The tool requires a clean matrix
    
    % --- Call the interactive mask function ---
    fprintf('  -> Please draw the specimen mask...\n');
    [mask_2d, ~] = interactive_specimen_mask(ci_matrix);
    
    % --- Map the 2D grid mask back to the original 1D EBSD data ---
    specimenlocation = false(length(ebsd_1d), 1);
    valid_ids_on_grid = ebsd_grid.id(~isnan(ebsd_grid.id));
    mask_values_on_grid = mask_2d(~isnan(ebsd_grid.id));
    specimenlocation(valid_ids_on_grid) = mask_values_on_grid;
    
    % --- Save the 1D mask to a .mat file ---
    fprintf('  -> Saving mask to %s\n', output_mat_path);
    save(output_mat_path, 'specimenlocation');
    
    fprintf('  -> Done.\n\n');
    
    % Close the figure created by the interactive tool
    close(gcf);
end

fprintf('--- All files processed. ---\n');