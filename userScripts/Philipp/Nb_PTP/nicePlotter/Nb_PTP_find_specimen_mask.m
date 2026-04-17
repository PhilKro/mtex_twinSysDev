%% Nb_PTP_find_specimen_mask.m
% Script to detect and isolate the true physical specimen from EBSD data.
% Uses a "Manual Edge Tracing" method with Magnetic Snapping: Computes a 
% gradient image, allows the user to trace rough boundaries, and snaps 
% those boundaries to the local gradient maximums.

%% 0. Initialization & Data Loading
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end

% Define increments (mapping TKD .ang files to concurrent DIC .hdf5 prefixes)
increments = struct();
increments(1).name = 'Increment 0';  increments(1).ang_file = 'insitu Specimen 1 0 Map Data 1_BW263.ang';
increments(2).name = 'Increment 2';  increments(2).ang_file = 'insitu Specimen 1 1 Map Data 2_BW263.ang';
increments(3).name = 'Increment 4';  increments(3).ang_file = 'insitu Specimen 1 2 Map Data 3_BW263.ang';
increments(4).name = 'Increment 6';  increments(4).ang_file = 'insitu Specimen 1 3 Map Data 4_BW263.ang';
increments(5).name = 'Increment 8';  increments(5).ang_file = 'insitu Specimen 1 4 Map Data 5_BW263.ang';
increments(6).name = 'Increment 10'; increments(6).ang_file = 'insitu Specimen 1 5 Map Data 6_BW263.ang';
increments(7).name = 'Increment 13'; increments(7).ang_file = 'insitu Specimen 1 7 Map Data 7_BW263.ang';

% Choose which increment to load
inc_idx = 7; % CHANGE THIS to the desired increment index
currentInc = increments(inc_idx);

if isempty(currentInc.ang_file)
    error('No .ang file specified for %s.', currentInc.name);
end

tkd_tilt = -10 * degree;
thresh_ci = 0.13;
CS = crystalSymmetry.load('Nb-Niobium.cif');

ang_path = fullfile(baseDir, currentInc.ang_file);
fprintf('Loading EBSD data from: %s\n', currentInc.ang_file);
ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');

%% 1. Grid the Data
fprintf('Gridding the EBSD data...\n');
ebsd_grid = ebsd.gridify;

%% 2. Extract the Map
fprintf('Extracting the CI map...\n');
if isfield(ebsd_grid.prop, 'ci')
    ci_matrix = ebsd_grid.prop.ci;
else
    ci_matrix = ebsd_grid.ci;
end

ci_matrix(isnan(ci_matrix)) = 0;

%% 3. Load Saved Parameters or Initialize Defaults
[ang_dir, ang_name, ~] = fileparts(ang_path);
param_file = fullfile(ang_dir, sprintf('%s_mask_params.mat', ang_name));

if exist(param_file, 'file')
    fprintf('Loading saved edge parameters...\n');
    load(param_file, 'p');
    % Handle structural updates if loading an old parameter file format
    if ~isfield(p, 'disp_min')
        p.disp_min = 0.0;
        p.disp_max = 0.1;
        p.left_edge_x = [size(ci_matrix,2)*0.2; size(ci_matrix,2)*0.2];
        p.left_edge_y = [1; size(ci_matrix,1)];
        p.right_edge_x = [size(ci_matrix,2)*0.8; size(ci_matrix,2)*0.8];
        p.right_edge_y = [1; size(ci_matrix,1)];
    end
    if ~isfield(p, 'snap_window')
        p.snap_window = 20;
    end
else
    % Initialize default parameters for the Edge Tracing method
    p.sigma = 2.0;       % Smoothing for gradient calculation
    p.disp_min = 0.0;    % Display contrast min
    p.disp_max = 0.1;    % Display contrast max
    p.snap_window = 20;  % Snapping window radius in pixels
    
    % Default boundary lines (straight down the sides)
    p.left_edge_x = [size(ci_matrix,2)*0.2; size(ci_matrix,2)*0.2];
    p.left_edge_y = [1; size(ci_matrix,1)];
    
    p.right_edge_x = [size(ci_matrix,2)*0.8; size(ci_matrix,2)*0.8];
    p.right_edge_y = [1; size(ci_matrix,1)];
end

%% 4. Interactive GUI for Mask Tuning
fprintf('Opening interactive mask tuning GUI...\n');
[mask_2d, p] = interactive_specimen_mask(ci_matrix, p);

% Save the parameters automatically if they were confirmed
save(param_file, 'p');
fprintf('Parameters saved to: %s\n', param_file);

%% 5. Assign the Property
fprintf('Assigning the mask to ebsd.prop.specimenlocation...\n');

% Flatten grid coordinates and the mask to map back to 1D EBSD
grid_x = ebsd_grid.x(:);
grid_y = ebsd_grid.y(:);
mask_flat = mask_2d(:);

% Remove padded NaNs that gridify introduces
valid_idx = ~isnan(grid_x) & ~isnan(grid_y);
grid_x = grid_x(valid_idx);
grid_y = grid_y(valid_idx);
mask_flat = mask_flat(valid_idx);

% Map back to the original 1D ebsd points using Nearest Neighbor
idx = knnsearch([grid_x, grid_y], [ebsd.x, ebsd.y]);
ebsd.prop.specimenlocation = mask_flat(idx);

% Save the mask to a .mat file for easy loading in fullAnalysis
mask_file = fullfile(ang_dir, sprintf('%s_specimen_mask.mat', ang_name));
specimenlocation = ebsd.prop.specimenlocation;
save(mask_file, 'specimenlocation');
fprintf('Mask successfully saved to: %s\n', mask_file);

fprintf('Done! Filter your data using: ebsd_clean = ebsd(ebsd.prop.specimenlocation);\n');

%% 6. Verify and Plot
figure('Name', 'Mask Verification');
plot(ebsd, ebsd.ci, 'micronbar', 'off'); 
colormap gray; hold on;
plot(ebsd(ebsd.prop.specimenlocation), 'FaceColor', 'r', 'FaceAlpha', 0.3);
title('1D EBSD Map with Specimen Mask Overlay (Red)', 'Interpreter', 'none');