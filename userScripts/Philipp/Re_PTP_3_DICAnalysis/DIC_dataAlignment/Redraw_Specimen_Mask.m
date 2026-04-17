%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Redraw Specimen Mask for Aligned EBSD Data                          %%%
%%%                                                                     %%%
%%% This script allows a user to interactively redraw the specimen mask %%%
%%% for a series of pre-processed .mat files containing EBSD data.      %%%
%%% It uses the 'interactive_specimen_mask' function to trace the       %%%
%%% specimen boundaries and then overwrites the original .mat file      %%%
%%% with the updated EBSD object containing the new mask.               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% 1. Setup
% Add path to helper functions if not already in the MATLAB path
% Ensure 'interactive_specimen_mask.m' is in this path.
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');

% MTEX preferences
plotx2east;
plotzIntoPlane;

%% 2. Define Input Files
% Define the base path for your data files
inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';

% List the base filenames (without the '_MTEX6.1.mat' suffix) that you
% want to process.
fnames_base = {
    fullfile(inputPath,'inc_0_newEBSD__withStrain_05-08-2025_10-27_MTEX6.1.mat'),...
    fullfile(inputPath,'inc_4_newEBSD__withStrain_05-08-2025_23-58_MTEX6.1.mat'),...
    fullfile(inputPath,'inc_6_newEBSD__withStrain_05-12-2025_15-07_MTEX6.1.mat'),...
    fullfile(inputPath,'inc_7_newEBSD__withStrain_05-12-2025_15-16_MTEX6.1.mat'),...
    fullfile(inputPath,'inc_8_newEBSD__withStrain_05-12-2025_15-23_MTEX6.1.mat')
    };

%% 3. Processing Loop
fprintf('--- Starting Specimen Mask Redrawing Process ---\n');

for i = 1:length(fnames_base)
    
    mat_filename = fnames_base{i};
    
    % Check if the file exists
    if ~exist(mat_filename, 'file')
        warning('File not found, skipping: %s', mat_filename);
        continue;
    end
    
    fprintf('Processing file: %s\n', mat_filename);
    
    % Load all variables from the .mat file into a struct
    loaded_data = load(mat_filename);
    
    % --- Find the EBSD object ---
    % The EBSD object might be named 'newEBSD' or 'ebsd'
    if isfield(loaded_data, 'newEBSD')
        ebsd = loaded_data.newEBSD;
        ebsd_var_name = 'newEBSD';
    elseif isfield(loaded_data, 'ebsd')
        ebsd = loaded_data.ebsd;
        ebsd_var_name = 'ebsd';
    else
        warning('Could not find an EBSD object named "newEBSD" or "ebsd" in the file. Skipping.');
        continue;
    end
    
    % Ensure data is gridified for matrix operations
    if ~isa(ebsd,'EBSDsquare')
        ebsd = ebsd.gridify;
    end
    
    % --- Get CI matrix for the mask function ---
    if isfield(ebsd.prop, 'ci')
        ci_matrix = ebsd.prop.ci;
    else
        ci_matrix = ebsd.ci;
    end
    
    % The mask function requires a clean matrix with no NaNs
    ci_matrix(isnan(ci_matrix)) = 0;
    
    % --- Initialize mask parameters ---
    % Check if old mask parameters exist to initialize the GUI
    init_p = struct();
    if isfield(ebsd.prop, 'specimen_mask_params')
        init_p = ebsd.prop.specimen_mask_params;
        fprintf('Found existing mask parameters. Initializing GUI with them.\n');
    end
    
    % --- Call the interactive mask function ---
    [final_mask, updated_p] = interactive_specimen_mask(ci_matrix);%, init_p);
    
    % --- Update the EBSD object and save ---
    ebsd.prop.specimenlocation = final_mask(:);
    ebsd.prop.specimen_mask_params = updated_p; % Store parameters for future edits
    loaded_data.(ebsd_var_name) = ebsd;
    
    fprintf('Saving updated mask to file...\n');
    save(mat_filename, '-struct', 'loaded_data');
    
    fprintf('Done with file %d of %d.\n\n', i, length(fnames_base));
end

fprintf('--- All files processed. ---\n');