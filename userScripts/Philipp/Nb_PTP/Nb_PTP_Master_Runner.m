%% Master Runner Script for PTP Nb Data Analysis
% This script sequentially runs the Importer and the Analysis script for 
% all specified increments, keeping the workspace clean and automated.

clear; close all; clc;
is_master_run = true; % Signals child scripts to skip their 'clear' commands

% Provide the name of the MAT file containing the fully merged EBSD & DIC 
% data from the Undeformed script (typically Increment 0)
% The analysis script relies on this!
undeformed_dic_mat = fullfile('onedriveData', 'Philipp', 'Nb_PTP', 'Processed_Output', 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');

% Choose which increments to run (e.g. 1 through 6)
increments_to_run = 1:6; 

for i = 1:length(increments_to_run)
    inc_idx = increments_to_run(i);
    
    fprintf('\n=======================================================\n');
    fprintf('STARTING PROCESSING FOR INCREMENT INDEX: %d\n', inc_idx);
    fprintf('=======================================================\n');
    
    % --- 1. Set Input struct for Script01 ---
    input = get_input_struct(inc_idx);
    
    % --- 2. Run the Importer Script ---
    if ~isempty(input)
        run('Script01_IMPORTER_Nb_Data_TKD_DIC_SE_alignment.m');
    else
        warning('No input defined for increment %d. Skipping Importer.', inc_idx);
    end
    
    % --- 3. Run the Analysis Script ---
    run('Nb_PTP_fullAnalysis.m');
    
    % Clean up variables from the scripts to avoid cross-contamination
    % but keep the variables needed for the main loop.
    clearvars -except is_master_run undeformed_dic_mat increments_to_run i inc_idx
end

fprintf('\nAll specified increments processed successfully!\n');

%% Local function to get input structure for each increment
function input = get_input_struct(idx)
    % Base defaults
    defaults.DIC_filename = 'Processing_subset=21_step=6';
    
    switch idx
        case 1 % Increment 0
            input = struct('fname', 'insitu Specimen 1 0 Map Data 1', 'SEImagefname', '0_10kv_DIC.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0001.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 0 Map Data 1_BW263_20251126_115536_Selection_Points_Backup_EBSD.mat');
        case 2 % Increment 1
            input = struct('fname', 'insitu Specimen 1 1 Map Data 2', 'SEImagefname', '1postTKD_10kv_DIC_0d.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0003.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 1 Map Data 2_BW263_20251126_132847_Selection_Points_Backup_EBSD.mat'); 
        case 3 % Increment 2
            input = struct('fname', 'insitu Specimen 1 2 Map Data 3', 'SEImagefname', '2postTKD_10kv_DIC_0d.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0005.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 2 Map Data3_BW263_20251126_134943_Selection_Points_Backup_%sEBSD.mat'); 
        case 4 % Increment 3
            input = struct('fname', 'insitu Specimen 1 3 Map Data 4', 'SEImagefname', '3postTKD_10kv_DIC_0.tif', ...
                'DICfoldername', fullfile('DIC_h5', 'Processing_subset=19_step=1'), 'DICfname', 'B0007.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 3 Map Data 4_BW263_20251127_154633_Selection_Points_Backup_%sEBSD.mat');
        case 5 % Increment 4
            input = struct('fname', 'insitu Specimen 1 4 Map Data 5', 'SEImagefname', '4postTKD_10kv_DIC_0.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0009.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 4 Map Data 5_BW263_20251201_092539_Selection_Points_Backup_%sEBSD.mat');
        case 6 % Increment 5
            input = struct('fname', 'insitu Specimen 1 5 Map Data 6', 'SEImagefname', '5postTKD_Take2_10kv_DIC_0.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0011.hdf5', ...
                'sel_points_fname', 'insitu Specimen 1 5 Map Data 6_BW263_20251201_100658_Selection_Points_Backup_%sEBSD.mat'); 
        case 7 % Increment 6
            input = struct('fname', 'insitu Specimen 1 7 Map Data 7', 'SEImagefname', '7b_afterDescend_10kv_DIC_0.tif', ...
                'DICfoldername', fullfile('DIC_h5', defaults.DIC_filename), 'DICfname', 'B0014.hdf5', ...
                'sel_points_fname', '');
        otherwise
            input = [];
    end
end