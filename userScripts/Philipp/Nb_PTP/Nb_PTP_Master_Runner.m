%% Master Runner Script for PTP Nb Data Analysis
% This script sequentially runs the Importer and the Analysis script for 
% all specified increments, keeping the workspace clean and automated.

clear; close all; clc;
is_master_run = true; % Signals child scripts to skip their 'clear' commands

% Provide the name of the MAT file containing the fully merged EBSD & DIC 
% data from the Undeformed script (typically Increment 0)
% The analysis script relies on this!
% Use pwd to create an absolute path and bypass OneDrive relative path bugs
master_pwd = pwd;
undeformed_dic_mat = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP', 'Processed_Output', 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');

% Choose which increments to run (e.g. 1 through 6)
increments_to_run = 1:14; 

try
    load(undeformed_dic_mat);
catch
    warning('No undeformed mat found running importer script.');
end


for i = 1:length(increments_to_run)
    inc_idx = increments_to_run(i);
    
    fprintf('\n=======================================================\n');
    fprintf('STARTING PROCESSING FOR INCREMENT INDEX: %d\n', inc_idx);
    fprintf('=======================================================\n');
    
    % --- 1. Set Input struct for Script01 ---
    currentInc = get_currentInc_struct(inc_idx);
    
    
    % --- 2. Run the Analysis Script ---
    run('C:\Users\phkr\MTEX_workEnv\userScripts\Philipp\Nb_PTP\Nb_PTP_fullAnalysis.m');
    
    % Clean up variables from the scripts to avoid cross-contamination
    % but keep the variables needed for the main loop.
    clearvars -except is_master_run undeformed_dic_mat newEBSD increments_to_run i inc_idx master_pwd
end

fprintf('\nAll specified increments processed successfully!\n');

%% Local function to get input structure for each increment
function currentInc = get_currentInc_struct(idx)
    increments = struct();
    
    increments(1).name = 'Increment 0';
    increments(1).ang_file = 'insitu Specimen 1 0 Map Data 1_BW263.ang';
    increments(1).dic_prefix = 'B0001';

    increments(2).name = 'Increment 1';
    increments(2).ang_file = '';
    increments(2).dic_prefix = 'B0002';
    
    increments(3).name = 'Increment 2';
    increments(3).ang_file = 'insitu Specimen 1 1 Map Data 2_BW263.ang';
    increments(3).dic_prefix = 'B0003';
    
    increments(4).name = 'Increment 3';
    increments(4).ang_file = '';
    increments(4).dic_prefix = 'B0004';
    
    increments(5).name = 'Increment 4';
    increments(5).ang_file = 'insitu Specimen 1 2 Map Data 3_BW263.ang';
    increments(5).dic_prefix = 'B0005';
    
    increments(6).name = 'Increment 5';
    increments(6).ang_file = '';
    increments(6).dic_prefix = 'B0006';
    
    increments(7).name = 'Increment 6';
    increments(7).ang_file = 'insitu Specimen 1 3 Map Data 4_BW263.ang';
    increments(7).dic_prefix = 'B0007';
    
    increments(8).name = 'Increment 7';
    increments(8).ang_file = '';
    increments(8).dic_prefix = 'B0008';
    
    increments(9).name = 'Increment 8';
    increments(9).ang_file = 'insitu Specimen 1 4 Map Data 5_BW263.ang';
    increments(9).dic_prefix = 'B0009';
    
    increments(10).name = 'Increment 9';
    increments(10).ang_file = '';
    increments(10).dic_prefix = 'B0010';
    
    increments(11).name = 'Increment 10';
    increments(11).ang_file = 'insitu Specimen 1 5 Map Data 6_BW263.ang';
    increments(11).dic_prefix = 'B0011';
    
    increments(12).name = 'Increment 11';
    increments(12).ang_file = '';
    increments(12).dic_prefix = 'B0012';
    
    increments(13).name = 'Increment 12';
    increments(13).ang_file = '';
    increments(13).dic_prefix = 'B0013';
    
    increments(14).name = 'Increment 13';
    increments(14).ang_file = 'insitu Specimen 1 7 Map Data 7_BW263.ang';
    increments(14).dic_prefix = 'B0014';
    
    currentInc = increments(idx);
end