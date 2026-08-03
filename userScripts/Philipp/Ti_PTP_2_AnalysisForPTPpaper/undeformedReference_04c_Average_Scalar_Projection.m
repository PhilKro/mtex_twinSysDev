%% undeformedReference_04c_Average_Scalar_Projection.m
% Standalone Evaluation of Scalar Projection Factor Statistics from Twin Analysis
% -------------------------------------------------------------------------
% Loads precomputed per-variant matrices (magnitudeX, shapeMismatch) from undeformedReference_04_TwinAnalysis
% and computes statistical metrics (mean, standard deviation, median) for the
% scalar projection factor, evaluated specifically for individual variants
% (default: Variant 6 where scalar projection factor > 0.80).
% -------------------------------------------------------------------------
clear; clc;

%% 0. Configuration
target_inc      = 'B00004'; % Increment prefix to evaluate
target_variant  = 6;        % Target twin variant ID to focus on (e.g., Variant 6)
min_mag_cutoff  = 0.80;     % Minimum scalar projection factor cutoff (e.g., > 0.80)
max_mismatch    = 0.50;     % Optional shape mismatch tolerance filter (default in undeformedReference_04_TwinAnalysis is 0.50)

% Resolve input .mat file path
try
    addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
    [registry, defaults] = a00_Dataset_Registry();
    resDir  = fullfile(defaults.dir, defaults.outDir_undeformed, 'Twin_Analysis_Results');
    matFile = fullfile(resDir, [target_inc, '_Twin_Analysis_Results.mat']);
catch
    % Hardcoded direct fallback path
    matFile = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2\Processed_Output_Undeformed\Twin_Analysis_Results\B00004_Twin_Analysis_Results.mat';
end

if ~exist(matFile, 'file')
    error('Target results file not found: %s\nPlease ensure undeformedReference_04_TwinAnalysis has been executed for increment %s.', matFile, target_inc);
end

%% 1. Load Precomputed Results
fprintf('========================================================================\n');
fprintf('  SCALAR PROJECTION FACTOR ANALYSIS (INCREMENT: %s)\n', target_inc);
fprintf('========================================================================\n');
fprintf('Loading data from:\n  %s\n\n', matFile);

data = load(matFile, 'magnitudeX', 'shapeMismatch');

% Ensure required per-variant matrices exist
if ~isfield(data, 'magnitudeX') || ~isfield(data, 'shapeMismatch')
    error('Required per-variant matrices (magnitudeX, shapeMismatch) not found in .mat file. Please re-run the updated undeformedReference_04_TwinAnalysis script first.');
end

N_vars    = size(data.magnitudeX, 1);
total_pts = size(data.magnitudeX, 2);

%% 2. Evaluate Target Variant Statistics (e.g., Variant 6 where x > 0.80)
v_mag  = data.magnitudeX(target_variant, :)';
v_mism = data.shapeMismatch(target_variant, :)';

% Condition A: Direct scalar projection factor threshold
mask_raw      = (v_mag > min_mag_cutoff) & ~isnan(v_mag);
vals_raw      = v_mag(mask_raw);

% Condition B: Including shape mismatch quality threshold
mask_filtered = mask_raw & (v_mism < max_mismatch);
vals_filtered = v_mag(mask_filtered);

fprintf('------------------------------------------------------------------------\n');
fprintf('  TARGET VARIANT #%d ANALYSIS (Scalar Projection Factor x > %.2f)\n', target_variant, min_mag_cutoff);
fprintf('------------------------------------------------------------------------\n');

if isempty(vals_raw)
    fprintf('  No datapoints found where Variant %d projection factor > %.2f.\n\n', target_variant, min_mag_cutoff);
else
    fprintf('1. Direct Threshold Criterion (Variant %d with x > %.2f):\n', target_variant, min_mag_cutoff);
    fprintf('   -> Matching Datapoints : %d (%.2f%% of all %d pixels)\n', numel(vals_raw), 100*numel(vals_raw)/total_pts, total_pts);
    fprintf('   -> MEAN Projection x   : %.4f\n', mean(vals_raw));
    fprintf('   -> Standard Deviation  : %.4f\n', std(vals_raw));
    fprintf('   -> Median              : %.4f\n', median(vals_raw));
    fprintf('   -> Min / Max           : %.4f / %.4f\n\n', min(vals_raw), max(vals_raw));

    fprintf('2. Quality Filtered Criterion (Variant %d with x > %.2f & Shape Mismatch < %.2f):\n', target_variant, min_mag_cutoff, max_mismatch);
    if isempty(vals_filtered)
        fprintf('   -> No datapoints satisfied both magnitude and shape mismatch criteria.\n\n');
    else
        fprintf('   -> Matching Datapoints : %d (%.2f%% of all %d pixels)\n', numel(vals_filtered), 100*numel(vals_filtered)/total_pts, total_pts);
        fprintf('   -> MEAN Projection x   : %.4f\n', mean(vals_filtered));
        fprintf('   -> Standard Deviation  : %.4f\n', std(vals_filtered));
        fprintf('   -> Median              : %.4f\n', median(vals_filtered));
        fprintf('   -> Min / Max           : %.4f / %.4f\n\n', min(vals_filtered), max(vals_filtered));
    end
end

%% 3. Summary Table Across All Variants (for comparison)
fprintf('========================================================================\n');
fprintf('  SUMMARY TABLE ACROSS ALL %d VARIANTS (Criterion: x > %.2f)\n', N_vars, min_mag_cutoff);
fprintf('========================================================================\n');
fprintf(' Var ID | Matching Pixels | Pixel %% | Mean (x) | Std Dev  | Median (x)\n');
fprintf(' -----------------------------------------------------------------------\n');

for v = 1:N_vars
    m_v = (data.magnitudeX(v, :)' > min_mag_cutoff) & ~isnan(data.magnitudeX(v, :)');
    v_vals = data.magnitudeX(v, m_v);
    
    if isempty(v_vals)
        fprintf('  # %3d | %15d | %6.2f%% | %8s | %8s | %10s\n', v, 0, 0.0, '-', '-', '-');
    else
        fprintf('  # %3d | %15d | %6.2f%% | %8.4f | %8.4f | %10.4f\n', ...
                v, numel(v_vals), 100*numel(v_vals)/total_pts, mean(v_vals), std(v_vals), median(v_vals));
    end
end
fprintf('========================================================================\n\n');
