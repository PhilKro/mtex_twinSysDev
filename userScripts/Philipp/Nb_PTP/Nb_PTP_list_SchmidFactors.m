%% Nb_PTP_list_SchmidFactors.m
% Lists all possible slip systems in Grain 1 and Grain 2, ranked by Schmid
% factor under uniaxial tension, and marks the systems identified as active.

clear; close all; clc;

%% 1. Initialization & Mtex Preferences
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Data Loading
baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
outDir  = fullfile(baseDir, 'Processed_Output');

% Load SSLIP data to get the reference slip systems and orientations
dic_prefix = 'B0005'; % Using B0005 as a representative increment to extract grains
sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', dic_prefix));

if ~exist(sslip_mat_path, 'file')
    error('SSLIP data not found: %s', sslip_mat_path);
end
fprintf('Loading SSLIP data from %s...\n', sslip_mat_path);
load(sslip_mat_path, 'ebsdID', 'sSLocalRef');

%% 3. Reconstruct grains to get orientations
grainThresh = 5 * degree;
minSize = 1000;
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

ori_g1 = grainsUndef(1).meanOrientation;
ori_g2 = grainsUndef(2).meanOrientation;

%% 4. Identify active slip systems from the struct setup
% Define major and minor active systems for Grain 1
% (Adjust these indices based on which ones you consider major vs minor)
major_g1_idx = [18, 30, 46, 37];
major_g1_mag = [0.0131,0.0027,0.0004, 0.006];
minor_g1_idx = [47, 11];
minor_g1_mag = [0.0005,0.0002];
active_g1_idx = unique([major_g1_idx, minor_g1_idx]);

% Define major and minor active systems for Grain 2
% (Adjust these indices based on which ones you consider major vs minor)
major_g2_idx = [27, 11, 47];
major_g2_mag = [0.0187,0.0036,0.0025];
minor_g2_idx = [36, 24, 33,  7, 29];
minor_g2_mag = [0.0007,0.0006,0.0006,0.0001,0.0001];
active_g2_idx = unique([major_g2_idx, minor_g2_idx]);

fprintf('Identified major active systems for Grain 1 (Indices): %s\n', num2str(major_g1_idx));
fprintf('Identified minor active systems for Grain 1 (Indices): %s\n', num2str(minor_g1_idx));
fprintf('Identified major active systems for Grain 2 (Indices): %s\n', num2str(major_g2_idx));
fprintf('Identified minor active systems for Grain 2 (Indices): %s\n\n', num2str(minor_g2_idx));

%% 5. Calculate Schmid Factors
stress = stressTensor.uniaxial(yvector);

% Get local slip systems for both grains
sS_g1 = ori_g1 * sSLocalRef;
sS_g2 = ori_g2 * sSLocalRef;

num_sys = length(sSLocalRef);

% Preallocate arrays for the table
Sys_Idx = (1:num_sys)';
b_ref = arrayfun(@(x) char(x.b.round), sSLocalRef(:), 'UniformOutput', false);
n_ref = arrayfun(@(x) char(x.n.round), sSLocalRef(:), 'UniformOutput', false);

SF_g1 = round(abs(sS_g1.SchmidFactor(stress)), 2);
SF_g1 = SF_g1(:);
SF_g2 = round(abs(sS_g2.SchmidFactor(stress)), 2);
SF_g2 = SF_g2(:);

% Map magnitudes
Mag_G1 = zeros(num_sys, 1);
for i = 1:length(major_g1_idx)
    Mag_G1(Sys_Idx == major_g1_idx(i)) = major_g1_mag(i);
end
for i = 1:length(minor_g1_idx)
    Mag_G1(Sys_Idx == minor_g1_idx(i)) = minor_g1_mag(i);
end

Mag_G2 = zeros(num_sys, 1);
for i = 1:length(major_g2_idx)
    Mag_G2(Sys_Idx == major_g2_idx(i)) = major_g2_mag(i);
end
for i = 1:length(minor_g2_idx)
    Mag_G2(Sys_Idx == minor_g2_idx(i)) = minor_g2_mag(i);
end

% Map activities to string categories
Activity_G1 = repmat({'-'}, num_sys, 1);
Activity_G1(ismember(Sys_Idx, major_g1_idx)) = {'Major'};
Activity_G1(ismember(Sys_Idx, minor_g1_idx)) = {'Minor'};

Activity_G2 = repmat({'-'}, num_sys, 1);
Activity_G2(ismember(Sys_Idx, major_g2_idx)) = {'Major'};
Activity_G2(ismember(Sys_Idx, minor_g2_idx)) = {'Minor'};

%% 6. Create and sort tables
T_G1 = table(Sys_Idx, b_ref, n_ref, SF_g1, Mag_G1, Activity_G1, 'VariableNames', {'Sys_Idx', 'b_ref', 'n_ref', 'Schmid_Factor', 'Magnitude', 'Activity'});
T_G1_sorted = sortrows(T_G1, 'Schmid_Factor', 'descend');
T_G1_sorted.SF_Rank = (1:height(T_G1_sorted))';
T_G1_sorted.Delta_SF = round(T_G1_sorted.Schmid_Factor(1) - T_G1_sorted.Schmid_Factor, 2);
T_G1_active = T_G1_sorted(~strcmp(T_G1_sorted.Activity, '-'), :);

T_G2 = table(Sys_Idx, b_ref, n_ref, SF_g2, Mag_G2, Activity_G2, 'VariableNames', {'Sys_Idx', 'b_ref', 'n_ref', 'Schmid_Factor', 'Magnitude', 'Activity'});
T_G2_sorted = sortrows(T_G2, 'Schmid_Factor', 'descend');
T_G2_sorted.SF_Rank = (1:height(T_G2_sorted))';
T_G2_sorted.Delta_SF = round(T_G2_sorted.Schmid_Factor(1) - T_G2_sorted.Schmid_Factor, 2);
T_G2_active = T_G2_sorted(~strcmp(T_G2_sorted.Activity, '-'), :);

%% 7. Display and Save Results
fprintf('==============================================================\nGRAIN 1: Active Slip Systems Ranked by Schmid Factor\n==============================================================\n');
disp(T_G1_active);
fprintf('\n==============================================================\nGRAIN 2: Active Slip Systems Ranked by Schmid Factor\n==============================================================\n');
disp(T_G2_active);

writetable(T_G1_active, fullfile(outDir, 'Grain1_SchmidFactors_Active.csv'));
writetable(T_G2_active, fullfile(outDir, 'Grain2_SchmidFactors_Active.csv'));
fprintf('Saved Grain 1 & Grain 2 tabular data to %s\n', outDir);

%% 8. Generate LaTeX Tables
fprintf('\n==============================================================\n');
fprintf('LaTeX Formatted Tables\n');
fprintf('==============================================================\n\n');

% --- Grain 1 LaTeX Table ---
fprintf('\\begin{table}[h!]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{cccccccc}\n');
fprintf('\\toprule\n');
fprintf('SF Rank & System Index & Slip Direction $\\mathbf{b}$ & Slip Plane $\\mathbf{n}$ & Schmid Factor & $\\Delta$ SF & Magnitude & Activity \\\\\\\\\n');
fprintf('\\midrule\n');
for i = 1:height(T_G1_active)
    sys_idx = T_G1_active.Sys_Idx(i);
    b_vals = round([sSLocalRef(sys_idx).b.u, sSLocalRef(sys_idx).b.v, sSLocalRef(sys_idx).b.w]);
    n_vals = round([sSLocalRef(sys_idx).n.h, sSLocalRef(sys_idx).n.k, sSLocalRef(sys_idx).n.l]);
    
    b_str = format_miller_latex(b_vals, true);
    n_str = format_miller_latex(n_vals, false);
    
    fprintf('  %d & %d & $%s$ & $%s$ & %.2f & %.2f & %.4f & %s \\\\\\\\\n', ...
        T_G1_active.SF_Rank(i), T_G1_active.Sys_Idx(i), b_str, n_str, ...
        T_G1_active.Schmid_Factor(i), T_G1_active.Delta_SF(i), T_G1_active.Magnitude(i), T_G1_active.Activity{i});
end
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\caption{Active slip systems ranked by Schmid factor under uniaxial tension for Grain 1.}\n');
fprintf('\\label{tab:sf_g1}\n');
fprintf('\\end{table}\n\n');

% --- Grain 2 LaTeX Table ---
fprintf('\\begin{table}[h!]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{cccccccc}\n');
fprintf('\\toprule\n');
fprintf('SF Rank & System Index & Slip Direction $\\mathbf{b}$ & Slip Plane $\\mathbf{n}$ & Schmid Factor & $\\Delta$ SF & Magnitude & Activity \\\\\\\\\n');
fprintf('\\midrule\n');
for i = 1:height(T_G2_active)
    sys_idx = T_G2_active.Sys_Idx(i);
    b_vals = round([sSLocalRef(sys_idx).b.u, sSLocalRef(sys_idx).b.v, sSLocalRef(sys_idx).b.w]);
    n_vals = round([sSLocalRef(sys_idx).n.h, sSLocalRef(sys_idx).n.k, sSLocalRef(sys_idx).n.l]);
    
    b_str = format_miller_latex(b_vals, true);
    n_str = format_miller_latex(n_vals, false);
    
    fprintf('  %d & %d & $%s$ & $%s$ & %.2f & %.2f & %.4f & %s \\\\\\\\\n', ...
        T_G2_active.SF_Rank(i), T_G2_active.Sys_Idx(i), b_str, n_str, ...
        T_G2_active.Schmid_Factor(i), T_G2_active.Delta_SF(i), T_G2_active.Magnitude(i), T_G2_active.Activity{i});
end
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\caption{Active slip systems ranked by Schmid factor under uniaxial tension for Grain 2.}\n');
fprintf('\\label{tab:sf_g2}\n');
fprintf('\\end{table}\n\n');

%% Helper Functions
function str = format_miller_latex(vals, is_dir)
    strs = cell(1, length(vals));
    for i = 1:length(vals)
        if vals(i) < 0
            strs{i} = sprintf('\\bar{%d}', abs(vals(i)));
        else
            strs{i} = sprintf('%d', vals(i));
        end
    end
    inner = strjoin(strs, '');
    if is_dir
        str = sprintf('[%s]', inner);
    else
        str = sprintf('(%s)', inner);
    end
end