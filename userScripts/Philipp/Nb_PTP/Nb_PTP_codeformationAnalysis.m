%% Nb_PTP_codeformationAnalysis.m
% Script to analyze the total deformation, co-deformation between grains,
% and the role of individual slip systems in accommodating compatibility.
%
% This addresses 3 main questions:
% 1. How good does the sum of identified slip systems account for the total deformation?
% 2. A measure of co-deformation of the two involved grains.
% 3. How do diffuse slip systems aid in the co-deformation?

% 1. Total Deformation vs Identified Slip: It directly compares the residualEeff (the strain unaccounted for by SSLIP) to the total Eeff extracted from DIC. It calculates the percentage of deformation unexplained by your slip systems for each grain and creates a visual map.
% 2. Measure of Co-Deformation: It creates 3D macroscopic strain tensors (averaging across the grain) for Grain 1 and Grain 2 based on both the original DIC data and the reconstructed SSLIP data. The macroscopic incompatibility is then quantified by the Frobenius norm of the difference between these tensors (||E2 - E1||).
% 3. Role of Diffuse Slip Systems: For each grain, the script isolates the "dominant" slip systems (those within 10% of the maximum magnitude) and the subsequent "diffuse/minor" systems. It computes the "missing" compatibility strain (the strain required to deform identically to the neighbor). By taking the tensor inner product of the diffuse slip's strain and the missing compatibility strain, we can tell if the diffuse slip is geometrically driven by the boundary mismatch. A positive inner product confirms the mismatch causes a driving force that triggers this diffuse slip to reduce incompatibility!

clear; close all; clc;

%% 1. Initialization & Directories
plotx2east; plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% --- Classification Parameters ---
dominant_threshold = 0.50; % Slip magnitude > 50% of max is considered "dominant"
diffuse_threshold  = 0.10; % Slip magnitude between 10% and 50% of max is "diffuse"
% ---------------------------------

if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir  = fullfile(baseDir, 'Processed_Output');

% Select which increment to analyze (can be adapted into a loop)
dic_prefix = 'B0005'; 

% Load SSLIP Data
sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', dic_prefix));
if ~exist(sslip_mat_path, 'file')
    error('SSLIP data not found for prefix %s.', dic_prefix);
end
fprintf('Loading SSLIP data for %s...\n', dic_prefix);
load(sslip_mat_path, 'ebsdID', 'sSLocalRef');

% Reconstruct grains
grainThresh = 5 * degree;
minSize = 1000;
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
[grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

% Assume the two largest grains are our bicrystal
[~, sorted_grain_idx] = sort(grainsUndef.numPixel, 'descend');
g1_id = grainsUndef(sorted_grain_idx(1)).id;
g2_id = grainsUndef(sorted_grain_idx(2)).id;

fprintf('\nAnalyzing Grain %d and Grain %d for prefix %s\n', g1_id, g2_id, dic_prefix);

% Plot Grain 1 and Grain 2 map
fig_grains = figure;
plot(ebsdID, ebsdID.orientations, 'micronbar', 'off');
hold on;
plot(grainsUndef.boundary, 'linecolor', 'k', 'linewidth', 2);
text(grainsUndef, grainsUndef.id, 'color', 'k', 'fontsize', 16, 'fontweight', 'bold', 'backgroundcolor', 'w', 'margin', 1);
title(sprintf('Grain Map (Prefix %s)', dic_prefix));
drawNow(gcm, 'figSize', 'large');

mask_g1 = (ebsdID.grainId == g1_id);
mask_g2 = (ebsdID.grainId == g2_id);

%% Q1: How good does the sum of identified slip systems account for total deformation?
% The SSLIP method inherently outputs this as the 'residualEeff'
warning('somehow residualEeff is not a property of EBSDId although SSLIP should return it???')
res_g1 = mean(ebsdID.prop.residualEeff(mask_g1), 'omitnan');
res_g2 = mean(ebsdID.prop.residualEeff(mask_g2), 'omitnan');

Eeff_DIC_g1 = mean(ebsdID.prop.Eeff(mask_g1), 'omitnan');
Eeff_DIC_g2 = mean(ebsdID.prop.Eeff(mask_g2), 'omitnan');

fprintf('\n=================================================================\n');
fprintf('Q1: SSLIP Reconstructed Deformation vs Total DIC Deformation\n');
fprintf('=================================================================\n');
fprintf('Grain %d: Mean DIC Eff Strain = %.4f | Mean SSLIP Residual = %.4f (%.1f%% unexplained)\n', ...
    g1_id, Eeff_DIC_g1, res_g1, (res_g1/Eeff_DIC_g1)*100);
fprintf('Grain %d: Mean DIC Eff Strain = %.4f | Mean SSLIP Residual = %.4f (%.1f%% unexplained)\n', ...
    g2_id, Eeff_DIC_g2, res_g2, (res_g2/Eeff_DIC_g2)*100);

% Visualization of Residuals
fig1 = figure;
plot(ebsdID, ebsdID.prop.residualEeff, 'micronbar', 'off');
hold on;
plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
text(grainsUndef, grainsUndef.id, 'color', 'w', 'fontsize', 14);
mtexColorMap viridis;
setColorRange([0, prctile(ebsdID.prop.residualEeff(:), 95)]);
mtexColorbar('title', 'SSLIP Residual Eff Strain');
title('Residual Effective Strain (Unaccounted Deformation)');

%% Helper variables for Q2 and Q3
% Calculate deformation gradients dynamically from the DIC displacement fields
% (Displacements are in pixels, so gradient over a 1-pixel grid directly gives unitless strain)
U_disp = ebsdID.prop.(sprintf('%s_Displacements_x_displacement', dic_prefix));
V_disp = ebsdID.prop.(sprintf('%s_Displacements_y_displacement', dic_prefix));

[H11, H12] = gradient(double(U_disp));
[H21, H22] = gradient(double(V_disp));

E_DIC_g1_mat = get_mean_DIC_strain_mat(mask_g1, H11, H12, H21, H22);
E_DIC_g2_mat = get_mean_DIC_strain_mat(mask_g2, H11, H12, H21, H22);

% Flatten the logical masks to match the 1D shape of slipIDcor
slip_mags_g1 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g1(:))), 2, 'omitnan');
slip_mags_g2 = mean(abs(ebsdID.prop.slipIDcor(:, mask_g2(:))), 2, 'omitnan');

ori1 = grainsUndef(g1_id).meanOrientation;
ori2 = grainsUndef(g2_id).meanOrientation;

sSLocal_1 = ori1 * sSLocalRef;
sSLocal_2 = ori2 * sSLocalRef;

H_slip_g1_mat = get_reconstruct_slip_strain_mat(sSLocal_1, slip_mags_g1);
E_slip_g1_mat = 0.5 * (H_slip_g1_mat + H_slip_g1_mat');

H_slip_g2_mat = get_reconstruct_slip_strain_mat(sSLocal_2, slip_mags_g2);
E_slip_g2_mat = 0.5 * (H_slip_g2_mat + H_slip_g2_mat');

%% Q2: Measure of Co-deformation
fprintf('\n=================================================================\n');
fprintf('Q2: Measure of Co-deformation (Macroscopic Compatibility)\n');
fprintf('=================================================================\n');
Delta_E_DIC = E_DIC_g2_mat - E_DIC_g1_mat;
incompat_DIC = norm(Delta_E_DIC, 'fro');
fprintf('Norm of DIC Strain Difference ||E_2 - E_1||:   %.4f\n', incompat_DIC);

Delta_E_slip = E_slip_g2_mat - E_slip_g1_mat;
incompat_slip = norm(Delta_E_slip, 'fro');
fprintf('Norm of SSLIP Strain Difference ||E_2 - E_1||: %.4f\n', incompat_slip);

%% Q3: How do diffuse slip systems aid in co-deformation?
fprintf('\n=================================================================\n');
fprintf('Q3: Role of Diffuse Slip in Co-deformation\n');
fprintf('=================================================================\n');
fprintf('If the inner product of the diffuse strain with the "missing" boundary strain is POSITIVE,\n');
fprintf('it means the boundary mismatch generates a driving force that favors this diffuse slip.\n\n');

% Analyze Grain 1 adapting to Grain 2
analyze_slip_role(g1_id, slip_mags_g1, sSLocal_1, E_slip_g1_mat, E_slip_g2_mat, sprintf('Grain %d', g2_id), ...
    dominant_threshold, diffuse_threshold, ebsdID, grainsUndef, mask_g1);

% Analyze Grain 2 adapting to Grain 1
analyze_slip_role(g2_id, slip_mags_g2, sSLocal_2, E_slip_g2_mat, E_slip_g1_mat, sprintf('Grain %d', g1_id), ...
    dominant_threshold, diffuse_threshold, ebsdID, grainsUndef, mask_g2);


%% Helper Functions

function E_mat = get_mean_DIC_strain_mat(mask, H11, H12, H21, H22)
    H = [mean(H11(mask), 'omitnan'), mean(H12(mask), 'omitnan'), 0; ...
         mean(H21(mask), 'omitnan'), mean(H22(mask), 'omitnan'), 0; ...
         0, 0, 0];
    E_mat = 0.5 * (H + H');
end

function H_slip_mat = get_reconstruct_slip_strain_mat(sS, slip_mags)
    H_slip_mat = zeros(3,3);
    for i = 1:length(sS)
        if slip_mags(i) > 0
                b_vec = sS(i).b.xyz';
                n_vec = sS(i).n.xyz';
            M = b_vec * n_vec';
            H_slip_mat = H_slip_mat + slip_mags(i) * M;
        end
    end
end

function analyze_slip_role(grain_id, slip_mags, sSLocal, E_grain_total_mat, E_neighbor_mat, neighbor_name, dom_thresh, diffuse_thresh, ebsdID, grainsUndef, mask)
    [sorted_mags, sort_idx] = sort(slip_mags, 'descend');
    
    % Find dominant systems
    max_mag = sorted_mags(1);
    dom_mask = sorted_mags > dom_thresh * max_mag;
    num_dom = sum(dom_mask);
    if num_dom == 0
        fprintf('No dominant slip systems found for Grain %d.\n', grain_id);
        return; 
    end
    
    dom_idx = sort_idx(1:num_dom);
    
    % Find diffuse systems
    diffuse_mask = (sorted_mags <= dom_thresh * max_mag) & (sorted_mags > diffuse_thresh * max_mag);
    diffuse_idx = sort_idx(diffuse_mask);
    num_diffuse = length(diffuse_idx);
    
    incompat_actual = norm(E_neighbor_mat - E_grain_total_mat, 'fro');
    
    fprintf('--- Analysis for Grain %d adapting to %s ---\n', grain_id, neighbor_name);
    fprintf('Number of dominant slip systems: %d\n', num_dom);
    fprintf('Number of diffuse slip systems: %d\n', num_diffuse);
    fprintf('Actual Incompatibility (All slips): %.4f\n', incompat_actual);
    
    % Set up plotting for this grain
    num_plots = num_dom + num_diffuse;
    if num_plots > 0
        fig_slips = figure;
        fig_slips.Name = sprintf('Slip Systems - Grain %d', grain_id);
        plot_cols = min(4, num_plots);
        plot_rows = ceil(num_plots / plot_cols);
        newMtexFigure('layout', [plot_rows, plot_cols]);
        
        plot_counter = 1;
        
        % Plot Dominant Systems
        for i = 1:num_dom
            nextAxis(ceil(plot_counter/plot_cols), mod(plot_counter-1, plot_cols)+1);
            sys_idx = dom_idx(i);
            
            % Isolate data for this grain to plot
            slip_data = abs(ebsdID.prop.slipIDcor(sys_idx, :));
            plot_data = nan(1, length(ebsdID));
            plot_data(mask(:)) = slip_data(mask(:));
            
            plot(ebsdID, plot_data, 'micronbar', 'off');
            mtexColorMap viridis;
            hold on;
            plot(grainsUndef.boundary, 'linewidth', 1, 'linecolor', 'w');
            
            sS_sys = sSLocal(sys_idx);
            title(sprintf('DOMINANT\nSys %d (Rank %d)\nMag: %.4f\nb=[%s] n=(%s)', ...
                sys_idx, i, slip_mags(sys_idx), char(sS_sys.b.round), char(sS_sys.n.round)), 'FontSize', 10);
            
            plot_counter = plot_counter + 1;
        end
    end
    
    for i = 1:length(diffuse_idx)
        sys_idx = diffuse_idx(i);
        sys_mag = slip_mags(sys_idx);
        if sys_mag == 0, continue; end
        
        mags_diffuse = zeros(size(slip_mags));
        mags_diffuse(sys_idx) = sys_mag;
        H_diffuse_mat = get_reconstruct_slip_strain_mat(sSLocal, mags_diffuse);
        E_diffuse_mat = 0.5 * (H_diffuse_mat + H_diffuse_mat');
        
        E_missing_hypothetical = E_neighbor_mat - (E_grain_total_mat - E_diffuse_mat);
        incompat_without_diffuse = norm(E_missing_hypothetical, 'fro');
        
        % Tensor inner product (double dot product)
        inner_prod = sum(E_missing_hypothetical(:) .* E_diffuse_mat(:));
        
        fprintf('  Diffuse Slip Rank %d (System Index %d, Magnitude=%.4f):\n', num_dom + i, sys_idx, sys_mag);
        
        fprintf('    Inner Product with Missing Strain: %.6f ', inner_prod);
        if inner_prod > 0
            fprintf('-> POSITIVE (Mismatch acts as driving force!)\n');
        else
            fprintf('-> NEGATIVE (Not driven by this boundary mismatch)\n');
        end
        
        fprintf('    Hypothetical Incompatibility (without diffuse): %.4f ', incompat_without_diffuse);
        if incompat_actual < incompat_without_diffuse
            fprintf('(REDUCED by %.2f%% - aids co-deformation)\n', 100*(incompat_without_diffuse - incompat_actual)/incompat_without_diffuse);
        else
            fprintf('(INCREASED by %.2f%% - does not aid co-deformation)\n', 100*(incompat_actual - incompat_without_diffuse)/incompat_without_diffuse);
        end
        
        % Plot Diffuse System
        if num_plots > 0
            nextAxis(ceil(plot_counter/plot_cols), mod(plot_counter-1, plot_cols)+1);
            
            slip_data = abs(ebsdID.prop.slipIDcor(sys_idx, :));
            plot_data = nan(1, length(ebsdID));
            plot_data(mask(:)) = slip_data(mask(:));
            
            plot(ebsdID, plot_data, 'micronbar', 'off');
            mtexColorMap parula; % Different colormap to distinguish diffuse
            hold on;
            plot(grainsUndef.boundary, 'linewidth', 1, 'linecolor', 'w');
            
            sS_sys = sSLocal(sys_idx);
            
            if inner_prod > 0
                title_color = [0 0.5 0]; role_str = '(Aids)';
            else
                title_color = [0.8 0 0]; role_str = '(Hinders)';
            end
            
            t = title(sprintf('DIFFUSE %s\nSys %d (Rank %d)\nMag: %.4f\nb=[%s] n=(%s)', ...
                role_str, sys_idx, num_dom + i, sys_mag, char(sS_sys.b.round), char(sS_sys.n.round)), 'FontSize', 10);
            t.Color = title_color;
            
            plot_counter = plot_counter + 1;
        end
    end
    
    if num_plots > 0
        drawNow(gcm, 'figSize', 'large');
    end
    fprintf('\n');
end