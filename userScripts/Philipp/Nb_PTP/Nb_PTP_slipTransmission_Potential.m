%% Nb_PTP_slipTransmission_Potential.m
% Script to calculate potential slip transmission from one active system
% to ALL possible systems in a neighboring grain.

clear; close all; clc;

%% 1. Initialization & Mtex Preferences
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Data Loading
baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
outDir  = fullfile(baseDir, 'Processed_Output');

%% 3. Input Structs Setup
% Define a list of configurations to analyze.
% For each analysis, specify the source grain/slip system and the target grain.
analyses = struct();

analyses(1).dic_prefix = 'B0005';
analyses(1).grain1_id  = 2; % Source grain
analyses(1).ss1_rank   = 2; % Rank of the active slip system in the source grain
analyses(1).grain2_id  = 1; % Target grain
analyses(1).ss2_rank   = 1; % Actual active target system

analyses(2).dic_prefix = 'B0005';
analyses(2).grain1_id  = 2;
analyses(2).ss1_rank   = 4;
analyses(2).grain2_id  = 1;
analyses(2).ss2_rank   = 1;

analyses(3).dic_prefix = 'B0005';
analyses(3).grain1_id  = 2;
analyses(3).ss1_rank   = 5;
analyses(3).grain2_id  = 1;
analyses(3).ss2_rank   = 1;

% analyses(4).dic_prefix = 'B0003';
% analyses(4).grain1_id  = 2;
% analyses(4).ss1_rank   = 1;
% analyses(4).grain2_id  = 1;
% analyses(4).ss2_rank   = 2;
% 
% analyses(5).dic_prefix = 'B0007';
% analyses(5).grain1_id  = 2;
% analyses(5).ss1_rank   = 2;
% analyses(5).grain2_id  = 1;
% analyses(5).ss2_rank   = 1;
% 
% analyses(6).dic_prefix = 'B0007';
% analyses(6).grain1_id  = 2;
% analyses(6).ss2_rank   = 3;
% analyses(6).ss1_rank   = 2;
% analyses(6).grain2_id  = 1;

%% 4. Analysis loop
for i_analysis = 1:length(analyses)
    current_analysis = analyses(i_analysis);
    dic_prefix = current_analysis.dic_prefix;
    grain1_id  = current_analysis.grain1_id;
    ss1_rank   = current_analysis.ss1_rank;
    grain2_id  = current_analysis.grain2_id;
    ss2_rank   = current_analysis.ss2_rank;
    
    % Load SSLIP data
    sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', dic_prefix));
    if ~exist(sslip_mat_path, 'file')
        warning('SSLIP data not found for prefix %s. Skipping.', dic_prefix);
        continue;
    end
    load(sslip_mat_path, 'ebsdID', 'sSLocalRef');

    % Reconstruct grains
    grainThresh = 5 * degree;
    minSize = 1000;
    [grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
    ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
    [grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

    % Get the active slip system in the source grain (Grain 1)
    sslip_activity_sum = sum(abs(ebsdID.prop.slipIDcor), 2, 'omitnan'); % Rank by magnitude
    [~, sorted_sslip_idx] = sort(sslip_activity_sum, 'descend');
    ss1_idx = sorted_sslip_idx(ss1_rank);
    ss2_idx = sorted_sslip_idx(ss2_rank); % The actual active system in G2

    ori1 = grainsUndef(grain1_id).meanOrientation;
    ori2 = grainsUndef(grain2_id).meanOrientation;
    
    sS_g1_active = ori1 * sSLocalRef(ss1_idx);
    sS_g2_active = ori2 * sSLocalRef(ss2_idx);
    
    % Determine effective slip direction from the sign of the mean slip activity
    % This is now possible because posConstr=0 in SSLIP allows for negative slip activity.
    slip_sign_1 = sign(mean(ebsdID.prop.slipIDcor(ss1_idx, ebsdID.grainId == grain1_id), 'omitnan'));
    if slip_sign_1 == 0, slip_sign_1 = 1; end % Default to positive if mean is zero
    slip_sign_2 = sign(mean(ebsdID.prop.slipIDcor(ss2_idx, ebsdID.grainId == grain2_id), 'omitnan'));
    if slip_sign_2 == 0, slip_sign_2 = 1; end

    stress = stressTensor.uniaxial(yvector); % As requested for Schmid factor calc
    sf_g1_active = abs(sS_g1_active.SchmidFactor(stress));

    fprintf('\n\n================================================================\n');
    fprintf('--- Processing Analysis %d: Potential transmission from G%d (Rank %d, SF=%.3f) to G%d (Actual Rank %d) for prefix %s ---\n', ...
            i_analysis, grain1_id, ss1_rank, sf_g1_active, grain2_id, ss2_rank, dic_prefix);
    fprintf('================================================================\n');

    % --- Calculate primary m' and define the "Residual Mismatch Tensor" ---
    b1_eff = sS_g1_active.b * slip_sign_1;
    b2_eff = sS_g2_active.b * slip_sign_2;
    m_prime_active = cos(angle(b1_eff, b2_eff, 'noSymmetry')) * cos(angle(sS_g1_active.n, sS_g2_active.n, 'noSymmetry'));
    
    % Tensors for simple shear: M = b (x) n
    M1_eff = (b1_eff.normalize.xyz)' * sS_g1_active.n.normalize.xyz;
    M2_active_eff = (b2_eff.normalize.xyz)' * sS_g2_active.n.normalize.xyz;
    
    % Boundary mismatch left over after the primary active system accommodates its share
    H_residual = M1_eff - m_prime_active * M2_active_eff;

    % Get all possible slip systems for the target grain (Grain 2)
    sS_g2_all = ori2 * sSLocalRef;

    % --- Calculate parameters for all potential systems in Grain 2 ---
    num_systems = length(sS_g2_all);
    sys_indices = (1:num_systems)';
    b_vectors_ref = arrayfun(@(x) char(x.b.round), sSLocalRef, 'UniformOutput', false);
    n_vectors_ref = arrayfun(@(x) char(x.n.round), sSLocalRef, 'UniformOutput', false);
    m_prime_vals = nan(num_systems, 1);
    schmid_factors = nan(num_systems, 1);
    activity_ranks = nan(num_systems, 1);
    combined_metric = nan(num_systems, 1);
    residual_accom = nan(num_systems, 1);
    outofplane_comps = nan(num_systems, 1);

    for j = 1:num_systems
        sS_g2_current = sS_g2_all(j);
        
        % Luster-Morris parameter (m')
        % Compare effective incoming shear to the potential target shear. Allow for negative m'.
        m_prime_vals(j) = cos(angle(b1_eff, sS_g2_current.b, 'noSymmetry')) * cos(angle(sS_g1_active.n, sS_g2_current.n, 'noSymmetry'));
        
        % Schmid Factor (for uniaxial tension along y-vector)
        schmid_factors(j) = abs(sS_g2_current.SchmidFactor(stress));
        
        % Combined Transmission Factor
        % Use abs(m') as a large negative m' is still a strong (but opposing) geometric interaction
        combined_metric(j) = abs(m_prime_vals(j)) * schmid_factors(j);
        
        % How well does this system geometrically accommodate the remaining boundary mismatch? (Tensor Double Dot Product)
        Mj = sS_g2_current.b.normalize.xyz' * sS_g2_current.n.normalize.xyz;
        residual_accom(j) = sum(H_residual(:) .* Mj(:));
        
        % Out-of-plane component of Burgers vector (dot product with z-vector)
        outofplane_comps(j) = abs(dot(sS_g2_current.b.normalize, zvector));
        
        % Global Activity Rank of the system
        activity_ranks(j) = find(sorted_sslip_idx == j, 1);
    end

    % --- Create, Sort, and Display the Results Table ---
    results_table = table(sys_indices, activity_ranks, b_vectors_ref, n_vectors_ref, m_prime_vals, schmid_factors, combined_metric, residual_accom, outofplane_comps, ...
        'VariableNames', {'Sys_Idx', 'Activity_Rank', 'b_ref', 'n_ref', 'm_prime', 'SchmidFactor_Y', 'm_prime_x_SF', 'Res_Accom', 'b_z_component'});

    sorted_table = sortrows(results_table, 'm_prime_x_SF', 'descend');

    fprintf('\n--- Potential Transmission Systems in Grain %d (sorted by m'' x SF) ---\n', grain2_id);
    disp(sorted_table);

    % Save the table to a CSV file
    csv_filename = fullfile(outDir, sprintf('PotentialTransmission_%s_G%d-R%d_to_G%d.csv', dic_prefix, grain1_id, ss1_rank, grain2_id));
    writetable(sorted_table, csv_filename);
    fprintf('Saved potential transmission data to: %s\n', csv_filename);

    % --- Scatter Plot Visualization (Heatmap of Likelihood) ---
    fig_scatter = figure;
    fig_scatter.Name = 'Transmission Likelihood Heatmap';
    
    [M_grid, SF_grid] = meshgrid(linspace(0, 1, 200), linspace(0, 0.5, 200));
    Metric_grid = M_grid .* SF_grid;
    contourf(M_grid, SF_grid, Metric_grid, 30, 'LineStyle', 'none');
    mtexColorMap viridis;
    c = colorbar;
    c.Label.String = 'Transmission Metric (m'' \times SF)';
    hold on;
    
    % Find the actual active system in the sorted table
    idx_active = find(sorted_table.Sys_Idx == ss2_idx);
    idx_others = setdiff(1:height(sorted_table), idx_active);
    
    % Plot all non-active target systems
    scatter(sorted_table.m_prime(idx_others), sorted_table.SchmidFactor_Y(idx_others), 60, 'w', 'filled', 'MarkerEdgeColor', 'k');
    % Highlight the actual active target system
    if ~isempty(idx_active)
        scatter(sorted_table.m_prime(idx_active), sorted_table.SchmidFactor_Y(idx_active), 250, 'r', 'p', 'filled', 'MarkerEdgeColor', 'k');
    end
    
    % Label all systems with their System Index
    for k = 1:height(sorted_table)
        if k == idx_active
            text(sorted_table.m_prime(k) + 0.02, sorted_table.SchmidFactor_Y(k), sprintf('Active: Sys %d', sorted_table.Sys_Idx(k)), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);
        else
            text(sorted_table.m_prime(k) + 0.01, sorted_table.SchmidFactor_Y(k), num2str(sorted_table.Sys_Idx(k)), 'Color', 'w', 'FontSize', 8);
        end
    end

    xlabel('Luster-Morris Parameter (m'')');
    ylabel('Target Schmid Factor (SF)');
    title(sprintf('Transmission Likelihood: G%d (Sys %d) -> G%d', grain1_id, ss1_idx, grain2_id));
    
    scatter_path = fullfile(outDir, sprintf('PotentialTransmission_Scatter_%s_G%d-R%d_to_G%d.jpg', dic_prefix, grain1_id, ss1_rank, grain2_id));
    exportScaledFigure(fig_scatter, scatter_path, 'target_px', 2500, 'noDateString');
    fprintf('Saved likelihood heatmap to: %s\n', scatter_path);
    close(fig_scatter);

    % --- Visualization ---
    % Plot the source system and the actually active target system
    
    % Create a combined activity map for visualization of the source activity
    activity_map = nan(size(ebsdID));
    in_grain1 = (ebsdID.grainId == grain1_id);
    slip_data_1 = abs(ebsdID.prop.slipIDcor(ss1_idx, :));
    activity_map(in_grain1) = slip_data_1(in_grain1);
    
    fig_potential = figure;
    plot(ebsdID, activity_map, 'micronbar', 'off'); 
    mtexColorMap viridis;
    hold on;

    plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
    text(grainsUndef, grainsUndef.id, 'color', 'w', 'fontsize', 14);

    % Plot slip traces (red) and Burgers vectors (blue) for source and actual active target
    quiver(grainsUndef(grain1_id), sS_g1_active.trace, 'color', 'r', 'linewidth', 3);
    quiver(grainsUndef(grain2_id), sS_g2_active.trace, 'color', 'r', 'linewidth', 3, 'linestyle', '--');
    
    b1 = sS_g1_active.b.normalize;
    b2 = sS_g2_active.b.normalize;
    quiver(grainsUndef(grain1_id), b1, 'color', 'b', 'linewidth', 3);
    quiver(grainsUndef(grain2_id), b2, 'color', 'b', 'linewidth', 3, 'linestyle', '--');

    % Create text for the info box
    info_str = sprintf(['Analysis for %s:\n' ...
        'Source: G%d, Active Sys Rank %d (Sys %d), SF = %.3f\n  b=[%s], n=(%s)\n\n' ...
        'Actual Active Target in G%d (Rank %d):\n  Sys %d with m''xSF = %.3f\n  m'' = %.3f, SF = %.3f\n  b=[%s], n=(%s)'], ...
        dic_prefix, grain1_id, ss1_rank, ss1_idx, sf_g1_active, ...
        char(sSLocalRef(ss1_idx).b.round), char(sSLocalRef(ss1_idx).n.round), ...
        grain2_id, ss2_rank, ...
        ss2_idx, sorted_table.m_prime_x_SF(idx_active), sorted_table.m_prime(idx_active), sorted_table.SchmidFactor_Y(idx_active), ...
        char(sSLocalRef(ss2_idx).b.round), char(sSLocalRef(ss2_idx).n.round));
        
    annotation('textbox', [0.02 0.7 0.3 0.25], 'String', info_str, ...
        'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'FitBoxToText', 'on', 'FontSize', 10, 'Margin', 5);
        
    title(sprintf('Potential Transmission from G%d to G%d', grain1_id, grain2_id));
    
    % Export the figure
    fig_path = fullfile(outDir, sprintf('PotentialTransmission_Plot_%s_G%d-R%d_to_G%d.jpg', dic_prefix, grain1_id, ss1_rank, grain2_id));
    exportScaledFigure(fig_potential, fig_path, 'target_px', 3500, 'noDateString');
    fprintf('Saved summary plot to: %s\n', fig_path);
    
    close(fig_potential);
end