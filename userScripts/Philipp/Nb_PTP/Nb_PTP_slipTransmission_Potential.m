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
analyses(1).sg_id      = 2; % Source grain ID
analyses(1).sg_ss_idx  = 18; % Index of the active slip system in the source grain
analyses(1).tg_id      = 1; % Target grain ID
analyses(1).tg_ss_rank = 27; % Actual active target system rank for comparison

analyses(2).dic_prefix = 'B0005';
analyses(2).sg_id      = 2;
analyses(2).sg_ss_idx  = 11;
analyses(2).tg_id      = 1; 
analyses(2).tg_ss_rank = 27;

analyses(3).dic_prefix = 'B0005';
analyses(3).sg_id      = 2;
analyses(3).sg_ss_idx  = 47;
analyses(3).tg_id      = 1; 
analyses(3).tg_ss_rank = 27;

analyses(4).dic_prefix = 'B0005';
analyses(4).sg_id      = 2;
analyses(4).sg_ss_idx  = 36;
analyses(4).tg_id      = 1; 
analyses(4).tg_ss_rank = 27;

analyses(5).dic_prefix = 'B0005';
analyses(5).sg_id      = 2;
analyses(5).sg_ss_idx  = 24;
analyses(5).tg_id      = 1; 
analyses(5).tg_ss_rank = 27;

% 
% analyses(5).dic_prefix = 'B0007';
% analyses(5).sg_id      = 2;
% analyses(5).sg_ss_idx  = 2;
% analyses(5).tg_id      = 1;
% analyses(5).tg_ss_rank = 1;
% 
% analyses(6).dic_prefix = 'B0007';
% analyses(6).sg_id      = 2;
% analyses(6).sg_ss_idx  = 2;
% analyses(6).tg_id      = 1;
% analyses(6).tg_ss_rank = 3;
%
% analyses(7).dic_prefix = 'B0003';
% analyses(7).sg_id      = 2;
% analyses(7).sg_ss_idx  = 1;
% analyses(7).tg_id      = 1;
% analyses(7).tg_ss_rank = 2;
%% 4. Analysis loop
for i_analysis = 1:length(analyses)
    current_analysis = analyses(i_analysis);
    dic_prefix = current_analysis.dic_prefix;
    sg_id      = current_analysis.sg_id;
    sg_ss_idx  = current_analysis.sg_ss_idx;
    tg_id      = current_analysis.tg_id;
    tg_ss_rank = current_analysis.tg_ss_rank;
    
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

    % Rank slip systems by overall activity to find the rank of the actual
    % active system in the target grain for comparison.
    sslip_activity_sum = sum(abs(ebsdID.prop.slipIDcor), 2, 'omitnan'); % Rank by magnitude
    [~, sorted_sslip_idx] = sort(sslip_activity_sum, 'descend');
    
    % The source system is now given by its index directly.
    % The actual active system in the target grain is found by its rank.
    tg_ss_idx = sorted_sslip_idx(tg_ss_rank);

    ori_sg = grainsUndef(sg_id).meanOrientation;
    ori_tg = grainsUndef(tg_id).meanOrientation;
    
    sS_sg_active = ori_sg * sSLocalRef(sg_ss_idx);
    sS_tg_active = ori_tg * sSLocalRef(tg_ss_idx);
    
    % Determine effective slip direction from the sign of the mean slip activity
    % This is now possible because posConstr=0 in SSLIP allows for negative slip activity.
    slip_sign_sg = sign(mean(ebsdID.prop.slipIDcor(sg_ss_idx, ebsdID.grainId == sg_id), 'omitnan'));
    if slip_sign_sg == 0, slip_sign_sg = 1; end % Default to positive if mean is zero
    slip_sign_tg = sign(mean(ebsdID.prop.slipIDcor(tg_ss_idx, ebsdID.grainId == tg_id), 'omitnan'));
    if slip_sign_tg == 0, slip_sign_tg = 1; end

    stress = stressTensor.uniaxial(yvector); % As requested for Schmid factor calc
    sf_sg_active = abs(sS_sg_active.SchmidFactor(stress));

    fprintf('\n\n================================================================\n');
    fprintf('--- Processing Analysis %d: Potential transmission from G%d (Sys Idx %d, SF=%.3f) to G%d (Actual Rank %d) for prefix %s ---\n', ...
            i_analysis, sg_id, sg_ss_idx, sf_sg_active, tg_id, tg_ss_rank, dic_prefix);
    fprintf('================================================================\n');

    % --- Calculate primary m' and define the "Residual Mismatch Tensor" ---
    b_sg_eff = sS_sg_active.b * slip_sign_sg;
    b_tg_eff = sS_tg_active.b * slip_sign_tg;
    m_prime_active = abs(cos(angle(b_sg_eff, b_tg_eff, 'noSymmetry')) * cos(angle(sS_sg_active.n, sS_tg_active.n, 'noSymmetry')));

    % Tensors for simple shear: M = b (x) n
    M_sg_eff = (b_sg_eff.normalize.xyz)' * sS_sg_active.n.normalize.xyz;
    M_tg_active_eff = (b_tg_eff.normalize.xyz)' * sS_tg_active.n.normalize.xyz;
    
    % Boundary mismatch left over after the primary active system accommodates its share
    H_residual = M_sg_eff - m_prime_active * M_tg_active_eff;

    % Get all possible slip systems for the target grain (Grain 2)
    sS_tg_all = ori_tg * sSLocalRef;

    % --- Cluster slip systems by plane and burger vector families ---
    cs = ebsdID.CS;
    num_systems_ref = length(sSLocalRef);
    plane_families = {Miller(1,1,0,cs,'hkl'), Miller(1,1,2,cs,'hkl'), Miller(1,2,3,cs,'hkl')};
    unique_burgers_ref = unique(sSLocalRef.b,'noSymmetry','antipodal');

    plane_family_id = zeros(num_systems_ref, 1);
    burger_family_id = zeros(num_systems_ref, 1);

    for j = 1:num_systems_ref
        current_n = sSLocalRef(j).n;
        current_b = sSLocalRef(j).b;
        
        % Assign plane family ID
        if angle(current_n, plane_families{1}) < 1e-5
            plane_family_id(j) = 1; % {110}
        elseif angle(current_n, plane_families{2}) < 1e-5
            plane_family_id(j) = 2; % {112}
        elseif angle(current_n, plane_families{3}) < 1e-5
            plane_family_id(j) = 3; % {123}
        end

        % Assign burger family ID
        for b_idx = 1:length(unique_burgers_ref)
            if angle(current_b, unique_burgers_ref(b_idx), 'noSymmetry','antipodal') < 1e-5
                burger_family_id(j) = b_idx;
                break;
            end
        end
    end

    % --- Calculate parameters for all potential systems in Grain 2 ---
    num_systems = length(sS_tg_all);
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
        sS_tg_current = sS_tg_all(j);
        
        % Luster-Morris parameter (m')
        % Compare effective incoming shear to the potential target shear.
        m_prime_vals(j) = abs(cos(angle(b_sg_eff, sS_tg_current.b, 'noSymmetry')) * cos(angle(sS_sg_active.n, sS_tg_current.n, 'noSymmetry')));
      
        % Schmid Factor (for uniaxial tension along y-vector)
        schmid_factors(j) = abs(sS_tg_current.SchmidFactor(stress));
        
        % Combined Transmission Factor
        combined_metric(j) = m_prime_vals(j) * schmid_factors(j);
        
        % How well does this system geometrically accommodate the remaining boundary mismatch? (Tensor Double Dot Product)
        Mj = sS_tg_current.b.normalize.xyz' * sS_tg_current.n.normalize.xyz;
        residual_accom(j) = sum(H_residual(:) .* Mj(:));
        
        % Out-of-plane component of Burgers vector (dot product with z-vector)
        outofplane_comps(j) = abs(dot(sS_tg_current.b.normalize, zvector));
        
        % Global Activity Rank of the system
        activity_ranks(j) = find(sorted_sslip_idx == j, 1);
    end

    % --- Create, Sort, and Display the Results Table ---
    results_table = table(sys_indices, activity_ranks, b_vectors_ref, n_vectors_ref, m_prime_vals, schmid_factors, combined_metric, residual_accom, outofplane_comps, plane_family_id, burger_family_id, ...
        'VariableNames', {'Sys_Idx', 'Activity_Rank', 'b_ref', 'n_ref', 'm_prime', 'SchmidFactor_Y', 'm_prime_x_SF', 'Res_Accom', 'b_z_component', 'plane_family_id', 'burger_family_id'});

    sorted_table = sortrows(results_table, 'm_prime_x_SF', 'descend');

    fprintf('\n--- Potential Transmission Systems in Grain %d (sorted by m'' x SF) ---\n', tg_id);
    disp(sorted_table);

    % Save the table to a CSV file
    csv_filename = fullfile(outDir, sprintf('PotentialTransmission_%s_G%d-Sys%d_to_G%d.csv', dic_prefix, sg_id, sg_ss_idx, tg_id));
    writetable(sorted_table, csv_filename);
    fprintf('Saved potential transmission data to: %s\n', csv_filename);

    % --- Scatter Plot Visualization (Clustered by system type) ---
    fig_scatter = figure;
    fig_scatter.Name = 'Transmission Likelihood Scatter';
    hold on;
    grid on;

    % Define markers and colors
    plane_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; % Blue, Orange, Yellow for {110}, {112}, {123}
    burger_markers = {'o', 's', 'd', '^'}; % For the 4 <111> directions

    % Plot non-active systems first
    for k = 1:height(sorted_table)
        sys_idx = sorted_table.Activity_Rank(k);
        % if sys_idx == tg_ss_idx, continue; end % Skip the active one for now

        m_prime = sorted_table.m_prime(k);
        sf = sorted_table.SchmidFactor_Y(k);
        
        plane_id = sorted_table.plane_family_id(k);
        burger_id = sorted_table.burger_family_id(k);
        
        text(sorted_table.m_prime(k) + 0.015, sorted_table.SchmidFactor_Y(k), num2str(sorted_table.Sys_Idx(k)), 'Color', 'k', 'FontSize', 8, 'BackgroundColor', 'none', 'Margin', 1);

        if plane_id > 0 && burger_id > 0
            scatter(m_prime, sf, 80, burger_markers{burger_id}, ...
                'MarkerFaceColor', plane_colors(plane_id, :), ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8);
        end
    end

    % % Now plot the active system on top
    % idx_active = find(sorted_table.Sys_Idx == tg_ss_idx);
    % if ~isempty(idx_active)
    %     m_prime = sorted_table.m_prime(idx_active);
    %     sf = sorted_table.SchmidFactor_Y(idx_active);
    % 
    %     plane_id = sorted_table.plane_family_id(idx_active);
    % 
    %     if plane_id > 0
    %         scatter(m_prime, sf, 300, 'p', ...
    %             'MarkerFaceColor', plane_colors(plane_id, :), ...
    %             'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    %     end
    % end

    % % Add labels
    % for k = 1:height(sorted_table)
    %     text(sorted_table.m_prime(k) + 0.015, sorted_table.SchmidFactor_Y(k), num2str(sorted_table.Sys_Idx(k)), ...
    %         'Color', 'k', 'FontSize', 8, 'BackgroundColor', 'none', 'Margin', 1);
    % end

    hold off;
    xlim([0, 1]);
    ylim([0, 0.5]);
    xlabel('Luster-Morris Parameter (m'')');
    ylabel('Target Schmid Factor (SF)');
    title(sprintf('Transmission Likelihood: G%d (Sys %d) -> G%d', sg_id, sg_ss_idx, tg_id));
    
    scatter_path = fullfile(outDir, sprintf('PotentialTransmission_Scatter_%s_G%d-Sys%d_to_G%d.jpg', dic_prefix, sg_id, sg_ss_idx, tg_id));
    exportScaledFigure(fig_scatter, scatter_path, 'target_px', 2500, 'noDateString');
    fprintf('Saved likelihood heatmap to: %s\n', scatter_path);

    % Create a custom legend in a new figure
    fig_legend = figure('Visible', 'off');
    ax_leg = axes(fig_legend);
    hold(ax_leg, 'on');
    legend_handles = [];
    legend_labels = {};
    plane_family_names = {'\{110\}', '\{112\}', '\{123\}'};
    for p_idx = 1:size(plane_colors, 1)
        h = plot(ax_leg, nan, nan, 's', 'MarkerFaceColor', plane_colors(p_idx,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
        legend_handles = [legend_handles, h];
        legend_labels = [legend_labels, sprintf('Plane %s', plane_family_names{p_idx})];
    end
    for b_idx = 1:length(unique_burgers_ref)
        h = plot(ax_leg, nan, nan, burger_markers{b_idx}, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
        legend_handles = [legend_handles, h];
        label_str = sprintf('b %s', char(unique_burgers_ref(b_idx).round));
        legend_labels = [legend_labels, label_str];
    end
    h = plot(ax_leg, nan, nan, 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5);
    legend_handles = [legend_handles, h];
    legend_labels = [legend_labels, 'Active System'];
    legend(legend_handles, legend_labels, 'Location', 'best');
    axis(ax_leg, 'off');
    legend_path = fullfile(outDir, sprintf('PotentialTransmission_Legend_%s_G%d-Sys%d_to_G%d.jpg', dic_prefix, sg_id, sg_ss_idx, tg_id));
    exportScaledFigure(fig_legend, legend_path, 'target_px', 800, 'noDateString');
    close(fig_legend);

    close(fig_scatter);

    % --- Visualization ---
    % Plot the source system and the actually active target system
    
    % Create a combined activity map for visualization of the source activity
    activity_map = nan(size(ebsdID));
    in_sg = (ebsdID.grainId == sg_id);
    slip_data_sg = abs(ebsdID.prop.slipIDcor(sg_ss_idx, :));
    activity_map(in_sg) = slip_data_sg(in_sg);
    
    fig_potential = figure;
    plot(ebsdID, activity_map, 'micronbar', 'off'); 
    mtexColorMap viridis;
    hold on;

    plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
    text(grainsUndef, grainsUndef.id, 'color', 'w', 'fontsize', 14);

    % Plot slip traces (red) and Burgers vectors (blue) for source and actual active target
    quiver(grainsUndef(sg_id), sS_sg_active.trace, 'color', 'r', 'linewidth', 3);
    quiver(grainsUndef(tg_id), sS_tg_active.trace, 'color', 'r', 'linewidth', 3, 'linestyle', '--');
    
    b_sg = sS_sg_active.b.normalize;
    b_sg.antipodal  =true;
    b_tg = sS_tg_active.b.normalize;
    b_tg.antipodal  =true;
    quiver(grainsUndef(sg_id), b_sg, 'color', 'b', 'linewidth', 3);
    quiver(grainsUndef(tg_id), b_tg, 'color', 'b', 'linewidth', 3, 'linestyle', '--');

    % Create text for the info box
    idx_active = find(sorted_table.Sys_Idx == tg_ss_idx);
    info_str = sprintf(['Analysis for %s:\n' ...
        'Source: G%d, Active Sys Idx %d, SF = %.3f\n  b=[%s], n=(%s)\n\n' ...
        'Actual Active Target in G%d (Rank %d):\n  Sys %d with m''xSF = %.3f\n  m'' = %.3f, SF = %.3f\n  b=[%s], n=(%s)'], ...
        dic_prefix, sg_id, sg_ss_idx, sf_sg_active, ...
        char(sSLocalRef(sg_ss_idx).b.round), char(sSLocalRef(sg_ss_idx).n.round), ...
        tg_id, tg_ss_rank, ...
        tg_ss_idx, sorted_table.m_prime_x_SF(idx_active), sorted_table.m_prime(idx_active), sorted_table.SchmidFactor_Y(idx_active), ...
        char(sSLocalRef(tg_ss_idx).b.round), char(sSLocalRef(tg_ss_idx).n.round));
        
    annotation('textbox', [0.02 0.7 0.3 0.25], 'String', info_str, ...
        'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'FitBoxToText', 'on', 'FontSize', 10, 'Margin', 5);
        
    title(sprintf('Potential Transmission from G%d to G%d', sg_id, tg_id));
    
    % Export the figure
    fig_path = fullfile(outDir, sprintf('PotentialTransmission_Plot_%s_G%d-Sys%d_to_G%d.jpg', dic_prefix, sg_id, sg_ss_idx, tg_id));
    exportScaledFigure(fig_potential, fig_path, 'target_px', 3500, 'noDateString');
    fprintf('Saved summary plot to: %s\n', fig_path);
    
    close(fig_potential);
end