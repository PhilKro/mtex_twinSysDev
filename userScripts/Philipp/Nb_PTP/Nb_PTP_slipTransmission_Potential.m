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
analyses(1).sg_ss_idx  = 27; % Index of the active slip system in the source grain
analyses(1).tg_id      = 1; % Target grain ID
analyses(1).tg_ss_rank = 18; % Actual active target system rank for comparison

analyses(2).dic_prefix = 'B0005';
analyses(2).sg_id      = 2;
analyses(2).sg_ss_idx  = 11;
analyses(2).tg_id      = 1; 
analyses(2).tg_ss_rank = 18;

analyses(3).dic_prefix = 'B0005';
analyses(3).sg_id      = 2;
analyses(3).sg_ss_idx  = 47;
analyses(3).tg_id      = 1; 
analyses(3).tg_ss_rank = 18;

analyses(4).dic_prefix = 'B0005';
analyses(4).sg_id      = 2;
analyses(4).sg_ss_idx  = [36,24];
analyses(4).tg_id      = 1; 
analyses(4).tg_ss_rank = 18;

% analyses(5).dic_prefix = 'B0005';
% analyses(5).sg_id      = 2;
% analyses(5).sg_ss_idx  = 24;
% analyses(5).tg_id      = 1; 
% analyses(5).tg_ss_rank = 18;

% Example: combine two source systems into one analysis step
% analyses(6).dic_prefix = 'B0005';
% analyses(6).sg_id      = 2;
% analyses(6).sg_ss_idx  = [24, 36];
% analyses(6).tg_id      = 1;
% analyses(6).tg_ss_rank = 18;
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
    
    if isnumeric(sg_ss_idx)
        sg_ss_idx = sg_ss_idx(:)';
    end

    num_sg_systems = numel(sg_ss_idx);
    if num_sg_systems > 1
        sg_ss_label = strjoin(arrayfun(@num2str, sg_ss_idx, 'UniformOutput', false), '-');
    else
        sg_ss_label = num2str(sg_ss_idx);
    end

    csv_filename = fullfile(outDir, sprintf('PotentialTransmission_%s_G%d-Sys%s_to_G%d.csv', dic_prefix, sg_id, sg_ss_label, tg_id));

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
    
    sS_sg_actives = ori_sg * sSLocalRef(sg_ss_idx);
    sS_tg_active = ori_tg * sSLocalRef(tg_ss_idx);
    
    sS_tg_all = ori_tg * sSLocalRef;

    % Check if CSV exists and load if possible
    if exist(csv_filename, 'file')
        fprintf('Loading existing potential transmission data from %s...\n', csv_filename);
        sorted_table = readtable(csv_filename);
        
        % Ensure the table has the necessary columns for m_prime ranges
        if ~ismember('m_prime_min', sorted_table.Properties.VariableNames)
            fprintf('Recomputing m_prime ranges for loaded data...\n');
            sorted_table.m_prime_min = nan(height(sorted_table), 1);
            sorted_table.m_prime_max = nan(height(sorted_table), 1);
            sorted_table.m_prime_mid = nan(height(sorted_table), 1);
            sorted_table.m_prime_avg = nan(height(sorted_table), 1);
            
            for j = 1:height(sorted_table)
                sys_idx = sorted_table.Sys_Idx(j);
                sS_tg_current = sS_tg_all(sys_idx);
                m_prime_source = nan(num_sg_systems, 1);
                
                for s = 1:num_sg_systems
                    sS_k = sS_sg_actives(s);
                    slip_sign = sign(mean(ebsdID.prop.slipIDcor(sg_ss_idx(s), ebsdID.grainId == sg_id), 'omitnan'));
                    if slip_sign == 0, slip_sign = 1; end
                    b_eff_k = sS_k.b * slip_sign;
                    m_prime_source(s) = abs(cos(angle(b_eff_k, sS_tg_current.b, 'noSymmetry')) * ...
                        cos(angle(sS_k.n, sS_tg_current.n, 'noSymmetry')));
                end
                
                sorted_table.m_prime_min(j) = min(m_prime_source);
                sorted_table.m_prime_max(j) = max(m_prime_source);
                sorted_table.m_prime_mid(j) = 0.5 * (min(m_prime_source) + max(m_prime_source));
                sorted_table.m_prime_avg(j) = mean(m_prime_source);
            end
        end
    else
        fprintf('Calculating potential transmission data...\n');
        
        sg_activity_by_system = zeros(num_sg_systems, 1);
        slip_sign_sg = zeros(num_sg_systems, 1);
        M_sg_eff = zeros(3, 3);
        m_prime_active_per_source = nan(num_sg_systems, 1);

        for k = 1:num_sg_systems
            idx = sg_ss_idx(k);
            slip_values = ebsdID.prop.slipIDcor(idx, ebsdID.grainId == sg_id);
            sg_activity_by_system(k) = mean(abs(slip_values), 'omitnan');
            slip_sign_sg(k) = sign(mean(slip_values, 'omitnan'));
            if slip_sign_sg(k) == 0
                slip_sign_sg(k) = 1;
            end

            sS_k = sS_sg_actives(k);
            b_eff_k = sS_k.b * slip_sign_sg(k);
            n_eff_k = sS_k.n;
            M_k = b_eff_k.normalize.xyz' * n_eff_k.normalize.xyz;
            M_sg_eff = M_sg_eff + M_k * sg_activity_by_system(k);

            % Primary m' for the selected target active system
            m_prime_active_per_source(k) = abs(cos(angle(b_eff_k, sS_tg_active.b, 'noSymmetry')) * ...
                cos(angle(n_eff_k, sS_tg_active.n, 'noSymmetry')));
        end

        if sum(sg_activity_by_system) > 0
            M_sg_eff = M_sg_eff / sum(sg_activity_by_system);
        end

        % --- Calculate primary m' and define the "Residual Mismatch Tensor" ---
        b_tg_eff = sS_tg_active.b * sign(mean(ebsdID.prop.slipIDcor(tg_ss_idx, ebsdID.grainId == tg_id), 'omitnan'));
        if sign(mean(ebsdID.prop.slipIDcor(tg_ss_idx, ebsdID.grainId == tg_id), 'omitnan')) == 0
            b_tg_eff = sS_tg_active.b;
        end
        M_tg_active_eff = (b_tg_eff.normalize.xyz)' * sS_tg_active.n.normalize.xyz;
        
        % Boundary mismatch left over after the primary active system accommodates its share
        H_residual = M_sg_eff - m_prime_active * M_tg_active_eff;

        % Get all possible slip systems for the target grain (Grain 2)
        % sS_tg_all = ori_tg * sSLocalRef;

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
        m_prime_min = nan(num_systems, 1);
        m_prime_max = nan(num_systems, 1);
        m_prime_mid = nan(num_systems, 1);
        m_prime_avg = nan(num_systems, 1);
        schmid_factors = nan(num_systems, 1);
        activity_ranks = nan(num_systems, 1);
        combined_metric = nan(num_systems, 1);
        residual_accom = nan(num_systems, 1);
        outofplane_comps = nan(num_systems, 1);

        for j = 1:num_systems
            sS_tg_current = sS_tg_all(j);
            m_prime_source = nan(num_sg_systems, 1);

            for k = 1:num_sg_systems
                sS_k = sS_sg_actives(k);
                b_eff_k = sS_k.b * slip_sign_sg(k);
                m_prime_source(k) = abs(cos(angle(b_eff_k, sS_tg_current.b, 'noSymmetry')) * ...
                    cos(angle(sS_k.n, sS_tg_current.n, 'noSymmetry')));
            end

            m_prime_min(j) = min(m_prime_source);
            m_prime_max(j) = max(m_prime_source);
            m_prime_mid(j) = 0.5 * (m_prime_min(j) + m_prime_max(j));
            m_prime_avg(j) = mean(m_prime_source);

            % Schmid Factor (for uniaxial tension along y-vector)
            schmid_factors(j) = abs(sS_tg_current.SchmidFactor(stress));
            
            % Combined Transmission Factor
            combined_metric(j) = m_prime_mid(j) * schmid_factors(j);
            
            % How well does this system geometrically accommodate the remaining boundary mismatch? (Tensor Double Dot Product)
            Mj = sS_tg_current.b.normalize.xyz' * sS_tg_current.n.normalize.xyz;
            residual_accom(j) = sum(H_residual(:) .* Mj(:));
            
            % Out-of-plane component of Burgers vector (dot product with z-vector)
            outofplane_comps(j) = abs(dot(sS_tg_current.b.normalize, zvector));
            
            % Global Activity Rank of the system
            activity_ranks(j) = find(sorted_sslip_idx == j, 1);
        end

        % --- Create, Sort, and Display the Results Table ---
        results_table = table(sys_indices, activity_ranks, b_vectors_ref, n_vectors_ref, m_prime_min, m_prime_max, m_prime_mid, m_prime_avg, schmid_factors, combined_metric, residual_accom, outofplane_comps, plane_family_id, burger_family_id, ...
            'VariableNames', {'Sys_Idx', 'Activity_Rank', 'b_ref', 'n_ref', 'm_prime_min', 'm_prime_max', 'm_prime', 'm_prime_avg', 'SchmidFactor_Y', 'm_prime_x_SF', 'Res_Accom', 'b_z_component', 'plane_family_id', 'burger_family_id'});

        sorted_table = sortrows(results_table, 'm_prime_x_SF', 'descend');

        fprintf('\n--- Potential Transmission Systems in Grain %d (sorted by m'' x SF) ---\n', tg_id);
        disp(sorted_table);

        % Save the table to a CSV file
        writetable(sorted_table, csv_filename);
        fprintf('Saved potential transmission data to: %s\n', csv_filename);
    end

    % Compute common variables needed for plotting
    stress = stressTensor.uniaxial(yvector);
    sf_sg_active = mean(abs([sS_sg_actives.SchmidFactor(stress)]));
    m_prime_active_per_source = nan(num_sg_systems, 1);
    for k = 1:num_sg_systems
        sS_k = sS_sg_actives(k);
        slip_sign = sign(mean(ebsdID.prop.slipIDcor(sg_ss_idx(k), ebsdID.grainId == sg_id), 'omitnan'));
        if slip_sign == 0, slip_sign = 1; end
        b_eff_k = sS_k.b * slip_sign;
        m_prime_active_per_source(k) = abs(cos(angle(b_eff_k, sS_tg_active.b, 'noSymmetry')) * ...
            cos(angle(sS_k.n, sS_tg_active.n, 'noSymmetry')));
    end
    m_prime_active = mean(m_prime_active_per_source);

    fprintf('\n\n================================================================\n');
    fprintf('--- Processing Analysis %d: Potential transmission from G%d (Sys Idx %s, SF=%.3f) to G%d (Actual Rank %d) for prefix %s ---\n', ...
            i_analysis, sg_id, sg_ss_label, sf_sg_active, tg_id, tg_ss_rank, dic_prefix);
    fprintf('================================================================\n');



    % --- Scatter Plot Visualization (Clustered by system type) ---
    fig_scatter = figure;
    fig_scatter.Name = 'Transmission Likelihood Scatter';
    hold on;
    grid on;

    % Define colors and markers for source systems
    source_colors = lines(num_sg_systems);  % Use MATLAB's lines colormap for distinct colors
    source_markers = {'o', 's', 'd', '^', 'v', '>', '<'};  % Cycle through markers if more than 7 sources
    if num_sg_systems > length(source_markers)
        source_markers = repmat(source_markers, 1, ceil(num_sg_systems / length(source_markers)));
    end
    source_markers = source_markers(1:num_sg_systems);

    % Plot each target system with points for each source system
    legend_handles = [];
    legend_labels = {};
    for k = 1:height(sorted_table)
        sys_idx = sorted_table.Sys_Idx(k);
        sS_tg_current = sS_tg_all(sys_idx);
        sf = sorted_table.SchmidFactor_Y(k);
        
        for s = 1:num_sg_systems
            sS_k = sS_sg_actives(s);
            slip_sign = sign(mean(ebsdID.prop.slipIDcor(sg_ss_idx(s), ebsdID.grainId == sg_id), 'omitnan'));
            if slip_sign == 0, slip_sign = 1; end
            b_eff_k = sS_k.b * slip_sign;
            m_prime_s = abs(cos(angle(b_eff_k, sS_tg_current.b, 'noSymmetry')) * ...
                cos(angle(sS_k.n, sS_tg_current.n, 'noSymmetry')));
            
            h = scatter(m_prime_s, sf, 100, source_markers{s}, ...
                'MarkerFaceColor', source_colors(s, :), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
            
            if k == 1  % Collect handles only once per source
                legend_handles = [legend_handles, h];
                legend_labels = [legend_labels, sprintf('Source Sys %d', sg_ss_idx(s))];
            end
        end
    end

    % Add legend
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'best');
    end

    hold off;
    xlim([0, 1]);
    ylim([0, 0.5]);
    xlabel('Luster-Morris Parameter (m'')');
    ylabel('Target Schmid Factor (SF)');
    title(sprintf('Transmission Likelihood: G%d (Sys %s) -> G%d', sg_id, sg_ss_label, tg_id));
    
    scatter_path = fullfile(outDir, sprintf('Nodiff_PotentialTransmission_Scatter_%s_G%d-Sys%s_to_G%d.jpg', dic_prefix, sg_id, sg_ss_label, tg_id));
    exportScaledFigure(fig_scatter, scatter_path, 'target_px', 2500, 'noDateString');
    fprintf('Saved likelihood heatmap to: %s\n', scatter_path);
    
    scatter_annot_path = fullfile(outDir, sprintf('PotentialTransmission_Scatter_Annotated_%s_G%d-Sys%s_to_G%d.jpg', dic_prefix, sg_id, sg_ss_label, tg_id));
    exportScaledFigure(fig_scatter, scatter_annot_path, 'target_px', 2500, 'noDateString');
    fprintf('Saved annotated likelihood heatmap to: %s\n', scatter_annot_path);

    close(fig_scatter);

    % --- Visualization ---
    % Plot the source system and the actually active target system
    
    % --- Load and apply specimen mask for background ---
    mask_file_undef = fullfile(baseDir, 'insitu Specimen 1 0 Map Data 1_BW263_specimen_mask.mat');
    mask_aligned = true(size(ebsdID)); % Default all true if no mask found
    if exist(mask_file_undef, 'file')
        load(mask_file_undef, 'specimenlocation');
        if isfield(ebsdID.prop, 'oldId')
            old_id = ebsdID.prop.oldId;
            valid_id = ~isnan(old_id);
            mask_aligned = false(size(old_id));
            mask_aligned(valid_id) = specimenlocation(old_id(valid_id));
        end
    end

    % Create a combined activity map for visualization of the source activity
    activity_map = nan(size(ebsdID));
    in_sg = (ebsdID.grainId == sg_id);
    slip_data_sg = sum(abs(ebsdID.prop.slipIDcor(sg_ss_idx, :)), 1);
    activity_map(in_sg) = slip_data_sg(in_sg);
    
    fig_potential = figure;
    % Plot black background for the whole specimen
    plot(ebsdID, zeros(length(ebsdID), 3), 'micronbar', 'off');
    hold on;
    
    % Plot SSLIP activity for the source grain
    plot(ebsdID, activity_map, 'micronbar', 'off'); 
    mtexColorMap viridis;
    setColorRange([0,0.15])

    plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
    %text(grainsUndef, grainsUndef.id, 'color', 'w', 'fontsize', 14);

    % Plot slip traces (red) and Burgers vectors (blue) ONLY for source grain
    for k = 1:num_sg_systems
        quiver(grainsUndef(sg_id), sS_sg_actives(k).trace, 'color', 'r', 'linewidth', 3);
    end

    for k = 1:num_sg_systems
        b_sg = sS_sg_actives(k).b.normalize;
        b_sg.antipodal = true;
        quiver(grainsUndef(sg_id), b_sg, 'color', 'b', 'linewidth', 3);
    end

    % Create text for the info box
    idx_active = find(sorted_table.Sys_Idx == tg_ss_idx);
    source_b_list = strjoin(arrayfun(@(x) char(sSLocalRef(x).b.round), sg_ss_idx, 'UniformOutput', false), ', ');
    source_n_list = strjoin(arrayfun(@(x) char(sSLocalRef(x).n.round), sg_ss_idx, 'UniformOutput', false), ', ');
    info_str = sprintf(['Analysis for %s:\n' ...
        'Source: G%d, Active Sys Idx %s, SF = %.3f\n  b=[%s], n=(%s)\n\n' ...
        'Actual Active Target in G%d (Rank %d):\n  Sys %d with m''xSF = %.3f\n  m'' = %.3f, SF = %.3f\n  b=[%s], n=(%s)'], ...
        dic_prefix, sg_id, sg_ss_label, sf_sg_active, ...
        source_b_list, source_n_list, ...
        tg_id, tg_ss_rank, ...
        tg_ss_idx, sorted_table.m_prime_x_SF(idx_active), sorted_table.m_prime(idx_active), sorted_table.SchmidFactor_Y(idx_active), ...
        char(sSLocalRef(tg_ss_idx).b.round), char(sSLocalRef(tg_ss_idx).n.round));

    annotation('textbox', [0.02 0.02 0.3 0.25], 'String', info_str, ...
        'Interpreter', 'none', 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', ...
        'FitBoxToText', 'on', 'FontSize', 10, 'Margin', 5, 'VerticalAlignment', 'bottom');
        
    title(sprintf('Potential Transmission from G%d (Sys %s) to G%d', sg_id, sg_ss_label, tg_id));
    
    % Export the figure
    fig_path = fullfile(outDir, sprintf('PotentialTransmission_Plot_%s_G%d-Sys%s_to_G%d.jpg', dic_prefix, sg_id, sg_ss_label, tg_id));
    exportScaledFigure(fig_potential, fig_path, 'target_px', 3500, 'noDateString');
    fprintf('Saved summary plot to: %s\n', fig_path);
    
    close(fig_potential);
end