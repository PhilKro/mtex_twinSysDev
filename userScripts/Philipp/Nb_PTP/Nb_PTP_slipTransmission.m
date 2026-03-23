%% Nb_PTP_slipTransmission.m
% Script to calculate slip transmission parameters (m') 
% between manually selected slip systems in neighboring grains.

clear; close all; clc;

%% 1. Initialization & Mtex Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Data Loading
baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
outDir  = fullfile(baseDir, 'Processed_Output');

%% 3. Input Structs Setup
% Define a list of configurations to analyze
analyses = struct();

analyses(1).dic_prefix = 'B0003';
analyses(1).grain1_id  = 2;
analyses(1).ss1_rank   = 1;
analyses(1).grain2_id  = 1;
analyses(1).ss2_rank   = 2;

analyses(2).dic_prefix = 'B0005';
analyses(2).grain1_id  = 2;
analyses(2).ss1_rank   = 2;
analyses(2).grain2_id  = 1;
analyses(2).ss2_rank   = 1;

analyses(3).dic_prefix = 'B0005';
analyses(3).grain1_id  = 2;
analyses(3).ss1_rank   = 4;
analyses(3).grain2_id  = 1;
analyses(3).ss2_rank   = 1;

analyses(4).dic_prefix = 'B0005';
analyses(4).grain1_id  = 2;
analyses(4).ss1_rank   = 5;
analyses(4).grain2_id  = 1;
analyses(4).ss2_rank   = 1;
%% 4. Analysis loop
for i_analysis = 1:length(analyses)
    current_analysis = analyses(i_analysis);
    dic_prefix = current_analysis.dic_prefix;
    grain1_id  = current_analysis.grain1_id;
    ss1_rank   = current_analysis.ss1_rank;
    grain2_id  = current_analysis.grain2_id;
    ss2_rank   = current_analysis.ss2_rank;
    
    fprintf('\n--- Processing Analysis %d: %s ---\n', i_analysis, dic_prefix);

    % Load SSLIP data (which contains the reconstructed undeformed ebsd map and the reference slip systems)
    sslip_mat_path = fullfile(outDir, sprintf('SSLIP_Data_%s.mat', dic_prefix));
    if exist(sslip_mat_path, 'file')
        load(sslip_mat_path, 'ebsdID', 'sSLocalRef');
    else
        warning('SSLIP data not found for prefix %s. Skipping.', dic_prefix);
        continue;
    end

    % Reconstruct grains from ebsdID to get the mean orientations for the boundary calculations
    grainThresh = 5 * degree;
    minSize = 1000;
    [grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);
    ebsdID(grainsUndef(grainsUndef.numPixel<=minSize)) = 'notIndexed';
    [grainsUndef, ebsdID.grainId] = calcGrains(ebsdID, 'minPixel', minSize, 'angle', grainThresh);

    % 4. Rank Slip Systems by Activity
    % Rank SSLIP systems globally across the map (as done in fullAnalysis)
    sslip_activity_sum = sum(abs(ebsdID.prop.slipIDcor), 2, 'omitnan');
    [~, sorted_sslip_idx] = sort(sslip_activity_sum, 'descend');

    % Get the actual indices based on the specified ranks
    ss1_idx = sorted_sslip_idx(ss1_rank);
    ss2_idx = sorted_sslip_idx(ss2_rank);

    % 5. Initial Plot and Validation
    % Retrieve the mean orientations for the two grains
    ori1 = grainsUndef(grain1_id).meanOrientation;
    ori2 = grainsUndef(grain2_id).meanOrientation;

    % Rotate the reference slip systems to the local grain coordinates
    sSLocal_1 = ori1 * sSLocalRef(ss1_idx);
    sSLocal_2 = ori2 * sSLocalRef(ss2_idx);

    % Extract string representations of the slip systems for the title
    b1_str = char(sSLocalRef(ss1_idx).b.round);
    n1_str = char(sSLocalRef(ss1_idx).n.round);
    b2_str = char(sSLocalRef(ss2_idx).b.round);
    n2_str = char(sSLocalRef(ss2_idx).n.round);

    % Create a combined activity map for visualization
    activity_map = nan(size(ebsdID));
    in_grain1 = (ebsdID.grainId == grain1_id);
    in_grain2 = (ebsdID.grainId == grain2_id);
    
    slip_data_1 = abs(ebsdID.prop.slipIDcor(ss1_idx, :));
    slip_data_2 = abs(ebsdID.prop.slipIDcor(ss2_idx, :));
    
    activity_map(in_grain1) = slip_data_1(in_grain1);
    activity_map(in_grain2) = slip_data_2(in_grain2);

    % Plot the EBSD map, boundary, and the specified slip systems
    fig_check = figure;
    plot(ebsdID, activity_map, 'micronbar', 'off'); 
    mtexColorMap viridis;
    hold on;

    plot(grainsUndef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
    text(grainsUndef, grainsUndef.id, 'color', 'w', 'fontsize', 14);

    % Plot slip traces (red) and slip directions/Burgers vectors (white)
    quiver(grainsUndef(grain1_id), sSLocal_1.trace, 'color', 'r', 'linewidth', 3);
    quiver(grainsUndef(grain2_id), sSLocal_2.trace, 'color', 'r', 'linewidth', 3);
    warning('Somehow make the bs vector3ds antipodal otherwise they dont plot.')
    quiver(grainsUndef(grain1_id), sSLocal_1.b.normalize, 'color', 'b', 'linewidth', 3);
    quiver(grainsUndef(grain2_id), sSLocal_2.b.normalize, 'color', 'b', 'linewidth', 3);

    % Create text for the info box instead of title to prevent cut-offs
    info_str = sprintf('Prefix %s\nG%d (Rank %d): b=[%s] n=(%s)\nG%d (Rank %d): b=[%s] n=(%s)', ...
        dic_prefix, grain1_id, ss1_rank, b1_str, n1_str, grain2_id, ss2_rank, b2_str, n2_str);
        
    % Add textbox inside the figure (top left corner)
    tb = annotation('textbox', [0.02 0.75 0.3 0.2], 'String', info_str, ...
        'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'FitBoxToText', 'on', 'FontSize', 12, 'Margin', 5);
        
    % Prompt user
    choice = questdlg('Does this selection look correct? (Check figure)', 'Confirm Selection', 'Yes', 'No, skip', 'Yes');
    if strcmp(choice, 'No, skip')
        close(fig_check);
        fprintf('User skipped analysis %d.\n', i_analysis);
        continue;
    end

    % 6. Calculate Slip Transmission Parameters
    % Calculate the Luster-Morris parameter (m')
    % m' = |cos(angle between b1 and b2) * cos(angle between n1 and n2)|
    m_prime_val = abs(cos(angle(sSLocal_1.b, sSLocal_2.b, 'noSymmetry')) * cos(angle(sSLocal_1.n, sSLocal_2.n, 'noSymmetry')));

    % To calculate the LRB (Lee-Robertson-Birnbaum) criterion, we estimate the boundary normal
    gB = grainsUndef.boundary(num2str(grain1_id), num2str(grain2_id));

    if isempty(gB)
        warning('No common boundary found between Grain %d and Grain %d.', grain1_id, grain2_id);
        avg_gB_Dir = vector3d.X; % Fallback direction
    else
        % Smooth the boundary for a better, more representative normal estimation
        smoothGrains = smooth(grainsUndef, 10);
        gB_smooth = smoothGrains.boundary(num2str(grain1_id), num2str(grain2_id));
        avg_gB_Dir = mean(gB_smooth.direction);
    end

    % The surface boundary normal is perpendicular to the boundary direction and the z-vector (out of plane)
    avg_gB_Normal = cross(avg_gB_Dir, zvector);

    % Calculate the LRB parameter
    lrb_val = abs(cos(angle(sSLocal_1.b, sSLocal_2.b)) * cos(angle(cross(avg_gB_Normal, sSLocal_1.n), cross(avg_gB_Normal, sSLocal_2.n))));

    % 7. Display and Plot Results
    fprintf('Luster-Morris Parameter (m''): %.3f\n', m_prime_val);
    fprintf('LRB Parameter:                %.3f\n\n', lrb_val);

    % Update the figure with final parameters
    figure(fig_check);
    if ~isempty(gB), plot(gB, 'linecolor', 'y', 'linewidth', 3); end

    % Update textbox with final calculation results
    final_info_str = sprintf('%s\nm'' = %.3f | LRB = %.3f', info_str, m_prime_val, lrb_val);
    tb.String = final_info_str;
    
    % Export the figure
    fig_path = fullfile(outDir, sprintf('SlipTransmission_%s_G%d-R%d_to_G%d-R%d.jpg', dic_prefix, grain1_id, ss1_rank, grain2_id, ss2_rank));
    exportScaledFigure(fig_check, fig_path, 'target_px', 3500,'noDateString');
    fprintf('Saved slip transmission figure to %s\n', fig_path);
        
end