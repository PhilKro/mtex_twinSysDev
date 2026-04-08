%% Nb_PTP_prettyIPF_loop.m
% Script to loop over all Nb .ang files and generate pretty IPF Y maps
% overlaid on a tunable greyscale Confidence Index (CI) background.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Tuning Parameters
if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir = fullfile(baseDir, 'Pretty_Plots');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% --- "PHOTOSHOP CURVES" TUNING PARAMETERS ---
% 1. Black Point (0.0 to 1.0): CI values below this will be pure black.
blackPoint = 0.05; 

% 2. White Point (0.0 to 1.0): CI values above this will be pure white.
whitePoint = 0.98; 

% 3. Gamma (Midtones): 
%    < 1.0 makes the mid-range lighter (good if the map is too dark)
%    > 1.0 makes the mid-range darker (good if the map looks washed out)
%    = 1.0 applies no gamma curve.
gammaVal = 1.5; 

% 4. IPF Overlay Opacity (0.0 to 1.0)
ipf_alpha = 0.55;
% --------------------------------------------

% --- GND TUNING PARAMETERS ---
plot_gnds = true;          % Set to false to skip GND calculation and plotting
plot_ci_background = true; % Toggle CI background for GNDs
gnd_alpha = 0.65;          % Opacity of the GND overlay when CI is active
gnd_limits = [14.5, 16.5]; % Limits for log10(GND Density m^{-2})
% --------------------------------------------

% --- DIC STRAIN TUNING PARAMETERS ---
plot_dic_strain = false;    % Toggle DIC Strain calculation and plotting
% --------------------------------------------

% --- IPF TUNING PARAMETERS ---
plot_ipf = false;           % Toggle IPF map calculation and plotting
% --------------------------------------------

pf_dir = yvector;

%% 3. Load Data & Loop over .ang files
tkd_tilt = -10 * degree;
thresh_ci = 0.13;
CS = crystalSymmetry.load('Nb-Niobium.cif');

% --- Load Undeformed DIC Data (for DIC strain plotting) ---
if plot_dic_strain
    undeformed_dic_mat = fullfile(baseDir, 'Processed_Output', 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');
    if exist(undeformed_dic_mat, 'file')
        fprintf('Loading undeformed EBSD with DIC data...\n');
        undefData = load(undeformed_dic_mat, 'newEBSD');
        ebsdUndef = undefData.newEBSD;
        ebsdUndef = ebsdUndef.gridify;
        
        % Compute tuned CI for the undeformed EBSD
        ci_raw_undef = ebsdUndef.ci;
        ci_min_undef = min(ci_raw_undef(:));
        ci_max_undef = max(ci_raw_undef(:));
        ci_norm_undef = (ci_raw_undef - ci_min_undef) ./ (ci_max_undef - ci_min_undef);
        ci_tuned_undef = (ci_norm_undef - blackPoint) ./ (whitePoint - blackPoint);
        ci_tuned_undef(ci_tuned_undef < 0) = 0;
        ci_tuned_undef(ci_tuned_undef > 1) = 1;
        ci_tuned_undef = ci_tuned_undef .^ gammaVal;
        cibw_undef = [ci_tuned_undef(:), ci_tuned_undef(:), ci_tuned_undef(:)];
    else
        warning('Undeformed DIC MAT file not found at %s. DIC Strain plotting will be skipped.', undeformed_dic_mat);
        ebsdUndef = [];
        cibw_undef = [];
    end
end

% Find all .ang files in the directory
angFiles = dir(fullfile(baseDir, '*.ang'));

% angFiles = angFiles(1);

for i = 1:length(angFiles)
% for i = 8:8
    currentFile = angFiles(i).name;
    ang_path = fullfile(baseDir, currentFile);
    [~, baseName, ~] = fileparts(currentFile);
    
    fprintf('Processing %d/%d: %s\n', i, length(angFiles), currentFile);
    
    % Load EBSD Data and perform identical rotations as in your analysis scripts
    ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
    ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
    
    ebsd_1D = ebsd; % Keep 1D version for coordinate matching
    
    % --- Apply Specimen Mask to base EBSD ---
    mask_file = fullfile(baseDir, sprintf('%s_specimen_mask.mat', baseName));
    if exist(mask_file, 'file')
        load(mask_file, 'specimenlocation');
        ebsd(~specimenlocation) = 'notIndexed';
    end
    
    ebsd = ebsd.gridify;
    
    % --- Calculate Tuned CI Greyscale (Curves Adjustment) ---
    ci_raw = ebsd.ci;
    ci_min = min(ci_raw(:));
    ci_max = max(ci_raw(:));
    
    % Step 1: Strict normalization between 0 and 1
    ci_norm = (ci_raw - ci_min) ./ (ci_max - ci_min);
    
    % Step 2: Levels Clipping (Black and White points)
    ci_tuned = (ci_norm - blackPoint) ./ (whitePoint - blackPoint);
    ci_tuned(ci_tuned < 0) = 0; % clip darks
    ci_tuned(ci_tuned > 1) = 1; % clip highlights
    
    % Step 3: Gamma Correction (Midtone stretching/compressing)
    ci_tuned = ci_tuned .^ gammaVal;
    
    % Create Nx3 array for RGB plotting
    cibw = [ci_tuned(:), ci_tuned(:), ci_tuned(:)];
    
    % --- Plotting IPF ---
    if plot_ipf
        % ipf = ipfColorKey(ebsd.CS);
        ipf = ipfTSLKey(ebsd.CS);
        ipf.inversePoleFigureDirection = pf_dir;
        
        fig = figure('Name', currentFile, 'Visible', 'off'); % Set 'on' if you want to see them pop up
        plot(ebsd, cibw, 'micronbar', 'off');
        hold on;
        colors = ipf.orientation2color(ebsd('indexed').orientations);
        plot(ebsd('indexed'), colors, 'FaceAlpha', ipf_alpha);
        
        % Add a title and create a dynamic filename based on the IPF direction
        if isequal(pf_dir, xvector)
            dirStr = 'X';
        elseif isequal(pf_dir, yvector)
            dirStr = 'Y';
        elseif isequal(pf_dir, zvector)
            dirStr = 'Z';
        else
            dirStr = 'custom';
        end
        title(['IPF ', dirStr]);
        
        % Corrected saveFigure call - it saves the *current* active figure.
        % The first argument must be the filename string, not a figure handle.
        exportScaledFigure(gcf, fullfile(outDir, sprintf('%s_%sTSL_ipf_TunedIQ_nomic.jpg', baseName, dirStr)), 'noDateString');
        close;
    end
    
    % Map ang file to DIC prefix
    dic_prefix = '';
    switch currentFile
        case 'insitu Specimen 1 0 Map Data 1_BW263.ang', dic_prefix = 'B0001';
        case 'insitu Specimen 1 1 Map Data 2_BW263.ang', dic_prefix = 'B0003';
        case 'insitu Specimen 1 2 Map Data 3_BW263.ang', dic_prefix = 'B0005';
        case 'insitu Specimen 1 3 Map Data 4_BW263.ang', dic_prefix = 'B0007';
        case 'insitu Specimen 1 4 Map Data 5_BW263.ang', dic_prefix = 'B0009';
        case 'insitu Specimen 1 5 Map Data 6_BW263.ang', dic_prefix = 'B0011';
        case 'insitu Specimen 1 7 Map Data 7_BW263.ang', dic_prefix = 'B0014';
    end

    % --- Calculate and Plot DIC Strain ---
    if plot_dic_strain
        if isempty(dic_prefix)
            % Skip if not in the target list
        elseif ~exist('ebsdUndef', 'var') || isempty(ebsdUndef)
            warning('Undeformed DIC data not available. Skipping DIC Strain plot for %s.', dic_prefix);
        else
            u_field = [dic_prefix, '_Displacements_x_displacement'];
            v_field = [dic_prefix, '_Displacements_y_displacement'];
            
            if isfield(ebsdUndef.prop, u_field) && isfield(ebsdUndef.prop, v_field)
                fprintf('  Calculating and plotting DIC strain for %s...\n', dic_prefix);
                
                U = double(ebsdUndef.prop.(u_field));
                V = double(ebsdUndef.prop.(v_field));
                
                % DIC un-correlated points often export as exactly 0.
                % Replacing 0 with NaN prevents artificial "disconnected lines" in the gradient.
                U(U == 0) = NaN;
                V(V == 0) = NaN;
                
                % Apply the SSLIP Gaussian filter to smooth displacements before gradient calculation
                filt_opt.filt_std = 4;
                try
                    sm_data = filterDisplacements(U, V, filt_opt);
                    U = sm_data.U;
                    V = sm_data.V;
                catch
                    warning('filterDisplacements not found in path. Falling back to raw displacements.');
                end
                
                % Calculate gradients.
                % (Displacements are in pixels, so gradient over a 1-pixel grid directly gives unitless strain)
                [Hxx, Hxy] = gradient(U);
                [Hyx, Hyy] = gradient(V);
                
                % Calculate Effective Shear Strain
                try
                    Eeff_grid = calcEffectiveE(Hxx, Hxy, Hyx, Hyy);
                catch
                    Eeff_grid = sqrt(0.5*(Hxx - Hyy).^2 + 0.5*(Hxy + Hyx).^2);
                end
                Eeff = Eeff_grid(:);
                
                % Plot
                fig_strain = figure('Name', ['DIC Strain ', dic_prefix], 'Visible', 'off');
                
                if plot_ci_background
                    plot(ebsdUndef, cibw_undef, 'micronbar', 'off');
                    hold on;
                    plot(ebsdUndef, Eeff, 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
                else
                    plot(ebsdUndef, Eeff, 'micronbar', 'off');
                end
                
                mtexColorMap viridis;
                % Eeff_valid = Eeff(~isnan(Eeff));
                % if ~isempty(Eeff_valid)
                %     setColorRange([0, prctile(Eeff_valid, 98)]);
                % end
                setColorRange([0, 0.1]);
                mtexColorbar('title', 'Effective Shear Strain');
                title(sprintf('DIC Effective Shear Strain - %s', dic_prefix));
                
                exportScaledFigure(fig_strain, fullfile(outDir, sprintf('%s_DIC_Strain.jpg', baseName)), 'noDateString');
                close(fig_strain);
            else
                fprintf('  DIC displacement fields for %s not found in undeformed data.\n', dic_prefix);
            end
        end
    end
    
    % --- Calculate and Plot GNDs ---
    if plot_gnds
        gnd_file = fullfile(baseDir, 'Processed_Output', sprintf('GND_Data_%s.mat', dic_prefix));
        
        if isempty(dic_prefix) || ~exist(gnd_file, 'file')
            warning('GND data not found for %s (looked for %s). Skipping GND plotting.', currentFile, gnd_file);
        else
            fprintf('  Loading GNDs for %s from %s...\n', currentFile, sprintf('GND_Data_%s.mat', dic_prefix));
            gnd_data = load(gnd_file, 'ebsdDeformed', 'rho', 'factor', 'dS');
            
            ebsdSGND = gnd_data.ebsdDeformed;
            rho = gnd_data.rho;
            factor = gnd_data.factor;
            dS = gnd_data.dS;
            
            % --- Apply Specimen Mask to GND EBSD and Recalculate Grains ---
            if exist(mask_file, 'file')
                % Map the 1D mask via coordinates since ebsdSGND might be gridified
                valid_sgnd = ~isnan(ebsdSGND.x) & ~isnan(ebsdSGND.y);
                sgnd_x = ebsdSGND.x(valid_sgnd);
                sgnd_y = ebsdSGND.y(valid_sgnd);
                
                idx = knnsearch([ebsd_1D.x, ebsd_1D.y], [sgnd_x, sgnd_y]);
                
                sgnd_mask = false(size(ebsdSGND));
                sgnd_mask(valid_sgnd) = specimenlocation(idx);
                
                ebsdSGND(~sgnd_mask) = 'notIndexed';
            end
            
            fprintf('  Recalculating grains for GNDs (to exclude vacuum artifacts)...\n');
            [grainsDef, ebsdSGND.grainId] = calcGrains(ebsdSGND, 'minPixel', 1000, 'angle', 5*degree);
        
        % Total Edge and Screw
        GND_E = factor * sum(abs(rho(:, dS.isEdge)), 2, 'omitnan');
        GND_S = factor * sum(abs(rho(:, dS.isScrew)), 2, 'omitnan');
        
        % Unique <111> Burgers vectors (there are 4 in BCC)
        b_ref = unique(dS.b, 'noSymmetry', 'antipodal');
        
        GND_E_b = zeros(length(ebsdSGND), length(b_ref));
        GND_S_b = zeros(length(ebsdSGND), length(b_ref));
        GND_Tot_b = zeros(length(ebsdSGND), length(b_ref));
        
        for b_idx = 1:length(b_ref)
            match_idx = angle(dS.b, b_ref(b_idx), 'noSymmetry', 'antipodal') < 1e-5;
            edge_idx = match_idx & dS.isEdge;
            screw_idx = match_idx & dS.isScrew;
            
            GND_E_b(:, b_idx) = factor * sum(abs(rho(:, edge_idx)), 2, 'omitnan');
            GND_S_b(:, b_idx) = factor * sum(abs(rho(:, screw_idx)), 2, 'omitnan');
            GND_Tot_b(:, b_idx) = factor * sum(abs(rho(:, match_idx)), 2, 'omitnan');
        end
        
        % Plot 0: Total GND Densities
        fig_es = figure('Name', 'Total GNDs', 'Visible', 'off');
        % newMtexFigure('layout', [1, 2]);
        
        if plot_ci_background
            plot(ebsd, cibw, 'micronbar', 'off'); hold on;
            plot(ebsdSGND, log10(GND_E+GND_S), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
        else
            plot(ebsdSGND, log10(GND_E+GND_S), 'micronbar', 'off');
        end
        mtexColorMap parula;
        setColorRange(gnd_limits);
        %title('Total GND Density');
               
        %mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
        exportScaledFigure(fig_es, fullfile(outDir, sprintf('%s_GND_Total_nocolorbar.jpg', baseName)), 'noDateString');
        close(fig_es);
        
        % warning('skipping individual screw/edge plotting etc')
        % continue

        % Plot 1: Total Edge and Screw Densities
        fig_es = figure('Name', 'Total Edge and Screw GNDs', 'Visible', 'off');
        newMtexFigure('layout', [1, 2]);
        
        nextAxis(1, 1);
        if plot_ci_background
            plot(ebsd, cibw, 'micronbar', 'off'); hold on;
            plot(ebsdSGND, log10(GND_E), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
        else
            plot(ebsdSGND, log10(GND_E), 'micronbar', 'off');
        end
        mtexColorMap parula;
        setColorRange(gnd_limits);
        title('Total Edge GND Density');
        
        nextAxis(1, 2);
        if plot_ci_background
            plot(ebsd, cibw, 'micronbar', 'off'); hold on;
            plot(ebsdSGND, log10(GND_S), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
        else
            plot(ebsdSGND, log10(GND_S), 'micronbar', 'off');
        end
        mtexColorMap parula;
        setColorRange(gnd_limits);
        title('Total Screw GND Density');
        
        mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
        exportScaledFigure(fig_es, fullfile(outDir, sprintf('%s_GND_Total_Edge_Screw.jpg', baseName)), 'noDateString');
        close(fig_es);
        
        % Plot 2: Edge dislocations for the 4 <111> vectors
        fig_edge_b = figure('Name', 'Edge GNDs per b', 'Visible', 'off');
        newMtexFigure('layout', [2, 2]);
        for b_idx = 1:length(b_ref)
            nextAxis;
            if plot_ci_background
                plot(ebsd, cibw, 'micronbar', 'off'); hold on;
                plot(ebsdSGND, log10(GND_E_b(:, b_idx)), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
            else
                plot(ebsdSGND, log10(GND_E_b(:, b_idx)), 'micronbar', 'off');
            end
            mtexColorMap parula;
            setColorRange(gnd_limits);
            title(sprintf('Edge GNDs b=[%s]', char(b_ref(b_idx).round)));
        end
        mtexColorbar('title', 'log_{10}(GND_E m^{-2})');
        exportScaledFigure(fig_edge_b, fullfile(outDir, sprintf('%s_GND_Edge_per_b.jpg', baseName)), 'noDateString');
        close(fig_edge_b);
        
        % Plot 3: Screw dislocations for the 4 <111> vectors
        fig_screw_b = figure('Name', 'Screw GNDs per b', 'Visible', 'off');
        newMtexFigure('layout', [2, 2]);
        for b_idx = 1:length(b_ref)
            nextAxis;
            if plot_ci_background
                plot(ebsd, cibw, 'micronbar', 'off'); hold on;
                plot(ebsdSGND, log10(GND_S_b(:, b_idx)), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
            else
                plot(ebsdSGND, log10(GND_S_b(:, b_idx)), 'micronbar', 'off');
            end
            mtexColorMap parula;
            setColorRange(gnd_limits);
            title(sprintf('Screw GNDs b=[%s]', char(b_ref(b_idx).round)));
        end
        mtexColorbar('title', 'log_{10}(GND_S m^{-2})');
        exportScaledFigure(fig_screw_b, fullfile(outDir, sprintf('%s_GND_Screw_per_b.jpg', baseName)), 'noDateString');
        close(fig_screw_b);
        
        % Plot 4: Total GNDs per b (Edge + Screw combined)
        fig_tot_b = figure('Name', 'Total GNDs per b', 'Visible', 'off');
        newMtexFigure('layout', [2, 2]);
        for b_idx = 1:length(b_ref)
            nextAxis;
            if plot_ci_background
                plot(ebsd, cibw, 'micronbar', 'off'); hold on;
                plot(ebsdSGND, log10(GND_Tot_b(:, b_idx)), 'micronbar', 'off', 'FaceAlpha', gnd_alpha);
            else
                plot(ebsdSGND, log10(GND_Tot_b(:, b_idx)), 'micronbar', 'off');
            end
            hold on;
            plot(grainsDef.boundary, 'linecolor', 'w', 'linewidth', 1.5);
            mtexColorMap parula;
            setColorRange(gnd_limits);
            title(sprintf('Total GNDs b=[%s]', char(b_ref(b_idx).round)));
        end
        mtexColorbar('title', 'log_{10}(GND_{Tot} m^{-2})');
        exportScaledFigure(fig_tot_b, fullfile(outDir, sprintf('%s_GND_Total_per_b.jpg', baseName)), 'noDateString');
        close(fig_tot_b);
        end
    end
end
fprintf('\nFinished processing all files.\n');





%% grain calc

minSize = 30;
grainThresh = 5*degree;

[grains, ebsd.grainId] = calcGrains(ebsd, 'minPixel', minSize, 'angle', grainThresh);
ebsd(grains(grains.numPixel<=minSize)) = 'notIndexed';
[grains, ebsd.grainId] = calcGrains(ebsd, 'minPixel', minSize, 'angle', grainThresh);
%% grain plot with burgers vector quiver
figure('Name', 'Grain Map with b traces');
plot(grains);
hold on;

% Extract the reference orientation from grain 2 (used in analysis below)
ori115693 = grains(2).meanOrientation;

% Define the 4 unique <111> Burgers vectors for BCC and corresponding colors
b_vecs_leg = [Miller(1,1,1,ebsd.CS), Miller(-1,1,1,ebsd.CS), Miller(1,-1,1,ebsd.CS), Miller(1,1,-1,ebsd.CS)];
b_labels_leg = {'[1 1 1]', '[-1 1 1]', '[1 -1 1]', '[1 1 -1]'};
b_colors = lines(4);

h_leg = gobjects(4, 1);

for i = 1:4
    % Transform to specimen coordinates and set antipodal to true for line traces
    b_spec = ori115693 * b_vecs_leg(i);
    b_spec.antipodal = true;
    
    % Quiver the Burgers vector trace at the centroid of grain 2
    quiver(grains(2), b_spec, 'Color', b_colors(i,:), 'LineWidth', 4);
    
    % Create an invisible dummy line for the legend (MTEX quiver doesn't easily map to legends)
    h_leg(i) = plot(nan, nan, 'Color', b_colors(i,:), 'LineWidth', 4, 'DisplayName', sprintf('b = %s', b_labels_leg{i}));
end

legend(h_leg, 'Location', 'best', 'FontSize', 11);
hold off;




%% is the screw [11-1] and edge [-111] jsut an artifact from the GND calculations (BUT SEEMS NOT TO BE THE CASE) CHECK FOR SOME BIAS IN THE SELECTION OF dISLOCATIONS BASED ON B VECTOR DIRECTION AND CURVATURE COMPONENTS OF SCRES/EDGES
ori115693 = grains(2).meanOrientation;
dS_rot = ori115693 * dS;

edges = dS.isEdge;
screws = dS.isScrew;

% Define the 4 unique <111> Burgers vectors for BCC
b_vecs = [Miller(1,1,1,dS.CS), Miller(-1,1,1,dS.CS), Miller(1,-1,1,dS.CS), Miller(1,1,-1,dS.CS)];
b_labels = {'[1 1 1]', '[-1 1 1]', '[1 -1 1]', '[1 1 -1]'};

% Plot 1: Textboxes with percentages (2x4 grid)
fig_text = figure('Name', 'Curvature Footprints (Text)', 'Position', [100 100 1400 600]);

for i = 1:4
    b_target = b_vecs(i);
    is_b = angle(dS.b, b_target, 'noSymmetry', 'antipodal') < 1e-5;

    % Extract all individual 3x3 curvature matrices
    M_screw_all = curvature(dS_rot(screws & is_b).tensor).M;
    M_edge_all  = curvature(dS_rot(edges & is_b).tensor).M;

    % Sum the absolute values of the footprints to avoid cancellation
    M_screw_abs = sum(abs(M_screw_all), 3);
    M_edge_abs  = sum(abs(M_edge_all), 3);

    % Calculate unmeasurable percentage (3rd column) individually for each system
    col3_screw = sum(abs(M_screw_all(:,3,:)), 1);
    tot_screw  = sum(sum(abs(M_screw_all), 1), 2);
    pct_screw_all = (col3_screw(:) ./ tot_screw(:)) * 100;
    
    col3_edge = sum(abs(M_edge_all(:,3,:)), 1);
    tot_edge  = sum(sum(abs(M_edge_all), 1), 2);
    pct_edge_all = (col3_edge(:) ./ tot_edge(:)) * 100;

    pct_screw_avg = mean(pct_screw_all);
    pct_screw_var = var(pct_screw_all, 1); % Using weight 1 avoids NaN if there is only 1 system
    
    pct_edge_avg = mean(pct_edge_all);
    pct_edge_var = var(pct_edge_all, 1);

    % Screw (Top Row)
    figure(fig_text);
    subplot(2, 4, i); axis off;
    str_s = sprintf('Screw %s (Abs Sum)\n\n%8.3f %8.3f %8.3f\n%8.3f %8.3f %8.3f\n%8.3f %8.3f %8.3f\n\nUnmeasurable (3rd col):\nAvg: %.1f%% | Var: %.1f', b_labels{i}, M_screw_abs(:,:,1)', pct_screw_avg, pct_screw_var);
    text(0.5, 0.5, str_s, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Courier', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 10);

    % Edge (Bottom Row)
    subplot(2, 4, i+4); axis off;
    str_e = sprintf('Edge %s (Abs Sum)\n\n%8.3f %8.3f %8.3f\n%8.3f %8.3f %8.3f\n%8.3f %8.3f %8.3f\n\nUnmeasurable (3rd col):\nAvg: %.1f%% | Var: %.1f', b_labels{i}, M_edge_abs(:,:,1)', pct_edge_avg, pct_edge_var);
    text(0.5, 0.5, str_e, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Courier', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 10);
end

% Plot 2: MTEX Tensor Plots (2x4 grid)
fig_plots = figure('Name', 'MTEX Tensor Plots');
newMtexFigure('layout', [2, 4]);

for i = 1:4
    b_target = b_vecs(i);
    is_b = angle(dS.b, b_target, 'noSymmetry', 'antipodal') < 1e-5;

    tensor_screw = sum(curvature(dS_rot(screws & is_b).tensor));
    tensor_edge  = sum(curvature(dS_rot(edges & is_b).tensor));

    % Screw Plot (Top Row)
    nextAxis(1, i);
    plot(tensor_screw);
    title(sprintf('Screw %s', b_labels{i}));

    % Edge Plot (Bottom Row)
    nextAxis(2, i);
    plot(tensor_edge);
    title(sprintf('Edge %s', b_labels{i}));

    % Console Output: Direction of the Burgers vector with respect to the sample surface normal
    b_spec = ori115693 * b_target;
    fprintf('\n--- Alignment of %s Burgers vector ---\n', b_labels{i});
    fprintf('Angle to X (RD): %.1f°\n', angle(b_spec, xvector, 'antipodal') / degree);
    fprintf('Angle to Y (TD): %.1f°\n', angle(b_spec, yvector, 'antipodal') / degree);
    fprintf('Angle to Z (ND - Out of plane): %.1f°\n', angle(b_spec, zvector, 'antipodal') / degree);
end

% Balance color ranges across the MTEX subplots so magnitudes are visually comparable
setColorRange equal;