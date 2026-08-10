%% forwardDeformed_05_SSLIP.m
% Performs multi-slip SSLIP optimization on all Deformed Forward datasets 
% (from forwardDeformed_03_IMPORTER).
% Features: Toggleable slip system filtering, Combined Error Index, lightweight JPG export,
% and comprehensive Markdown reporting.
% Outputs results directly to: Processed_Output_Deformed_Forward/SSLIP_Results.

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\SSLIP\src');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% ========================================================================
%  USER CONFIGURATION & PRESETS (Modify as desired)
% =========================================================================
cfg = struct();

% --- PRESET A: All Slip Systems with Low SF Threshold ---
cfg.analysis_name    = 'Standard_AllSlipSystems_LowSF_rotationEnabled';
cfg.slip_system_mode = 'all';       % 'all', 'a_only', 'ca_only', or 'basal_pyr'
cfg.threshold_SF     = 0.01;        % Minimum Schmid Factor threshold
cfg.stress_vector    = yvector;     % Uniaxial stress loading axis
cfg.restrict2GrainId = 0;           % 0 = all, [1,2,3] grainIds
cfg.stressSign       = 1;           % -1 compression, 1 extension
cfg.enableRotation   = 1;           % enable rotation fields

cfg.sslip_method     = 1;           % IDMethod (1 = multi-slip optimization)
cfg.posConstr        = 1;           % Enforce non-negative shear activity
cfg.threshResidual   = 0.02;        % Residual threshold for correlation rate calculation
cfg.threshResidualFraction = 0.50;  % Threshold residual as fraction of effective shear
cfg.minEeff          = 0.01;        % Minimum effective shear strain required for metric evaluation
cfg.minSize          = 300;         % Minimum grain size in pixels to analyze
cfg.logscale = 0; 

% Lightweight JPG export settings
cfg.jpg_dpi          = 150;
cfg.jpg_target_px    = 1200;
cfg.clim_percentile  = 99;
cfg.plotBoundariesOverDisps = 0;
% =========================================================================

%% 1. Load Registry & Configure Step Selection
[registry, defaults] = a00_Dataset_Registry_DifferentialDIC();

% Filter for non-excluded deformed steps with valid TKD maps and DIC groups
defSteps = registry(~[registry.isReference] & ~[registry.exclude] & ...
                    ~cellfun(@isempty, {registry.tkd_fname}) & ~cellfun(@isempty, {registry.dic_group}));

if isempty(defSteps)
    fprintf('========================================================================\n');
    fprintf('NOTICE: No valid deformed step entries found in registry.\n');
    fprintf('========================================================================\n');
    return;
end

inDir  = fullfile(defaults.dir, defaults.outDir_deformed);
outDir = fullfile(defaults.dir, defaults.outDir_deformed, ['SSLIP_Results_', cfg.analysis_name]);
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========================================================================\n');
fprintf('--- STARTING DEFORMED FORWARD SSLIP BATCH ANALYSIS (%d Steps) ---\n', numel(defSteps));
fprintf('========================================================================\n');

%% 2. Iterate Over Each Deformed Step
for stepIdx = 2 % 1:numel(defSteps)
    defStep = defSteps(stepIdx);
    prefix = defStep.stepName;
    fprintf('\n------------------------------------------------------------------------\n');
    fprintf('Processing Step [%s] (%d/%d) - DIC Group: %s\n', prefix, stepIdx, numel(defSteps), defStep.dic_group);
    fprintf('------------------------------------------------------------------------\n');
    
    matPath = fullfile(inDir, [prefix, '_deformedEBSD_withForwardDIC.mat']);
    if ~exist(matPath, 'file')
        warning('Aligned dataset not found at:\n  %s\nRun forwardDeformed_03_IMPORTER first for Step [%s]. Skipping.', matPath, prefix);
        continue;
    end
    
    fprintf('Loading aligned dataset: %s...\n', matPath);
    dataLoaded = load(matPath);
    newEBSD = dataLoaded.newEBSD;
    
    %% 3. Retrieve or Reconstruct Grains
    grainThresh = 5*degree; 
    alpha = 1;

    if isfield(dataLoaded, 'grains')
        grains = dataLoaded.grains;
        fprintf('  -> Successfully loaded precomputed grains directly from .mat dataset.\n');
    else
        fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
        thresh_ci = defaults.thresh_ci;
        newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
        [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', cfg.minSize, 'angle', grainThresh, 'alpha', alpha);
    end

    % Save reference Grain Map figure for this step
    figG = figure('Visible', 'off');
    plot(grains, grains.meanOrientation); hold on;
    text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
    title(sprintf('%s: Analyzed Grains (>%d px)', prefix, cfg.minSize), 'Interpreter', 'none');
    hold off;
    exportScaledFigure(figG, fullfile(outDir, sprintf('%s_GrainMap.jpg', prefix)), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
    close(figG);

    %% 4. Locate the Forward DIC Displacement Fields
    all_fields = fields(newEBSD.prop);
    ux_idx = find(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true), 1);
    uy_idx = find(contains(all_fields, '_Displacements_y_displacement', 'IgnoreCase', true), 1);
    
    if isempty(ux_idx) || isempty(uy_idx)
        warning('Could not find x and y displacement fields in newEBSD.prop for step %s. Skipping.', prefix);
        continue;
    end
    ux_field = all_fields{ux_idx};
    uy_field = all_fields{uy_idx};
    
    hxx_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hxy_idx = find(contains(all_fields, 'du', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    hyx_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dx', 'IgnoreCase', true) & ~contains(all_fields, 'dy', 'IgnoreCase', true), 1);
    hyy_idx = find(contains(all_fields, 'dv', 'IgnoreCase', true) & contains(all_fields, 'dy', 'IgnoreCase', true) & ~contains(all_fields, 'dx', 'IgnoreCase', true), 1);
    
    hasGradients = ~isempty(hxx_idx) && ~isempty(hxy_idx) && ~isempty(hyx_idx) && ~isempty(hyy_idx);

    %% 5. Define Candidate Slip Systems based on Toggled Mode
    CS = newEBSD.CS;
    sS_basal_all = slipSystem.basal(CS);
    sS_prism_all = slipSystem.prismaticA(CS);
    sS_pyra_all  = slipSystem.pyramidalA(CS);
    sS_ca_all    = slipSystem.pyramidalCA(CS);

    switch lower(cfg.slip_system_mode)
        case {'all', 'everything'}
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all; sS_ca_all];
            mode_str = 'All Systems (Basal, Prism, Pyr <a>, Pyr <c+a>)';
        case {'a', 'a_only', 'only_a'}
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all];
            mode_str = 'Only <a> Type Systems (Basal, Prism, Pyr <a>)';
        case {'ca', 'c+a', 'ca_only', 'only_ca'}
            sS_raw = sS_ca_all;
            mode_str = 'Only <c+a> Type Systems (Pyr <c+a>)';
        case {'basal_pyramidal', 'basal_pyr', 'only_basal_pyr'}
            sS_raw = [sS_basal_all; sS_pyra_all; sS_ca_all];
            mode_str = 'Basal & Pyramidal Systems (Excluding Prismatic)';
        otherwise
            sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all; sS_ca_all];
            mode_str = 'Default: All Systems';
    end

    sS = symmetrise(sS_raw, 'antipodal');
    fprintf('  -> Slip System Mode [%s]: %d candidate systems initialized.\n', cfg.slip_system_mode, length(sS));

    %% 6. Evaluate Significant Grains
    valid_grains = grains(grains.numPixel > cfg.minSize);
    if isempty(valid_grains)
        warning('No grains larger than %d pixels in step %s. Skipping.', cfg.minSize, prefix);
        continue;
    end

    grains_ids = valid_grains.id; 
    if isfield(cfg, 'restrict2GrainId') && ~isequal(cfg.restrict2GrainId, 0) && ~isempty(cfg.restrict2GrainId) && all(cfg.restrict2GrainId > 0)
        grains_ids = intersect(grains_ids, cfg.restrict2GrainId);
        if isempty(grains_ids)
            warning('None of the restricted grain IDs [%s] were found among significant grains. Skipping.', num2str(cfg.restrict2GrainId));
            continue;
        end
        fprintf('  -> Restricting analysis to %d specified grain(s).\n', length(grains_ids));
    else
        fprintf('  -> Analyzing all %d significant grains...\n', length(grains_ids));
    end
    
    summary_results = struct('grainId', {}, 'numPixel', {}, 'comb_index', {}, 'corr_rate', {}, 'total_gamma', {}, 'top_sys', {}, 'baseName', {});

    combinedMatPath = fullfile(outDir, sprintf('%s_%s_CombinedSSLIP.mat', prefix, cfg.slip_system_mode));
    if exist(combinedMatPath, 'file')
        fprintf('  -> Loading combined SSLIP results from: %s\n', combinedMatPath);
        load(combinedMatPath, 'stepResults');
    else
        stepResults = struct();
    end

    for g_i = 1:length(grains_ids)
        clear IDoptions;
        close all;
        grainId = grains_ids(g_i);

        fprintf('\n   [Grain %d (%d px)] Initializing slip geometry...\n', grainId, grains(grainId).numPixel);

        % Project global slip systems into local grain lattice orientation
        sSLocal = grains(grainId).meanOrientation * sS;

        stressState = cfg.stressSign * stressTensor.uniaxial(cfg.stress_vector);
        
        % Reconfigure slip systems to ensure positive Schmid factor
        sf_vals = sSLocal.SchmidFactor(stressState);
        sf_signs = sign(sf_vals);
        sf_signs(sf_signs == 0) = 1;
        sSLocal.b = sf_signs .* sSLocal.b;
        
        SF_crit = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
        sSLocal_filt = sSLocal(SF_crit);
        sS_filt = sS(SF_crit);

        if isempty(sS_filt)
            warning('   No slip systems with SF >= %.2f for Step [%s] Grain %d. Skipping grain.', cfg.threshold_SF, prefix, grainId);
            continue;
        end

        % Extract EBSD data slice specific to this grain and gridify
        ebsd_grain = newEBSD(grains(grainId));
        ebsd_grain = ebsd_grain.gridify;

        % Convert displacement field to physical units via dPos
        U = ebsd_grain.prop.(ux_field) * newEBSD.dPos; 
        V = ebsd_grain.prop.(uy_field) * newEBSD.dPos;

        plotBaseName = sprintf('%s_Grain%d_%s', prefix, grainId, cfg.slip_system_mode);

        IDoptions.IDMethod    = cfg.sslip_method;
        IDoptions.filterSize  = 4;
        IDoptions.coarsegrain = 4;
        IDoptions.threshResidual = cfg.threshResidual;
        if isfield(cfg, 'threshResidualFraction')
            IDoptions.threshResidualFraction = cfg.threshResidualFraction;
        end
        IDoptions.minEeff     = cfg.minEeff;
        IDoptions.posConstr   = cfg.posConstr;
        IDoptions.plotSSLIP   = 0;
        IDoptions.layout      = [];
        IDoptions.plotTraces  = 1;
        IDoptions.stress      = stressState;
        IDoptions.casename    = plotBaseName;
        IDoptions.plotDefGrad = 1;
        IDoptions.saveFig     = 0;
        IDoptions.plotname    = fullfile(outDir, plotBaseName);
        IDoptions.sizeAdjust  = 0.5;
        IDoptions.plotResidual = 1;
        IDoptions.logscale    = cfg.logscale;
        IDoptions.enableRotation = cfg.enableRotation;

        if cfg.logscale, IDoptions.logmin = 0.01; end
        try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
        IDoptions.comment = sprintf('Grain %d [%s Mode]', grainId, cfg.slip_system_mode);
        IDoptions.NoSs = 1:length(sS_filt);

        if hasGradients
            IDoptions.Hxx = ebsd_grain.prop.(all_fields{hxx_idx});
            IDoptions.Hxy = ebsd_grain.prop.(all_fields{hxy_idx});
            IDoptions.Hyx = ebsd_grain.prop.(all_fields{hyx_idx});
            IDoptions.Hyy = ebsd_grain.prop.(all_fields{hyy_idx});
        end
        
        grainKey = sprintf('grain_%d', grainId);

        try
            if isfield(stepResults, grainKey)
                fprintf('   Loading existing SSLIP optimization results for Grain %d from combined mat...\n', grainId);
                ebsdID       = stepResults.(grainKey).ebsdID;
                sSLocal_filt = stepResults.(grainKey).sSLocal_filt;
                optOut       = stepResults.(grainKey).optOut;
                m_res        = stepResults.(grainKey).m_res;
                corr_rate    = stepResults.(grainKey).corr_rate;
            else
                fprintf('   Running SSLIP optimization (%d candidate systems)...\n', length(sS_filt));
                [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);
    
                % Save deformation gradient summary figure using new modular function
                if isfield(cfg, 'plotDefGrad') && cfg.plotDefGrad
                    plotSSLIP_DeformationFields(ebsdID, IDoptions);
                    defGrad_fig = gcf;
                    exportScaledFigure(defGrad_fig, fullfile(outDir, [plotBaseName, '_DefGrad.jpg']), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    close(defGrad_fig);
                end
    
                % Calculate robust quality metrics
                Eeff_eval = ebsdID.prop.Eeff(:);
                eval_mask = (Eeff_eval >= cfg.minEeff) & ~isnan(Eeff_eval);
                
                res_map  = ebsdID.prop.residualEeff(:);
                res_eval = res_map(eval_mask);
                total_attempted = sum(eval_mask);
                
                corr_mask = ~isnan(res_eval) & (res_eval < cfg.threshResidual);
                num_corr  = sum(corr_mask);
                corr_rate = (num_corr / max(1, total_attempted)) * 100;
                corr_frac = num_corr / max(1, total_attempted);
                
                if num_corr > 0
                    r_corr = mean(res_eval(corr_mask));
                else
                    r_corr = cfg.threshResidual;
                end
                m_res = r_corr / max(1e-3, corr_frac);
    
                % Save into the combined struct immediately after solving
                stepResults.(grainKey).ebsdID = ebsdID;
                stepResults.(grainKey).sSLocal_filt = sSLocal_filt;
                stepResults.(grainKey).optOut = optOut;
                stepResults.(grainKey).m_res = m_res;
                stepResults.(grainKey).corr_rate = corr_rate;
                save(combinedMatPath, 'stepResults', 'cfg', '-v7.3');
            end
            
            mean_acts = mean(ebsdID.prop.slipIDcor, 2, 'omitnan');
            [sort_acts, sort_idx] = sort(mean_acts, 'descend');
            
            top_str = '';
            for t_i = 1:min(3, length(sort_idx))
                sys_idx = sort_idx(t_i);
                sys_obj = sS_filt(sys_idx);
                if sys_obj.b.norm > 4.0, fam_str = 'Pyr <c+a>';
                elseif abs(sys_obj.n.z) > 0.95, fam_str = 'Basal <a>';
                elseif abs(sys_obj.n.z) < 0.05, fam_str = 'Prism <a>';
                else, fam_str = 'Pyr <a>'; end
                top_str = [top_str, sprintf('%s (%s) [SF=%.2f, γ=%.4f]; ', fam_str, sys_obj.b.char, abs(sSLocal_filt(sys_idx).SchmidFactor(stressState)), sort_acts(t_i))];
            end
            
            % Export CSV Table of Active Systems
            SF_col   = abs(sSLocal_filt.SchmidFactor(IDoptions.stress));
            SF_col   = round(SF_col(:)*1000)/1000;
            acts_col = mean_acts(:);
            id_col   = (1:length(acts_col))';
            b_col    = arrayfun(@(x) x.b.char, sS_filt(:), 'UniformOutput', false);
            n_col    = arrayfun(@(x) x.n.char, sS_filt(:), 'UniformOutput', false);
            
            num_acts = length(acts_col);
            num_sys  = length(sS_filt);
            if num_acts > num_sys
                for i_ex = num_sys+1:num_acts
                    b_col{i_ex, 1} = 'Rotation';
                    n_col{i_ex, 1} = 'N/A';
                    SF_col(i_ex, 1) = NaN;
                end
            end
            
            T_sys = table(id_col, b_col, n_col, SF_col, acts_col, ...
                          'VariableNames', {'SlipSystem_ID', 'Burgers_Vector_b', 'Plane_Normal_n', 'Schmid_Factor', 'Mean_Shear_Activity'});
            writetable(T_sys, fullfile(outDir, [plotBaseName, '_SlipSystems.csv']));
            
            % --- NEW MODULAR PLOTTING SECTION ---
            close all;
            optOut.grains = grains;
            optOut.plotGrainContext = 1;
            optOut.plotDominant = 1;
            optOut.num_top_systems = 3;
            optOut.saveFig = 0;
            
            % Call Dominant Activities Plotter
            plotSSLIP_DominantActivities(ebsdID, sSLocal_filt, optOut);
            exportScaledFigure(gcf, [optOut.plotname, '_SlipActivities.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
            close(gcf);
            
            % Call Residual Plotter with context
            plotSSLIP_Residual(ebsdID, optOut);
            exportScaledFigure(gcf, [optOut.plotname, '_Residual.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
            close(gcf);
            
            % Call Rotation Plotter if rotation was calculated
            if isfield(ebsdID.prop, 'rotationIDcor')
                plotSSLIP_Rotation(ebsdID, optOut);
                exportScaledFigure(gcf, [optOut.plotname, '_Rotation.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                close(gcf);
            end
            close all;

            summary_results(end+1).grainId   = grainId; %#ok<AGROW>
            summary_results(end).numPixel    = grains(grainId).numPixel;
            summary_results(end).comb_index  = m_res;
            summary_results(end).corr_rate   = corr_rate;
            summary_results(end).total_gamma = sum(mean_acts);
            summary_results(end).top_sys     = top_str;
            summary_results(end).baseName    = plotBaseName;

            fprintf('   -> Complete! Combined Error Index: %.6f | Correlation Rate: %.1f%%\n', m_res, corr_rate);
        catch ME
            fprintf(2, '   [ERROR] SSLIP processing failed for Grain %d in step %s:\n%s\n', grainId, prefix, getReport(ME));
        end
    end
    
    if ~isempty(summary_results) && isfield(summary_results, 'grainId') && isempty(summary_results(1).grainId)
        summary_results(1) = [];
    end

    if ~isempty(summary_results)
        % 1. Print Leaderboard
        fprintf('\n     =================================================================================================================================\n');
        fprintf('     STEP [%s] SSLIP SUMMARY LEADERBOARD (Ranked by Combined Error Index)\n', prefix);
        fprintf('     =================================================================================================================================\n');
        fprintf('     Rank | Grain ID | Pixels | Combined Error Index | Correlation %% | Total Shear γ | Top Active Slip Systems \n');
        fprintf('     ---------------------------------------------------------------------------------------------------------------------------------\n');
        all_err = [summary_results.comb_index];
        [~, sort_r] = sort(all_err, 'ascend');
        ranked_summary = summary_results(sort_r);
        
        for r_i = 1:length(ranked_summary)
            pr = ranked_summary(r_i);
            fprintf('     #%02d  | Grain %02d | %6d |       %.6f       |    %5.1f%%     |    %.4f     | %s \n', ...
                    r_i, pr.grainId, pr.numPixel, pr.comb_index, pr.corr_rate, pr.total_gamma, pr.top_sys);
        end
        fprintf('     =================================================================================================================================\n\n');
        
        % 2. Save Summary Table CSV
        T_sum = table((1:length(ranked_summary))', [ranked_summary.grainId]', [ranked_summary.numPixel]', ...
                      [ranked_summary.comb_index]', [ranked_summary.corr_rate]', [ranked_summary.total_gamma]', {ranked_summary.top_sys}', ...
                      'VariableNames', {'Rank', 'Grain_ID', 'Num_Pixels', 'Combined_Error_Index', 'Correlation_Rate_Pct', 'Total_Mean_Shear', 'Dominant_Slip_Systems'});
        writetable(T_sum, fullfile(outDir, sprintf('%s_Summary_Leaderboard.csv', prefix)));
        
        % 3. Generate Interactive Markdown Report
        md_path = fullfile(outDir, sprintf('%s_SSLIP_Report.md', prefix));
        fprintf('     -> Writing Markdown Report with embedded JPG maps to:\n        %s\n', md_path);
        
        fid = fopen(md_path, 'w');
        fprintf(fid, '# Standalone SSLIP Analysis Report\n\n');
        fprintf(fid, '- **DIC Step:** `%s`\n', prefix);
        fprintf(fid, '- **Slip System Mode:** `%s` (%s)\n', cfg.slip_system_mode, mode_str);
        fprintf(fid, '- **Schmid Factor Threshold:** `%.2f`\n', cfg.threshold_SF);
        fprintf(fid, '- **Residual Threshold:** `%.3f`\n\n', cfg.threshResidual);
        
        fprintf(fid, '## Analyzed Grain Map\n');
        fprintf(fid, '![Analyzed Grains](%s_GrainMap.jpg)\n\n---\n\n', prefix);
        
        fprintf(fid, '## Grain Analysis Leaderboard\n');
        fprintf(fid, 'Ranked by lowest **Combined Error Index** ($r_{\\text{comb}} = \\frac{r_{\\text{corr}}}{\\eta_{\\text{corr}}}$).\n\n');
        fprintf(fid, '| Rank | Grain ID | Pixels | Combined Error Index | Correlation %% | Total $\\gamma$ | Dominant Active Systems |\n');
        fprintf(fid, '| :---: | :---: | :---: | :---: | :---: | :---: | :--- |\n');
        for r_i = 1:length(ranked_summary)
            pr = ranked_summary(r_i);
            fprintf(fid, '| **#%02d** | Grain %d | `%d` | `%.6f` | `%5.1f%%` | `%.4f` | %s |\n', ...
                    r_i, pr.grainId, pr.numPixel, pr.comb_index, pr.corr_rate, pr.total_gamma, pr.top_sys);
        end
        fprintf(fid, '\n---\n\n');
        
        fprintf(fid, '## Detailed Slip Activity & Residual Maps per Grain\n');
        for r_i = 1:length(ranked_summary)
            pr = ranked_summary(r_i);
            b_name = pr.baseName;
            
            fprintf(fid, '### Rank #%02d: Grain %d\n', r_i, pr.grainId);
            fprintf(fid, '- **Combined Error Index:** `%.6f`\n', pr.comb_index);
            fprintf(fid, '- **Correlation Rate:** `%.1f%%`\n', pr.corr_rate);
            fprintf(fid, '- **Total Mean Shear Activity:** `%.4f`\n', pr.total_gamma);
            fprintf(fid, '- **Dominant Systems:** %s\n\n', pr.top_sys);
            
            if exist(fullfile(outDir, [b_name, '_DefGrad.jpg']), 'file')
                fprintf(fid, '#### SSLIP Deformation & Strain Fields ($U_x$, $E_{\\text{eff}}$, $H_{ij}$)\n');
                fprintf(fid, '![Deformation Fields %s](%s_DefGrad.jpg)\n\n', b_name, b_name);
            end
            
            fprintf(fid, '#### Slip System Activity Maps\n');
            fprintf(fid, '![Activity Maps %s](%s_SlipActivities.jpg)\n\n', b_name, b_name);
            
            fprintf(fid, '#### Residual Error Map\n');
            fprintf(fid, '![Residual Map %s](%s_Residual.jpg)\n\n', b_name, b_name);
            fprintf(fid, '---\n\n');
        end
        fclose(fid);
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL DEFORMED FORWARD SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
