%% step3to4_02_SSLIP_Analysis.m
% Performs multi-slip SSLIP optimization specifically on the Step 3 to 4 DIC correlation
% mapped directly onto the Step 03 (BW158) TKD reference grid.
% Features: Toggleable slip system filtering, Combined Error Index, lightweight JPG export,
% and comprehensive Markdown reporting.
% Outputs to: Processed_Output_Step3to4/SSLIP_Results

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
cfg.analysis_name    = 'Step3to4_AllSlipSystems_LowSF';
cfg.slip_system_mode = 'all';       % 'all', 'a_only', 'ca_only', or 'basal_pyr'
cfg.threshold_SF     = 0.01;        % Minimum Schmid Factor threshold
cfg.stress_vector    = yvector;     % Uniaxial stress loading axis
cfg.restrict2GrainId = 5;            % 0 = all, [1,2,3] grainIds
cfg.stressSign = 1;                %-1 compression, 1 extension

cfg.sslip_method     = 1;           % IDMethod (1 = multi-slip optimization)
cfg.posConstr        = 1;           % Enforce non-negative shear activity
cfg.threshResidual   = 0.02;        % Residual threshold for correlation rate calculation
cfg.threshResidualFraction = 0.50; % Threshold residual as fraction of effective shear
cfg.minEeff          = 0.01;        % Minimum effective shear strain required for metric evaluation
cfg.minSize          = 300;         % Minimum grain size in pixels to analyze
cfg.logscale = 0; 

% Lightweight JPG export settings
cfg.jpg_dpi          = 150;
cfg.jpg_target_px    = 1200;
cfg.clim_percentile  = 99;
cfg.plotBoundariesOverDisps = 0;

% % --- PRESET B: Only <a> Type Slip Systems ---
% cfg.analysis_name    = 'Step3to4_aType_Only';
% cfg.slip_system_mode = 'a_only';
% cfg.threshold_SF     = 0.01;
% cfg.stress_vector    = yvector;
% cfg.sslip_method     = 1;
% cfg.posConstr        = 1;
% cfg.threshResidual   = 0.02;
% cfg.minEeff          = 0.01;
% cfg.minSize          = 300;
% cfg.jpg_dpi          = 150;
% cfg.jpg_target_px    = 1200;
% =========================================================================

%% 1. Locate & Load Step 3 to 4 Aligned Dataset
[registry, defaults] = a00_Dataset_Registry();
inDir  = fullfile(defaults.dir, 'Processed_Output_Step3to4');
outDir = fullfile(inDir, ['SSLIP_Results_', cfg.analysis_name]);
if ~exist(outDir, 'dir'), mkdir(outDir); end

matFile = fullfile(inDir, 'step3to4_EBSD_withDIC.mat');
if ~exist(matFile, 'file')
    error('Aligned dataset not found at:\n  %s\nRun step3to4_01_IMPORTER first.', matFile);
end

fprintf('Loading aligned Step 3 to 4 dataset:\n  %s\n', matFile);
dataLoaded = load(matFile);
newEBSD = dataLoaded.newEBSD;

%% 2. Retrieve or Reconstruct Grains
grainThresh = 5*degree; 
alpha = 1;

if isfield(dataLoaded, 'grains')
    grains = dataLoaded.grains;
    fprintf('  -> Successfully loaded precomputed grains directly from dataset.\n');
else
    fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
    newEBSD(newEBSD.ci < defaults.thresh_ci) = 'notIndexed';
    [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', cfg.minSize, 'angle', grainThresh, 'alpha', alpha);
end

% Save overview grain map
figG = figure('Visible', 'off');
plot(grains, grains.meanOrientation); hold on;
text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold');
title('Step 03 (BW158) TKD: Analyzed Grains', 'Interpreter', 'none'); hold off;
exportScaledFigure(figG, fullfile(outDir, 'Step3to4_GrainMap.jpg'), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
close(figG);

%% 3. Extract DIC Loading Increments
all_fields = fields(newEBSD.prop);
disp_logical = contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true);
disp_fields  = all_fields(disp_logical);

prefix_list = cellfun(@(x) x(1:find(x=='_', 1)-1), disp_fields, 'UniformOutput', false);
prefix_list = unique(prefix_list);
fprintf('Found %d DIC loading group(s) mapped to grid.\n', length(prefix_list));

%% 3b. Overlay Grain Boundaries on Displacement Gradient Fields
if ~cfg.plotBoundariesOverDisps
    fprintf('\nSkip plotting Displacement Gradient & Grain Boundary Overlays...\n');
else
    fprintf('\nGenerating Displacement Gradient & Grain Boundary Overlays...\n');
    all_prop = fieldnames(newEBSD.prop);
    
    for p_idx = 1:length(prefix_list)
        prefix_ov = prefix_list{p_idx};
        
        dudx_idx = find(contains(all_prop, prefix_ov, 'IgnoreCase', true) & contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true));
        dudy_idx = find(contains(all_prop, prefix_ov, 'IgnoreCase', true) & contains(all_prop, 'du', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true));
        dvdx_idx = find(contains(all_prop, prefix_ov, 'IgnoreCase', true) & contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dx', 'IgnoreCase', true));
        dvdy_idx = find(contains(all_prop, prefix_ov, 'IgnoreCase', true) & contains(all_prop, 'dv', 'IgnoreCase', true) & contains(all_prop, 'dy', 'IgnoreCase', true));
        
        if ~isempty(dudx_idx) && ~isempty(dudy_idx) && ~isempty(dvdx_idx) && ~isempty(dvdy_idx)
            dudx_f = all_prop{dudx_idx(1)};
            dudy_f = all_prop{dudy_idx(1)};
            dvdx_f = all_prop{dvdx_idx(1)};
            dvdy_f = all_prop{dvdy_idx(1)};
            
            dudx_data = newEBSD.prop.(dudx_f);
            dudy_data = newEBSD.prop.(dudy_f);
            dvdx_data = newEBSD.prop.(dvdx_f);
            dvdy_data = newEBSD.prop.(dvdy_f);
            
            all_vals = [dudx_data(:); dudy_data(:); dvdx_data(:); dvdy_data(:)];
            all_vals = all_vals(~isnan(all_vals));
            % cmin = 0.2 * min(all_vals);
            % cmax = 0.2 * max(all_vals);
    
            cmin = prctile(all_vals, 100 - cfg.clim_percentile);
            cmax = prctile(all_vals, cfg.clim_percentile);
            
            grad_names = {'Hxx (du/dx)', 'Hxy (du/dy)', 'Hyx (dv/dx)', 'Hyy (dv/dy)'};
            grad_files = {'Hxx', 'Hxy', 'Hyx', 'Hyy'};
            grad_data  = {dudx_data, dudy_data, dvdx_data, dvdy_data};
            
            for g_i = 1:4
                fig_g = figure('Name', [grad_files{g_i} ' Overlay - ' prefix_ov], 'Position', [100, 100, 950, 900], 'Visible', 'off');
                t_lay = tiledlayout(fig_g, 1, 1, 'Padding', 'loose', 'TileSpacing', 'loose');
                nexttile(t_lay);
                plot(newEBSD, grad_data{g_i}, 'micronbar', 'off');
                colormap(gca, 'parula'); colorbar; caxis([cmin cmax]);
                hold on; plot(grains.boundary, 'lineColor', 'white', 'lineWidth', 1.5, 'micronbar', 'off');
                text(grains, grains.id, 'Color', 'white', 'FontWeight', 'bold'); hold off;
                
                % Assign title directly to the tiled layout to guarantee top margin reservation and prevent cropping during export
                title(t_lay, [grad_names{g_i} ' - ' strrep(prefix_ov, '_', '\_')], 'FontSize', 15, 'FontWeight', 'bold');
                
                drawnow; % Realize text bounds before exporting
                out_path = fullfile(outDir, sprintf('%s_DispGradient_%s.jpg', prefix_ov, grad_files{g_i}));
                exportScaledFigure(fig_g, out_path, 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                close(fig_g);
                fprintf('  -> Saved displacement gradient map: %s\n', out_path);
            end
        else
            warning('Displacement gradient fields incomplete for increment [%s]. Skipping overlay check figure.', prefix_ov);
        end
    end
end

%% 4. Define Candidate Slip Systems based on Toggled Mode
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
fprintf('Slip System Mode [%s]: %d candidate systems initialized.\n', cfg.slip_system_mode, length(sS));

%% 5. Evaluate Significant Grains Across DIC Increments
valid_grains = grains(grains.numPixel > cfg.minSize); 
if isempty(valid_grains)
    error('No valid grains (>%d px) found for analysis.', cfg.minSize);
end

grains_ids = valid_grains.id; 
if isfield(cfg, 'restrict2GrainId') && ~isequal(cfg.restrict2GrainId, 0) && ~isempty(cfg.restrict2GrainId) && all(cfg.restrict2GrainId > 0)
    grains_ids = intersect(grains_ids, cfg.restrict2GrainId);
    if isempty(grains_ids)
        error('None of the restricted grain IDs [%s] were found among significant grains (>%d px).', num2str(cfg.restrict2GrainId), cfg.minSize);
    end
    fprintf('Restricting analysis to %d specified grain(s): [%s]\n', length(grains_ids), num2str(grains_ids(:)'));
else
    fprintf('Analyzing all %d significant grains...\n', length(grains_ids));
end

summary_results = struct('grainId', {}, 'numPixel', {}, 'comb_index', {}, 'corr_rate', {}, 'total_gamma', {}, 'top_sys', {}, 'baseName', {});

for i_i = 1:length(prefix_list)
    prefix = prefix_list{i_i};
    fprintf('\n========================================================================\n');
    fprintf('Evaluating SSLIP for DIC Increment [%s] (%d/%d)\n', prefix, i_i, length(prefix_list));
    fprintf('========================================================================\n');
    
    for g_i = 1:length(grains_ids)
        clear IDoptions;
        close all;
        grainId = grains_ids(g_i);
        
        fprintf('   [Grain %d (%d px)] Initializing slip geometry...\n', grainId, grains(grainId).numPixel);
        sSLocal = grains(grainId).meanOrientation * sS;
        
        stressState = cfg.stressSign * stressTensor.uniaxial(cfg.stress_vector);
        
        % Reconfigure slip systems to ensure positive Schmid factor
        sf_vals = sSLocal.SchmidFactor(stressState);
        sf_signs = sign(sf_vals);
        sf_signs(sf_signs == 0) = 1;
        sSLocal.b = sf_signs .* sSLocal.b;
        
        SF_crit = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
        sSLocal_filt = sSLocal(SF_crit);
        sS_filt      = sS(SF_crit);
        
        if isempty(sS_filt)
            warning('   No slip systems with SF >= %.2f for Grain %d. Skipping.', cfg.threshold_SF, grainId);
            continue;
        end
        
        ebsd_grain = newEBSD(grains(grainId));
        ebsd_grain = ebsd_grain.gridify;
        
        ux_fname = [prefix, '_Displacements_x_displacement'];
        uy_fname = [prefix, '_Displacements_y_displacement'];
        
        if ~isfield(ebsd_grain.prop, ux_fname) || ~isfield(ebsd_grain.prop, uy_fname)
            warning('   Missing displacement data %s for grain %d. Skipping.', ux_fname, grainId);
            continue;
        end
        
        U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos; 
        V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
        
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
        IDoptions.plotSSLIP   = 0; % Disable automatic slow figure saving
        IDoptions.layout      = [];
        IDoptions.plotTraces  = 1;
        IDoptions.stress      = stressState;
        IDoptions.casename    = plotBaseName;
        IDoptions.plotDefGrad = 1;
        IDoptions.saveFig     = 0; % Avoid heavy PNG generation; export lightweight JPGs manually
        IDoptions.plotname    = fullfile(outDir, plotBaseName);
        IDoptions.sizeAdjust  = 0.5;
        IDoptions.plotResidual = 1;
        IDoptions.logscale    = cfg.logscale;

        IDoptions.Hxx = ebsd_grain.prop.B00001_Gradients_du_dx_RC;
        IDoptions.Hxy = ebsd_grain.prop.B00001_Gradients_du_dy_RC;
        IDoptions.Hyx = ebsd_grain.prop.B00001_Gradients_dv_dx_RC;
        IDoptions.Hyy = ebsd_grain.prop.B00001_Gradients_dv_dy_RC;

        if cfg.logscale
            IDoptions.logmin      = 0.01;
        end
        try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
        IDoptions.comment     = sprintf('Grain %d [%s Mode]', grainId, cfg.slip_system_mode);
        IDoptions.NoSs        = 1:length(sS_filt);
        resMatPath = fullfile(outDir, [plotBaseName, '_SSLIP_results.mat']);
        
        try
            %if isfile(resMatPath)
            if false
                fprintf('   [Cache Hit] Loading pre-computed SSLIP optimization from MAT file...\n');
                cached_data = load(resMatPath);
                ebsdID       = cached_data.ebsdID;
                optOut       = cached_data.optOut;
                sSLocal_filt = cached_data.sSLocal_filt;
                m_res        = cached_data.m_res;
                corr_rate    = cached_data.corr_rate;
            else
                fprintf('   Running SSLIP optimization (%d systems)...\n', length(sS_filt));
                [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);
                
                % Save deformation gradient summary figure produced by SSLIP
                if IDoptions.plotDefGrad
                    defGrad_fig = gcf;
                    exportScaledFigure(defGrad_fig, fullfile(outDir, [plotBaseName, '_DefGrad.jpg']), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    close(defGrad_fig);
                    fprintf('   -> Saved SSLIP plotDefGrad figure: %s\n', [plotBaseName, '_DefGrad.jpg']);
                end
                
                % Calculate robust quality metrics (Combined Error Index & Correlation Rate)
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
                
                % Save numerical MAT dataset immediately after solving
                save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId', 'cfg', 'm_res', 'corr_rate');
            end
            
            % Identify top active slip systems
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
            
            % Export CSV Table of Active Systems (Enforce column vectors)
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
            
            % Generate and export lightweight JPG Figures
            close all;
            optOut.saveFig = 0;
            plotSSLIP(ebsdID.prop.slipIDcor, ebsdID.prop.residualEeff, ebsdID, sSLocal_filt, optOut);
            
            figs = findall(0, 'Type', 'figure');
            [~, s_idx_f] = sort([figs.Number]);
            figs = figs(s_idx_f);
            
            if ~isempty(figs)
                exportScaledFigure(figs(1), [optOut.plotname, '_SlipActivities.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
            end
            if length(figs) >= 2
                exportScaledFigure(figs(2), [optOut.plotname, '_Residual.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
            end
            close all;
            
            % Record summary statistics
            summary_results(end+1).grainId   = grainId; %#ok<AGROW>
            summary_results(end).numPixel    = grains(grainId).numPixel;
            summary_results(end).comb_index  = m_res;
            summary_results(end).corr_rate   = corr_rate;
            summary_results(end).total_gamma = sum(mean_acts);
            summary_results(end).top_sys     = top_str;
            summary_results(end).baseName    = plotBaseName;
            
            fprintf('   -> Complete! Combined Error Index: %.6f | Correlation Rate: %.1f%%\n', m_res, corr_rate);
        catch ME
            fprintf(2, '   [ERROR] SSLIP processing failed for Grain %d in increment %s:\n%s\n', grainId, prefix, getReport(ME));
            for stk_i = 1:length(ME.stack)
                fprintf(2, '      in %s (line %d)\n', ME.stack(stk_i).name, ME.stack(stk_i).line);
            end
        end
    end

    % Remove empty root initialization if present
    if ~isempty(summary_results) && isfield(summary_results, 'grainId') && isempty(summary_results(1).grainId)
        summary_results(1) = [];
    end

    if ~isempty(summary_results)
        % 1. Print Leaderboard to Terminal
        fprintf('\n     =================================================================================================================================\n');
        fprintf('     STEP 3 TO 4 SSLIP SUMMARY LEADERBOARD - INCREMENT [%s] (Ranked by Combined Error Index)\n', prefix);
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
        fprintf(fid, '# Step 3 to 4 Standalone SSLIP Analysis Report\n\n');
        fprintf(fid, '- **DIC Increment:** `%s`\n', prefix);
        fprintf(fid, '- **Slip System Mode:** `%s` (%s)\n', cfg.slip_system_mode, mode_str);
        fprintf(fid, '- **Schmid Factor Threshold:** `%.2f`\n', cfg.threshold_SF);
        fprintf(fid, '- **Residual Threshold:** `%.3f`\n\n', cfg.threshResidual);
        
        fprintf(fid, '## Analyzed Grain Map & Displacement Gradient Overlays\n');
        fprintf(fid, '### Grain Boundary Hierarchy\n');
        fprintf(fid, '![Analyzed Grains](Step3to4_GrainMap.jpg)\n\n');
        fprintf(fid, '### Displacement Gradients vs. Grain Boundaries\n');
        fprintf(fid, '#### $H_{xx}$ (du/dx)\n');
        fprintf(fid, '![Hxx](%s_DispGradient_Hxx.jpg)\n\n', prefix);
        fprintf(fid, '#### $H_{xy}$ (du/dy)\n');
        fprintf(fid, '![Hxy](%s_DispGradient_Hxy.jpg)\n\n', prefix);
        fprintf(fid, '#### $H_{yx}$ (dv/dx)\n');
        fprintf(fid, '![Hyx](%s_DispGradient_Hyx.jpg)\n\n', prefix);
        fprintf(fid, '#### $H_{yy}$ (dv/dy)\n');
        fprintf(fid, '![Hyy](%s_DispGradient_Hyy.jpg)\n\n---\n\n', prefix);
        
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
fprintf('ALL STEP 3 TO 4 SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
