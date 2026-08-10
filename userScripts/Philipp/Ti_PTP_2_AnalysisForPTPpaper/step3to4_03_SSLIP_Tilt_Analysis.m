% =========================================================================
% step3to4_03_SSLIP_Tilt_Analysis.m
% =========================================================================
% SYSTEMATIC CRYSTAL ORIENTATION TILT SENSITIVITY FOR STEP 3 TO 4 SSLIP
%
% This script evaluates sensitivity of the multi-slip optimization to minor
% experimental variations in crystal mounting angle. By applying systematic
% 3D lattice rotations (-15° to +15° along both xvector and yvector) to
% selected grains, it determines the orientation alignment that minimizes
% the Combined Error Index.
%
% Features:
% - Grain ID filtering via cfg.restrict2GrainId (e.g., focus on Grain 5)
% - Instant MAT caching of solved optimization states for fast re-runs
% - Robust 2D contour heatmap visualization of quality vs. tilt angle
% - Automated Markdown diagnostic report generation with embedded activity maps
%
% Author: Philipp / Antigravity
% Date:   2026-08-05
% =========================================================================

clear all; clc; close all;
warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
warning('off', 'all');

% Ensure consistent MTEX plotting conventions
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 0. User Configuration & Toggles
cfg = struct();

% --- ANALYSIS NAMING & TILT GRID CONFIGURATION ---
cfg.analysis_name    = 'Step3to4_Tilt_Sensitivity_Grid';
cfg.tilt_mode        = 'grid';        % 'grid', 'x_series' (along xvector), 'y_series' (along yvector), or 'cross_series'
cfg.tilt_x_range     = -7:7:7;        % Angles in degrees (-15, -5, +5, +15)
cfg.tilt_y_range     = -7:7:7;        % Angles in degrees (-15, -5, +5, +15)
cfg.restrict2GrainId = 5;             % 0 = all, or array of specific ID(s) like [5]

% --- SSLIP CONTROLS & EXPORT PARAMETERS ---
cfg.slip_system_mode = 'all';         % 'all', 'a_only', 'ca_only', or 'basal_pyr'
cfg.threshold_SF     = 0.01;          % Minimum Schmid Factor threshold
cfg.stress_vector    = yvector;       % Uniaxial stress loading axis
cfg.sslip_method     = 1;             % IDMethod (1 = multi-slip optimization)
cfg.posConstr        = 1;             % Enforce non-negative shear activity
cfg.threshResidual   = 0.015;         % Residual threshold for correlation rate
cfg.minEeff          = 0.01;          % Minimum deformation threshold for evaluation
cfg.minSize          = 300;           % Minimum grain size in pixels to process
cfg.clim_percentile  = 99;            % Robust percentile limits for color scaling
cfg.jpg_target_px    = 1200;          % Target resolution width for exports
cfg.jpg_dpi          = 150;           % Target DPI for lightweight image outputs

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

%% 2. Retrieve or Reconstruct Grains & Apply Filtering
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
    fprintf('Restricting tilt grid analysis to %d specified grain(s): [%s]\n', length(grains_ids), num2str(grains_ids(:)'));
else
    fprintf('Analyzing all %d significant grains...\n', length(grains_ids));
end

%% 3. Extract DIC Loading Increments
all_fields = fields(newEBSD.prop);
disp_logical = contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true);
disp_fields  = all_fields(disp_logical);

if isempty(disp_fields)
    error('No DIC displacement fields found in dataset.');
end

prefix_list = cell(length(disp_fields), 1);
for p_i = 1:length(disp_fields)
    idx_delim = strfind(disp_fields{p_i}, '_Displacements');
    prefix_list{p_i} = disp_fields{p_i}(1:idx_delim-1);
end
prefix_list = unique(prefix_list);

%% 4. Define Candidate Slip Systems based on Mode
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
        mode_str = 'Only <a> Type Systems';
    case {'ca', 'c+a', 'ca_only', 'only_ca'}
        sS_raw = sS_ca_all;
        mode_str = 'Only <c+a> Type Systems';
    case {'basal_pyramidal', 'basal_pyr', 'only_basal_pyr'}
        sS_raw = [sS_basal_all; sS_pyra_all; sS_ca_all];
        mode_str = 'Basal and Pyramidal Only (No Prism)';
    otherwise
        sS_raw = [sS_basal_all; sS_prism_all; sS_pyra_all; sS_ca_all];
        mode_str = 'Default: All Systems';
end

sS = symmetrise(sS_raw, 'antipodal');
fprintf('Slip System Mode [%s]: %d candidate systems initialized.\n', cfg.slip_system_mode, length(sS));

%% 5. Execute Tilt Sensitivity Grid / Series
if ~isfield(cfg, 'tilt_mode'), cfg.tilt_mode = 'grid'; end
switch lower(cfg.tilt_mode)
    case {'x', 'x_series', 'xvector'}
        TX_vec = cfg.tilt_x_range(:); TY_vec = zeros(size(TX_vec));
        TX = TX_vec'; TY = TY_vec';
    case {'y', 'y_series', 'yvector'}
        TY_vec = cfg.tilt_y_range(:); TX_vec = zeros(size(TY_vec));
        TX = TX_vec'; TY = TY_vec';
    case {'cross', 'cross_series'}
        tx_part = cfg.tilt_x_range(:); ty_zero = zeros(size(tx_part));
        ty_part = cfg.tilt_y_range(cfg.tilt_y_range ~= 0); ty_part = ty_part(:); tx_zero = zeros(size(ty_part));
        TX = [tx_part; tx_zero]'; TY = [ty_zero; ty_part]';
    otherwise % 'grid', '2d'
        [TX, TY] = meshgrid(cfg.tilt_x_range, cfg.tilt_y_range);
end
num_runs = numel(TX);

for i_i = 1:length(prefix_list)
    prefix = prefix_list{i_i};
    ux_fname = [prefix, '_Displacements_x_displacement'];
    uy_fname = [prefix, '_Displacements_y_displacement'];
    
    for g_i = 1:length(grains_ids)
        grainId = grains_ids(g_i);
        originalGrainOri = grains(grainId).meanOrientation;
        
        fprintf('\n========================================================================\n');
        fprintf('EVALUATING TILT GRID: INCREMENT [%s] | GRAIN %d (%d px)\n', prefix, grainId, grains(grainId).numPixel);
        fprintf('Executing %d-Run Crystal Tilt Grid (X: %d° to %d°, Y: %d° to %d°)\n', ...
                num_runs, min(cfg.tilt_x_range), max(cfg.tilt_x_range), min(cfg.tilt_y_range), max(cfg.tilt_y_range));
        fprintf('========================================================================\n');
        
        tilt_results = repmat(struct(), num_runs, 1);
        res_grid     = zeros(size(TX));
        
        for run_i = 1:num_runs
            tx = TX(run_i);
            ty = TY(run_i);
            
            fprintf('   ---> [Run %02d/%02d] Evaluating Tilt X = %+03d°, Tilt Y = %+03d°...\n', run_i, num_runs, tx, ty);
            
            % 1. Apply 3D crystal orientation rotation
            rot_x = rotation.byAxisAngle(xvector, tx * degree);
            rot_y = rotation.byAxisAngle(yvector, ty * degree);
            rot_tilt = rot_y * rot_x;
            tiltedGrainOri = rot_tilt * originalGrainOri;
            
            % 2. Setup grain data slice with tilted orientation
            ebsd_grain = newEBSD(grains(grainId)).gridify;
            ebsd_grain.orientations = rot_tilt * ebsd_grain.orientations;
            
            if ~isfield(ebsd_grain.prop, ux_fname) || ~isfield(ebsd_grain.prop, uy_fname)
                warning('   Missing displacement data %s for grain %d. Skipping.', ux_fname, grainId);
                continue;
            end
            
            U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos; 
            V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
            
            % 3. Rotate slip systems and apply Schmid factor filter
            sSLocal = tiltedGrainOri * sS;
            stressState = stressTensor.uniaxial(cfg.stress_vector);
            
            sf_signs = sign(sSLocal.SchmidFactor(stressState));
            sf_signs(sf_signs == 0) = 1;
            sSLocal.n = sf_signs .* sSLocal.n;
            
            SF_crit = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
            sSLocal_filt = sSLocal(SF_crit);
            sS_filt      = sS(SF_crit);
            
            if isempty(sS_filt)
                warning('   No slip systems with SF >= %.2f for Run %02d. Skipping.', cfg.threshold_SF, run_i);
                continue;
            end
            
            % 4. Configure SSLIP optimization
            tag_str   = sprintf('TiltX_%+03d_TiltY_%+03d', tx, ty);
            clean_tag = strrep(strrep(tag_str, '+', 'p'), '-', 'm');
            case_name = sprintf('%s_Grain%02d_Run%02d_%s', prefix, grainId, run_i, clean_tag);
            
            clear IDoptions;
            IDoptions.IDMethod    = cfg.sslip_method;
            IDoptions.filterSize  = 4;
            IDoptions.coarsegrain = 4;
            IDoptions.threshResidual = cfg.threshResidual;
            IDoptions.minEeff     = cfg.minEeff;
            IDoptions.posConstr   = cfg.posConstr;
            IDoptions.plotSSLIP   = 0;
            IDoptions.layout      = [];
            IDoptions.plotTraces  = 1;
            IDoptions.stress      = stressState;
            IDoptions.casename    = case_name;
            IDoptions.plotDefGrad = 1;
            IDoptions.saveFig     = 0;
            IDoptions.plotname    = fullfile(outDir, case_name);
            IDoptions.sizeAdjust  = 0.5;
            IDoptions.plotResidual = 1;
            IDoptions.logscale    = 1;
            IDoptions.logmin      = 0.01;
            try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
            IDoptions.comment     = sprintf('Grain %d [%s | %s]', grainId, cfg.slip_system_mode, tag_str);
            IDoptions.NoSs        = 1:length(sS_filt);
            
            resMatPath = fullfile(outDir, [case_name, '_SSLIP_results.mat']);
            
            try
                if isfile(resMatPath)
                    fprintf('      [Cache Hit] Loading pre-computed optimization from MAT file...\n');
                    cached_data  = load(resMatPath);
                    ebsdID       = cached_data.ebsdID;
                    optOut       = cached_data.optOut;
                    sSLocal_filt = cached_data.sSLocal_filt;
                    m_res        = cached_data.m_res;
                    corr_rate    = cached_data.corr_rate;
                else
                    fprintf('      Running multi-slip optimization (%d systems)...\n', length(sS_filt));
                    [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);
                    
                    if IDoptions.plotDefGrad
                        defGrad_fig = gcf;
                        exportScaledFigure(defGrad_fig, fullfile(outDir, [case_name, '_DefGrad.jpg']), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                        close(defGrad_fig);
                        fprintf('      -> Saved SSLIP plotDefGrad summary figure: %s\n', [case_name, '_DefGrad.jpg']);
                    end
                    
                    % Quality metrics
                    Eeff_eval = ebsdID.prop.Eeff(:);
                    eval_mask = (Eeff_eval >= cfg.minEeff) & ~isnan(Eeff_eval);
                    res_map   = ebsdID.prop.residualEeff(:);
                    res_eval  = res_map(eval_mask);
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
                    
                    save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId', 'tx', 'ty', 'cfg', 'm_res', 'corr_rate');
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
                writetable(T_sys, fullfile(outDir, [case_name, '_SlipSystems.csv']));
                
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
                
                % Store diagnostics
                res_grid(run_i)           = m_res;
                tilt_results(run_i).run_id    = run_i;
                tilt_results(run_i).tx        = tx;
                tilt_results(run_i).ty        = ty;
                tilt_results(run_i).mean_res  = m_res;
                tilt_results(run_i).corr_rate = corr_rate;
                tilt_results(run_i).total_gam = sum(mean_acts);
                tilt_results(run_i).top_sys   = top_str;
                tilt_results(run_i).case_name = case_name;
                
                fprintf('      -> Complete! Combined Error Index: %.6f | Correlation Rate: %.1f%%\n', m_res, corr_rate);
            catch ME
                fprintf(2, '      [ERROR] Tilt evaluation failed for Run %02d (Grain %d):\n         %s\n', run_i, grainId, ME.message);
                for stk_i = 1:length(ME.stack)
                    fprintf(2, '         in %s (line %d)\n', ME.stack(stk_i).name, ME.stack(stk_i).line);
                end
            end
        end
        
        % --- SORT AND SAVE DIAGNOSTICS FOR GRAIN ---
        all_res = [tilt_results.mean_res];
        [~, sort_rank] = sort(all_res, 'ascend');
        ranked_results = tilt_results(sort_rank);
        
        % Print Leaderboard to Terminal
        fprintf('\n     =================================================================================================================================\n');
        fprintf('     CRYSTAL TILT SENSITIVITY LEADERBOARD - GRAIN %02d (Ranked by Combined Error Index)\n', grainId);
        fprintf('     =================================================================================================================================\n');
        fprintf('     Rank | Run ID | Tilt X | Tilt Y | Combined Error Index | Correlation %% | Total Shear γ | Top Active Slip Systems \n');
        fprintf('     ---------------------------------------------------------------------------------------------------------------------------------\n');
        for r_i = 1:length(ranked_results)
            pr = ranked_results(r_i);
            fprintf('     #%02d  | Run %02d | %+03d°   | %+03d°   |       %.6f       |    %5.1f%%     |    %.4f     | %s \n', ...
                    r_i, pr.run_id, pr.tx, pr.ty, pr.mean_res, pr.corr_rate, pr.total_gam, pr.top_sys);
        end
        fprintf('     =================================================================================================================================\n\n');
        
        % Save to CSV Table
        csv_path = fullfile(outDir, sprintf('%s_Grain%02d_Tilt_Leaderboard.csv', prefix, grainId));
        T_csv = table((1:length(ranked_results))', [ranked_results.run_id]', [ranked_results.tx]', [ranked_results.ty]', ...
                      [ranked_results.mean_res]', [ranked_results.corr_rate]', [ranked_results.total_gam]', {ranked_results.top_sys}', ...
                      'VariableNames', {'Rank', 'Run_ID', 'Tilt_X_deg', 'Tilt_Y_deg', 'Combined_Error_Index', 'Correlation_Rate_Pct', 'Total_Mean_Shear', 'Dominant_Slip_Systems'});
        writetable(T_csv, csv_path);
        
        % Save MAT dataset
        mat_path = fullfile(outDir, sprintf('%s_Grain%02d_Tilt_Data.mat', prefix, grainId));
        save(mat_path, 'tilt_results', 'ranked_results', 'TX', 'TY', 'res_grid', 'cfg', 'grainId');
        
        % --- GENERATE SENSITIVITY VISUALIZATION (2D HEATMAP OR 1D CURVE) ---
        fprintf('     -> Generating Tilt Sensitivity Residual Error Visualization...\n');
        figHeat = figure('Name', 'Tilt Sensitivity Visualization', 'Position', [200, 200, 800, 650], 'Visible', 'off');
        if size(TX, 1) > 1 && size(TX, 2) > 1 && ~isvector(TX)
            contourf(TX, TY, res_grid, 20, 'LineColor', 'k', 'LineWidth', 0.5);
            colormap(parula(256));
            cb = colorbar;
            cb.Label.String = 'Combined Error Index (r_{corr} / \eta_{corr})';
            xlabel('X-Vector Tilt (\Delta\theta_x) [deg]', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Y-Vector Tilt (\Delta\theta_y) [deg]', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('SSLIP Combined Quality Index vs Crystal Tilt Grid - Grain %d', grainId), 'FontSize', 13, 'FontWeight', 'bold');
            grid on; hold on;
            [min_val, min_lin] = min(res_grid(:));
            plot(TX(min_lin), TY(min_lin), 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            text(TX(min_lin), TY(min_lin), sprintf('  Min (X=%+d°, Y=%+d°: %.6f)', TX(min_lin), TY(min_lin), min_val), ...
                 'Color', 'w', 'FontWeight', 'bold', 'BackgroundColor', 'k');
            hold off;
        else
            if all(TY(:) == 0)
                plot(TX(:), res_grid(:), '-ob', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'c');
                xlabel('X-Vector Tilt (\Delta\theta_x) [deg]', 'FontSize', 12, 'FontWeight', 'bold');
                title_str = sprintf('SSLIP Quality Index vs X-Vector Tilt Series - Grain %d', grainId);
            elseif all(TX(:) == 0)
                plot(TY(:), res_grid(:), '-or', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'y');
                xlabel('Y-Vector Tilt (\Delta\theta_y) [deg]', 'FontSize', 12, 'FontWeight', 'bold');
                title_str = sprintf('SSLIP Quality Index vs Y-Vector Tilt Series - Grain %d', grainId);
            else
                plot(1:numel(res_grid), res_grid(:), '-og', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
                xlabel('Evaluation Run Index', 'FontSize', 12, 'FontWeight', 'bold');
                title_str = sprintf('SSLIP Quality Index across Cross-Series Tilts - Grain %d', grainId);
            end
            ylabel('Combined Error Index (r_{corr} / \eta_{corr})', 'FontSize', 12, 'FontWeight', 'bold');
            title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
            grid on; hold on;
            [min_val, min_lin] = min(res_grid(:));
            plot(min_lin, min_val, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            hold off;
        end
        
        heatmap_fname = sprintf('%s_Grain%02d_Tilt_Heatmap.jpg', prefix, grainId);
        exportScaledFigure(figHeat, fullfile(outDir, heatmap_fname), 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
        close(figHeat);
        
        % --- GENERATE INTERACTIVE MARKDOWN REPORT ---
        md_path = fullfile(outDir, sprintf('%s_Grain%02d_Tilt_Report.md', prefix, grainId));
        fprintf('     -> Writing interactive Markdown Report to:\n        %s\n', md_path);
        
        fid = fopen(md_path, 'w');
        fprintf(fid, '# Crystal Orientation Tilt Sensitivity Report (Step 3 to 4 SSLIP)\n\n');
        fprintf(fid, '- **DIC Increment:** `%s`\n', prefix);
        fprintf(fid, '- **Evaluated Grain ID:** `%d` (`%d` pixels)\n', grainId, grains(grainId).numPixel);
        fprintf(fid, '- **Slip System Mode:** `%s` (%s)\n', cfg.slip_system_mode, mode_str);
        fprintf(fid, '- **Tilt Sweep Bounds:** X $\\in [%d^\\circ, %d^\\circ]$, Y $\\in [%d^\\circ, %d^\\circ]$ ($%d^\\circ$ steps)\n\n', ...
                min(cfg.tilt_x_range), max(cfg.tilt_x_range), min(cfg.tilt_y_range), max(cfg.tilt_y_range), cfg.tilt_x_range(2)-cfg.tilt_x_range(1));
        
        fprintf(fid, '## 2D Combined Error Index Sensitivity Surface\n');
        fprintf(fid, 'This heatmap displays how experimental variations in crystal alignment alter the **Combined Error Index** ($r_{\\text{comb}} = \\frac{r_{\\text{corr}}}{\\eta_{\\text{corr}}}$) of the multi-slip optimization.\n\n');
        fprintf(fid, '![Tilt Sensitivity Heatmap](%s)\n\n---\n\n', heatmap_fname);
        
        fprintf(fid, '## Tilt Configuration Leaderboard\n');
        fprintf(fid, 'Ranked by lowest Combined Error Index and reporting overall Correlation Rate.\n\n');
        fprintf(fid, '| Rank | Run ID | Tilt X | Tilt Y | Combined Error Index | Correlation %% | Total $\\gamma$ | Dominant Active Systems |\n');
        fprintf(fid, '| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |\n');
        for r_i = 1:length(ranked_results)
            pr = ranked_results(r_i);
            fprintf(fid, '| **#%02d** | Run %02d | `%+03d°` | `%+03d°` | `%.6f` | `%5.1f%%` | `%.4f` | %s |\n', ...
                    r_i, pr.run_id, pr.tx, pr.ty, pr.mean_res, pr.corr_rate, pr.total_gam, pr.top_sys);
        end
        fprintf(fid, '\n---\n\n');
        
        fprintf(fid, '## Complete Visual Activity Maps Across All Tilt Configurations\n');
        for r_i = 1:length(ranked_results)
            pr = ranked_results(r_i);
            c_name = pr.case_name;
            
            fprintf(fid, '### Rank #%02d: Tilt X = %+03d°, Tilt Y = %+03d° (Run %02d)\n', r_i, pr.tx, pr.ty, pr.run_id);
            fprintf(fid, '- **Combined Error Index:** `%.6f`\n', pr.mean_res);
            fprintf(fid, '- **Correlation Rate:** `%.1f%%`\n', pr.corr_rate);
            fprintf(fid, '- **Total Mean Shear Activity:** `%.4f`\n', pr.total_gam);
            fprintf(fid, '- **Dominant Systems:** %s\n\n', pr.top_sys);
            fprintf(fid, '#### Slip System Activity Maps\n');
            fprintf(fid, '![Slip Activities](%s_SlipActivities.jpg)\n\n', c_name);
            fprintf(fid, '#### Residual Error Field\n');
            fprintf(fid, '![Residual Error](%s_Residual.jpg)\n\n---\n\n', c_name);
        end
        fclose(fid);
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL STEP 3 TO 4 TILT SENSITIVITY ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
