% =========================================================================
% undeformedReference_06d_SSLIP_Tilt_Analysis.m
% =========================================================================
% SYSTEMATIC CRYSTAL ORIENTATION TILT SENSITIVITY FOR TWIN DECOUPLED SSLIP
%
% This tool applies systematic angular tilts (-15° to +15° in 10° increments
% along both xvector and yvector, yielding a 16-run parameter grid) to the 
% primary parent crystal orientation before performing multiplicative twin 
% shear decoupling and secondary SSLIP multi-slip optimization.
%
% Designed to address ambiguities where residual error alone is insufficient 
% to assess optimization quality, this script outputs high-clarity, lightweight
% JPG activity maps for ALL active slip systems across every tilt combination,
% along with a 2D error contour heatmap and an interactive Markdown Report.
%
% Author: Philipp / Antigravity
% Date:   2026-08-04
% =========================================================================

clear all; clc; close all;
warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
warning('off', 'all');

% Ensure consistent MTEX plotting conventions
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 0. User Configuration Toggles
cfg = struct();

% --- ANALYSIS NAMING & OUTPUT ROUTING ---
cfg.analysis_name        = 'Tilt_Sensitivity_Grid';

% --- TILT GRID & SERIES PARAMETERS ---
cfg.tilt_mode            = 'grid';     % 'grid', 'x_series' (along xvector), 'y_series' (along yvector), or 'cross_series'
cfg.tilt_x_range         = -15:10:15; % Angles in degrees (-15, -5, +5, +15)
cfg.tilt_y_range         = -15:10:15; % Angles in degrees (-15, -5, +5, +15)

% --- TWIN DECOUPLING PARAMETERS ---
cfg.target_variants      = [6];     % Specify twin variant indices to evaluate
cfg.shear_scaling_factor = 1.0;     % Twin shear scaling (1.0 = theoretical)
cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
cfg.process_increments   = [4];     % Indices of DIC increments to process

% --- SSLIP CONTROLS & EXPORT OPTIMIZATION ---
cfg.slip_system_mode     = 'a_only';   % Options: 'all', 'a_only', 'ca_only', 'basal_prism'
cfg.threshold_SF         = 0.01;    % Minimum Schmid factor threshold
cfg.stress_vector        = yvector; % Uniaxial stress direction for SF calculation
cfg.sslip_method         = 1;       % Mode 1: Combined multi-slip minimization
cfg.posConstr            = 1;       % Enforce positive shear constraint
cfg.threshResidual       = 0.01;    % Residual threshold
cfg.minEeff              = 0.01;    % Minimum effective shear strain required for evaluation
cfg.clim_percentile      = 99;      % Robust upper limit percentile for plotting
cfg.jpg_target_px        = 1200;    % Reduced resolution for crisp yet lightweight JPG exports
cfg.jpg_dpi              = 150;     % Reduced DPI for fast file saving and compact storage

%% 1. Load Undeformed Reference EBSD Dataset & Setup Output Routing
try
    [~, defaults] = a00_Dataset_Registry();
    inDir = fullfile(defaults.dir, defaults.outDir_undeformed);
    baseDir = defaults.dir;
catch
    baseDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2\';
    inDir = fullfile(baseDir, 'Processed_Output_Undeformed');
    defaults.thresh_ci = 0.10;
end

matFiles = dir(fullfile(inDir, '*_undeformedEBSD_withAllDIC.mat'));
if isempty(matFiles)
    error('No aligned reference .mat file found in %s. Run undeformedReference_03_IMPORTER first.', inDir);
end
[~, idx] = max([matFiles.datenum]);
matname = matFiles(idx).name;
fprintf('Loading undeformed reference dataset:\n  %s\n', fullfile(inDir, matname));
dataLoaded = load(fullfile(inDir, matname));
newEBSD = dataLoaded.newEBSD;
if isfield(dataLoaded, 'defaults')
    defaults = dataLoaded.defaults;
end

% Setup dedicated analysis output directory
outDir = fullfile(inDir, ['SSLIP_', cfg.analysis_name]);
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created dedicated output directory:\n  %s\n', outDir);
end

%% 2. Retrieve Grains & Identify Primary Parent Grain
grainThresh = 5*degree;
if isfield(dataLoaded, 'grains')
    grains = dataLoaded.grains;
else
    newEBSD(newEBSD.ci < defaults.thresh_ci) = 'notIndexed';
    [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', cfg.minSize, 'angle', grainThresh);
end

valid_grains = grains(grains.numPixel > cfg.minSize);
[~, maxIdx] = max(valid_grains.numPixel);
parentGrain = valid_grains(maxIdx);
originalParentOri = parentGrain.meanOrientation;

fprintf('Primary Parent Matrix identified: Grain ID %d (%d pixels)\n', parentGrain.id, parentGrain.numPixel);

%% 3. Define {10-12} Tensile Twinning and HCP Slip Systems
CS = newEBSD.CS;
tS_1012 = twinSystem.hexagonal_1012(CS);
tS_sym  = tS_1012.symmetrise;

% Instantiate candidate slip system families for HCP Titanium
sS_basal = slipSystem.basal(CS);
sS_prism = slipSystem.prismaticA(CS);
sS_pyrA  = slipSystem.pyramidalA(CS);
sS_pyrCA = slipSystem.pyramidalCA(CS);

switch lower(cfg.slip_system_mode)
    case 'a_only'
        raw_sS = [sS_basal; sS_prism; sS_pyrA];
        fprintf('  -> Slip Mode [%s]: Only <a> type systems (Basal, Prism, Pyr <a>)\n', cfg.slip_system_mode);
    case 'ca_only'
        raw_sS = sS_pyrCA;
        fprintf('  -> Slip Mode [%s]: Only <c+a> type systems (Pyramidal <c+a>)\n', cfg.slip_system_mode);
    case 'basal_prism'
        raw_sS = [sS_basal; sS_prism];
        fprintf('  -> Slip Mode [%s]: Basal <a> and Prismatic <a> systems only\n', cfg.slip_system_mode);
    otherwise
        raw_sS = [sS_basal; sS_prism; sS_pyrA; sS_pyrCA];
        fprintf('  -> Slip Mode [all]: ALL systems (Basal <a>, Prismatic <a>, Pyramidal <a>, Pyramidal <c+a>)\n');
end
sS = symmetrise(raw_sS, 'antipodal');

%% 4. Discover DIC Increments
all_fields = fields(newEBSD.prop);
displacement_fields = all_fields(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true));
prefix_list = unique(cellfun(@(x) x(1:find(x=='_', 1)-1), displacement_fields, 'UniformOutput', false));

if isempty(cfg.process_increments)
    inc_indices = 1:length(prefix_list);
else
    inc_indices = intersect(cfg.process_increments, 1:length(prefix_list));
end

%% 5. EXECUTION LOOP: 2D TILT SENSITIVITY SWEEP
for inc_idx = inc_indices
    prefix = prefix_list{inc_idx};
    fprintf('\n========================================================================\n');
    fprintf('EVALUATING DIC INCREMENT [%s] FOR TILT SENSITIVITY\n', prefix);
    fprintf('========================================================================\n');
    
    ux_fname = [prefix, '_Displacements_x_displacement'];
    uy_fname = [prefix, '_Displacements_y_displacement'];
    if ~isfield(newEBSD.prop, ux_fname) || ~isfield(newEBSD.prop, uy_fname)
        warning('Missing displacement fields for %s. Skipping.', prefix);
        continue;
    end
    
    hxx_field = [prefix, '_Gradients_du_dx_RC'];
    hxy_field = [prefix, '_Gradients_du_dy_RC'];
    hyx_field = [prefix, '_Gradients_dv_dx_RC'];
    hyy_field = [prefix, '_Gradients_dv_dy_RC'];
    
    if isfield(newEBSD.prop, hxx_field)
        Hxx_raw = double(newEBSD.prop.(hxx_field));
        Hxy_raw = double(newEBSD.prop.(hxy_field));
        Hyx_raw = double(newEBSD.prop.(hyx_field));
        Hyy_raw = double(newEBSD.prop.(hyy_field));
    else
        error('No displacement gradients found in the imported ebsd object for increment %s.', prefix);
    end
    
    F_tot_11 = 1 + Hxx_raw; F_tot_12 = Hxy_raw;
    F_tot_21 = Hyx_raw;     F_tot_22 = 1 + Hyy_raw;
    
    for v = cfg.target_variants
        twinMis = tS_sym(v).parentTwinMisorientation;
        
        for factor_val = cfg.shear_scaling_factor(:)'
            % Setup 2D Tilt Grid or 1D Tilt Series
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
            tilt_results = repmat(struct(), num_runs, 1);
            res_grid     = zeros(size(TX));
            
            fprintf('\n   ---------------------------------------------------------------------\n');
            fprintf('   Executing %d-Run Crystal Tilt Grid (-15° to +15°) | Variant %d | Scale %.2f\n', num_runs, v, factor_val);
            fprintf('   ---------------------------------------------------------------------\n');
            
            for run_i = 1:num_runs
                tx = TX(run_i);
                ty = TY(run_i);
                
                fprintf('   ---> [Run %02d/%02d] Evaluating Tilt X = %+03d°, Tilt Y = %+03d°...\n', run_i, num_runs, tx, ty);
                
                % 1. Apply 3D crystal orientation rotation
                rot_x = rotation.byAxisAngle(xvector, tx * degree);
                rot_y = rotation.byAxisAngle(yvector, ty * degree);
                rot_tilt = rot_y * rot_x;
                
                tiltedParentOri   = rot_tilt * originalParentOri;
                parentOri_twinned = tiltedParentOri * twinMis;
                
                % 2. Recalculate twin deformation gradients under tilted lattice frame
                H_twin_global = tiltedParentOri * tS_sym.displacementGradient;
                H_tv          = H_twin_global(v).M(1:2, 1:2);
                F_twin_2d     = eye(2) + factor_val * H_tv;
                F_twin_inv    = inv(F_twin_2d);
                
                a11 = F_twin_inv(1,1); a12 = F_twin_inv(1,2);
                a21 = F_twin_inv(2,1); a22 = F_twin_inv(2,2);
                
                F_slip_11 = F_tot_11 .* a11 + F_tot_12 .* a21;
                F_slip_12 = F_tot_11 .* a12 + F_tot_12 .* a22;
                F_slip_21 = F_tot_21 .* a11 + F_tot_22 .* a21;
                F_slip_22 = F_tot_21 .* a12 + F_tot_22 .* a22;
                
                newEBSD.prop.Hxx_sub = F_slip_11 - 1;
                newEBSD.prop.Hxy_sub = F_slip_12;
                newEBSD.prop.Hyx_sub = F_slip_21;
                newEBSD.prop.Hyy_sub = F_slip_22 - 1;
                
                % 3. Setup grain data slice with tilted orientation
                grainId = parentGrain.id;
                ebsd_grain = newEBSD(grains(grainId)).gridify;
                ebsd_grain.orientations = repmat(parentOri_twinned, size(ebsd_grain.orientations));
                
                U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos;
                V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
                
                % 4. Rotate slip systems and apply Schmid factor filter
                sSLocal = parentOri_twinned * sS;
                stressState = stressTensor.uniaxial(cfg.stress_vector);
                
                sf_signs = sign(sSLocal.SchmidFactor(stressState));
                sf_signs(sf_signs == 0) = 1;
                sSLocal.n = sf_signs .* sSLocal.n;
                
                if cfg.threshold_SF > 0
                    SF_crit = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
                    sSLocal_filt = sSLocal(SF_crit);
                    sS_filt      = sS(SF_crit);
                else
                    sSLocal_filt = sSLocal;
                    sS_filt      = sS;
                end
                
                if isempty(sS_filt)
                    warning('No slip systems meeting SF threshold for Run %02d. Skipping.', run_i);
                    continue;
                end
                
                % 5. Configure SSLIP optimization (disable slow internal plotting)
                tag_str   = sprintf('TiltX_%+03d_TiltY_%+03d', tx, ty);
                clean_tag = strrep(strrep(tag_str, '+', 'p'), '-', 'm');
                case_name = sprintf('Run%02d_Var%d_Factor%03d_%s', run_i, v, round(factor_val*100), clean_tag);
                
                clear IDoptions;
                IDoptions.IDMethod       = cfg.sslip_method;
                IDoptions.filterSize     = 4;
                IDoptions.coarsegrain    = 4;
                IDoptions.threshResidual = cfg.threshResidual;
                IDoptions.minEeff        = cfg.minEeff;
                IDoptions.posConstr      = cfg.posConstr;
                IDoptions.plotSSLIP      = 0; % Disabled to execute fast
                IDoptions.saveFig        = 0; % Disabled to output lightweight JPGs directly
                IDoptions.layout         = [];
                IDoptions.plotTraces     = 1;
                IDoptions.stress         = stressState;
                IDoptions.casename       = case_name;
                IDoptions.plotDefGrad    = 0;
                IDoptions.sizeAdjust     = 0.5;
                IDoptions.plotResidual   = 1;
                IDoptions.logscale       = 0;
                IDoptions.clim_percentile = cfg.clim_percentile;
                try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
                IDoptions.comment        = sprintf('Tilt: X=%+03d°, Y=%+03d° | Var %d | Scale %.2f', tx, ty, v, factor_val);
                IDoptions.NoSs           = 1:length(sS_filt);
                IDoptions.plotname       = fullfile(outDir, case_name);
                
                IDoptions.Hxx = ebsd_grain.prop.Hxx_sub;
                IDoptions.Hxy = ebsd_grain.prop.Hxy_sub;
                IDoptions.Hyx = ebsd_grain.prop.Hyx_sub;
                IDoptions.Hyy = ebsd_grain.prop.Hyy_sub;
                
                % Execute SSLIP solver
                [ebsdID, optOut] = SSLIP(ebsd_grain, U, V, sSLocal_filt, IDoptions);
                
                % 6. Manually trigger plotSSLIP and save crisp lightweight JPGs
                close all;
                optOut.saveFig = 0;
                plotSSLIP(ebsdID, sSLocal_filt, optOut);
                
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
                
                % 7. Record robust quality metrics (Combined Error Index & Correlation Rate)
                % Scales actual mean residual by inverse correlation fraction to avoid threshold domination
                % while remaining sensitive to variations in fit quality and unexplained map area.
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
                res_grid(run_i) = m_res;
                
                mean_acts = mean(ebsdID.prop.slipIDcor, 1, 'omitnan');
                [sort_acts, sort_idx] = sort(mean_acts, 'descend');
                
                top_str = '';
                for t_i = 1:min(3, length(sort_idx))
                    sys_idx = sort_idx(t_i);
                    sys_obj = sS_filt(sys_idx);
                    
                    % Classify family name
                    if sys_obj.b.norm > 4.0
                        fam_str = 'Pyr <c+a>';
                    elseif abs(sys_obj.n.z) > 0.95
                        fam_str = 'Basal <a>';
                    elseif abs(sys_obj.n.z) < 0.05
                        fam_str = 'Prism <a>';
                    else
                        fam_str = 'Pyr <a>';
                    end
                    
                    top_str = [top_str, sprintf('%s (%s) [SF=%.2f, γ=%.4f]; ', fam_str, sys_obj.b.char, abs(sSLocal_filt(sys_idx).SchmidFactor(stressState)), sort_acts(t_i))];
                end
                
                tilt_results(run_i).run_id    = run_i;
                tilt_results(run_i).tx        = tx;
                tilt_results(run_i).ty        = ty;
                tilt_results(run_i).mean_res  = m_res;
                tilt_results(run_i).corr_rate = corr_rate;
                tilt_results(run_i).total_gam = sum(mean_acts);
                tilt_results(run_i).top_sys   = top_str;
                tilt_results(run_i).case_name = case_name;
            end
            
            % --- SORT AND SAVE DIAGNOSTICS ---
            all_res = [tilt_results.mean_res];
            [~, sort_rank] = sort(all_res, 'ascend');
            ranked_results = tilt_results(sort_rank);
            
            % Print Leaderboard to Terminal
            fprintf('\n     =================================================================================================================================\n');
            fprintf('     CRYSTAL TILT SENSITIVITY LEADERBOARD (Ranked by Combined Error Index)\n');
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
            csv_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_Tilt_Leaderboard.csv', prefix, parentGrain.id, v));
            T_csv = table((1:length(ranked_results))', [ranked_results.run_id]', [ranked_results.tx]', [ranked_results.ty]', ...
                          [ranked_results.mean_res]', [ranked_results.corr_rate]', [ranked_results.total_gam]', {ranked_results.top_sys}', ...
                          'VariableNames', {'Rank', 'Run_ID', 'Tilt_X_deg', 'Tilt_Y_deg', 'Combined_Error_Index', 'Correlation_Rate_Pct', 'Total_Mean_Shear', 'Dominant_Slip_Systems'});
            writetable(T_csv, csv_path);
            
            % Save MAT dataset
            mat_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_Tilt_Data.mat', prefix, parentGrain.id, v));
            save(mat_path, 'tilt_results', 'ranked_results', 'TX', 'TY', 'res_grid', 'cfg', 'v');
            
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
                title(sprintf('SSLIP Combined Quality Index vs Crystal Tilt - Grain %d (Var %d)', parentGrain.id, v), 'FontSize', 13, 'FontWeight', 'bold');
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
                    title_str = sprintf('SSLIP Quality Index vs X-Vector Tilt Series - Grain %d (Var %d)', parentGrain.id, v);
                elseif all(TX(:) == 0)
                    plot(TY(:), res_grid(:), '-or', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'y');
                    xlabel('Y-Vector Tilt (\Delta\theta_y) [deg]', 'FontSize', 12, 'FontWeight', 'bold');
                    title_str = sprintf('SSLIP Quality Index vs Y-Vector Tilt Series - Grain %d (Var %d)', parentGrain.id, v);
                else
                    plot(1:numel(res_grid), res_grid(:), '-og', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
                    xlabel('Evaluation Run Index', 'FontSize', 12, 'FontWeight', 'bold');
                    title_str = sprintf('SSLIP Quality Index across Cross-Series Tilts - Grain %d (Var %d)', parentGrain.id, v);
                end
                ylabel('Combined Error Index (r_{corr} / \eta_{corr})', 'FontSize', 12, 'FontWeight', 'bold');
                title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
                grid on; hold on;
                [min_val, min_lin] = min(res_grid(:));
                plot(min_lin, min_val, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
                hold off;
            end
            
            heatmap_fname = sprintf('%s_Grain%d_Var%d_Tilt_Heatmap.png', prefix, parentGrain.id, v);
            exportScaledFigure(figHeat, fullfile(outDir, heatmap_fname), 'nodatestring');
            close(figHeat);
            
            % --- GENERATE INTERACTIVE MARKDOWN REPORT ---
            md_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_Tilt_Report.md', prefix, parentGrain.id, v));
            fprintf('     -> Writing interactive Markdown Report with embedded activity maps to:\n        %s\n', md_path);
            
            fid = fopen(md_path, 'w');
            fprintf(fid, '# Crystal Orientation Tilt Sensitivity Report (Twin Decoupled SSLIP)\n\n');
            fprintf(fid, '- **DIC Increment:** `%s`\n', prefix);
            fprintf(fid, '- **Parent Grain ID:** `%d`\n', parentGrain.id);
            fprintf(fid, '- **Twin Variant:** `%d`\n', v);
            fprintf(fid, '- **Shear Scaling Factor:** `%.2f`\n', factor_val);
            fprintf(fid, '- **Tilt Sweep Bounds:** X $\\in [%d^\\circ, %d^\\circ]$, Y $\\in [%d^\\circ, %d^\\circ]$ ($%d^\\circ$ steps)\n\n', ...
                    min(cfg.tilt_x_range), max(cfg.tilt_x_range), min(cfg.tilt_y_range), max(cfg.tilt_y_range), cfg.tilt_x_range(2)-cfg.tilt_x_range(1));
            
            fprintf(fid, '## 2D Combined Error Index Sensitivity Surface\n');
            fprintf(fid, 'This heatmap illustrates how variations in crystal mounting angle alter the **Combined Error Index** ($r_{\\text{comb}} = \\frac{\\text{mean correlated residual}}{\\text{correlation fraction}}$) of the multi-slip optimization, ensuring sensitive quality discrimination across different mounting alignments.\n\n');
            fprintf(fid, '![Tilt Sensitivity Heatmap](%s)\n\n---\n\n', heatmap_fname);
            
            fprintf(fid, '## Tilt Configuration Leaderboard\n');
            fprintf(fid, 'Ranked by lowest Combined Error Index and reporting overall Correlation Rate (%% of evaluated pixels below residual threshold).\n\n');
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
                fprintf(fid, '![Activity Maps %s](%s_SlipActivities.jpg)\n\n', c_name, c_name);
                
                fprintf(fid, '#### Residual Error Map\n');
                fprintf(fid, '![Residual Map %s](%s_Residual.jpg)\n\n', c_name, c_name);
                fprintf(fid, '---\n\n');
            end
            fclose(fid);
            
            fprintf('     -> Complete! Saved all %d tilt maps, CSV leaderboard, Heatmap, and Markdown Report to:\n        %s\n', num_runs, outDir);
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL CRYSTAL TILT SENSITIVITY RUNS COMPLETE!\n');
fprintf('========================================================================\n');
