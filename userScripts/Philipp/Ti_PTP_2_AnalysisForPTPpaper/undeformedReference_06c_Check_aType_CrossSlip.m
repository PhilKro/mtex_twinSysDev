%% undeformedReference_06c_Check_aType_CrossSlip.m
% Exhaustive <a>-Type Colinear Burgers Vector Cross-Slip Diagnostics
% -------------------------------------------------------------------------
% This script tests the physical hypothesis of localized <a>-type cross-slip
% by evaluating SSLIP in Mode 1 (combined multi-slip optimization), but strictly
% restricting the candidate slip systems to colinear-b pairs!
%
% Workflow & Physics:
%   1. Applies multiplicative twin shear decoupling (F_slip = F_tot * F_twin^-1)
%      in the twinned crystal lattice orientation.
%   2. Instantiates all 12 symmetrized <a> type slip systems (3 Basal, 3 Prism,
%      6 Pyramidal <a>).
%   3. Automatically clusters systems by identical Burgers vector line direction
%      (3 unique <a> Burgers vectors, each possessing exactly 4 glide planes:
%      1 Basal + 1 Prism + 2 Pyramidal).
%   4. Generates all 18 possible colinear-b cross-slip pairs (nchoosek(1:4,2) * 3).
%   5. Evaluates SSLIP Mode 1 separately for each pair, logs the residual error,
%      ranks all cross-slip combinations from best to worst fit, and outputs
%      detailed diagnostic comparison plots for the winning combinations!
% -------------------------------------------------------------------------

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\geometry');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\SSLIP\src');

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 0. User Configuration Toggles
cfg = struct();

% --- ANALYSIS NAMING & OUTPUT ROUTING ---
cfg.analysis_name        = 'Check_aType_CrossSlip_Pairs';

% --- TWIN DECOUPLING PARAMETERS ---
cfg.target_variants      = [6];     % Specify twin variant indices to evaluate
cfg.shear_scaling_factor = 1.0;     % Twin shear scaling (1.0 = theoretical)
cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
cfg.process_increments   = [4];     % Indices of DIC increments to process

% --- SSLIP CONTROLS & EXPORT OPTIMIZATION ---
cfg.stress_vector        = yvector; % Uniaxial stress direction for Schmid factor calculation
cfg.sslip_method         = 2;       %
cfg.posConstr            = 0;       % Enforce positive shear constraint
cfg.threshResidual       = 0.02;    % Residual threshold
cfg.threshResidualFraction = 0.50; % Threshold residual as fraction of effective shear
cfg.minEeff              = 0.02;    % Minimum effective shear strain required for evaluation

cfg.clim_percentile      = 99;      % Robust upper limit percentile for plotting
cfg.jpg_target_px        = 1000;    % Reduced pixel resolution for lightweight JPG exports
cfg.jpg_dpi              = 150;     % Reduced DPI for faster file saving and compact disk footprint

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
CS = newEBSD.CS;

%% 3. Instantiate ALL <a> Type Slip Systems
sS_basal = slipSystem.basal(CS);
sS_prism = slipSystem.prismaticA(CS);
sS_pyrA  = slipSystem.pyramidalA(CS);
sS_a_all = symmetrise([sS_basal; sS_prism; sS_pyrA], 'antipodal');
fprintf('Instantiated %d total <a> type slip systems (3 Basal, 3 Prism, 6 Pyr <a>).\n', length(sS_a_all));

% Define tensile twin variants
tS_1012 = twinSystem.hexagonal_1012(CS);
tS_sym  = tS_1012.symmetrise;
H_twin_global = originalParentOri * tS_sym.displacementGradient;

%% 4. Discover DIC Increments
all_fields = fields(newEBSD.prop);
displacement_fields = all_fields(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true));
prefix_list = unique(cellfun(@(x) x(1:find(x=='_', 1)-1), displacement_fields, 'UniformOutput', false));

if isempty(cfg.process_increments)
    inc_indices = 1:length(prefix_list);
else
    inc_indices = intersect(cfg.process_increments, 1:length(prefix_list));
end

%% 5. EXECUTION LOOP: COLINEAR-B CROSS-SLIP PAIR EVALUATION
for inc_idx = inc_indices
    prefix = prefix_list{inc_idx};
    fprintf('\n========================================================================\n');
    fprintf('EVALUATING DIC INCREMENT [%s] FOR CROSS-SLIP DIAGNOSTICS\n', prefix);
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
        parentOri_twinned = originalParentOri * twinMis;
        
        for factor_val = cfg.shear_scaling_factor
            factor_tag = sprintf('Factor_%04d', round(factor_val * 100));
            
            % Multiplicative decoupling
            H_tv = H_twin_global(v).M(1:2, 1:2);
            F_twin_2d = eye(2) + factor_val * H_tv;
            F_twin_inv = inv(F_twin_2d);
            
            newEBSD.prop.Hxx_sub = (F_tot_11 * F_twin_inv(1,1) + F_tot_12 * F_twin_inv(2,1)) - 1;
            newEBSD.prop.Hxy_sub = (F_tot_11 * F_twin_inv(1,2) + F_tot_12 * F_twin_inv(2,2));
            newEBSD.prop.Hyx_sub = (F_tot_21 * F_twin_inv(1,1) + F_tot_22 * F_twin_inv(2,1));
            newEBSD.prop.Hyy_sub = (F_tot_21 * F_twin_inv(1,2) + F_tot_22 * F_twin_inv(2,2)) - 1;
            
            for g_i = 1:length(valid_grains)
                grainId = valid_grains(g_i).id;
                fprintf('  -> [Var %d | Grain %d | %s] Clustering <a> systems by colinear Burgers vector...\n', v, grainId, factor_tag);
                
                sSLocal = parentOri_twinned * sS_a_all;
                stressState = stressTensor.uniaxial(cfg.stress_vector);
                
                % Flip normal vectors to guarantee positive Schmid Factor
                sf_signs = sign(sSLocal.SchmidFactor(stressState));
                sf_signs(sf_signs == 0) = 1;
                sSLocal.n = sf_signs .* sSLocal.n;
                
                % --- GROUP SYSTEMS BY UNIQUE BURGERS VECTOR DIRECTION ---
                b_groups = {};
                for s_idx = 1:length(sSLocal)
                    matched = false;
                    for grp_i = 1:length(b_groups)
                        ref_b = sSLocal(b_groups{grp_i}(1)).b;
                        % Check if parallel or antipodal within tolerance
                        if angle(sSLocal(s_idx).b, ref_b) < 1e-3 || angle(sSLocal(s_idx).b, -ref_b) < 1e-3
                            b_groups{grp_i} = [b_groups{grp_i}, s_idx];
                            matched = true;
                            break;
                        end
                    end
                    if ~matched
                        b_groups{end+1} = [s_idx];
                    end
                end
                
                fprintf('     Identified %d unique <a> Burgers vector line directions.\n', length(b_groups));
                
                % Slice EBSD to grain
                ebsd_grain = newEBSD(grains(grainId)).gridify;
                ebsd_grain.orientations = repmat(parentOri_twinned, size(ebsd_grain.orientations));
                U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos;
                V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
                
                % Setup baseline IDoptions (Plotting OFF for automated pair screening)
                clear IDoptions;
                IDoptions.IDMethod       = 1; % Mode 1: Combined multi-slip
                IDoptions.filterSize     = 4;
                IDoptions.coarsegrain    = 4;
                IDoptions.threshResidual = cfg.threshResidual;
                IDoptions.minEeff        = cfg.minEeff;
                IDoptions.posConstr      = cfg.posConstr;
                IDoptions.plotSSLIP      = 0; % Disable interactive figure popup during automated loop
                IDoptions.plotTraces     = 1;
                IDoptions.stress         = stressState;
                IDoptions.Hxx            = ebsd_grain.prop.Hxx_sub;
                IDoptions.Hxy            = ebsd_grain.prop.Hxy_sub;
                IDoptions.Hyx            = ebsd_grain.prop.Hyx_sub;
                IDoptions.Hyy            = ebsd_grain.prop.Hyy_sub;
                IDoptions.clim_percentile = cfg.clim_percentile;
                try, IDoptions.cmap = viridis(512); catch, IDoptions.cmap = parula(256); end
                
                % --- EXHAUSTIVE PAIRWISE EVALUATION ---
                pair_results = struct([]);
                p_cnt = 0;
                
                for grp_i = 1:length(b_groups)
                    sys_indices = b_groups{grp_i};
                    pairs = nchoosek(sys_indices, 2);
                    
                    for p_i = 1:size(pairs, 1)
                        p_cnt = p_cnt + 1;
                        idx1 = pairs(p_i, 1);
                        idx2 = pairs(p_i, 2);
                        sS_pair = sSLocal([idx1, idx2]);
                        
                        % Label systems by family type
                        fam_names = cell(1, 2);
                        for k = 1:2
                            idx_orig = pairs(p_i, k);
                            if idx_orig <= 3, fam_names{k} = 'Basal <a>';
                            elseif idx_orig <= 6, fam_names{k} = 'Prism <a>';
                            else, fam_names{k} = 'Pyr <a>'; end
                        end
                        
                        pair_str = sprintf('Group %d | %s (Sys %d) + %s (Sys %d)', grp_i, fam_names{1}, idx1, fam_names{2}, idx2);
                        
                        IDoptions.NoSs     = [1, 2];
                        IDoptions.casename = sprintf('Pair_%02d_Grp%d', p_cnt, grp_i);
                        IDoptions.plotname = fullfile(outDir, IDoptions.casename);
                        
                        % Run SSLIP solver
                        try
                            [ebsdID_out, optOut] = SSLIP(ebsd_grain, U, V, sS_pair, IDoptions);
                            res_map = ebsdID_out.prop.residualEeff(:);
                            
                            % Calculate robust quality metrics (Combined Error Index & Correlation Rate)
                            Eeff_eval = ebsdID_out.prop.Eeff(:);
                            eval_mask = (Eeff_eval >= cfg.minEeff) & ~isnan(Eeff_eval);
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
                            mean_res = r_corr / max(1e-3, corr_frac);
                            
                            act_map  = ebsdID_out.prop.slipIDcor;
                            mean_act1 = mean(act_map(1, ~isnan(act_map(1,:))), 'all', 'omitnan');
                            mean_act2 = mean(act_map(2, ~isnan(act_map(2,:))), 'all', 'omitnan');
                        catch ME
                            mean_res = NaN; corr_rate = 0; mean_act1 = NaN; mean_act2 = NaN; ebsdID_out = []; optOut = [];
                        end
                        
                        pair_results(p_cnt).pair_id    = p_cnt;
                        pair_results(p_cnt).group_id   = grp_i;
                        pair_results(p_cnt).indices    = [idx1, idx2];
                        pair_results(p_cnt).families   = fam_names;
                        pair_results(p_cnt).pair_name  = pair_str;
                        pair_results(p_cnt).mean_res   = mean_res;
                        pair_results(p_cnt).corr_rate  = corr_rate;
                        pair_results(p_cnt).sf         = abs(sS_pair.SchmidFactor(stressState))';
                        pair_results(p_cnt).mean_act   = [mean_act1, mean_act2];
                        pair_results(p_cnt).ebsdID     = ebsdID_out;
                        pair_results(p_cnt).sS_pair    = sS_pair;
                        pair_results(p_cnt).optOut     = optOut;
                    end
                end
                
                % --- RANK PAIRS BY COMBINED ERROR INDEX ---
                all_res = [pair_results.mean_res];
                [~, sort_idx] = sort(all_res, 'ascend');
                ranked_results = pair_results(sort_idx);
                
                % Print clean comparison leaderboard to terminal
                fprintf('\n     ===========================================================================================================\n');
                fprintf('     COLINEAR-B CROSS-SLIP PAIR LEADERBOARD (Ranked by Combined Error Index)\n');
                fprintf('     ===========================================================================================================\n');
                fprintf('     Rank | Pair ID | Cross-Slip Combination                     | Combined Error Index | Correlation %% | SF 1 / SF 2 \n');
                fprintf('     -----------------------------------------------------------------------------------------------------------\n');
                for r_i = 1:length(ranked_results)
                    pr = ranked_results(r_i);
                    fprintf('     #%02d  | Pair %02d | %-42s |       %.6f       |    %5.1f%%     | %.2f / %.2f \n', ...
                            r_i, pr.pair_id, pr.pair_name, pr.mean_res, pr.corr_rate, pr.sf(1), pr.sf(2));
                end
                fprintf('     ===========================================================================================================\n\n');
                
                % Save full diagnostic table to CSV
                csv_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_CrossSlip_Leaderboard.csv', prefix, grainId, v));
                sf_mat  = reshape([ranked_results.sf], 2, [])';
                act_mat = reshape([ranked_results.mean_act], 2, [])';
                T_csv = table((1:length(ranked_results))', [ranked_results.pair_id]', {ranked_results.pair_name}', ...
                              [ranked_results.mean_res]', [ranked_results.corr_rate]', sf_mat(:,1), sf_mat(:,2), act_mat(:,1), act_mat(:,2), ...
                              'VariableNames', {'Rank', 'Pair_ID', 'Combination_Name', 'Combined_Error_Index', 'Correlation_Rate_Pct', 'Schmid_Factor_1', 'Schmid_Factor_2', 'Mean_Shear_Activity_1', 'Mean_Shear_Activity_2'});
                writetable(T_csv, csv_path);
                
                % Save MAT file for fast re-plotting
                mat_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_CrossSlip_Data.mat', prefix, grainId, v));
                save(mat_path, 'ranked_results', 'cfg', 'grainId', 'v');
                
                % --- AUTOMATED PLOTTING & LIGHTWEIGHT JPG EXPORT FOR ALL COMBINATIONS ---
                fprintf('     -> Generating lightweight JPG figures for ALL %d cross-slip combinations (%dx%d px / %d DPI)...\n', ...
                        length(ranked_results), cfg.jpg_target_px, cfg.jpg_target_px, cfg.jpg_dpi);
                
                for r_i = 1:length(ranked_results)
                    win_pr = ranked_results(r_i);
                    if isempty(win_pr.ebsdID), continue; end
                    
                    win_opt = win_pr.optOut;
                    win_opt.saveFig = 0; % Disable slow internal saveFigure (avoids heavy PNGs)
                    win_opt.casename = sprintf('Rank%02d_Pair%02d_Grp%d_Var%d', r_i, win_pr.pair_id, win_pr.group_id, v);
                    win_opt.plotname = fullfile(outDir, win_opt.casename);
                    win_opt.comment  = sprintf('Rank #%02d Cross-Slip: %s (SF: %.2f / %.2f)', r_i, win_pr.pair_name, win_pr.sf(1), win_pr.sf(2));
                    
                    close all;
                    % Re-trigger plotSSLIP directly on precomputed fields
                    plotSSLIP(win_pr.ebsdID.prop.slipIDcor, win_pr.ebsdID.prop.residualEeff, win_pr.ebsdID, win_pr.sS_pair, win_opt);
                    
                    % Export all open figures instantly as lightweight JPGs via exportScaledFigure
                    figs = findall(0, 'Type', 'figure');
                    [~, s_idx_f] = sort([figs.Number]);
                    figs = figs(s_idx_f);
                    
                    if ~isempty(figs)
                        exportScaledFigure(figs(1), [win_opt.plotname, '_SlipActivities.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    end
                    if length(figs) >= 2
                        exportScaledFigure(figs(2), [win_opt.plotname, '_Residual.jpg'], 'target_px', cfg.jpg_target_px, 'dpi', cfg.jpg_dpi, 'nodatestring');
                    end
                    close all;
                end
                
                % --- GENERATE COMPREHENSIVE MARKDOWN REPORT ---
                md_path = fullfile(outDir, sprintf('%s_Grain%d_Var%d_CrossSlip_Report.md', prefix, grainId, v));
                fprintf('     -> Writing interactive Markdown Report with embedded JPG figures to:\n        %s\n', md_path);
                
                fid = fopen(md_path, 'w');
                fprintf(fid, '# <a>-Type Colinear Burgers Vector Cross-Slip Diagnostics Report\n\n');
                fprintf(fid, '- **DIC Increment:** `%s`\n', prefix);
                fprintf(fid, '- **Grain ID:** `%d`\n', grainId);
                fprintf(fid, '- **Twin Variant:** `%d`\n', v);
                fprintf(fid, '- **Shear Scaling Factor:** `%.2f`\n\n', factor_val);
                
                fprintf(fid, '## Physical Hypothesis & Method Summary\n');
                fprintf(fid, 'In HCP Titanium, screw dislocations with identical $\\langle a \\rangle$ Burgers vectors ($\\mathbf{b} = \\frac{1}{3}\\langle 11\\overline{2}0 \\rangle$) can freely cross-slip between Basal, Prismatic, and Pyramidal $\\langle a \\rangle$ planes.\n');
                fprintf(fid, 'This analysis evaluates SSLIP in Mode 1 (`IDMethod = 1`), restricting candidate systems exclusively to pairwise colinear-$\\mathbf{b}$ combinations to identify which primary + cross-slip pair best explains the observed experimental surface deformation.\n\n');
                
                fprintf(fid, '## Colinear-b Cross-Slip Leaderboard\n');
                fprintf(fid, 'Ranked by lowest **Combined Error Index** ($r_{\\text{comb}} = \\frac{\\text{mean correlated residual}}{\\text{correlation fraction}}$) to ensure sharp quality differentiation without constant threshold flattening.\n\n');
                fprintf(fid, '| Rank | Pair ID | Cross-Slip Combination | Combined Error Index | Correlation %% | SF 1 | SF 2 | Mean Shear $\\gamma_1$ | Mean Shear $\\gamma_2$ |\n');
                fprintf(fid, '| :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n');
                for r_i = 1:length(ranked_results)
                    pr = ranked_results(r_i);
                    fprintf(fid, '| **#%02d** | Pair %02d | %s | `%.6f` | `%5.1f%%` | `%.2f` | `%.2f` | `%.4f` | `%.4f` |\n', ...
                            r_i, pr.pair_id, pr.pair_name, pr.mean_res, pr.corr_rate, pr.sf(1), pr.sf(2), pr.mean_act(1), pr.mean_act(2));
                end
                fprintf(fid, '\n---\n\n');
                
                fprintf(fid, '## Detailed Pairwise Maps & Residual Distributions\n\n');
                for r_i = 1:length(ranked_results)
                    pr = ranked_results(r_i);
                    case_str = sprintf('Rank%02d_Pair%02d_Grp%d_Var%d', r_i, pr.pair_id, pr.group_id, v);
                    
                    fprintf(fid, '### Rank #%02d: %s (Pair %02d)\n', r_i, pr.pair_name, pr.pair_id);
                    fprintf(fid, '- **Combined Error Index:** `%.6f`\n', pr.mean_res);
                    fprintf(fid, '- **Correlation Rate:** `%.1f%%`\n', pr.corr_rate);
                    fprintf(fid, '- **Schmid Factors:** System 1 = `%.2f`, System 2 = `%.2f`\n', pr.sf(1), pr.sf(2));
                    fprintf(fid, '- **Mean Shear Activity:** $\\gamma_1$ = `%.4f`, $\\gamma_2$ = `%.4f`\n\n', pr.mean_act(1), pr.mean_act(2));
                    
                    fprintf(fid, '#### Slip Activities (System 1 & System 2)\n');
                    fprintf(fid, '![Slip Activities %s](%s_SlipActivities.jpg)\n\n', case_str, case_str);
                    
                    fprintf(fid, '#### Residual Error Map\n');
                    fprintf(fid, '![Residual Map %s](%s_Residual.jpg)\n\n', case_str, case_str);
                    fprintf(fid, '---\n\n');
                end
                fclose(fid);
                
                fprintf('     -> Complete! Saved all %d combination maps, CSV leaderboard, and Markdown Report to:\n        %s\n', length(ranked_results), outDir);
            end
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL COLINEAR-B CROSS-SLIP DIAGNOSTIC RUNS COMPLETE!\n');
fprintf('========================================================================\n');
