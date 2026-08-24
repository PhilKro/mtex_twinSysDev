%% undeformedReference_06_SSLIP_subtractedTwinShear.m
% Decoupled Twin Shear SSLIP Analysis for Undeformed Reference Dataset
% -------------------------------------------------------------------------
% Implements rigorous multiplicative sequence-of-events decomposition:
%    F_total = F_slip * F_twin   =>   F_slip = F_total * (F_twin)^-1
% where F = I + H.
%
% Additionally transforms the parent crystal reference lattice by the twin 
% variant misorientation, evaluating secondary slip activity across the 
% entire field within the actual active twinned crystal orientation.
% -------------------------------------------------------------------------

clear; close all; clc;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\geometry');
addpath('c:\Users\phkr\MTEX_workEnv\userScripts\SSLIP'); initSSLIP();

plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 0. User Configuration Toggles
cfg = struct();

% --- ANALYSIS NAMING & OUTPUT ROUTING ---
% Assign a distinct name to your analysis run. This name creates a dedicated 
% sub-directory inside 'Processed_Output_Undeformed' and tags all saved figures,
% CSVs, and MAT files to prevent overwriting results from different runs.
% Examples: 'Standard_AllSystems_SF0p10', 'aTypeOnly_NoSF', 'SingleSlipMode_AllSystems'
% cfg.analysis_name        = 'Standard_AllSystems_SF0p10';

% --- TWIN DECOUPLING PARAMETERS ---
% cfg.target_variants      = [6];     % Specify twin variant indices (e.g. [1] or 1:6) to decouple and evaluate
% cfg.shear_scaling_factor = 0:0.5:2; % Array of scaling factors for parameter series analysis (0 = no correction; 1.0 = theoretical; 1.91 = apparent DIC)
% cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
% cfg.process_increments   = [4];     % Indices of DIC increments to process (empty [] processes all found)

% --- SLIP SYSTEM SELECTION & SCHMID FACTOR THRESHOLDING ---
% Choose active slip systems:
%   'all'        : All <a> and <c+a> systems (Basal <a>, Prismatic <a>, Pyramidal <a>, Pyramidal <c+a>)
%   'a_only'     : Only <a> type systems (Basal <a>, Prismatic <a>, Pyramidal <a>)
%   'ca_only'    : Only <c+a> type systems (Pyramidal <c+a>)
%   'basal_prism': Only Basal <a> and Prismatic <a>
% cfg.slip_system_mode     = 'all';   
% cfg.threshold_SF         = 0.10;    % Minimum Schmid Factor (set <= 0 to disregard SF threshold and evaluate all selected systems)
% cfg.stress_vector        = yvector; % Uniaxial stress direction for Schmid factor calculation

% --- SSLIP OPTIMIZATION & THRESHOLD CONTROLS ---
%   IDMethod 1 : Combined & minimized slip ID (Standard multi-slip optimization)
%   IDMethod 2 : Constrained slip ID (Only constraint, no minimization)
%   IDMethod 3 : Single slip ID (Evaluates each slip system individually against experimental gradients)
% cfg.sslip_method         = 1;       

% Adjusting thresholds is critical, particularly when running single slip mode (IDMethod = 3):
%   - If threshResidual is too high -> accepts almost everything (false positive activity)
%   - If threshResidual is too low  -> rejects valid slip activity (identifies nothing)
% cfg.threshResidual       = 0.02;    % Residual threshold || H^exp - H^th || (default: 0.02)
% cfg.threshResidualFraction = 0.2;   % Optional: Threshold residual as fraction of total effective shear (primarily used in IDMethod = 3)
% cfg.minEeff              = 0.01;    % Minimum effective shear strain required for evaluation (default: 0.01)
% cfg.maxE                 = 0.05;    % Optional: Lock maximum colorbar limit for slip activity plots
% cfg.clim_percentile      = 99;      % Optional: Percentile scaling threshold for automatic color limits (default: 99)


% =========================================================================
% PRESET EXAMPLES (Copy & paste over the fields above to quickly switch runs)
% =========================================================================
% --- PRESET A: Standard Multi-Slip with All Systems & SF Threshold ---
cfg.analysis_name    = 'Standard_MultiSlip_AllSystems';

cfg.target_variants      = [6];     % Specify twin variant indices (e.g. [1] or 1:6) to decouple and evaluate
cfg.shear_scaling_factor = 1; % Array of scaling factors for parameter series analysis (0 = no correction; 1.0 = theoretical; 1.91 = apparent DIC)
cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
cfg.process_increments   = [5];     % Indices of DIC increments to process (empty [] processes all found)

cfg.slip_system_mode = 'all';
cfg.threshold_SF     = 0.10;  % Very low SF threshold 0.01
cfg.stress_vector    = yvector; % Uniaxial stress direction for Schmid factor calculation

cfg.sslip_method     = 1;
cfg.posConstr        = 1;
cfg.threshResidual   = 0.01;
cfg.threshResidualFraction = 0.25; % Threshold residual as fraction of effective shear
cfg.minEeff          = 0.01;    % Minimum effective shear strain required for evaluation, just some value

% % --- PRESET B: Only <a> Type Slip Systems, different SF Threshold ---
% cfg.analysis_name    = 'aType_Only_VeryLowSFThreshold_DifferentShearFacs';
% 
% cfg.target_variants      = [6];     % Specify twin variant indices (e.g. [1] or 1:6) to decouple and evaluate
% cfg.shear_scaling_factor = 1:0.2:2; % Array of scaling factors for parameter series analysis (0 = no correction; 1.0 = theoretical; 1.91 = apparent DIC)
% cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
% cfg.process_increments   = [4];     % Indices of DIC increments to process (empty [] processes all found)
% 
% cfg.slip_system_mode = 'a_only';
% cfg.threshold_SF     = 0.01;  % Very low SF threshold 0.01
% cfg.stress_vector    = yvector; % Uniaxial stress direction for Schmid factor calculation
% 
% cfg.sslip_method     = 1;
% cfg.posConstr        = 1;
% cfg.threshResidual   = 0.005;
% cfg.threshResidualFraction = 0.50; % Threshold residual as fraction of effective shear
% cfg.minEeff          = 0.01;    % Minimum effective shear strain required for evaluation, just some value

% % --- PRESET C: Single Slip Mode (IDMethod 3) with All Systems ---
% cfg.analysis_name    = 'SingleSlipMode_AllSystems';
% 
% cfg.target_variants      = [6];     % Specify twin variant indices (e.g. [1] or 1:6) to decouple and evaluate
% cfg.shear_scaling_factor = 1; % Array of scaling factors for parameter series analysis (0 = no correction; 1.0 = theoretical; 1.91 = apparent DIC)
% cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
% cfg.process_increments   = [4];     % Indices of DIC increments to process (empty [] processes all found)
% 
% cfg.slip_system_mode = 'all';
% cfg.threshold_SF     = 0.00;  % Disregards SF threshold entirely
% cfg.stress_vector    = yvector; % Uniaxial stress direction for Schmid factor calculation
% 
% cfg.sslip_method     = 3;     % Single slip identification mode
% cfg.threshResidual   = 0.02;  % Tune this carefully: higher = more permissive, lower = stricter
% cfg.threshResidualFraction = 0.2; % Threshold residual as fraction of effective shear
% cfg.minEeff          = 0.01;    % Minimum effective shear strain required for evaluation

% % --- PRESET D: All Slip Systems, different SF Threshold, different shear scaling ---
% cfg.analysis_name    = 'allSlip_VeryLowSFThreshold_DifferentShearFacs';
% 
% cfg.target_variants      = [6];     % Specify twin variant indices (e.g. [1] or 1:6) to decouple and evaluate
% cfg.shear_scaling_factor = 1.2:0.2:2.2; % Array of scaling factors for parameter series analysis (0 = no correction; 1.0 = theoretical; 1.91 = apparent DIC)
% cfg.minSize              = 300;     % Minimum grain size in pixels to evaluate
% cfg.process_increments   = [4];     % Indices of DIC increments to process (empty [] processes all found)
% 
% cfg.slip_system_mode = 'all';
% cfg.threshold_SF     = 0.01;  % Very low SF threshold 0.01
% cfg.stress_vector    = yvector; % Uniaxial stress direction for Schmid factor calculation
% 
% cfg.sslip_method     = 1;
% cfg.posConstr        = 1;
% cfg.threshResidual   = 0.005;
% cfg.threshResidualFraction = 0.50; % Threshold residual as fraction of effective shear
% cfg.minEeff          = 0.01;    % Minimum effective shear strain required for evaluation, just some value
% 
% % =========================================================================

%% 1. Load Registry & Locate Undeformed Reference Dataset
[registry, defaults] = a00_Dataset_Registry();
inDir = fullfile(defaults.dir, defaults.outDir_undeformed);

% Configure dedicated analysis output directory based on analysis_name
if ~isfield(cfg, 'analysis_name') || isempty(strtrim(cfg.analysis_name))
    outSubdirName = 'SSLIP_DecoupledTwinShear';
    file_tag      = '';
else
    clean_name    = regexprep(cfg.analysis_name, '[\\/:*?"<>| ]', '_');
    outSubdirName = ['SSLIP_DecoupledTwinShear_', clean_name];
    file_tag      = ['_', clean_name];
end
outDir = fullfile(defaults.dir, defaults.outDir_undeformed, outSubdirName);
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('Analysis Run [%s]:\n  -> Output directory created at: %s\n', cfg.analysis_name, outDir);

matFiles = dir(fullfile(inDir, '*_undeformedEBSD_withAllDIC.mat'));
if isempty(matFiles)
    error('No aligned reference .mat file found in %s. Run undeformedReference_03_IMPORTER first.', inDir);
end
[~, idx] = max([matFiles.datenum]);
matname = matFiles(idx).name;
fprintf('Loading undeformed reference dataset:\n  %s\n', fullfile(inDir, matname));
dataLoaded = load(fullfile(inDir, matname));
newEBSD = dataLoaded.newEBSD;

%% 2. Retrieve or Reconstruct Grains & Identify Primary Parent Grain
grainThresh = 5*degree;
alpha = 1;

if isfield(dataLoaded, 'grains')
    grains = dataLoaded.grains;
    fprintf('  -> Successfully loaded precomputed grains directly from reference dataset.\n');
else
    fprintf('  -> Precomputed grains not found. Calculating grains dynamically...\n');
    thresh_ci = defaults.thresh_ci;
    newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
    [grains, newEBSD.grainId, newEBSD.mis2mean] = calcGrains(newEBSD, 'minPixel', cfg.minSize, 'angle', grainThresh, 'alpha', alpha);
end

% Identify primary parent grain matrix (largest grain by pixel area)
valid_grains = grains(grains.numPixel > cfg.minSize);
if isempty(valid_grains)
    error('No valid grains (>%d px) found in dataset.', cfg.minSize);
end
[~, maxIdx] = max(valid_grains.numPixel);
parentGrain = valid_grains(maxIdx);
originalParentOri = parentGrain.meanOrientation;

fprintf('  -> Primary Parent Matrix identified: Grain ID %d (%d pixels)\n', parentGrain.id, parentGrain.numPixel);

%% 3. Define {10-12} Tensile Twinning and HCP Slip Systems
CS = newEBSD.CS;
tS_1012 = twinSystem.hexagonal_1012(CS);
tS_sym  = tS_1012.symmetrise; % All 6 primary {10-12} variants

% Transform theoretical twin displacement gradients into global specimen frame using original parent orientation
H_twin_global = originalParentOri * tS_sym.displacementGradient;

% Define available slip system types for HCP Titanium
sS_basal = slipSystem.basal(CS);
sS_prism = slipSystem.prismaticA(CS);
sS_pyrA  = slipSystem.pyramidalA(CS);
sS_pyrCA = slipSystem.pyramidalCA(CS);

switch lower(cfg.slip_system_mode)
    case 'a_only'
        raw_sS = [sS_basal; sS_prism; sS_pyrA];
        fprintf('  -> Slip System Mode [%s]: Only <a> type systems (Basal <a>, Prismatic <a>, Pyramidal <a>)\n', cfg.slip_system_mode);
    case 'ca_only'
        raw_sS = sS_pyrCA;
        fprintf('  -> Slip System Mode [%s]: Only <c+a> type systems (Pyramidal <c+a>)\n', cfg.slip_system_mode);
    case 'basal_prism'
        raw_sS = [sS_basal; sS_prism];
        fprintf('  -> Slip System Mode [%s]: Basal <a> and Prismatic <a> systems only\n', cfg.slip_system_mode);
    case 'all'
        raw_sS = [sS_basal; sS_prism; sS_pyrA; sS_pyrCA];
        fprintf('  -> Slip System Mode [%s]: ALL systems (Basal <a>, Prismatic <a>, Pyramidal <a>, Pyramidal <c+a>)\n', cfg.slip_system_mode);
    otherwise
        warning('Unknown slip_system_mode "%s". Defaulting to "all" (<a> and <c+a> systems).', cfg.slip_system_mode);
        raw_sS = [sS_basal; sS_prism; sS_pyrA; sS_pyrCA];
end

sS = symmetrise(raw_sS, 'antipodal');
fprintf('  -> Total symmetrized slip systems instantiated: %d\n', length(sS));

%% 4. Discover DIC Increments in Reference Dataset
all_fields = fields(newEBSD.prop);
displacement_fields = all_fields(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true));
prefix_list = unique(cellfun(@(x) x(1:find(x=='_', 1)-1), displacement_fields, 'UniformOutput', false));

if isempty(cfg.process_increments)
    inc_indices = 1:length(prefix_list);
else
    inc_indices = intersect(cfg.process_increments, 1:length(prefix_list));
end
fprintf('Evaluating %d DIC increment(s): %s\n', length(inc_indices), mat2str(inc_indices));

%% 5. Evaluate Multiplicative Twin Decoupled SSLIP across Increments
for idx_inc = inc_indices
    prefix = prefix_list{idx_inc};
    fprintf('\n========================================================================\n');
    fprintf('Processing Increment [%s] (%d/%d)\n', prefix, idx_inc, length(prefix_list));
    fprintf('========================================================================\n');
    
    % --- Retrieve or Calculate Raw Experimental Displacement Gradients ---
    ux_fname = [prefix, '_Displacements_x_displacement'];
    uy_fname = [prefix, '_Displacements_y_displacement'];
    if ~isfield(newEBSD.prop, ux_fname) || ~isfield(newEBSD.prop, uy_fname)
        warning('Missing displacement fields for %s. Skipping.', prefix);
        continue;
    end
    
    hxx_field = [prefix, '_Gradients_du_dx_RC']; hxy_field = [prefix, '_Gradients_du_dy_RC'];
    hyx_field = [prefix, '_Gradients_dv_dx_RC']; hyy_field = [prefix, '_Gradients_dv_dy_RC'];
    
    if isfield(newEBSD.prop, hxx_field)
        Hxx_raw = double(newEBSD.prop.(hxx_field));
        Hxy_raw = double(newEBSD.prop.(hxy_field));
        Hyx_raw = double(newEBSD.prop.(hyx_field));
        Hyy_raw = double(newEBSD.prop.(hyy_field));
    else
        error('No displacement gradients found in the imported ebsd object.')
    end
    
    % Total experimental deformation gradient: F_total = I + H_raw
    F_tot_11 = 1 + Hxx_raw;
    F_tot_12 = Hxy_raw;
    F_tot_21 = Hyx_raw;
    F_tot_22 = 1 + Hyy_raw;
    
    % --- Plot 1: Displacement Gradient Field BEFORE Subtraction ---
    all_vals_raw = [Hxx_raw(:); Hxy_raw(:); Hyx_raw(:); Hyy_raw(:)];
    all_vals_raw = all_vals_raw(~isnan(all_vals_raw));
    cmin_raw = 0.05 * min(all_vals_raw);
    cmax_raw = 0.05 * max(all_vals_raw);
    
    figRaw = figure('Name', 'Displacement Gradient Verification BEFORE Subtraction', 'Visible', 'off', 'Position', [150, 150, 1200, 900]);
    newMtexFigure('layout', [2, 2]);
    raw_vals   = {Hxx_raw, Hxy_raw, Hyx_raw, Hyy_raw};
    raw_titles = {'du/dx (H_{11} Raw)', 'du/dy (H_{12} Raw)', 'dv/dx (H_{21} Raw)', 'dv/dy (H_{22} Raw)'};
    raw_cb_lbl = {'du/dx', 'du/dy', 'dv/dx', 'dv/dy'};
    for k = 1:4
        if k > 1, nextAxis; end
        plot(newEBSD, raw_vals{k}); hold on;
        plot(grains.boundary('indexed', 'notIndexed'), 'lineW', 1, 'color', 'k');
        %quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
        %quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
        title(raw_titles{k});
        colormap(gca, 'parula');
        caxis([cmin_raw cmax_raw]);
        hold off;
    end
    mtexColorbar
    sgtitle(sprintf('[%s] Experimental Displacement Gradients BEFORE Subtraction', prefix), 'Interpreter', 'none');
    exportScaledFigure(figRaw, fullfile(outDir, sprintf('%s_BeforeSubtraction_Gradients.png', prefix)), 'nodatestring');
    close(figRaw);
    
    % --- Loop Over Selected Twin Variants ---
    for v = cfg.target_variants
        if v > numel(tS_sym)
            warning('Variant index %d out of bounds (max %d). Skipping.', v, numel(tS_sym));
            continue;
        end
        
        fprintf('\n   ---------------------------------------------------------------------\n');
        fprintf('   Multiplicative Decoupling & Twinned Lattice SSLIP for Variant %d\n', v);
        fprintf('   ---------------------------------------------------------------------\n');
        
        % Step A: Transform reference orientation into twinned crystal lattice
        twinMis = tS_sym(v).parentTwinMisorientation;

        parentOri_twinned = originalParentOri * twinMis;
        num_factors = length(cfg.shear_scaling_factor);
        
        % --- Loop Over Twinning Shear Scaling Factor Series ---
        for f_idx = 1:num_factors
            factor_val = cfg.shear_scaling_factor(f_idx);
            is_no_correction = abs(factor_val) < 1e-6;
            factor_tag = strrep(sprintf('Factor_%0.2f', factor_val), '.', 'p');
            
            if is_no_correction
                fprintf('\n     =================================================================\n');
                fprintf('     Evaluating Variant %d | Factor 0.00: NO TWIN SHEAR CORRECTION APPLIED\n', v);
                fprintf('     =================================================================\n');
                fig_header_str = sprintf('[%s] Secondary Slip Gradients (Var %d | NO TWIN SHEAR CORRECTION APPLIED, Factor = 0)', prefix, v);
            else
                fprintf('\n     =================================================================\n');
                fprintf('     Evaluating Variant %d | Twinning Shear Scaling Factor: %.2f (F_twin = I + %.2f * H_twin)\n', v, factor_val, factor_val);
                fprintf('     =================================================================\n');
                fig_header_str = sprintf('[%s] Multiplicatively Decoupled Slip Gradients (Variant %d | Shear Scale = %.2f)', prefix, v, factor_val);
            end
            
            % Step B: Multiplicative sequence-of-events decoupling: F_slip = F_tot * F_twin^-1
            H_tv = H_twin_global(v).M(1:2, 1:2);
            H_tv_scaled = factor_val * H_tv;
            
            F_twin_2d  = eye(2) + H_tv_scaled;
            F_twin_inv = inv(F_twin_2d);
            
            a11 = F_twin_inv(1,1); a12 = F_twin_inv(1,2);
            a21 = F_twin_inv(2,1); a22 = F_twin_inv(2,2);
            
            F_slip_11 = F_tot_11 .* a11 + F_tot_12 .* a21;
            F_slip_12 = F_tot_11 .* a12 + F_tot_12 .* a22;
            F_slip_21 = F_tot_21 .* a11 + F_tot_22 .* a21;
            F_slip_22 = F_tot_21 .* a12 + F_tot_22 .* a22;
            
            % Residual secondary slip displacement gradient: H_slip = F_slip - I
            newEBSD.prop.Hxx_sub = F_slip_11 - 1;
            newEBSD.prop.Hxy_sub = F_slip_12;
            newEBSD.prop.Hyx_sub = F_slip_21;
            newEBSD.prop.Hyy_sub = F_slip_22 - 1;
            
            % Step C: Export Visual Check Figure AFTER Subtraction (Decoupled Gradients)
            Hxx_s = newEBSD.prop.Hxx_sub; Hxy_s = newEBSD.prop.Hxy_sub;
            Hyx_s = newEBSD.prop.Hyx_sub; Hyy_s = newEBSD.prop.Hyy_sub;
            
            figSub = figure('Name', 'Decoupled Secondary Slip Gradients AFTER Subtraction', 'Visible', 'off', 'Position', [150, 150, 1200, 900]);
            newMtexFigure('layout', [2, 2]);
            sub_fields = {Hxx_s, Hxy_s, Hyx_s, Hyy_s};
            sub_titles = {'H_{11} Decoupled Slip', 'H_{12} Decoupled Slip', 'H_{21} Decoupled Slip', 'H_{22} Decoupled Slip'};
            for k = 1:4
                if k > 1, nextAxis; end
                plot(newEBSD, sub_fields{k}); hold on;
                plot(grains.boundary('indexed', 'notIndexed'), 'lineW', 1, 'color', 'k');
                quiver(grains, xvector, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'xvector');
                quiver(grains, yvector, 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'yvector');
                title(sprintf('%s (Var %d | %s)', sub_titles{k}, v, factor_tag), 'Interpreter', 'none');
                colormap(gca, 'parula');
                caxis([cmin_raw cmax_raw]);
                hold off;
            end
            mtexColorbar
            sgtitle(fig_header_str, 'Interpreter', 'none');
            exportScaledFigure(figSub, fullfile(outDir, sprintf('%s_DecoupledGradients_Var%d_%s%s.png', prefix, v, factor_tag, file_tag)), 'nodatestring');
            close(figSub);
            
            % Step D: Run SSLIP optimization using Twinned Orientation and Decoupled Shear
            for g_i = 1:length(valid_grains)
                clear IDoptions
                close all
                grainId = valid_grains(g_i).id;
                
                fprintf('     [Var %d | Grain %d (%d px) | %s] Evaluating secondary slip in twinned crystal lattice...\n', v, grainId, valid_grains(g_i).numPixel, factor_tag);
                
                % Project slip systems using the twinned parent orientation
                sSLocal = parentOri_twinned * sS;
                stressState = stressTensor.uniaxial(cfg.stress_vector);
                
                % Reconfigure slip systems to ensure positive Schmid factor under assumed loading
                sf_signs = sign(sSLocal.SchmidFactor(stressState));
                sf_signs(sf_signs == 0) = 1;
                sSLocal.n = sf_signs .* sSLocal.n;
                
                if cfg.threshold_SF <= 0
                    fprintf('        -> SF filter disabled (cfg.threshold_SF <= 0): Retaining all %d selected slip systems without SF thresholding.\n', length(sSLocal));
                    sSLocal_filt = sSLocal;
                    sS_filt      = sS;
                else
                    SF_crit      = abs(sSLocal.SchmidFactor(stressState)) >= cfg.threshold_SF;
                    sSLocal_filt = sSLocal(SF_crit);
                    sS_filt      = sS(SF_crit);
                    fprintf('        -> Applied SF >= %.2f filter: Retaining %d of %d slip systems.\n', cfg.threshold_SF, length(sS_filt), length(sSLocal));
                end
                
                if isempty(sS_filt)
                    warning('     No slip systems satisfying criteria for Grain %d under twinned lattice. Skipping.', grainId);
                    continue;
                end
                
                % Slice EBSD to grain, gridify, and overwrite orientations with twinned lattice orientation
                ebsd_grain = newEBSD(grains(grainId)).gridify;
                ebsd_grain.orientations = repmat(parentOri_twinned, size(ebsd_grain.orientations));
                
                U = ebsd_grain.prop.(ux_fname) * newEBSD.dPos;
                V = ebsd_grain.prop.(uy_fname) * newEBSD.dPos;
                
                plotBaseName = sprintf('%s_SSLIP_Grain%d_Var%d_%s%s', prefix, grainId, v, factor_tag, file_tag);
                
                clear ssl_cfg
                close all
                
                % -------------------------------------------------------------
                % Configure SSLIP run (Modern Architecture)
                % -------------------------------------------------------------
                ssl_cfg = struct();
                ssl_cfg.preprocess.filterSize = 4;
                ssl_cfg.preprocess.coarsegrain = 4;
                
                ssl_cfg.solver.IDMethod = cfg.sslip_method;
                ssl_cfg.solver.threshResidual = cfg.threshResidual;
                if isfield(cfg, 'threshResidualFraction')
                    ssl_cfg.solver.threshResidualFraction = cfg.threshResidualFraction;
                end
                ssl_cfg.solver.minEeff = cfg.minEeff;
                ssl_cfg.solver.posConstr = cfg.posConstr; 
                ssl_cfg.solver.NoSs = 1:length(sS_filt);
                
                % Turn on Vectorized ADMM solver
                ssl_cfg.solver.ADMM = 1;
                ssl_cfg.solver.ADMM_iters = 150;
                
                ssl_cfg.plot.stress = stressState;
                ssl_cfg.plot.casename = plotBaseName;
                ssl_cfg.plot.plotDefGrad = 1;
                ssl_cfg.plot.saveFig = 1;
                ssl_cfg.plot.plotname = fullfile(outDir, plotBaseName);
                ssl_cfg.plot.ori = parentOri_twinned;
                ssl_cfg.plot.outputReport = fullfile(outDir, [plotBaseName, '_TwinnedLattice_SlipSystems.csv']);
                ssl_cfg.plot.sizeAdjust = 0.5;
                ssl_cfg.plot.plotResidual = 1;
                ssl_cfg.plot.plotTraces = 1;
                if isfield(cfg, 'logscale'), ssl_cfg.plot.logscale = cfg.logscale; else, ssl_cfg.plot.logscale = 0; end
                if isfield(cfg, 'maxE'), ssl_cfg.plot.maxE = cfg.maxE; end
                if isfield(cfg, 'clim_percentile'), ssl_cfg.plot.clim_percentile = cfg.clim_percentile; end
                ssl_cfg.plot.logmin = 0.01;
                try, ssl_cfg.plot.cmap = viridis(512); catch, ssl_cfg.plot.cmap = parula(256); end
                
                if is_no_correction
                    ssl_cfg.plot.comment = sprintf('Inc_%s_Grain_%d_TwinnedVar_%d_NoCorrection%s', prefix, grainId, v, file_tag);
                else
                    ssl_cfg.plot.comment = sprintf('Inc_%s_Grain_%d_TwinnedVar_%d_%s%s', prefix, grainId, v, factor_tag, file_tag);
                end
                
                % Feed decoupled secondary slip gradients into optimizer
                DefData = struct();
                DefData.Hxx = ebsd_grain.prop.Hxx_sub;
                DefData.Hxy = ebsd_grain.prop.Hxy_sub;
                DefData.Hyx = ebsd_grain.prop.Hyx_sub;
                DefData.Hyy = ebsd_grain.prop.Hyy_sub;
                
                try
                    fprintf('     Running SSLIP optimization (Vectorized)...\n');
                    [ebsdID, prepData, optOut] = SSLIP(ebsd_grain, DefData, sSLocal_filt, ssl_cfg);
                    
                    % External plotting decoupling
                    plotSSLIP(ebsdID, sSLocal_filt, optOut);
                    
                    resMatPath = fullfile(outDir, [plotBaseName, '_SSLIP_results.mat']);
                    save(resMatPath, 'ebsdID', 'sSLocal_filt', 'optOut', 'grainId', 'parentOri_twinned', 'v', 'factor_val', 'is_no_correction', 'cfg');
                    
                    fprintf('     -> Complete! Saved figures and CSV table to %s\n', outDir);
                catch ME
                    warning('     SSLIP optimization failed for Grain %d (Var %d | %s) in increment %s:\n        %s', grainId, v, factor_tag, prefix, ME.message);
                end
            end
        end
    end
end

fprintf('\n========================================================================\n');
fprintf('ALL DECOUPLED TWINNED LATTICE SSLIP ANALYSES COMPLETE!\n');
fprintf('Results saved in: %s\n', outDir);
fprintf('========================================================================\n');
