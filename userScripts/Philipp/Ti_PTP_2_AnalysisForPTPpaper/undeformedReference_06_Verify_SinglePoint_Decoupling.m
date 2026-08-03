%% undeformedReference_06_Verify_SinglePoint_Decoupling.m
% Standalone Diagnostic Verification of Multiplicative Twin Shear Decoupling
% -------------------------------------------------------------------------
% Evaluates the exact mathematical transformations at a specified target 
% coordinate (default: [x=1280/2, y=900]) and logs step-by-step matrices to
% the Command Window for interactive validation and review.
% -------------------------------------------------------------------------
clear; clc; close all;

addpath('c:\Users\phkr\MTEX_workEnv\userScripts\HelperFunctions');
addpath('c:\Users\phkr\MTEX_workEnv\geometry');

%% 0. Configuration
target_xy       = [1280 / 2, 900]; % Target EBSD datapoint [X, Y] coordinate
target_variant  = 6;               % Twin variant index to analyze in detail (1-6)
stress_vec      = yvector;         % Tensile loading direction for Schmid Factor comparison
inc             = 4;               % DIC increment

%% 1. Load Reference Dataset & Identify Primary Parent
[registry, defaults] = a00_Dataset_Registry();
inDir = fullfile(defaults.dir, defaults.outDir_undeformed);
matFiles = dir(fullfile(inDir, '*_undeformedEBSD_withAllDIC.mat'));
if isempty(matFiles)
    error('No reference .mat file found in %s. Run undeformedReference_03_IMPORTER first.', inDir);
end
[~, idx] = max([matFiles.datenum]);
matname = matFiles(idx).name;

fprintf('========================================================================\n');
fprintf('  SINGLE-POINT TWIN SHEAR DECOUPLING VERIFICATION\n');
fprintf('========================================================================\n');
fprintf('Dataset: %s\n\n', matname);
dataLoaded = load(fullfile(inDir, matname));
newEBSD    = dataLoaded.newEBSD;

% Identify grains and largest parent matrix
if isfield(dataLoaded, 'grains')
    grains = dataLoaded.grains;
else
    [grains, newEBSD.grainId, ~] = calcGrains(newEBSD, 'minPixel', 300, 'angle', 5*degree);
end
valid_grains = grains(grains.numPixel > 300);
[~, maxIdx] = max(valid_grains.numPixel);
parentGrain = valid_grains(maxIdx);

%% 2. Locate Target Datapoint in EBSD Field
dist = sqrt((newEBSD.x - target_xy(1)).^2 + (newEBSD.y - target_xy(2)).^2);
[min_dist, pt_idx] = min(dist);
ebsd_pt = newEBSD(pt_idx);

fprintf('--> Target Search Coordinates : [X = %.2f, Y = %.2f]\n', target_xy(1), target_xy(2));
fprintf('--> Matched EBSD Point Index  : #%d (distance: %.3f units)\n', pt_idx, min_dist);
fprintf('    Matched Exact Coordinates : [X = %.2f, Y = %.2f]\n', ebsd_pt.x, ebsd_pt.y);
fprintf('    Phase: %s | Confidence Index (CI): %.3f\n\n', ebsd_pt.mineral, ebsd_pt.ci);

%% 3. Theoretical Twin Displacement Gradients
CS = newEBSD.CS;
tS_1012 = twinSystem.hexagonal_1012(CS);
tS_sym  = tS_1012.symmetrise; 

% Global transformation using original parent matrix orientation
H_twin_global = parentGrain.meanOrientation * tS_sym.displacementGradient;
H_tv = H_twin_global(target_variant).M(1:2, 1:2);
F_tv = eye(2) + H_tv;
F_tv_inv = inv(F_tv);

fprintf('========================================================================\n');
fprintf('  THEORETICAL TWIN TENSORS (VARIANT %d)\n', target_variant);
fprintf('========================================================================\n');
fprintf('1. Theoretical Twin Displacement Gradient H_twin (2x2 In-Plane):\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', H_tv(1,1), H_tv(1,2), H_tv(2,1), H_tv(2,2));
fprintf('2. Theoretical Twin Deformation Gradient F_twin = I + H_twin:\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', F_tv(1,1), F_tv(1,2), F_tv(2,1), F_tv(2,2));
fprintf('3. Inverse Twin Deformation Gradient F_twin^(-1):\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', F_tv_inv(1,1), F_tv_inv(1,2), F_tv_inv(2,1), F_tv_inv(2,2));

%% 4. Evaluate First Available DIC Increment at Point
all_fields = fields(newEBSD.prop);
disp_fields = all_fields(contains(all_fields, '_Displacements_x_displacement', 'IgnoreCase', true));
prefixes = unique(cellfun(@(x) x(1:find(x=='_', 1)-1), disp_fields, 'UniformOutput', false));

if isempty(prefixes)
    error('No DIC increments found in dataset.');
end

% Pick the latest/last increment to observe significant deformation, or first available
eval_inc = prefixes{inc}; 
fprintf('========================================================================\n');
fprintf('  EXPERIMENTAL DIC DATA DECOUPLING AT POINT (INCREMENT: %s)\n', eval_inc);
fprintf('========================================================================\n');

% Retrieve raw gradients at point #pt_idx
hxx_f = [eval_inc, '_Gradients_du_dx_RC']; hxy_f = [eval_inc, '_Gradients_du_dy_RC'];
hyx_f = [eval_inc, '_Gradients_dv_dx_RC']; hyy_f = [eval_inc, '_Gradients_dv_dy_RC'];

if isfield(newEBSD.prop, hxx_f)
    H_raw = [newEBSD.prop.(hxx_f)(pt_idx), newEBSD.prop.(hxy_f)(pt_idx); ...
             newEBSD.prop.(hyx_f)(pt_idx), newEBSD.prop.(hyy_f)(pt_idx)];
else
    error('Pre-computed gradient field not found on dataset for %s.', eval_inc);
end

F_total = eye(2) + H_raw;

fprintf('1. Raw Experimental Displacement Gradient (H_raw) at Point #%d:\n', pt_idx);
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', H_raw(1,1), H_raw(1,2), H_raw(2,1), H_raw(2,2));
fprintf('2. Total Experimental Deformation Gradient F_total = I + H_raw:\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', F_total(1,1), F_total(1,2), F_total(2,1), F_total(2,2));

% --- Rigorous Multiplicative Decoupling ---
F_slip_mult = F_total * F_tv_inv;
H_slip_mult = F_slip_mult - eye(2);

% --- Old Simple Additive Subtraction ---
H_slip_add = H_raw - H_tv;

% --- Difference (Nonlinear Cross-Term Error) ---
diff_matrix = H_slip_add - H_slip_mult;

fprintf('------------------------------------------------------------------------\n');
fprintf('  DECOUPLED SECONDARY SLIP RESIDUALS\n');
fprintf('------------------------------------------------------------------------\n');
fprintf('[RIGOROUS] Multiplicative Decoupled Slip Gradient H_slip = (F_total * F_twin^-1) - I:\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', H_slip_mult(1,1), H_slip_mult(1,2), H_slip_mult(2,1), H_slip_mult(2,2));

fprintf('[APPROXIMATE] Simple Additive Subtracted Slip Gradient H_add = H_raw - H_twin:\n');
fprintf('       [ %9.6f   %9.6f ]\n       [ %9.6f   %9.6f ]\n\n', H_slip_add(1,1), H_slip_add(1,2), H_slip_add(2,1), H_slip_add(2,2));

fprintf('[DIFFERENCE] Nonlinear Cross-Term Distortion (H_add - H_slip):\n');
fprintf('       [ %+9.6f   %+9.6f ]\n       [ %+9.6f   %+9.6f ]\n', diff_matrix(1,1), diff_matrix(1,2), diff_matrix(2,1), diff_matrix(2,2));
fprintf('    Max absolute deviation at this point: %.6f\n\n', max(abs(diff_matrix(:))));

%% 5. Crystallographic Lattice Reorientation & Schmid Factor Comparison
twinMis = tS_sym(target_variant).parentTwinMisorientation;
parentOri_twinned = parentGrain.meanOrientation * twinMis;

fprintf('========================================================================\n');
fprintf('  CRYSTALLOGRAPHIC REORIENTATION & SCHMID FACTORS (TENSILE Y-VECTOR)\n');
fprintf('========================================================================\n');
[aE_old,bE_old,cE_old] = parentGrain.meanOrientation.Euler;
eu_old = [aE_old,bE_old,cE_old] ./ degree;
[aE_new,bE_new,cE_new] = parentOri_twinned.Euler;
eu_new = [aE_new,bE_new,cE_new] ./ degree;
fprintf('Original Parent Euler Angles (phi1, Phi, phi2) : (%.1f°, %.1f°, %.1f°)\n', eu_old(1), eu_old(2), eu_old(3));
fprintf('Twinned Lattice Euler Angles (phi1, Phi, phi2) : (%.1f°, %.1f°, %.1f°)\n\n', eu_new(1), eu_new(2), eu_new(3));

% Define slip systems
sS = symmetrise([slipSystem.basal(CS); slipSystem.prismaticA(CS)], 'antipodal');
stressState = stressTensor.uniaxial(stress_vec);

sS_orig_proj = parentGrain.meanOrientation * sS;
SF_orig      = abs(SchmidFactor(sS_orig_proj, stressState));

sS_twin_proj = parentOri_twinned * sS;
SF_twin      = abs(SchmidFactor(sS_twin_proj, stressState));

[~, idx_sort] = sort(SF_twin, 'descend');
fprintf('Top 5 Secondary Slip Systems inside Twinned Lattice (vs Original Parent):\n');
fprintf('  ID | Family      | Burgers Vector | Plane Normal | SF (Twinned) | SF (Original)\n');
fprintf('  -----------------------------------------------------------------------------\n');
for j = 1:min(5, length(idx_sort))
    s_id = idx_sort(j);
    fam_name = 'Prismatic';
    if contains(lower(char(sS(s_id))), '0001') || contains(lower(char(sS(s_id).n)), '0001'), fam_name = 'Basal    '; end
    fprintf('  #%2d | %-11s | %-14s | %-12s | %12.4f | %13.4f\n', ...
            s_id, fam_name, sS(s_id).b.char, sS(s_id).n.char, SF_twin(s_id), SF_orig(s_id));
end
fprintf('========================================================================\n\n');
