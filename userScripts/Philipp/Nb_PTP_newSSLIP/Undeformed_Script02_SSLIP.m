%% Undeformed_Script02_SSLIP.m
% Runs the modern modular SSLIP solver over the full EBSD map for all DIC increments.
% Solves the map grain-by-grain and merges the results back into a global ebsdID.
% Saves SSLIP_Data_BXXXX.mat files to the output directory.

clear; close all; clc;
plotx2east; plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 1. Load Registry and Aligned Data
[registry, defaults] = a00_Dataset_Registry();
outPath = fullfile(defaults.dir, defaults.outDir_undeformed);

% Search for the newest alignment file
search_pattern = fullfile(outPath, '*_undeformedEBSD_withAllDIC.mat');
files = dir(search_pattern);
if isempty(files)
    error('No alignment file found matching %s. Run Undeformed_Script01_Alignment.m first.', search_pattern);
end
[~, idx] = sort([files.datenum], 'descend');
master_alignment_file = fullfile(outPath, files(idx(1)).name);

if length(files) > 1
    warning('Multiple alignment files found. Loading the most recent one:\n%s', master_alignment_file);
end

fprintf('Loading master aligned EBSD data...\n');
load(master_alignment_file, 'newEBSD', 'grains');

%% 2. Setup Map and References
% Get the actual EBSD grid we will solve on
ebsdMap = newEBSD.gridify;

% Setup Global Reference Slip Systems
CS = ebsdMap.CS;
sSLocalRef = slipSystem.bcc(CS);
sSLocalRef = symmetrise(sSLocalRef(:), 'antipodal');
stress = stressTensor.uniaxial(yvector);

cfg = struct();
cfg.preprocess = struct('filterSize', 4, 'coarsegrain', 4);
cfg.solver = struct('IDMethod', 1, 'threshResidual', 0.02, 'minEeff', 0.01, 'posConstr', 1, 'enableRotation', 1, 'stress', stress);
cfg.solver.NoSs = 1:length(sSLocalRef);

minSize = 1000;
valid_grains = grains(grains.numPixel > minSize);
if isempty(valid_grains)
    error('No grains found larger than %d pixels.', minSize);
end

%% 3. Loop over all increments in the registry
for i = [3,5,7] %1:length(registry)
    step = registry(i);
    
    if (isfield(step, 'exclude') && step.exclude) || isempty(step.dic_fname)
        continue;
    end
    
    prefix = strrep(step.dic_fname, '.hdf5', '');
    
    % Generate valid field names exactly as Script 01 does
    field_Hxx = matlab.lang.makeValidName([prefix, '_Gradients_du÷dx - RC']);
    field_Hxy = matlab.lang.makeValidName([prefix, '_Gradients_du÷dy - RC']);
    field_Hyx = matlab.lang.makeValidName([prefix, '_Gradients_dv÷dx - RC']);
    field_Hyy = matlab.lang.makeValidName([prefix, '_Gradients_dv÷dy - RC']);
    
    if ~isfield(ebsdMap.prop, field_Hxx)
        warning('Gradient field %s not found in newEBSD. Skipping %s.', field_Hxx, prefix);
        continue;
    end
    
    fprintf('\n======================================================\n');
    fprintf('Running SSLIP for %s (%s)...\n', step.stepName, prefix);
    
    % Initialize global tracking arrays for this increment
    num_systems = length(sSLocalRef);
    global_slipIDcor = zeros(numel(ebsdMap), num_systems);
    global_residualEeff = zeros(numel(ebsdMap), 1);
    
    % Prepare global DefData (though we will subset it per grain)
    global_DefData = struct();
    global_DefData.Hxx = ebsdMap.prop.(field_Hxx);
    global_DefData.Hxy = ebsdMap.prop.(field_Hxy);
    global_DefData.Hyx = ebsdMap.prop.(field_Hyx);
    global_DefData.Hyy = ebsdMap.prop.(field_Hyy);
    
    % --- 4. Loop over grains ---
    for g_i = 1:length(valid_grains)
        grainId = valid_grains(g_i).id;
        fprintf('  -> Processing Grain %d (%d px)...\n', grainId, valid_grains(g_i).numPixel);
        
        % 4.1 Project global slip systems into local grain orientation
        sSLocal = valid_grains(g_i).meanOrientation * sSLocalRef;
        
        % 4.2 Extract sub-EBSD and gridify
        % CRITICAL: MTEX .gridify destroys the global .id values and assigns local bounding-box IDs!
        % We MUST stash the global IDs in a custom property before gridifying so they survive!
        ebsd_subset = ebsdMap(grains(grainId));
        ebsd_subset.prop.global_id = ebsd_subset.id;
        ebsd_grain = ebsd_subset.gridify;
        
        % 4.3 Subset DefData for this grain
        % We extract the properties directly from the gridified ebsd_grain
        DefData_grain = struct();
        DefData_grain.Hxx = ebsd_grain.prop.(field_Hxx);
        DefData_grain.Hxy = ebsd_grain.prop.(field_Hxy);
        DefData_grain.Hyx = ebsd_grain.prop.(field_Hyx);
        DefData_grain.Hyy = ebsd_grain.prop.(field_Hyy);
        
        % 4.4 Run SSLIP
        [grain_ebsdID, prepData, optOut, ~] = SSLIP(ebsd_grain, DefData_grain, sSLocal, cfg);
        
        % 4.4.1 Plot Raw Activities
        opt_plot = struct();
        plotDir = fullfile(outPath, 'Plots', 'ActivityFields');
        if ~exist(plotDir, 'dir'), mkdir(plotDir); end
        opt_plot.plotname = fullfile(plotDir, sprintf('%s_Grain%d_RawActivities', prefix, grainId));
        opt_plot.saveFig = 1;
        opt_plot.bundleSimilar = 0; % Force raw plotting
        opt_plot.plotGrainContext = 0;
        opt_plot.ori = valid_grains(g_i).meanOrientation;

        opt_plot.plotSlipActivities = 1;
        opt_plot.plotResidual = 1;
        opt_plot.plotRotation = 1;
        opt_plot.plotDefGrad = 1;

        figs_before = findobj('type','figure');
        plotSSLIP(grain_ebsdID, sSLocal, opt_plot);
        figs_after = findobj('type','figure');
        close(setdiff(figs_after, figs_before));
        
        % 4.5 Map the solved subset back to the global arrays
        % Mathematically upscale the coarse-grained solver results (if applicable) 
        % back onto the original high-res grain grid and flatten for native MTEX .id assignment.
        [slip_hr, res_hr, valid_idx] = upscaleSSLIP(grain_ebsdID, ebsd_grain);
        
        % Inject directly into the global array using the native linear memory indices
        global_slipIDcor(valid_idx, :) = slip_hr;
        global_residualEeff(valid_idx) = res_hr;
    end
    
    % --- 5. Package and Save ---
    fprintf('Packaging global SSLIP results...\n');
    ebsdID = ebsdMap;
    ebsdID.prop.slipIDcor = global_slipIDcor;
    ebsdID.prop.residualEeff = global_residualEeff;
    
    % The output must contain ebsdID and sSLocalRef for downstream scripts
    out_mat = fullfile(outPath, sprintf('SSLIP_Data_%s.mat', prefix));
    save(out_mat, 'ebsdID', 'sSLocalRef', 'prepData', 'optOut');
    fprintf('Saved %s\n', out_mat);
    close all
end

fprintf('\nSSLIP Processing complete for all increments!\n');
