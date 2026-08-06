%% Re_PTP3_Launch_Interactive_DIC_Plotter.m
% Orchestration script to ingest Re_PTP_3 DIC deformation series, calculate
% spatial strain tensor fields, and launch the interactive visualization studio
% with scalar property / colormap mapping and live colorbars.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% Add paths for required alignment routines and interactive GUI studio
addpath('userScripts\NanoMech_Aignment_Matlab\src');
addpath('userScripts\Philipp\Ti_PTP_2_AnalysisForPTPpaper');

%% 2. Directories and Tuning Parameters
inputPath  = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
outputPath = fullfile(inputPath, 'PrettyPlots_Interactive_DIC');

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    fprintf('Created output directory: %s\n', outputPath);
end

DIC_file = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20241218_Re_7dc_tkdPTP_3\DaVis\10kV_Images\export\Processing_subset=19_step=6_relativeToFirst_median3\1st_OrderShapeFun.hdf5';
sel_points_fname = char(fullfile(inputPath, 'inc_0_Selection_Points_Backup_EBSD_05-08-2025_10-23.mat'));

% EBSD and SE inputs for Undeformed Specimen (inc_0)
ang_file = fullfile(inputPath, 'EBSD_angs', 'Re_7deg_insitu Specimen 1 0 Map Data 1_BW158.ang');
se_file  = fullfile(inputPath, 'in_Lens_SE', '00_DIC_0.tif');

CS = {'notIndexed', crystalSymmetry('6/mmm', [2.761 2.761 4.458], 'x||a', 'mineral', 'Re')};
new_pixel_size = 5.4687 * 10^-3;
filterSize     = 3; % Filter window size for calcStrainsFromRawDisp

%% 3. Load and Align Undeformed EBSD and SE Data
fprintf('Loading and aligning undeformed domain dataset...\n');

% Setup EBSD entry
DataName = 'EBSD0';
rawData.(DataName).fname = ang_file;
rawData.(DataName).dtype = 'EBSD';
rawData.(DataName).CS    = CS;
rawData.(DataName).data  = [];

% Setup SE entry
DataName = 'SE_raw';
rawData.(DataName).fname = se_file;
rawData.(DataName).psize = new_pixel_size;
rawData.(DataName).dtype = 'SEM';
rawData.(DataName).imregion = [1281 1 1279 1279]; % startx starty numx-1 numy-1

rawData.dataNames = fieldnames(rawData);
rawData = loadRawData(rawData);

% Alignment options
options.Ndof = 6;
options.Nsel = 4;
options.case = 'EBSD';
options.selectionIPF = 1;
options.EBSD_selection = 'ci';
options.colormap = 'gray';
options.sel_points_fname = sel_points_fname;

alignment_data_EBSD = pointSelectionForAlignment(rawData.SE_raw, rawData.EBSD0, options);

% Warp to consistent grid domain
[Xnew, Ynew] = meshgrid(0:new_pixel_size:max(rawData.SE_raw.X(:)), 0:new_pixel_size:max(rawData.SE_raw.Y(:))); 
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data_EBSD);
newEBSD = newEBSD.gridify;
newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');  

reforientation = orientation.byEuler(0*degree, 0*degree, 0*degree, newEBSD.CS);
newEBSD.orientations = newEBSD.orientations.project2FundamentalRegion(reforientation);

% Precompute Grain Boundary Topology
fprintf('  -> Precomputing grain boundary networks (>150 px, 5°)...\n');
[grains, ~] = calcGrains(newEBSD, 'minPixel', 150, 'angle', 5*degree);

%% 4. Initialize Interactive GUI Options & Batch Progression
opts = struct('exportDir', outputPath, ...
              'defaultIpfDir', 'X', 'defaultIpfKey', 'ipfColorKey', ...
              'defaultBgField', 'ci', 'invert_bg', false, ...
              'blackPoint', 0.05, 'whitePoint', 0.95, 'gammaVal', 1.5, ...
              'alpha_val', 0.65, 'defaultGbStyle', 'Black Line', ...
              'defaultIpfThreshProp', 'ci', ...
              'defaultOverlayMode', 'Scalar Property / Strain Map', ...
              'defaultPropOverlay', 'geff', ...
              'defaultColormap', 'viridis', ...
              'cmin_val', 0.0, 'cmax_val', 0.4);

% Extract DIC interpolation query grid
region = rawData.SE_raw.imregion;
startX = region(1); startY = region(2);
X_query_px = (Xnew / new_pixel_size) + startX - 1;
Y_query_px = (Ynew / new_pixel_size) + startY - 1;

inc_labels = {'0', '4', '6', '7', '8'};
dic_fields = {'/B00001', '/B00005', '/B00008', '/B00009', '/B00010'};
N = length(inc_labels);

%% 5. Loop Over Deformed Increments & Launch Interactive Studio
for i = 1:N
    inc_str = inc_labels{i};
    datasetName = sprintf('Re_PTP3_Inc_%s', inc_str);
    fprintf('\n========================================================================\n');
    fprintf('PREPARING INTERACTIVE STUDIO [%d/%d]: %s\n', i, N, datasetName);
    
    field = dic_fields{i};
    
    % Ingest raw HDF5 DIC displacement & coordinate arrays
    X_dic = double(h5read(DIC_file, [field,'/Coordinates/x']));
    Y_dic = double(h5read(DIC_file, [field,'/Coordinates/y']));
    u_dic = double(h5read(DIC_file, [field,'/Displacements/x-displacement']));
    v_dic = double(h5read(DIC_file, [field,'/Displacements/y-displacement']));
    
    Fxx = double(h5read(DIC_file, [field,'/Gradients/du÷dx - RC']));
    Fxy = double(h5read(DIC_file, [field,'/Gradients/du÷dy - RC']));
    Fyx = double(h5read(DIC_file, [field,'/Gradients/dv÷dx - RC']));
    Fyy = double(h5read(DIC_file, [field,'/Gradients/dv÷dy - RC']));
    
    % Interpolate DIC fields directly to undeformed EBSD spatial coordinate domain
    interpGrid = @(dataMap) griddata(X_dic(:), Y_dic(:), dataMap(:), X_query_px, Y_query_px, 'linear');
    
    U_interp = interpGrid(u_dic);
    V_interp = interpGrid(v_dic);
    Fxx_interp = interpGrid(Fxx);
    Fyy_interp = interpGrid(Fyy);
    Fxy_interp = interpGrid(Fxy);
    Fyx_interp = interpGrid(Fyx);
    
    % Compute smoothed strain tensors using windowed kernel
    data_disp.X = newEBSD.x;
    data_disp.Y = newEBSD.y;
    data_disp.U = U_interp(:) * new_pixel_size;
    data_disp.V = V_interp(:) * new_pixel_size;
    
    try
        [~, ~, straindata] = calcStrainsFromRawDisp(data_disp, filterSize);
        geff = sqrt(((straindata.Exx - straindata.Eyy)/2).^2 + straindata.Exy.^2);
        has_tensor = true;
    catch ME
        warning(ME.identifier, 'calcStrainsFromRawDisp fallback triggered (%s). Using raw Fxx/Fyy gradients.', ME.message);
        geff = sqrt(((Fxx_interp(:) - Fyy_interp(:))/2).^2 + ((Fxy_interp(:) - Fyx_interp(:))/2).^2);
        has_tensor = false;
    end
    
    % Clone EBSD object and attach computed scalar domain properties to .prop struct
    ebsd_inc = newEBSD;
    ebsd_inc.prop.geff   = reshape(geff, size(newEBSD));
    ebsd_inc.prop.U_disp = reshape(U_interp, size(newEBSD));
    ebsd_inc.prop.V_disp = reshape(V_interp, size(newEBSD));
    
    if has_tensor
        ebsd_inc.prop.Exx = reshape(straindata.Exx, size(newEBSD));
        ebsd_inc.prop.Eyy = reshape(straindata.Eyy, size(newEBSD));
        ebsd_inc.prop.Exy = reshape(straindata.Exy, size(newEBSD));
    else
        ebsd_inc.prop.Fxx = reshape(Fxx_interp, size(newEBSD));
        ebsd_inc.prop.Fyy = reshape(Fyy_interp, size(newEBSD));
        ebsd_inc.prop.Fxy = reshape(Fxy_interp, size(newEBSD));
    end
    
    % Configure studio execution parameters
    opts.datasetName = datasetName;
    opts.currentStep = i;
    opts.totalSteps  = N;
    
    % Launch Interactive Plotter Studio
    [status, opts] = Interactive_PrettyPlotter(ebsd_inc, grains, opts);
    
    if strcmpi(status, 'exit')
        fprintf('\n[Batch Aborted] User terminated DIC batch analysis at increment %s.\n', inc_str);
        break;
    elseif strcmpi(status, 'skip')
        fprintf('\n[Batch Skipped] Advancing directly to next increment without saving.\n');
        continue;
    else
        fprintf('\n[Batch Complete] Advancing to increment %d/%d while preserving interactive settings.\n', min(i+1, N), N);
    end
end

fprintf('========================================================================\n');
fprintf('ALL DIC FORWARD-DEFORMATION DATASETS PROCESSED SUCCESSFULLY!\n');
fprintf('========================================================================\n');
