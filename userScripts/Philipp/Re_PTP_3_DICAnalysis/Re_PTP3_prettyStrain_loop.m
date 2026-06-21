%% Re_PTP3_prettyStrain_loop.m
% Script to loop over Re_PTP_3 increments and generate pretty plots 
% of the effective shear strain overlaid on a tunable grayscale 
% Confidence Index (CI) background of the undeformed specimen.

clear; close all; clc;

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% add path of src scripts
addpath('userScripts\NanoMech_Aignment_Matlab\src');

%% 2. Directories and Tuning Parameters
inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
outputPath = fullfile(inputPath, 'PrettyPlots_Strain_Undeformed');

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    fprintf('Created output directory: %s\n', outputPath);
end

DIC_file = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20241218_Re_7dc_tkdPTP_3\DaVis\10kV_Images\export\Processing_subset=19_step=6_relativeToFirst_median3\1st_OrderShapeFun.hdf5';
sel_points_fname = char(fullfile(inputPath, 'inc_0_Selection_Points_Backup_EBSD_05-08-2025_10-23.mat'));

% EBSD and SE inputs for Undeformed Specimen (inc_0)
ang_file = fullfile(inputPath, 'EBSD_angs', 'Re_7deg_insitu Specimen 1 0 Map Data 1_BW158.ang');
se_file = fullfile(inputPath, 'in_Lens_SE', '00_DIC_0.tif');

CS = {'notIndexed', crystalSymmetry('6/mmm', [2.761 2.761 4.458], 'x||a', 'mineral', 'Re')};
new_pixel_size = 5.4687*10^-3;

% --- "PHOTOSHOP CURVES" TUNING PARAMETERS (for CI) ---
% 1. Black Point (0.0 to 1.0): CI values below this will be pure black.
blackPoint = 0.05; 

% 2. White Point (0.0 to 1.0): CI values above this will be pure white.
whitePoint = 0.98; 

% 3. Gamma (Midtones): 
gammaVal = 1.5; 

% --- DIC STRAIN TUNING PARAMETERS ---
strain_alpha = 0.65;          % Opacity of the strain overlay
strainColorRange = [0, 0.4];  % Color range for Effective Shear Strain
filterSize = 3;               % Filter size for calcStrainsFromRawDisp
strainthresh = 1;             % Cutoff for bad strain pixels

%% 3. Load and Align Undeformed EBSD and SE Data
fprintf('Loading and aligning undeformed data...\n');

% Setup EBSD entry
DataName = 'EBSD0';
rawData.(DataName).fname = ang_file;
rawData.(DataName).dtype = 'EBSD';
rawData.(DataName).CS = CS;
rawData.(DataName).data = [];

% Setup SE entry
DataName = 'SE_raw';
rawData.(DataName).fname = se_file;
rawData.(DataName).psize = new_pixel_size;
rawData.(DataName).dtype = 'SEM';
rawData.(DataName).imregion = [1281 1 1279 1279]; % startx starty numx-1 numy-1

rawData.dataNames = fieldnames(rawData);
rawData = loadRawData(rawData);

% Alignment
options.Ndof = 6;
options.Nsel = 4;
options.case = 'EBSD';
options.selectionIPF = 1;
options.EBSD_selection = 'ci';
options.colormap = 'gray';
options.sel_points_fname = sel_points_fname;

alignment_data_EBSD = pointSelectionForAlignment(rawData.SE_raw, rawData.EBSD0, options);

% Warp to new grid
[Xnew, Ynew] = meshgrid(0:new_pixel_size:max(rawData.SE_raw.X(:)), 0:new_pixel_size:max(rawData.SE_raw.Y(:))); 
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data_EBSD);
newEBSD = newEBSD.gridify;
newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');  

reforientation = orientation.byEuler(0*degree, 0*degree, 0*degree, newEBSD.CS);
newEBSD.orientations = newEBSD.orientations.project2FundamentalRegion(reforientation);

%% 4. Prepare Undeformed CI map
fprintf('\nPreparing undeformed CI map...\n');

ci_raw = newEBSD.ci;
ci_raw(isnan(ci_raw)) = 0;
ci_min = min(ci_raw(:));
ci_max = max(ci_raw(:));

% Normalize and apply curve adjustments
ci_norm = (ci_raw - ci_min) ./ (ci_max - ci_min);
ci_tuned = (ci_norm - blackPoint) ./ (whitePoint - blackPoint);
ci_tuned(ci_tuned < 0) = 0; 
ci_tuned(ci_tuned > 1) = 1; 
ci_tuned = ci_tuned .^ gammaVal;

% Create RGB Nx3 array
cibw_undef = [ci_tuned(:), ci_tuned(:), ci_tuned(:)];
cibw_undef(~newEBSD.isIndexed, :) = 1; % Set non-indexed background to white

%% 4b. Extract Colormap for Manual RGB Blending
fig_tmp = figure('Visible', 'off');
mtexColorMap viridis;
cmap_strain = colormap(fig_tmp);
close(fig_tmp);

%% 5. Loop over Increments, Interpolate DIC, Calculate Strains and Plot
inc_labels = {'0', '4', '6', '7', '8'};
dic_fields = {'/B00001', '/B00005', '/B00008', '/B00009', '/B00010'};
N = length(inc_labels);

% Extract DIC interpolation setup
region = rawData.SE_raw.imregion;  % [startX startY width height]
startX = region(1);
startY = region(2);

X_query_px = (Xnew / new_pixel_size) + startX - 1;
Y_query_px = (Ynew / new_pixel_size) + startY - 1;

for i = 1:N
    fprintf('Processing DIC data for increment %s...\n', inc_labels{i});
    field = dic_fields{i};
    
    % Read DIC data
    X_dic = double(h5read(DIC_file, [field,'/Coordinates/x']));
    Y_dic = double(h5read(DIC_file, [field,'/Coordinates/y']));
    
    u_dic = double(h5read(DIC_file, [field,'/Displacements/x-displacement']));
    v_dic = double(h5read(DIC_file, [field,'/Displacements/y-displacement']));
    
    % Get displacement gradients calculated by DaVis (fallback)
    Fxx = double(h5read(DIC_file, [field,'/Gradients/du÷dx - RC']));
    Fxy = double(h5read(DIC_file, [field,'/Gradients/du÷dy - RC']));
    Fyx = double(h5read(DIC_file, [field,'/Gradients/dv÷dx - RC']));
    Fyy = double(h5read(DIC_file, [field,'/Gradients/dv÷dy - RC']));
    
    % Interpolate undeformed DIC fields directly to the undeformed EBSD grid
    interpGrid = @(dataMap) griddata(X_dic(:), Y_dic(:), dataMap(:), X_query_px, Y_query_px, 'linear');
    
    U_interp = interpGrid(u_dic);
    V_interp = interpGrid(v_dic);
    
    Fxx_interp = interpGrid(Fxx);
    Fyy_interp = interpGrid(Fyy);
    Fxy_interp = interpGrid(Fxy);
    Fyx_interp = interpGrid(Fyx);
    
    % Recalculate strains using filter
    data_disp.X = newEBSD.x;
    data_disp.Y = newEBSD.y;
    data_disp.U = U_interp(:) * new_pixel_size; % Convert from pixels to physical units
    data_disp.V = V_interp(:) * new_pixel_size;

    try
        [~, ~, straindata] = calcStrainsFromRawDisp(data_disp, filterSize);
        % Effective shear strain formula
        geff = sqrt(((straindata.Exx - straindata.Eyy)/2).^2 + straindata.Exy.^2);
    catch ME
        warning(ME.identifier, 'calcStrainsFromRawDisp failed (%s). Falling back to interpolated Fxx/Fyy/Fxy/Fyx fields.', ME.message);
        geff = sqrt(((Fxx_interp(:) - Fyy_interp(:))/2).^2 + ((Fxy_interp(:) - Fyx_interp(:))/2).^2);
    end
    
    plotData = geff;
    
    % Clean up bad strain data
    % plotData(plotData > strainthresh) = NaN;
    % plotData(plotData == 0) = NaN;
    % plotData(~newEBSD.isIndexed) = NaN; % Mask out non-indexed strain regions
    
    % --- Manual Alpha Blending ---
    % Normalize strain to [0, 1] based on strainColorRange
    normStrain = (plotData(:) - strainColorRange(1)) / (strainColorRange(2) - strainColorRange(1));
    normStrain(normStrain < 0) = 0;
    normStrain(normStrain > 1) = 1;
    
    % Map normalized strain to colormap indices
    numColors = size(cmap_strain, 1);
    colorIndices = round(normStrain * (numColors - 1)) + 1;
    
    % Handle NaNs in strain data
    validMask = ~isnan(plotData(:));
    colorIndices(~validMask) = 1; % Dummy index for NaNs to prevent errors
    
    % Extract RGB colors for the strain
    strainRGB = cmap_strain(colorIndices, :);
    
    % Combine CI background and Strain RGB via alpha blending
    combinedRGB = cibw_undef;
    combinedRGB(validMask, :) = (1 - strain_alpha) * cibw_undef(validMask, :) + strain_alpha * strainRGB(validMask, :);
    
    % --- Create Figure ---
    fig_strain = figure('Name', ['Effective Shear Strain - Inc ', inc_labels{i}], 'Visible', 'off');
    
    plot(newEBSD, combinedRGB, 'micronbar', 'on'); % Plot combined RGB
    
    mtexColorMap viridis;
    setColorRange(strainColorRange);
    mtexColorbar('title', 'Effective Shear Strain');
    title(sprintf('Effective Shear Strain - Increment %s', inc_labels{i}));
    
    saveName = fullfile(outputPath, sprintf('Re_PTP3_inc_%s_Strain_on_UndeformedCI.jpg', inc_labels{i}));
    exportScaledFigure(fig_strain, saveName, 'noDateString');
    close(fig_strain);
end

fprintf('\nFinished generating pretty strain plots for Re_PTP_3 (Undeformed Mapping)!\n');