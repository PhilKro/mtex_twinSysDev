%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data alignment and warping script for Re_PTP_cAxis_2                %%%
%%% Layers undeformed TKD data onto SE information and DIC info.        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

%% 0. Setup and Options

% Add path of src scripts (assuming they are reachable in your MTEX env)
addpath('src');

% MTEX Preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% Crystal symmetry for HCP Rhenium
CS = {'notIndexed', crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re')};

% Base output directory
outDir = 'C:\Users\phkr\MTEX_workEnv\onedriveData\Philipp\Re_cAxis_PTP2';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

options.in_situIncrementName = 'Re_cAxis_PTP2';

%% 1. Define Input Files
% TKD Map (Undeformed Reference)
ang_file = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0003_Re_cAxis_2\Experiment\Re_c_Axis_PTP_2\SphericalIndexing\Re_c_Axis_PTP_2 Specimen 2 Site 1 Map Data 2_BW123.ang';

% RAW SE Image
se_file = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240913_Re_c_tkdPTP_2\DIC\original\01.tif';

% DaVis DIC File
dic_file = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240913_Re_c_tkdPTP_2\DIC\DaVis\Processing_subset=63_step=20_01\B0001.hdf5';

% Note: Using the first DIC step (/B00001). Change to /B00002 if you wish 
% to project the deformed state strains onto the reference microstructure.
dic_field = '/B00002';

% Adjust this if your actual SE pixel size differs from the typical 5.4687 nm
new_pixel_size = 5.4687*10^-3; 

%% 2. Data Import
fprintf('Setting up data import structure...\n');

% Setup EBSD entry
DataName = 'EBSD0';
rawData.(DataName).fname = ang_file;
rawData.(DataName).dtype = 'EBSD';
rawData.(DataName).CS = CS;
rawData.(DataName).data = [];
% rawData.(DataName).rot = rotation.byAxisAngle(zvector, -90*degree);
clear DataName

% Setup SE entry
DataName = 'SE_raw';
rawData.(DataName).fname = se_file;
rawData.(DataName).psize = new_pixel_size;
rawData.(DataName).dtype = 'SEM';
% Define imregion if you need to crop the SE databar (e.g. [startx starty width height])
rawData.(DataName).imregion = [1 1 2048 2048]; 
clear DataName

rawData.dataNames = fieldnames(rawData);

% Execute generic loader
rawData = loadRawData(rawData);

% Rotate SE image 90 degrees clockwise (MTEX rotation objects don't apply to pure arrays)
rawData.SE_raw.data = rot90(rawData.SE_raw.data, -1);

if size(rawData.SE_raw.data, 1) ~= size(rawData.EBSD0.data, 2)
    warning('Not a square image adjusting grid') 
    % Re-calculate the grid based on the newly rotated image dimensions
    [rawData.SE_raw.m, rawData.SE_raw.n] = size(rawData.SE_raw.data);
    [rawData.SE_raw.Xp, rawData.SE_raw.Yp] = meshgrid(1:rawData.SE_raw.n, 1:rawData.SE_raw.m);
    rawData.SE_raw.X = (rawData.SE_raw.Xp - 1) * rawData.SE_raw.psize;
    rawData.SE_raw.Y = (rawData.SE_raw.Yp - 1) * rawData.SE_raw.psize;
end

%% 3. Alignment of EBSD to SE Data
fprintf('--- POINT SELECTION FOR ALIGNMENT ---\n');
REF = 'SE_raw';     % Undistorted SE image as reference
DEF = 'EBSD0';      % EBSD needs to be mapped to the SE

options.Ndof = 6;               % 1st order polynomial warping
options.Nsel = 4;               % Number of points to select
options.case = 'Re_cAxis';
options.selectionIPF = 1;       % Overlay IPF orientation for visual aid
options.EBSD_selection = 'ci';  % Use Confidence Index
options.colormap = 'gray';

% Note: Uncomment the following line if you have pre-selected points saved
options.sel_points_fname = 'C:\Users\phkr\MTEX_workEnv\onedriveData\Philipp\Re_cAxis_PTP2\Re_cAxis_PTP2_alignment_data_04-13-2026_13-37.mat';

alignment_data_EBSD = pointSelectionForAlignment(rawData.(REF), rawData.(DEF), options);

% Save the alignment selection
% save(fullfile(outDir, sprintf('%s_alignment_data_%s.mat', options.in_situIncrementName, datestr(now,'mm-dd-yyyy_HH-MM'))), 'alignment_data_EBSD');
clear options

%% 4. Warp EBSD onto the New SE Grid
fprintf('Warping EBSD onto SE grid...\n');
[Xnew, Ynew] = meshgrid(0:new_pixel_size:max(rawData.SE_raw.X(:)), 0:new_pixel_size:max(rawData.SE_raw.Y(:))); 

newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid, Xnew, Ynew, alignment_data_EBSD);
newEBSD = newEBSD.gridify;

% Merge SE into the EBSD variable properties
newEBSD = warpIncludeDataInMtex(newEBSD, rawData.SE_raw, 'SE_raw');
newEBSD.prop.SE_raw = reshape(newEBSD.prop.SE_raw, size(newEBSD));
%% 5. Extract DIC Data from DaVis HDF5
fprintf('Extracting DIC Data from %s...\n', dic_field);

% 5.1 Coordinates
X = double(h5read(dic_file, [dic_field, '/Coordinates/x']));
Y = double(h5read(dic_file, [dic_field, '/Coordinates/y']));

% 5.2 Displacements
displacements.u = double(h5read(dic_file, [dic_field, '/Displacements/x-displacement']));
displacements.v = double(h5read(dic_file, [dic_field, '/Displacements/y-displacement']));

% 5.3 Strains
strains.exx = double(h5read(dic_file, [dic_field, '/Strains/Transverse strain Ɛxx - RC']));
strains.exy = double(h5read(dic_file, [dic_field, '/Strains/Shear strain Ɛxy - RC']));
strains.eyy = double(h5read(dic_file, [dic_field, '/Strains/Axial strain Ɛyy - RC']));

% 5.4 Gradients
gradients.Fxx = double(h5read(dic_file, [dic_field, '/Gradients/du÷dx - RC']));
gradients.Fxy = double(h5read(dic_file, [dic_field, '/Gradients/du÷dy - RC']));
gradients.Fyx = double(h5read(dic_file, [dic_field, '/Gradients/dv÷dx - RC']));
gradients.Fyy = double(h5read(dic_file, [dic_field, '/Gradients/dv÷dy - RC']));

% 5.5 Other
try
    other.correlationValue = double(h5read(dic_file, [dic_field, '/Others/Correlation value']));
catch
    warning('Correlation value could not be loaded, skipping...');
end

% Calculate forward deformed positions
X_f = X + displacements.u;
Y_f = Y + displacements.v;

%% 6. Interpolate DIC Data onto EBSD Grid
fprintf('Interpolating DIC onto the combined EBSD/SE grid...\n');

if isfield(rawData.SE_raw, 'imregion') && ~isempty(rawData.SE_raw.imregion)
    startX = rawData.SE_raw.imregion(1);
    startY = rawData.SE_raw.imregion(2);
    endX = startX + rawData.SE_raw.imregion(3);
    endY = startY + rawData.SE_raw.imregion(4);
else
    startX = 1;
    startY = 1;
    endX = size(rawData.SE_raw.data, 2);
    endY = size(rawData.SE_raw.data, 1);
end

% Restrict bounds securely to the new MTEX EBSD matrix size
cropX = startX:endX;
cropY = startY:endY;

cropX = cropX(cropX <= size(Xnew, 2));
cropY = cropY(cropY <= size(Xnew, 1));

% Physical coordinates to store in the struct
X_crop = Xnew(cropY, cropX);
Y_crop = Ynew(cropY, cropX);

% The DIC data from DaVis is calculated in PIXELS of the reference image.
% Since the MTEX grid (X_crop, Y_crop) is in physical units (microns), we must 
% convert it to pixel coordinates to correctly map the subset step size.
X_crop_px = X_crop / new_pixel_size;
Y_crop_px = Y_crop / new_pixel_size;

% Helper function for safe, linear interpolation (flatten inputs to ignore orientation mismatches)
interpGrid = @(dataMap) griddata(X_f(:), Y_f(:), dataMap(:), X_crop_px, Y_crop_px, 'linear');

% Execute interpolations
DICdata_interp = struct('X', X_crop', 'Y', Y_crop');

fieldGroups = {displacements, strains, gradients};
for g = 1:length(fieldGroups)
    fields = fieldnames(fieldGroups{g});
    for i = 1:length(fields)
        f = fields{i};
        DICdata_interp.(f) = interpGrid(fieldGroups{g}.(f))';
    end
end

if exist('other', 'var') && isfield(other, 'correlationValue')
    DICdata_interp.CorrValue = interpGrid(other.correlationValue)';
end

% Calculate Effective Shear Strain
if isfield(DICdata_interp, 'Fxx')
    DICdata_interp.geff = sqrt(((DICdata_interp.Fxx - DICdata_interp.Fyy)/2).^2 + ...
                               ((DICdata_interp.Fxy - DICdata_interp.Fyx)/2).^2);
end

% Inject into newEBSD properties
new_fields = fieldnames(DICdata_interp);

% Loop through and append each one to the EBSD object
for i = 1:length(new_fields)
    fname = new_fields{i};
    newEBSD.prop.(fname) = DICdata_interp.(fname);
end

%% 7. Verification Plotting
fprintf('Generating verification plots...\n');
% Basic CI Filter and grain generation
newEBSD(newEBSD.ci < 0.13) = 'notIndexed';
[newGrains] = calcGrains(newEBSD, 'angle', 10*degree); 

% Ensure SE_raw image is suitable for rgb overlay
SE_gray = mat2gray(medfilt2(double(newEBSD.prop.SE_raw), [3 3])); 
SE_rgb = repmat(SE_gray, 1, 3);

%%% Plot 1: Strain Overlay on SE
fig1 = figure('Name', 'Effective Strain');
plot(newEBSD, SE_rgb); 
axis image off; hold on;

if isfield(newEBSD.prop, 'geff')
    gefTest = newEBSD.prop.geff;
    % gefTest(newEBSD.prop.CorrValue < 0.8) = NaN; % Example threshold filter
    
    plot(newEBSD, gefTest);%, 'FaceAlpha', 0.65);
    colormap('viridis');
    mtexColorbar('title', 'Effective Shear Strain');
else
    title('No Gradient Data to display strains');
end
plot(newGrains.boundary, 'lineColor', 'w', 'linewidth', 0.5);

%%% Plot 2: IPF Overlay on SE
fig2 = figure('Name', 'IPF Map Overlay');
plot(newEBSD, SE_rgb); 
axis image off; hold on;

plot(newEBSD, newEBSD.orientations, 'FaceAlpha', 0.6);
plot(newGrains.boundary, 'lineColor', 'k', 'linewidth', 1);
title('IPF map projected over SE Micrograph');

% Export results visually to the outDir
exportgraphics(fig1, fullfile(outDir, sprintf('StrainMap_SE_Overlay_%s.png', datestr(now,'mm-dd-yyyy_HH-MM'))), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, sprintf('IPF_SE_Overlay_%s.png', datestr(now,'mm-dd-yyyy_HH-MM'))), 'Resolution', 300);

%% 8. Save 
fprintf('Saving results ...\n');

% Save combined dataset
dic_increment = strrep(dic_field, '/', '');
saveFile = fullfile(outDir, sprintf('%s_AlignedData_%s_%s.mat', alignment_data_EBSD.options.in_situIncrementName, dic_increment, datestr(now,'mm-dd-yyyy_HH-MM')));
save(saveFile, 'newEBSD');
fprintf('Data successfully saved to %s\n', saveFile);


fprintf('Done.\n');
