%% clear workspace

clear 
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;
%% setup input parameters and paths
% set up defaults for any missing input fields

% defaults
defaults.dir            = fullfile('onedriveData','Philipp','Nb_PTP');
defaults.SEImagefolder  = fullfile('DIC_SE_Images');          % folder under input.dir
defaults.startx         = 1280;                             % crop region start x
defaults.starty         = 1;                            % crop region start y
defaults.numx           = 1280;                         % crop region size x
defaults.numy           = 1280;                         % crop region size y
defaults.bandWidth      = 263;                          % BW of the spherical Indexing
defaults.tkd_tilt = -10 * degree;                       % define tilt for TKD
defaults.thresh_ci = 0.01;                              % define threshold for CI, 0.15 is usually safe
defaults.DIC_filename = 'Processing_subset=21_step=6';


% ----- Input structs -----

%     % Increment 0
% input = struct( ...
%     'fname',             'insitu Specimen 1 0 Map Data 1', ...    % base EBSD filename (no extension)                   
%     'SEImagefname',      '0_10kv_DIC.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0001.hdf5',... % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname', 'insitu Specimen 1 0 Map Data 1_BW263_20251126_115536_Selection_Points_Backup_EBSD.mat'); % mat file with preselected points for alignment                     

%     % Increment 1
% input = struct( ...
%     'fname',             'insitu Specimen 1 1 Map Data 2', ...    % base EBSD filename (no extension)                   
%     'SEImagefname',      '1postTKD_10kv_DIC_0d.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0003.hdf5',...                       % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname', 'insitu Specimen 1 1 Map Data 2_BW263_20251126_132847_Selection_Points_Backup_EBSD.mat'); % mat file with preselected points for alignment 


%      % Increment 2
% input = struct( ...
%     'fname',             'insitu Specimen 1 2 Map Data 3', ...    % base EBSD filename (no extension)                   
%     'SEImagefname',      '2postTKD_10kv_DIC_0d.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0005.hdf5'                       % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname', 'insitu Specimen 1 2 Map Data3_BW263_20251126_134943_Selection_Points_Backup_%sEBSD.mat'); % mat file with preselected points for alignment 

%     % Increment 3
% input = struct( ...
%     'fname',             'insitu Specimen 1 3 Map Data 4', ...    % base EBSD filename (no extension)                   
%     'SEImagefname',      '3postTKD_10kv_DIC_0.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', 'Processing_subset=19_step=1'),... % DIC HDF5 filename
%     'DICfname',          'B0007.hdf5',...                       % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname',  'insitu Specimen 1 3 Map Data 4_BW263_20251127_154633_Selection_Points_Backup_%sEBSD'); % mat file with preselected points for alignment 

%     % Increment 4
% input = struct( ...
%     'fname',             'insitu Specimen 1 4 Map Data 5', ...    % base EBSD filename (no extension)                   
%        'SEImagefname',     '4postTKD_10kv_DIC_0.tif',...%'SEImagefname',      '4postTKD_10kv_DIC_0.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0009.hdf5',... % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname', 'insitu Specimen 1 4 Map Data 5_BW263_20251201_092539_Selection_Points_Backup_%sEBSD.mat'); % mat file with preselected points for alignment                    

%     % Increment 5
% input = struct( ...
%     'fname',             'insitu Specimen 1 5 Map Data 6', ...    % base EBSD filename (no extension)                   
%     'SEImagefname',      '5postTKD_Take2_10kv_DIC_0.tif', ...          % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0011.hdf5',...                       % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname',  'insitu Specimen 1 5 Map Data 6_BW263_20251201_100658_Selection_Points_Backup_%sEBSD.mat'); % mat file with preselected points for alignment 

%         % Increment 6
% input = struct( ...
%     'fname', 'insitu Specimen 1 7 Map Data 7', ...                  % base EBSD filename (no extension)                   
%     'SEImagefname', '7b_afterDescend_10kv_DIC_0.tif',...                            % SE image filename
%     'DICfoldername' ,    fullfile('DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
%     'DICfname',          'B0014.hdf5',...                       % folder under defaults.dir containing DIC .h5 file
%     'sel_points_fname',  ''); % mat file with preselected points for alignment 


%% import setup: merge input with defaults and save to JSON
% set up defaults for any missing input fields

% merge provided input with defaults
fn = fieldnames(defaults);
for k = 1:numel(fn)
    if ~isfield(input, fn{k}) || isempty(input.(fn{k}))
        input.(fn{k}) = defaults.(fn{k});
    end
end

% construct full paths
input.SEImageFullPath = fullfile(input.dir, input.SEImagefolder, input.SEImagefname);
input.DIC_path  = fullfile(input.dir, input.DICfoldername, input.DICfname);
input.ang_path        = fullfile(input.dir, [input.fname,'_BW',num2str(defaults.bandWidth),'.ang']);
input.ctf_path = fullfile(input.dir, [input.fname,'.ctf']);

options.dateString = datestr(now,'yyyymmdd_HHMMSS');
options.in_situIncrementName = [input.fname,'_BW',num2str(defaults.bandWidth),'_',options.dateString];

% (Optional) show summary
fprintf('Using IO config:\n  dir: %s\n  fname: %s\n  SE image: %s\n  crop: [%d %d %d %d]\n  DIC: %s (path: %s)\n', ...
    input.dir, input.fname, input.SEImageFullPath, input.startx, input.starty, input.numx, input.numy, input.DICfname, input.DIC_path);

% save ioConfig to a JSON file with datestring in the filename
outname = sprintf([options.in_situIncrementName, '_ioConfig.json']);

jsonText = jsonencode(input);
% write JSON to file
fid = fopen(outname,'w');
if fid < 0
    error('Could not open %s for writing.', outname);
end
fprintf(fid, '%s', jsonText);
fclose(fid);
fprintf('Saved input to %s\n', outname);


% Note:
% - Replace the hardcoded occurrences later in the script with the variables:
%     SEImageFullPath  -> rawData.<name>.fname
%     startx,starty,numx,numy -> rawData.<name>.imregion
%     DIC_fname and DICFullFolder -> use when reading DIC HDF5
% The script already uses variables named 'dir' and 'fname', so those are set.

%% import EBSD data and perform orientation corrections
% crystal symmetry, make sure to use correct space group and c/a ratio

%Ti
% CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

%Re
% CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');
%CS = crystalSymmetry.load("Rhenium.cif");

%nice to see bending and warping
%CS = crystalSymmetry('622',[2.761 2.761 4.458],'x||a','mineral','Re'); %



%Mg
%CS = crystalSymmetry('6/mmm',[1.570796370506287 1.570796370506287 2.094395102393195],'x||a','mineral','Magnesium');
%CS = crystalSymmetry.load("Mg-Magnesium.cif");

% Nb https://www.crystallography.net/cod/9008546.html
% CS = crystalSymmetry.load('m3m', [3.3004,3.3004,3.3004], 'mineral','Niobium');
CS = crystalSymmetry.load("Nb-Niobium.cif");

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(input.ang_path,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + defaults.tkd_tilt),'keepXY');
ebsd(ebsd.ci < defaults.thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

% Project to fundamental 
targetFile = fullfile(input.dir, 'refOri.mat');
try
    % Load only the variable 'refOrientation' from the specific file
    data = load(targetFile, 'refOrientation');
    refOrientation = data.refOrientation;

    % Apply the orientation projection using the loaded data
    ebsd.orientations = ebsd.orientations.project2FundamentalRegion(refOrientation);

catch
    % If the file is missing or 'refOrientation' is not inside it, trigger this error
    error('Could not find the file at %s. Please use the same reference orientation for all imports.', targetFile);
end

% import the ctfs to extract BC and MAD



ctf_ebsd = EBSD.load(input.ctf_path,CS,'convertEuler2SpatialReferenceFrame','setting 2');

props = {'bands','bc','bs','error','mad'};

for i = 1:numel(props)
    propName = props{i};
    ebsd.prop.(propName) = ctf_ebsd.prop.(propName);
end

clear ctf_ebsd

% ebsd0=ebsd;
% ebsd=ebsd.gridify;
% ebsd0=ebsd0.gridify;

%% test plots of ebsd data
ipf=ipfColorKey(ebsd.CS);

figure;

ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors)
title('IPF Y')


nextAxis
plot(ebsd,ebsd.prop.ci)
title('Spherical Indexing CI')
colormap gray

nextAxis
plot(ebsd,ebsd.prop.bc)
title('Band Contrast')

nextAxis
plot(ebsd,ebsd.prop.bs)
title('Band Slope')

nextAxis
plot(ebsd,ebsd.prop.error)
title('Error')

nextAxis
plot(ebsd,ebsd.prop.mad)
title('Mean Angular Deviation')

nextAxis
plot(ebsd,ebsd.prop.bands)
title('Bands')

exportScaledFigure(gcf,[input.fname,options.dateString,'ebsd_.png'],'nodatestring', 'extrawidth')

%%
% -------------------------------------------------------------------------
% --- SE / DIC loading, EBSD to SE alignment, verification plotting -----
% -------------------------------------------------------------------------
% The following section is adapted from `Point_Alignment_Example_Zn.m`.
% Edit the `rawData` filenames and pixel sizes below to match your Nb dataset.

% Define rawData entries for SEM / SE and EBSD (use existing `ebsd0` as grid)
DataName = 'SE_raw';
% EDIT: set the correct SE filename inside your `data\Nb_PTP` folder
rawData.(DataName).fname = fullfile(input.SEImageFullPath);   
rawData.(DataName).psize = 5.4687e-3;                               % pixel size in um (edit if needed)
rawData.(DataName).dtype = 'SEM';
% Optional: define imregion if only a subregion of the tif should be used
rawData.(DataName).imregion = [input.startx input.starty input.numx-1 input.numy-1];
clear DataName

% Add the EBSD grid as an entry so the alignment helper functions can use it
DataName = 'EBSD0';
rawData.(DataName).fname = input.fname;
rawData.(DataName).data = ebsd;
% rawData.(DataName).data_grid = ebsd0;   % `ebsd0` was created earlier in this script
rawData.(DataName).dtype = 'EBSD';
clear DataName

% store list of data names
rawData.dataNames = fieldnames(rawData);

% Load raw data using existing helper (will populate X,Y and image fields)
rawData = loadRawData(rawData);

%% ----------------------
% Align EBSD to SE (reference = SE_raw, deformed = EBSD0)
REF = 'SE_raw';
DEF = 'EBSD0';

options.Ndof = 6;            % polynomial degrees of freedom (6 -> 1st order)
options.Nsel = 4;            % number of selection points (>= Ndof*0.5 recommended)
options.case = 'EBSD';
options.selectionIPF = 1;    % show IPF during point selection
options.EBSD_selection = 'ci';
options.colormap = 'gray';

% If you have a pre-saved selection point file, set options.sel_points_fname

if ~isempty(input.sel_points_fname)
    options.sel_points_fname = input.sel_points_fname;
end


alignment_data_EBSD = pointSelectionForAlignment(rawData.(REF),rawData.(DEF),options);

if ~exist('SelectionFnameEBSD','var')
    save([options.in_situIncrementName, '_', options.case,'_alignment_data_EBSD'],'alignment_data_EBSD')
end

%% ----------------------
% Create a new grid matched to the SE reference and warp EBSD onto it
psize_new = rawData.SE_raw.psize;  % use SE pixel size (edit as needed)

[Xnew,Ynew] = meshgrid(0:psize_new:max(rawData.SE_raw.X(:)),0:psize_new:max(rawData.SE_raw.Y(:))); 

newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid,Xnew,Ynew,alignment_data_EBSD);

%save the id state so the ebsd objects can be refenced with each other
newEBSD.prop.originalId = newEBSD.prop.oldId;
originalEBSD = ebsd;
originalEBSD.prop.originalId = ebsd.id;

% Include the SE image into the new EBSD MTEX variable for plotting/verification
newEBSD = warpIncludeDataInMtex(newEBSD,rawData.SE_raw,'SE_raw');


newEBSD = newEBSD.gridify;

% If DIC data was provided earlier and loaded into rawData.DIC, include it too
if isfield(rawData,'DIC')
    newEBSD = warpIncludeDataInMtex(newEBSD,rawData.DIC,'DIC');
end

%% ----------------------
% Verification plotting: overlay EBSD-derived grain boundaries on SE image
newGrains = calcGrains(newEBSD,'angle',15*degree, 'minPixel',5);

figure;
h = newMtexFigure('nrows',1,'ncols',1,'visible','on','Position',get(0,'ScreenSize'));
plot(newEBSD,newEBSD.prop.SE_raw,'micronbar','off')
hold on
plot(newGrains.boundary,'lineColor','red','linewidth',1.5)
title('Aligned EBSD boundaries over SE image')
exportScaledFigure(gcf,[input.fname,options.dateString,'AlignedEBSD_SE_.png'],'extraheight','nodatestring')
% %% ---------------------
% % Save aligned newEBSD for later use
% save([options.in_situIncrementName,'_aligned_newEBSD'], 'newEBSD', 'originalEBSD')

%% ---------------------
% DIC DATA IMPORT SECTION
%
% This block reads DIC data from the .h5 file, detects orientation,
% replicates low-resolution DIC fields to full SE resolution,
% and crops the region of interest.
% The goal is to map DIC data onto the same geometry as the SE image,
% so it can later be forward-deformed and interpolated to the EBSD grid.
% ---------------------------------------------------------------

% Read HDF5 structure
info = h5info(input.DIC_path);

% Target size defined by SE image
targetM = rawData.SE_raw.origM;
targetN = rawData.SE_raw.origN;
targetSize = [targetM, targetN];

% If the SE image is square, orientation matching cannot be inferred
if targetSize(1) == targetSize(2)
    warning('Can not determine orientation of SE image to DIC data, because dataset is square. May be transposed!')
end

% Recursively collect all dataset paths in the HDF5 file
datasetPaths = collectDatasets(info);

% Detect coordinate datasets for determining orientation
xcoordPaths = datasetPaths(cellfun(@(x) endsWith(x, '/Coordinates/x'), datasetPaths));
ycoordPaths = datasetPaths(cellfun(@(x) endsWith(x, '/Coordinates/y'), datasetPaths));
if isempty(xcoordPaths) || isempty(ycoordPaths)
    error('Could not find coordinate field ending with /Coordinates/x /Coordinates/y in the HDF5 file');
end

% Load coordinate grids
xcoord = h5read(input.DIC_path, xcoordPaths{1});
ycoord = h5read(input.DIC_path, ycoordPaths{1});
dataSize = size(xcoord);

% Transpose coordinate grid if its aspect ratio mismatches the SE image
if xor(dataSize(1) > dataSize(2), targetSize(1) > targetSize(2))
    ycoord = ycoord';
    xcoord = xcoord';
    dataSize = dataSize';
end

% Determine the DIC grid spacing (pixel pitch)
stepsize = abs(xcoord(2,1)-xcoord(1,1));
if stepsize == 0, stepsize = abs(xcoord(1,2)-xcoord(1,1)); end
if stepsize == 0
    error('Could not determine stepsize for coordinate replication. Step size is zero.');
end

% Prepare container
rawData.DIC = struct();
rawData.DIC.names = {};

% Region of interest in SE/DIC images (1-based MATLAB indexing)
x0 = input.startx;
y0 = input.starty;
x1 = input.startx + input.numx - 1;
y1 = input.starty + input.numy - 1;

% -------------------------------------------------------------
% Read each dataset and decide:
%  • If it matches SE resolution → crop to ROI
%  • If smaller (i.e. DIC stepsize > 1) → replicate (nearest neighbor) using kron → crop to ROI
%  • Coordinate datasets get synthetic index-maps → crop to ROI
% -------------------------------------------------------------
for k = 1:numel(datasetPaths)
    dsPath = datasetPaths{k};

    % Attempt to read dataset
    try
        data = h5read(input.DIC_path, dsPath);
        dataSize = size(data);

        % Transpose if orientation mismatches SE
        if xor(dataSize(1) > dataSize(2), targetSize(1) > targetSize(2))
            data = data';
            dataSize = dataSize';
        elseif and(dataSize(1) == dataSize(2), targetSize(1) == targetSize(2))
            warning('Data and target size are square. The DIC data may be transposed.')
        end
        
    catch ME
        warning('Could not read dataset %s: %s', dsPath, ME.message);
        continue
    end

    % Only process 2D numeric data
    if ~ismatrix(data) || ~isnumeric(data)
        continue
    end

    sz = size(data);

    if isequal(sz, targetSize)
        % Full-size → direct crop of ROI
        subdata = extractRegion(data, x0, y0, x1, y1, dsPath);

    else
        % Smaller dataset → need to upscale to SE resolution

        % Coordinate fields get synthetic grids based on 1...N index maps
        if contains(dsPath, '/Coordinates/x') || contains(dsPath, '/Coordinates/y')

            % Detect orientation of the grid sweep
            if data(1,1) < data(2,1)
                % Column-wise sweep
                subdata = repmat((1:targetM), targetN,1)';
            elseif data(1,1) < data(1,2)
                % Row-wise sweep
                subdata = repmat((1:targetN)', 1,targetM)';
            else
                error('SUMTINGWONG, your coordinates are not increasing in any direction of the 2D array');
            end

            % Crop synthetic grid to ROI
            if size(subdata,1) >= y1 && size(subdata,2) >= x1
                subdata = subdata(y0:y1, x0:x1);
            end

        else
            % Displacement fields or similar small data → upscale via kron
            % stepsize tells how many SE pixels correspond to one DIC point
            subdata = kron(data, ones(stepsize));

            % Crop to ROI
            if size(subdata,1) >= y1 && size(subdata,2) >= x1
                subdata = subdata(y0:y1, x0:x1);
            end

            % Status message for debugging
            fprintf('Replicated %s from size [%d %d] with stepsize %d to [%d %d]\n', ...
                dsPath, size(data,1), size(data,2), stepsize, size(subdata,1), size(subdata,2));
        end
    end

    % Generate MATLAB-friendly field name
    fname = matlab.lang.makeValidName(strrep(dsPath, '/', '_'));

    % Attempt to strip unnecessary prefixes
    tok = regexp(fname, '_B\d+_(.+)$', 'tokens', 'once');
    if ~isempty(tok)
        tok = tok{:};
    end
    if ~isempty(tok) && ischar(tok)
        shortName = tok;
    else
        parts = strsplit(fname, '_');
        if numel(parts) >= 3
            shortName = strjoin(parts(3:end), '_');
        else
            shortName = fname;
        end
    end

    % Store dataset
    rawData.DIC.(shortName) = subdata;
    rawData.DIC.names{end+1} = shortName;
end

% -------------------------------------------------------------------
% FORWARD DEFORMATION OF DIC DATA
% Move undeformed DIC pixels into their deformed positions:
%
%     (x_def, y_def) = (x + Ux, y + Uy)
%
% Then interpolate irregular forward-mapped DIC field onto the regular
% EBSD grid. A distance-based mask excludes unsupported areas.
% -------------------------------------------------------------------

requiredFields = {'Coordinates_x','Coordinates_y','Displacements_x_displacement','Displacements_y_displacement'};
if ~all(isfield(rawData.DIC, requiredFields))
    error('Missing required DIC fields: %s', strjoin(requiredFields, ', '));
end

% Output container
rawData.forwardInterpolatedDIC = struct();
rawData.forwardInterpolatedDIC.names = {};

% Extract undeformed grid
x = double(rawData.DIC.Coordinates_x);
y = double(rawData.DIC.Coordinates_y);

% Compute forward-deformed coordinates
fx = x + double(rawData.DIC.Displacements_x_displacement);
fy = y + double(rawData.DIC.Displacements_y_displacement);

rawData.forwardInterpolatedDIC.Coordinates_x = fx;
rawData.forwardInterpolatedDIC.Coordinates_y = fy;

% Skip coordinate fields when interpolating
skipFields = {'Coordinates_x','Coordinates_y'};

% ---------------------------------------------------------------
% Create a distance threshold mask so scatteredInterpolant doesn't
% produce garbage far away from real data.
%
% dmax = maximum allowed distance from nearest DIC data point
% ---------------------------------------------------------------
dmax = 5;

% Points where real displacement exists (ignore zeros)
mask = rawData.DIC.Displacements_x_displacement ~= 0;

% Data points used by the interpolant
xy_data = [fx(mask), fy(mask)];

% Query points = entire EBSD grid (before forward deformation)
xy_query = [x(:),  y(:)];

% Compute nearest-neighbor distance for each EBSD gridpoint
[~, dist] = knnsearch(xy_data, xy_query);

% Distance map in EBSD grid shape
distGrid = reshape(dist, size(x));

% Chugg it into a EBSD prop
newEBSD.prop.DIC_distance_from_actual_data = distGrid;

fprintf('Interpolation always takes some time. GO GET A COFFEE!\n')

% ---------------------------------------------------------------
% Forward interpolation for each DIC field
% ---------------------------------------------------------------
for i = 1:numel(rawData.DIC.names)
    fname = rawData.DIC.names{i};
    if ismember(fname, skipFields)
        continue
    end

    dataMat = rawData.DIC.(fname);
    if ~ismatrix(dataMat) || ~isnumeric(dataMat)
        warning('Skipping %s: not a numeric 2D array.', fname);
        continue
    end

    try
        % Flatten data
        V = double(dataMat(:));

        % Build scattered interpolant in forward-deformed coordinates
        F = scatteredInterpolant( ...
            fx(mask), fy(mask), V(mask), ...
            'linear', 'none');

        % Sample interpolant on EBSD grid
        V_deformed = F(x, y);

        % Remove values farther than dmax from any known DIC point
        V_deformed(distGrid > dmax) = NaN;

        rawData.forwardInterpolatedDIC.(fname) = V_deformed;
        rawData.forwardInterpolatedDIC.names{end+1} = fname;

    catch ME
        warning('Interpolation failed for %s: %s', fname, ME.message);
    end
end

% -----------------------------------------------------------------
% INSERT FORWARD INTERPOLATED DIC FIELDS INTO newEBSD.prop
%
% This matches the DIC grid to the EBSD grid:
%   - direct match
%   - transpose match
%   - reshaping if the total number of elements fits
% -----------------------------------------------------------------

dicNames = rawData.forwardInterpolatedDIC.names;
if isempty(dicNames)
    warning('No DIC datasets found to import into newEBSD.prop.');
else
    % Determine target EBSD property grid size
    propFields = fieldnames(newEBSD.prop);
    if isempty(propFields)
        error('newEBSD.prop has no properties to infer grid size.');
    end
    refSize = size(newEBSD.prop.(propFields{1}));

    for i = 1:numel(dicNames)
        origName = dicNames{i};
        data = rawData.forwardInterpolatedDIC.(origName);

        if ~ismatrix(data) || ~isnumeric(data)
            warning('Skipping %s: not a numeric 2D array.', origName);
            continue
        end

        % Match shape
        if isequal(size(data), refSize)
            fprintf('Dataset %s is inserted into the EBSD variable\n', origName);
            newEBSD.prop.(origName) = data;

        elseif isequal(size(data), [refSize(2) refSize(1)])
            fprintf('Transposed dataset %s is inserted into the EBSD variable\n', origName);
            newEBSD.prop.(origName) = data.';

        elseif numel(data) == prod(refSize)
            fprintf('Reshaped dataset %s is inserted into the EBSD variable\n', origName);
            newEBSD.prop.(origName) = reshape(data, refSize);

        else
            warning('Skipping %s: incompatible size [%d %d], expected [%d %d].', ...
                origName, size(data,1), size(data,2), refSize(1), refSize(2));
        end
    end
end

%% Check if EBSD and DIC coordinate systems are increase in same direction

% Extract 2×2 corner blocks
DIC_x = rawData.DIC.Coordinates_x(1:2,1:2);
DIC_y = rawData.DIC.Coordinates_y(1:2,1:2);

ebsd_x = newEBSD(1:2,1:2).x;
ebsd_y = newEBSD(1:2,1:2).y;

% Determine dominant x-direction increase for DIC
DIC_x_inc = (DIC_x(2,1) > DIC_x(1,1)) || (DIC_x(1,2) > DIC_x(1,1));
% Determine dominant x-direction increase for EBSD
EBSD_x_inc = (ebsd_x(2,1) > ebsd_x(1,1)) || (ebsd_x(1,2) > ebsd_x(1,1));

% Determine dominant y-direction increase for DIC
DIC_y_inc = (DIC_y(2,1) > DIC_y(1,1)) || (DIC_y(1,2) > DIC_y(1,1));
% Determine dominant y-direction increase for EBSD
EBSD_y_inc = (ebsd_y(2,1) > ebsd_y(1,1)) || (ebsd_y(1,2) > ebsd_y(1,1));

% Check consistency
if DIC_x_inc ~= EBSD_x_inc
    warning('DIC and EBSD coordinate systems disagree in the x-direction.');
end

if DIC_y_inc ~= EBSD_y_inc
    warning('DIC and EBSD coordinate systems disagree in the y-direction.');
end

%% ----------------------
% compute edge mask and plot using plot(newEBSD,...)
% if ~isfloat(newEBSD.prop.SE_raw)
%     img = double(newEBSD.prop.SE_raw);
% end
% img = newEBSD.prop.SE_raw;
% img = mat2gray(img); % normalize to [0,1]
% % use MATLAB's edge which accepts double images
% edgeMask = edge(img, 'Canny');
% 
% figure();
% plot(newEBSD,newEBSD.prop.Strains_MaximumShearStrain_RC)
% hold on
% plot(newEBSD, double(edgeMask), 'FaceAlpha',0.5,'micronbar','off')
% title('SE edges (Canny)')
figure;
% Step 1: Convert SE_raw to grayscale RGB
SE_gray = mat2gray(newEBSD.prop.SE_raw); % Normalize to [0,1]
SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image

% Step 2: Create a new MTEX figure and plot SE_raw as background
plot(newEBSD, SE_rgb); 
axis image off
hold on

maxStrainTest = newEBSD.prop.Displacements_y_displacement;
%maxStrainTest(newEBSD.prop.Others_CorrelationValue>0.9) = 0;

% Step 3: Plot Fyx with jet colormap and transparency
plot(newEBSD, maxStrainTest, 'FaceAlpha', 0.7);
title('Forward-deformed data')
colormap(viridis);
%caxis([0 0.2]);
mtexColorbar('title', 'Displacement y');
% exportScaledFigure(gcf,[input.fname,options.dateString,'AlignedDIC_SE_MaximumShearStrain.png'],'extraheight','nodatestring')

figure;
% Step 1: Convert SE_raw to grayscale RGB
SE_gray = mat2gray(newEBSD.prop.SE_raw); % Normalize to [0,1]
SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image

% Step 2: Create a new MTEX figure and plot SE_raw as background
plot(newEBSD, SE_rgb); 
axis image off
hold on

maxStrainTest = newEBSD.prop.Others_CorrelationValue;
%maxStrainTest(newEBSD.prop.Others_CorrelationValue>0.9) = 0;

% Step 3: Plot Fyx with jet colormap and transparency
plot(newEBSD, maxStrainTest, 'FaceAlpha', 0.7);
title('Forward-deformed data')
colormap(viridis);
%caxis([0 0.2]);
mtexColorbar('title', 'Displacement y');
% exportScaledFigure(gcf,[input.fname,options.dateString,'AlignedDIC_SE_MaximumShearStrain.png'],'extraheight','nodatestring')


figure;
% Step 1: Convert SE_raw to grayscale RGB
SE_gray = mat2gray(newEBSD.prop.SE_raw); % Normalize to [0,1]
SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image

% Step 2: Create a new MTEX figure and plot SE_raw as background
plot(newEBSD, SE_rgb); 
axis image off
hold on

maxStrainTest = rawData.DIC.Displacements_y_displacement;
%maxStrainTest(newEBSD.prop.Others_CorrelationValue>0.9) = 0;

% Step 3: Plot Fyx with jet colormap and transparency
plot(newEBSD, maxStrainTest, 'FaceAlpha', 0.7);
colormap(viridis);
title('Non-deformed data')
%caxis([0 0.2]);
mtexColorbar('title', 'Displacement y');

% 
% figure();
% h = newMtexFigure('nrows',1,'ncols',1,'visible','on','Position',get(0,'ScreenSize'));
% plot(newEBSD,newEBSD.Strains_MaximumShearStrain_RC,'micronbar','off')
% 
% caxis([0,0.2])
% hold on
% plot(newGrains.boundary,'lineColor','red','linewidth',1.5)
% title('Aligned EBSD boundaries over maximum shear strain')
% mtexColorbar
% exportScaledFigure(gcf,[input.fname,options.dateString,'AlignedDIC_EBSDboundaries_MaximumShearStrain.png'],'extraheight','nodatestring')

%% ----------------------
% Save newEBSD with DIC for later use
save([options.in_situIncrementName,'_aligned_newEBSD_withDIC'], 'newEBSD', 'originalEBSD')
