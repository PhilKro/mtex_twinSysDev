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
defaults.DIC_filename = 'Processing_subset=25_step=8';
defaults.outDir         = 'Processed_Output';           % folder under input.dir to save output files


input = struct( ...
    'fname',             'insitu Specimen 1 0 Map Data 1', ...    % base EBSD filename (no extension)                   
    'SEImagefname',      '0_10kv_DIC.tif', ...          % SE image filename
    'DICfoldername' ,    fullfile('onedriveData','Philipp','Nb_PTP','DIC_h5', defaults.DIC_filename),... % DIC HDF5 filename
    'sel_points_fname', 'insitu Specimen 1 0 Map Data 1_BW263_20251126_115536_Selection_Points_Backup_EBSD.mat'); % mat file with preselected points for alignment                     

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
input.outPath         = fullfile(input.dir, input.outDir);
if ~exist(input.outPath, 'dir')
    mkdir(input.outPath);
end

input.SEImageFullPath = fullfile(input.dir, input.SEImagefolder, input.SEImagefname);
% input.DIC_path  = fullfile(input.dir, input.DICfoldername, input.DICfname);
input.ang_path        = fullfile(input.dir, [input.fname,'_BW',num2str(defaults.bandWidth),'.ang']);
input.ctf_path = fullfile(input.dir, [input.fname,'.ctf']);

options.dateString = datestr(now,'yyyymmdd_HHMMSS');
options.in_situIncrementName = [input.fname,'_BW',num2str(defaults.bandWidth),'_',options.dateString];

% (Optional) show summary
fprintf('Using IO config:\n  dir: %s\n  fname: %s\n  SE image: %s\n  crop: [%d %d %d %d]\n  DIC: %s (path: %s)\n', ...
    input.dir, input.fname, input.SEImageFullPath, input.startx, input.starty, input.numx, input.numy, input.DICfoldername);

% save ioConfig to a JSON file with datestring in the filename
outname = fullfile(input.outPath, sprintf([options.in_situIncrementName, '_undeformed_ioConfig.json']));

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

exportScaledFigure(gcf,fullfile(input.outPath, [input.fname,options.dateString,'ebsd_.png']),'nodatestring', 'extrawidth')

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
    save(fullfile(input.outPath, [options.in_situIncrementName, '_', options.case,'_alignment_data_EBSD']),'alignment_data_EBSD')
end

%% ----------------------
% Create a new grid matched to the SE reference and warp EBSD onto it
psize_new = rawData.SE_raw.psize;  % use SE pixel size (edit as needed)

[Xnew,Ynew] = meshgrid(0:psize_new:max(rawData.SE_raw.X(:)),0:psize_new:max(rawData.SE_raw.Y(:))); 

newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid,Xnew,Ynew,alignment_data_EBSD);

% save the id state so the ebsd objects can be refenced with each other
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
exportScaledFigure(gcf,fullfile(input.outPath, [input.fname,options.dateString,'AlignedEBSD_SE_.png']),'extraheight','nodatestring')
% %% ---------------------
% % Save aligned newEBSD for later use
% save([options.in_situIncrementName,'_aligned_newEBSD'], 'newEBSD', 'originalEBSD')

%% --- DIC DATA LOADING AND PROCESSING ---

% 1. Setup file iteration
dicDir = input.DICfoldername; % Assume input.DIC_filename is the directory path
filePattern = fullfile(dicDir, 'B*.hdf5');
dicFiles = dir(filePattern);

% Sort files numerically by name to ensure correct time sequence
fileNames = {dicFiles.name};
numStrings = regexp(fileNames, '(\d+)\.hdf5$', 'tokens', 'once');
numbers = cellfun(@(x) str2double(x{1}), numStrings, 'UniformOutput', true);
[~, sortedIdx] = sort(numbers);
dicFiles = dicFiles(sortedIdx);

if isempty(dicFiles)
    error('No HDF5 files matching pattern 0*.h5 found in the directory: %s', dicDir);
end

% Target size defined by SE image (assuming rawData.SE_raw is available)
targetM = rawData.SE_raw.origM;
targetN = rawData.SE_raw.origN;
targetSize = [targetM, targetN];

% Required HDF5 dataset paths to extract for each frame
requiredDatasetPaths = {
    'Coordinates/x',
    'Coordinates/y',
    'Displacements/x-displacement',
    'Displacements/y-displacement'
    };

% Prepare container
rawData.DIC = struct();
rawData.DIC.names = {};

% Region of interest in SE/DIC images (1-based MATLAB indexing)
x0 = input.startx;
y0 = input.starty;
x1 = input.startx + input.numx - 1;
y1 = input.starty + input.numy - 1;

% === Support functions to resolve dataset paths ===

% function allPaths = collectH5Paths(infoStruct, currentPath)
% if nargin < 2, currentPath = ''; end
% allPaths = {};
% for i = 1:numel(infoStruct.Datasets)
% p = fullfile(currentPath, infoStruct.Datasets(i).Name);
% allPaths{end+1} = ['/' strrep(p,'\','/')]; %#ok<AGROW>
% end
% for i = 1:numel(infoStruct.Groups)
% grp = infoStruct.Groups(i);
% sub = fullfile(currentPath, grp.Name);
% allPaths = [allPaths, collectH5Paths(grp, strrep(sub,'\','/'))]; %#ok<AGROW>
% end
% end

% function fullPath = findDatasetPath(allPaths, requiredRelPath)
% req = strsplit(requiredRelPath, '/');
% reqTail = req(max(end-1,1):end);
% fullPath = '';
% for i = 1:numel(allPaths)
% p = allPaths{i};
% parts = strsplit(p, '/');
% if length(parts) >= length(reqTail)
% if isequal(parts(end-numel(reqTail)+1:end), reqTail)
% fullPath = p;
% return;
% end
% end
% end
% end


% === INITIAL DETECTION: Use first file to resolve paths ===
firstFilePath = fullfile(dicDir, dicFiles(1).name);
info = h5info(firstFilePath);
allPaths = collectH5Paths(info);


resolvedPaths = struct();
for i = 1:numel(requiredDatasetPaths)
    rel = requiredDatasetPaths{i};
    full = findDatasetPath(allPaths, rel);
    if isempty(full)
        error('Dataset %s not found in %s', rel, dicFiles(1).name);
    end
    resolvedPaths.(matlab.lang.makeValidName(rel)) = full;
end


% Load coordinates
xcoord = h5read(firstFilePath, resolvedPaths.Coordinates_x);
ycoord = h5read(firstFilePath, resolvedPaths.Coordinates_y);


dataSize = size(xcoord);
if xor(dataSize(1) > dataSize(2), targetM > targetN)
    xcoord = xcoord';
    ycoord = ycoord';
    dataSize = size(xcoord);
    isTransposed = true;
else
    isTransposed = false;
end


stepsize = abs(xcoord(2,1)-xcoord(1,1));
if stepsize == 0, stepsize = abs(xcoord(1,2)-xcoord(1,1)); end
if stepsize == 0, error('Stepsize could not be determined.'); end


% === MAIN LOOP ===
for fileIdx = 1:numel(dicFiles)
    currentFile = dicFiles(fileIdx).name;
    currentFilePath = fullfile(dicDir, currentFile);


    frameNum = numbers(sortedIdx(fileIdx));
    frameStr = sprintf('%05d', frameNum);


    for k = 1:numel(requiredDatasetPaths)
        rel = requiredDatasetPaths{k};
        % Build full dataset path dynamically for this file
        info = h5info(currentFilePath);
        allPaths = collectH5Paths(info);      % same function as before
        fullPath = findDatasetPath(allPaths, rel);

        if isempty(fullPath)
            error('Dataset %s not found in %s', rel, currentFile);
        end

        data = h5read(currentFilePath, fullPath);

        if isTransposed
            data = data';
        end


        sz = size(data);
        if isequal(sz, targetSize)
            subdata = extractRegion(data, x0, y0, x1, y1, rel);
        else
            if contains(rel, 'Coordinates')
                if data(1,1) < data(2,1)
                    subdata = repmat(1:targetN, targetM, 1);
                elseif data(1,1) < data(1,2)
                    subdata = repmat((1:targetM)', 1, targetN);
                else
                    subdata = kron(data, ones(stepsize));
                end
                subdata = extractRegion(subdata, x0, y0, x1, y1, rel);
            else
                subdata = kron(data, ones(stepsize));
                subdata = extractRegion(subdata, x0, y0, x1, y1, rel);
            end
        end

        % Generate MATLAB-friendly field name
        fname = matlab.lang.makeValidName(strrep(currentFile, '/', '_'));

        % Extract clean frame name (e.g., B0001)
        tok = regexp(fname, '(B\d+)', 'tokens', 'once');
        if ~isempty(tok)
            baseName = tok{1};
        else
            baseName = fname;
        end

        % Remove trailing _hdf5 if present
        baseName = regexprep(baseName, '_hdf5$', '');

        % Extract dataset suffix (Coordinates_x, Displacements_y_displacement, etc.)
        relParts = strsplit(rel, '/');
        dsSuffix = strjoin(relParts, '_');

        % Build raw name and then sanitize it
        rawName = [baseName '_' dsSuffix];
        shortName = matlab.lang.makeValidName(rawName);

        rawData.DIC.(shortName) = subdata;
        rawData.DIC.names{end+1} = shortName;
    end
end

% -----------------------------------------------------------------
% INSERT DIC FIELDS INTO newEBSD.prop
%
% This matches the DIC grid to the EBSD grid:
%   - direct match
%   - transpose match
%   - reshaping if the total number of elements fits
% -----------------------------------------------------------------

dicNames = rawData.DIC.names;
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
        data = rawData.DIC.(origName);

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
DIC_x = rawData.DIC.B0001_Coordinates_x(1:2,1:2);
DIC_y = rawData.DIC.B0001_Coordinates_y(1:2,1:2);

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

figure;
% Step 1: Convert SE_raw to grayscale RGB
SE_gray = mat2gray(newEBSD.prop.SE_raw); % Normalize to [0,1]
SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image

% Step 2: Create a new MTEX figure and plot SE_raw as background
plot(newEBSD, SE_rgb); 
axis image off
hold on

maxStrainTest = newEBSD.prop.B0005_Displacements_y_displacement;
%maxStrainTest(newEBSD.prop.Others_CorrelationValue>0.9) = 0;

% Step 3: Plot Fyx with jet colormap and transparency
plot(newEBSD, maxStrainTest(:), 'FaceAlpha', 0.7);
title('First step data')
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


maxStrainTest = newEBSD.prop.(dicNames{end});
%maxStrainTest(newEBSD.prop.Others_CorrelationValue>0.9) = 0;

% Step 3: Plot Fyx with jet colormap and transparency
plot(newEBSD, maxStrainTest(:), 'FaceAlpha', 0.7);
title('Last step data')
colormap(viridis);
%caxis([0 0.2]);
mtexColorbar('title', 'Displacement y');
% exportScaledFigure(gcf,[input.fname,options.dateString,'AlignedDIC_SE_MaximumShearStrain.png'],'extraheight','nodatestring')
%% ----------------------
% Save newEBSD with DIC for later use
finalOutName = fullfile(input.outPath, [options.in_situIncrementName, input.DIC_filename, '_undeformedEBSD_withAllDIC.mat']);
save(finalOutName, 'newEBSD', 'originalEBSD')
fprintf('Saved as: %s\n', finalOutName)
