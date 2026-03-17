clear
close all
load('CS_cS_andPlanes.mat')

%% Set MTEX plotting conventions
plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% Define file directories
ang_dir = 'data\Ti_PTP\Ti_cAxis_PTP_1\';
h5_dir  = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0007_Ti_cAxis_PTP_1\h5oina';  % change if needed

% ang_fname = 'Insitu Specimen 1 Site 8 Map Data 7_BW123.ang';
% ang_fname = 'Insitu Specimen 1 Site 7 Map Data 6_BW158.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 6 Map Data 5_BW123.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 5 Map Data 4_BW123.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 4 Map Data 3_BW123.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 3 Map Data 2_BW123.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 2 Map Data 1_BW123.ang'; 
% ang_fname = 'Insitu Specimen 1 Site 2 Map Data 1_BW263.ang';
ang_file  = fullfile(ang_dir, ang_fname);

% Try to find matching .h5 file automatically
h5_fname = regexprep(ang_fname, '_BW\d+\.ang$', '.h5oina');
h5_file = fullfile(h5_dir, h5_fname);

if ~isfile(h5_file)
    error('H5 file not found for %s. Stopping execution.', ang_fname);
else
    fprintf('Found matching H5 file:\n%s\n', h5_file);
end

%% Read TKD tilt from H5 (if exists)

info = h5info(h5_file);
try
    tkd_tilt = h5read(h5_file, '/1/EBSD/Header/Tilt Angle');
    fprintf('Read in tkd tilt from H5 file: %s\n', tkd_tilt);
catch
    error('Could not read tilt angle from H5.');
end
%% Read iu EBSD data from .ang file
thresh_ci = 0.12;

% Define crystal symmetry
% Titanium
CS = crystalSymmetry('6/mmm', [2.95 2.95 4.68], 'x||a', 'mineral', 'Ti');

% Load EBSD data
ebsd = EBSD.load(ang_file, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');

% Apply TKD tilt correction and rotation
ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector, 180*degree);
ebsd = rotate(ebsd, rot, 'keepXY');

%% Read existing Euler data from H5
try
    originalEuler = h5read(h5_file, '/1/EBSD/Data/Euler');
catch
    error('Could not read /1/EBSD/Data/Euler from H5 file. Stopping.');
end

% Backup original Euler as /1/EBSD/Data/HoughEuler
data_info = h5info(h5_file, '/1/EBSD/Data');
dataset_names = {data_info.Datasets.Name};

if any(strcmp(dataset_names, 'HoughEuler'))
    % Prompt user before overwriting
    choice = questdlg('Dataset /1/EBSD/Data/HoughEuler already exists. Overwrite?', ...
                      'Overwrite Confirmation', ...
                      'Yes', 'No', 'No');
    if strcmp(choice, 'Yes')
        fprintf('Overwriting existing /1/EBSD/Data/HoughEuler...\n');
        h5write(h5_file, '/1/EBSD/Data/HoughEuler', originalEuler);
        fprintf('Backup updated.\n');
    else
        fprintf('Keeping existing /1/EBSD/Data/HoughEuler. No overwrite performed.\n');
    end
else
    % Create and write new dataset
    h5create(h5_file, '/1/EBSD/Data/HoughEuler', size(originalEuler), 'Datatype', 'single');
    h5write(h5_file, '/1/EBSD/Data/HoughEuler', originalEuler);
    fprintf('Created and wrote new backup /1/EBSD/Data/HoughEuler.\n');
end
%% Find the symetrically equivalent orientations with the lowest overall sum of euler angles
Symmiorientations = [];
for i = 1:length(ebsd)
    cur = ebsd(i);
    % check if the phase is empty or the orientation is NaN, then assign a
    % NaN orientation BUT OF THE correct CS otherwise (if the first entry is NaN)
    % everything is triclinic in the assignment ebsd.orientations = Symmiorientations
    if isempty(cur.indexedPhasesId) | isnan(cur.orientations.phi1)
        Symmiorientations = [Symmiorientations;orientation.byEuler(NaN,NaN,NaN,ebsd.CS)];
        continue
    end
    symmi = symmetrise(cur.orientations,'proper');
    %chooses the symmetrized orientation where phi1 Phi and phi2 are below
    %pi which corresponds to the format chosen in the h5 files (probably
    %redundant with the minimization afterwards (checks after the for loop)
    symmi = symmi(symmi.phi1<pi & symmi.Phi<pi & symmi.phi2<pi);
    if isempty(symmi)
        warning('symmi is empty')
        symmi
    end
    %chooses the symmetrized orientation with the lowest overall sum of
    %angles
    [~, min_idx] = min(sum([symmi.phi1,symmi.Phi,symmi.phi2], 2));
    % if i == 6213 || cur.oldId == 6213
    %     symmi
    %     symmi(min_idx)
    % end
    Symmiorientations = [Symmiorientations;symmi(min_idx)];
end
%check if all Eulerangles are below pi (or NaN)
assert(all(Symmiorientations.phi1(:) < pi | isnan(Symmiorientations.phi1(:)))&all(Symmiorientations.Phi(:) < pi | isnan(Symmiorientations.Phi(:)))&all(Symmiorientations.phi2(:) < pi | isnan(Symmiorientations.phi2(:))), 'Error: Not all values are below pi.');

%% Write new Euler data if dimensions match
newEuler = Symmiorientations.Euler;
newEuler = newEuler';

if isequal(size(originalEuler), size(newEuler))
    fprintf('Writing new Euler data to %s...\n', h5_file);
    h5write(h5_file, '/1/EBSD/Data/Euler', newEuler);
    fprintf('Successfully updated Euler dataset in H5.\n');
else
    warning('Size mismatch: H5 Euler data is %s, EBSD Euler data is %s. Skipping write.', ...
        mat2str(size(originalEuler)), mat2str(size(newEuler)));
end