%% Nb_PTP_line_profile.m
% Script to interactively draw a line across an EBSD map and plot the
% misorientation profile and orientation evolution along that line.

clear; close all; clc;

error('The line extraction is currently not working. Try plotting plot(ebsd_line, ebsd_line.orientations) to see if the datapoints are correctly extracted.');

%% 1. Initialization & MTEX Preferences
plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Increment Selection
if exist('master_pwd', 'var')
    baseDir = fullfile(master_pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
else
    baseDir = fullfile(pwd, 'onedriveData', 'Philipp', 'Nb_PTP');
end
outDir = fullfile(baseDir, 'Line_Profiles');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% --- SELECT INCREMENT TO ANALYZE ---
% This should correspond to a deformed TKD scan from your experiments.
ang_file_to_analyze = 'insitu Specimen 1 3 Map Data 4_BW263.ang';
% -------------------------------------

%% 3. Load and Prepare Data
fprintf('Loading data from: %s\n', ang_file_to_analyze);
ang_path = fullfile(baseDir, ang_file_to_analyze);
if ~exist(ang_path, 'file')
    error('Specified .ang file not found: %s', ang_path);
end

% Standard parameters from your other scripts
tkd_tilt = -10 * degree;
thresh_ci = 0.13;
grainThresh = 5 * degree;
minSize = 10;
CS = crystalSymmetry.load('Nb-Niobium.cif');

% Load EBSD Data and perform identical rotations as in your analysis scripts
ebsd = EBSD.load(ang_path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
ebsd = ebsd.gridify;

% Reconstruct grains for context
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'), 'minPixel', minSize, 'angle', grainThresh);

%% 4. Interactive Line Selection
fprintf('Please select a line segment on the map by clicking two points.\n');

% Plot the map for selection
fig_map = figure;
ipf = ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection = xvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'), colors, 'micronbar', 'off');
hold on;
plot(grains.boundary, 'lineWidth', 1.5, 'lineColor', 'k');
title('Click two points to define a line profile');

% Get line from user input and draw it on the map
lineSec = ginput(2);
line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'color', 'w', 'linestyle', '--');
line(lineSec(:,1), lineSec(:,2), 'linewidth', 1, 'color', 'k', 'linestyle', '-');
title('Profile line selected');
drawnow;

% Extract the profile data using the MTEX spatialProfile function
ebsd_line = spatialProfile(ebsd, lineSec);

% Save the selection map for reference
[~, baseName, ~] = fileparts(ang_file_to_analyze);
profile_name = sprintf('%s_profile_%s', baseName, datestr(now,'yyyymmdd_HHMMSS'));
exportScaledFigure(fig_map, fullfile(outDir, [profile_name, '_map.png']), 'noDateString');

%% 5. Plot Misorientation Profile
fig_profile = figure;

% Calculate distance along the profile from the start point
dist = sqrt((ebsd_line.x - ebsd_line.x(1)).^2 + (ebsd_line.y - ebsd_line.y(1)).^2);

% Plot 1: Misorientation angle relative to the first point on the line (cumulative)
plot(dist, angle(ebsd_line(1).orientations, ebsd_line.orientations) / degree, 'LineWidth', 2, 'DisplayName', 'Cumulative (from start)');
hold on;

% Plot 2: Point-to-point misorientation gradient (local)
dist_grad = 0.5 * (dist(1:end-1) + dist(2:end));
grad_angle = angle(ebsd_line(1:end-1).orientations, ebsd_line(2:end).orientations) / degree;
plot(dist_grad, grad_angle, 'LineWidth', 2, 'DisplayName', 'Local Gradient (point-to-point)');

% Formatting
grid on;
xlabel('Distance along profile (\mum)');
ylabel('Misorientation Angle (degree)');
title({'Misorientation along selected profile', baseName}, 'Interpreter', 'none');
legend('show', 'Location', 'best');

% Save the profile plot
exportScaledFigure(fig_profile, fullfile(outDir, [profile_name, '_plot.png']), 'noDateString');

%% 6. Plot Orientations in IPF
fig_ipf = figure;

% Plot orientations along the line in an inverse pole figure
% Colorize the points by their distance along the line
plotIPDF(ebsd_line.orientations, [xvector, yvector, zvector], ...
  'property', dist, 'markersize', 15, 'antipodal');

mtexColorbar('title', 'Distance along profile (\mum)');
mtexTitle('Orientation Evolution along Profile');

% Save the IPF plot
exportScaledFigure(fig_ipf, fullfile(outDir, [profile_name, '_IPF.png']), 'noDateString');

fprintf('\nFinished. Profile plots saved in:\n%s\n', outDir);