%% Initialization
clear; close all; clc;

% Setup MTEX plotting conventions
plotx2east;
plotzOutOfPlane; % Standard projection (change to plotzIntoPlane if preferred)

%% 1. Load Data
% Ensure the file path is correct
filePath = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0010_Ti_Polycrystal\Project 1 01_Specimen_WedgeShape Site 1 Map Data 1.h5oina';

% Load Aztec .h5oina EBSD data. 
% 'convertEuler2SpatialReferenceFrame' usually handles the standard Aztec 70 degree setup perfectly.
ebsd = loadEBSD_h5oina(filePath, 'convertEuler2SpatialReferenceFrame');

% Isolate indexed data and extract the crystal symmetry (assuming Alpha-Ti)
ebsd_indexed = ebsd('indexed');
if iscell(ebsd_indexed.CS)
    CS = ebsd_indexed.CS{1}; % Pick dominant indexed phase if multiple exist
else
    CS = ebsd_indexed.CS;
end

%% 2. Grain Reconstruction
% Define threshold for grain boundary and minimum pixel size
grainThresh = 5 * degree;
minPixelSize = 5;

% Reconstruct grains
[grains, ebsd.grainId] = calcGrains(ebsd_indexed, 'angle', grainThresh, 'minPixel', minPixelSize);

% Smooth grain boundaries (crucial to get reliable/straight boundary trace directions)
grains = smooth(grains, 5);

%% 3. Plot Grains (IPF Z Map)
figure;
ipfKey = ipfColorKey(CS);
ipfKey.inversePoleFigureDirection = zvector;
grainColors = ipfKey.orientation2color(grains.meanOrientation);

% Plot grains colored by IPF-Z
plot(grains, grainColors, 'micronbar', 'on');
hold on;

% Overlay all grain boundaries in thin, semi-transparent black
plot(grains.boundary, 'linecolor', [0.3 0.3 0.3], 'linewidth', 1, 'DisplayName', 'All Grain Boundaries');

%% 4. Find Specific Grain Boundaries
% Extract internal boundaries between indexed grains
gB = grains.boundary('indexed', 'indexed');

% Get the mean orientations of the adjacent grains for each boundary segment
ori1 = grains(gB.grainId(:,1)).meanOrientation;
ori2 = grains(gB.grainId(:,2)).meanOrientation;

% Define the c-axis direction for Titanium (HCP)
c_axis = Miller(0, 0, 0, 1, CS);

% Calculate c-axis in specimen coordinates for both adjacent grains
c1 = ori1 * c_axis;
c2 = ori2 * c_axis;

% Set antipodal flag to true (since the c-axis direction represents a line without polarity)
c1.antipodal = true;
c2.antipodal = true;

% Extract boundary trace directions and set to antipodal
gB_dir = gB.direction;
gB_dir.antipodal = true;

% --- Define Tolerances ---
tol_xy = 15 * degree;    % Tolerance for c-axis being in the XY plane
tol_trace = 15 * degree; % Tolerance for c-axis being parallel to the boundary trace
tol_far = 30 * degree;   % Tolerance for the adjacent grain being "randomly" oriented (far from alignment)

% Condition A: c-axis is roughly in the XY plane 
% (Angle with zvector should be close to 90 degrees; because c1 is antipodal, max angle is 90)
cond1_xy = (90*degree - angle(c1, zvector)) < tol_xy;
cond2_xy = (90*degree - angle(c2, zvector)) < tol_xy;

% Condition B: c-axis is roughly parallel to the grain boundary trace
cond1_trace = angle(c1, gB_dir) < tol_trace;
cond2_trace = angle(c2, gB_dir) < tol_trace;

% Condition C: The neighboring grain should be far from this alignment
% (We define "far" as its c-axis misaligning with the boundary trace by at least `tol_far`)
far1 = angle(c1, gB_dir) > tol_far;
far2 = angle(c2, gB_dir) > tol_far;

% Combine conditions to identify aligned grains
is_g1_aligned = cond1_xy & cond1_trace;
is_g2_aligned = cond2_xy & cond2_trace;

% Final logical filter: ONE grain perfectly matches the condition, the OTHER is completely far off
target_condition = (is_g1_aligned & far2) | (is_g2_aligned & far1);

% Extract the specific boundaries
target_gB = gB(target_condition);

%% 5. Plot the Specific Boundaries
% Overlay the matching boundaries in thick red
plot(target_gB, 'linecolor', 'g', 'linewidth', 3, 'DisplayName', 'C-Axis Aligned Boundaries');
title('Grains and Boundaries with Specific C-Axis Alignment');

% Optional: Quiver plot to visually verify the matching c-axes
% (Draws a magenta line representing the aligned c-axis across the boundary)
% aligned_g1_gB = gB(is_g1_aligned & far2);
% aligned_g2_gB = gB(is_g2_aligned & far1);
% quiver(aligned_g1_gB, grains(aligned_g1_gB.grainId(:,1)).meanOrientation * c_axis, 'color', 'm', 'autoScale', 'off');
% quiver(aligned_g2_gB, grains(aligned_g2_gB.grainId(:,2)).meanOrientation * c_axis, 'color', 'm', 'autoScale', 'off');

legend('Location', 'bestout');
hold off;

fprintf('Found %d boundary segments matching the specific criteria.\n', length(target_gB));

%% 6. Plot Pattern Quality (Band Contrast / CI) Map
% Aztec .h5oina files typically store pattern quality as Band Contrast (bc)
figure;

if isfield(ebsd.prop, 'bc')
    pq_data = ebsd.bc;
    title_str = 'Band Contrast Map';
elseif isfield(ebsd.prop, 'ci')
    pq_data = ebsd.ci;
    title_str = 'Confidence Index (CI) Map';
else
    pq_data = ebsd.mad; % Mean Angular Deviation fallback
    title_str = 'Mean Angular Deviation (MAD) Map';
end

% Normalize to [0, 1] to create a proper grayscale RGB array
pq_norm = (pq_data - min(pq_data(:))) ./ (max(pq_data(:)) - min(pq_data(:)));
pq_rgb = [pq_norm(:), pq_norm(:), pq_norm(:)];

plot(ebsd, pq_rgb, 'micronbar', 'on');
hold on;
plot(target_gB, 'linecolor', 'g', 'linewidth', 3, 'DisplayName', 'C-Axis Aligned Boundaries');
title(title_str);
hold off;
