%% Nb_PTP_3D_GB_Slip_Viewer.m
% Interactive 3D visualization of the Grain Boundary and Slip Planes
%
% This script loads undeformed EBSD data, identifies the left/right grains
% forming the bicrystal, calculates their active slip systems, and plots a 
% 3D representation of the slip planes, slip directions, and the grain 
% boundary. The grain boundary can be rotated +/- 10 degrees around the 
% Y-axis via a slider, which dynamically updates the slip plane traces on it.

clear; close all; clc;

%% 1. Initialization & Mtex Preferences
plotx2east;
plotzIntoPlane;
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% 2. Directories and Data Loading
% Determine paths (assuming this script is in Pretty_Plotter directory)
% scriptDir = fileparts(mfilename('fullpath'));
% baseDir = fullfile(scriptDir, '..');
outDir  = fullfile('C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Nb_PTP\', 'Processed_Output');

undeformed_dic_mat = fullfile(outDir, 'insitu Specimen 1 0 Map Data 1_BW263_20260320_101109Processing_subset=25_step=8_undeformedEBSD_withAllDIC.mat');

if ~exist(undeformed_dic_mat, 'file')
    error('Cannot find undeformed EBSD data at %s. Please check the path.', undeformed_dic_mat);
end

fprintf('Loading EBSD data...\n');
load(undeformed_dic_mat, 'newEBSD');
ebsd = newEBSD;

%% 3. Grain Reconstruction
fprintf('Reconstructing grains...\n');
thresh_ci = 0.13;
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
[grains, ebsd.grainId] = calcGrains(ebsd, 'minPixel', 1000, 'angle', 5*degree);

% Find the two largest grains (the bicrystal)
[~, sortedIdx] = sort(grains.numPixel, 'descend');
bicrystalGrains = grains(sortedIdx(1:min(2, length(grains))));

% Sort left to right based on centroid X coordinate
[~, left_to_right] = sort(bicrystalGrains.centroid.x);
bicrystalGrains = bicrystalGrains(left_to_right);
leftGrain = bicrystalGrains(1);
rightGrain = bicrystalGrains(2);

fprintf('Identified Left Grain ID: %d, Right Grain ID: %d\n', leftGrain.id, rightGrain.id);

%% 4. Define Active Slip Systems
% Selected indices for the active slip systems. Change these based on your findings!
left_ss_idx  = 18; 
right_ss_idx = 27; 

% Set up BCC slip systems as done in the full analysis
sS_BCC_plain = slipSystem.bcc(ebsd.CS);
sS_BCC_select = sS_BCC_plain(1:3); % Restricted BCC systems
sS_BCC = symmetrise(sS_BCC_select, 'antipodal');

% Transform into specimen reference frame using the grain's mean orientation
sS_left = leftGrain.meanOrientation * sS_BCC(left_ss_idx);
sS_right = rightGrain.meanOrientation * sS_BCC(right_ss_idx);

fprintf('\nLeft Slip System (Idx %d):\n', left_ss_idx);
disp(sS_left);
fprintf('Right Slip System (Idx %d):\n', right_ss_idx);
disp(sS_right);

%% 5. Setup the 3D Viewer GUI
fig = figure('Name', '3D GB & Slip Plane Viewer', 'Position', [100, 100, 950, 700]);
ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.8, 0.7]);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel(ax, 'X (East)'); ylabel(ax, 'Y (South)'); zlabel(ax, 'Z (Into Plane)');

% Set initial view corresponding to MTEX defaults (viewing from -Z towards +Z)
view(ax, 0, -90);
rotate3d(ax, 'on');

L = 2; % Scale factor for plotting plane geometries

% Plot Left Slip Plane (Red)
[V_L, F_L] = createPlanePatch(sS_left.n, L);
patch('Faces', F_L, 'Vertices', V_L, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'Parent', ax, 'DisplayName', 'Left Slip Plane');

% Plot Right Slip Plane (Green)
[V_R, F_R] = createPlanePatch(sS_right.n, L);
patch('Faces', F_R, 'Vertices', V_R, 'FaceColor', 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'Parent', ax, 'DisplayName', 'Right Slip Plane');

% Plot Slip Directions (Burgers Vectors)
b_L = sS_left.b.normalize;
b_R = sS_right.b.normalize;
quiver3(ax, 0, 0, 0, b_L.x, b_L.y, b_L.z, L*0.9, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'DisplayName', 'Left Slip Dir');
quiver3(ax, 0, 0, 0, b_R.x, b_R.y, b_R.z, L*0.9, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'DisplayName', 'Right Slip Dir');

% Initialize Grain Boundary Plane (Blue)
gb_patch = patch('Faces', [1 2 3 4], 'Vertices', zeros(4,3), 'FaceColor', 'b', 'FaceAlpha', 0.35, 'EdgeColor', 'b', 'Parent', ax, 'DisplayName', 'Grain Boundary');

% Initialize Trace Lines
trace_L_line = plot3(ax, [0 0], [0 0], [0 0], 'r-', 'LineWidth', 4, 'DisplayName', 'Left Trace on GB');
trace_R_line = plot3(ax, [0 0], [0 0], [0 0], 'g-', 'LineWidth', 4, 'DisplayName', 'Right Trace on GB');

legend(ax, 'Location', 'bestoutside');

%% 6. GUI Controls
% Slider for GB Rotation (+- 10 degrees)
slider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Min', -20, 'Max', 20, 'Value', 0, ...
    'Position', [250, 40, 400, 20]);

% Text label for the Slider
txt = uicontrol('Parent', fig, 'Style', 'text', ...
    'Position', [250, 65, 400, 20], ...
    'String', 'GB Rotation (around Y): 0.0 deg', ...
    'FontSize', 11, 'FontWeight', 'bold');

% Text label for the calculated angles
txt_angles = uicontrol('Parent', fig, 'Style', 'text', ...
    'Position', [150, 15, 600, 20], ...
    'String', 'Angles - Traces: -- | Burgers: -- | Normals: --', ...
    'FontSize', 11, 'FontWeight', 'bold');

% Button to return to initial MTEX view
btn = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Reset View', ...
    'Position', [680, 40, 100, 30], ...
    'Callback', @(src, event) view(ax, 0, -90));

% Attach Callback to Slider
slider.Callback = @(src, event) updateGB(src.Value, txt, txt_angles, gb_patch, trace_L_line, trace_R_line, sS_left, sS_right, L);

% Run the update once to initialize the GB plane and traces at 0 degrees
updateGB(0, txt, txt_angles, gb_patch, trace_L_line, trace_R_line, sS_left, sS_right, L);


%% Helper Functions

function [V, F] = createPlanePatch(n_vec, size)
    % Create a square patch perpendicular to the normal vector
    n = [n_vec.x; n_vec.y; n_vec.z];
    n = n / norm(n);
    
    % Find an orthogonal vector
    if abs(n(1)) < 0.9
        u = [1; 0; 0];
    else
        u = [0; 1; 0];
    end
    u = u - dot(u, n)*n;
    u = u / norm(u);
    v = cross(n, u); % ensures v is also unit length
    
    V = [u + v, -u + v, -u - v, u - v]' * size;
    F = [1 2 3 4];
end

function updateGB(angle_deg, txt_obj, txt_angles_obj, gb_patch, trace_L_line, trace_R_line, sS_left, sS_right, L)
    % Update text label
    txt_obj.String = sprintf('GB Rotation (around Y): %0.1f deg', angle_deg);
    
    % Calculate new GB normal by rotating xvector around yvector
    rot = rotation.byAxisAngle(yvector, angle_deg * degree);
    n_gb_mtex = rot * xvector;
    n_gb = [n_gb_mtex.x; n_gb_mtex.y; n_gb_mtex.z];
    n_gb = n_gb / norm(n_gb);
    
    % Update Grain Boundary Patch Vertices
    u_gb = [0; 1; 0]; % Y vector is always inside the grain boundary
    v_gb = cross(n_gb, u_gb);
    V = [u_gb + v_gb, -u_gb + v_gb, -u_gb - v_gb, u_gb - v_gb]' * L;
    gb_patch.Vertices = V;
    
    % Subfunction to update trace line efficiently
    t_L = updateTrace(trace_L_line, n_gb, sS_left.n, L);
    t_R = updateTrace(trace_R_line, n_gb, sS_right.n, L);
    
    % Calculate acute angles
    if norm(t_L) > 0 && norm(t_R) > 0
        trace_angle = acosd(abs(dot(t_L, t_R)));
    else
        trace_angle = NaN;
    end
    
    b1 = [sS_left.b.x; sS_left.b.y; sS_left.b.z]; b1 = b1 / norm(b1);
    b2 = [sS_right.b.x; sS_right.b.y; sS_right.b.z]; b2 = b2 / norm(b2);
    b_angle = acosd(abs(dot(b1, b2)));
    
    n1 = [sS_left.n.x; sS_left.n.y; sS_left.n.z]; n1 = n1 / norm(n1);
    n2 = [sS_right.n.x; sS_right.n.y; sS_right.n.z]; n2 = n2 / norm(n2);
    n_angle = acosd(abs(dot(n1, n2)));
    
    txt_angles_obj.String = sprintf('Angles - Traces: %0.1f deg | Burgers: %0.1f deg | Normals: %0.1f deg', trace_angle, b_angle, n_angle);
end

function t = updateTrace(line_obj, n_gb, n_sS_mtex, L)
    n_sS = [n_sS_mtex.x; n_sS_mtex.y; n_sS_mtex.z];
    t = cross(n_gb, n_sS); % Trace line is the cross product of the normals
    if norm(t) > 1e-5
        t = t / norm(t);
        t_plot = t * L * 1.05; % Extend slightly out of the boundaries
        set(line_obj, 'XData', [-t_plot(1) t_plot(1)], 'YData', [-t_plot(2) t_plot(2)], 'ZData', [-t_plot(3) t_plot(3)]);
    else
        t = [0; 0; 0];
    end
end