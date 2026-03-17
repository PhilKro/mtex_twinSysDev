%% MTEX Batch Processing: Grain and Twin Boundary Plots
clear close all

% --- Setup ---
dirPath = 'data\Ti_PTP\Ti_cAxis_PTP_2\'; % The folder containing your .ang files
filePattern = fullfile(dirPath, '*.ang');
files = dir(filePattern);

% Parameters
tkd_tilt = -11 * degree;
thresh_ci = 0.15;
grainThresh = 5 * degree;
minSize = 3;

% Crystal Symmetry (Ti)
CS = crystalSymmetry('6/mmm', [2.95 2.95 4.68], 'x||a', 'mineral', 'Ti');

% Define Twinning System ({10-12})
tS = twinningSystem.hexagonal_1012_1012(CS);

% Global plotting settings
plotx2east
plotzIntoPlane

% --- Loop through each file ---
for i = 1:length(files)
    baseFileName = files(i).name;
    fullFileName = fullfile(dirPath, baseFileName);
    fprintf('Processing: %s\n', baseFileName);
    
    try
        % 1. Load and Clean Data
        ebsd = EBSD.load(fullFileName, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
        
        % TKD Rotations
        ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
        ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
        ebsd = rotate(ebsd, rotation.byAxisAngle(zvector, 180*degree), 'keepXY');
        
        % Gridify for boundary calculation
        ebsd = ebsd.gridify;
        
        % 2. Grain Reconstruction
        [grains, ebsd.grainId] = calcGrains(ebsd, 'minPixel', minSize, 'angle', grainThresh);
        grains = smooth(grains, 2); % Optional: smoother boundaries
        
        % Extract boundaries
        gB = grains.boundary('Ti', 'Ti');
        % Identify twin boundaries within 5 degrees of the theoretical misorientation
        twinninggB = gB(angle(gB.misorientation, tS.parentTwinMisorientation) < 5*degree);
        
        % 3. Prepare Plotting Colors (IPF + CI)
        % Calculate CI grayscale
        ci_norm = (ebsd.ci - min(ebsd.ci, [], 'all')) ./ (max(ebsd.ci, [], 'all') - min(ebsd.ci, [], 'all'));
        cibw = repmat(ci_norm, 1, 3); 
        
        % IPF Key for Y-direction
        ipfKey = ipfColorKey(CS);
        ipfKey.inversePoleFigureDirection = yvector;
        colors = ipfKey.orientation2color(ebsd('indexed').orientations);
        
        % 4. Generate and Save Plot
        f1 = figure('Visible', 'off'); % Run in background
        
        % Plot CI (Base layer)
        plot(ebsd, cibw, 'micronbar', 'off');
        hold on
        
        % Plot IPF orientations with transparency
        plot(ebsd('indexed'), colors, 'FaceAlpha', 0.5);
        
        % Plot General Boundaries (Dark Gray)
        plot(gB, 'linecolor', [0.2, 0.2, 0.2], 'linewidth', 1.0);
        
        % Plot Twin Boundaries (Blue)
        plot(twinninggB, 'linecolor', 'b', 'linewidth', 2.0);
        
        title(['Twin Boundaries: ', baseFileName], 'Interpreter', 'none');
        
        % Save Figure
        saveName = fullfile(dirPath, [baseFileName, '_twinPlot.png']);
        exportScaledFigure(f1, saveName, 'extraheight');
        
        close(f1);
        fprintf('Saved: %s\n', saveName);
        
    catch ME
        fprintf('Error processing %s: %s\n', baseFileName, ME.message);
    end
end

fprintf('Batch processing complete.\n');