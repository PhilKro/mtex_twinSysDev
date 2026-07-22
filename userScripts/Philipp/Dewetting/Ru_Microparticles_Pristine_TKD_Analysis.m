%% Ru Pristine TKD Analysis
clear
close all

plotx2east
plotzIntoPlane
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% Directories and Input Setup
inDir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0011_Ru_Microparticles_pristine_TKD';
outDir = fullfile(inDir, 'Analysis');

if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created missing directory: %s\n', outDir);
end

% Specify your ang files here. Example:
ang_files = {'Project 1 Specimen 1 Site 3 Map Data 6_Map_Data_6_BW263.ang','Project 1 Specimen 1 Site 4 Map Data 9_Map_Data_9_BW263.ang','Project 1 Specimen 1 Site 5 Map Data 11_Map_Data_11_BW263.ang'}; 
% Populate this array with the actual ang files you want to analyze later

%% Universal Parameters
tkd_tilt = -20 * degree; % Adjust based on your TKD setup (from analysis_P48.m)
thresh_ci = 0.10;
grainThresh = 5 * degree;
minSize = 10;
alpha = 2.2;

%% Crystal Symmetry for Ru (hcp)
CS = crystalSymmetry('6/mmm', [2.706 2.706 4.282], 'x||a', 'mineral', 'Ru');

%% Process Files
for i = 1:length(ang_files)
    file = ang_files{i};
    path = fullfile(inDir, file);
    
    if ~exist(path, 'file')
        warning('File not found: %s. Skipping...', path);
        continue;
    end
    
    [~, baseName, ~] = fileparts(file);
    fprintf('Processing %s...\n', baseName);
    
    %% Load and clean EBSD data
    ebsd = EBSD.load(path, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    
    % Perform corrections for tkd tilt and misaligned axes
    ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, -70*degree + tkd_tilt), 'keepXY');
    ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
    ebsd.phase(ebsd.isnan) = -1;
    
    rot = rotation.byAxisAngle(zvector, 180*degree);
    ebsd = rotate(ebsd, rot, 'keepXY');
    ebsd = rotate(ebsd, rot);
    
    %% TKD Plot with CI background
    fprintf('  -> Generating IPF Y with CI background...\n');
    ipf = ipfColorKey(ebsd.CS);
    ipf.inversePoleFigureDirection = yvector;
    
    % CI background
    cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
    cibw = cibw.xyz;
    
    figure;
    plot(ebsd, cibw, 'micronbar', 'off');
    hold on;
    colors = ipf.orientation2color(ebsd('indexed').orientations);
    plot(ebsd('indexed'), colors, 'FaceAlpha', 0.5);
    title('IPF Y with CI overlay');
    
    exportScaledFigure(gcf, fullfile(outDir, [baseName, '_IPF_Y_CI_overlay.jpg']), 'target_px', 4000, 'noDateString');

    figure;
    ipf.inversePoleFigureDirection = zvector;

    plot(ebsd, cibw, 'micronbar', 'off');
    hold on;
    colors = ipf.orientation2color(ebsd('indexed').orientations);
    plot(ebsd('indexed'), colors, 'FaceAlpha', 0.5);
    title('IPF Z with CI overlay');
    
    exportScaledFigure(gcf, fullfile(outDir, [baseName, '_IPF_Z_CI_overlay.jpg']), 'target_px', 4000, 'noDateString');
    
    %% Grain Reconstruction
    fprintf('  -> Reconstructing grains...\n');
    [grains, ebsd.grainId] = calcGrains(ebsd, 'minPixel', minSize, 'angle', grainThresh, 'alpha', alpha, 'boundary', 'tight');
    
    %% GND Analysis
    fprintf('  -> Computing curvature and fitting dislocations for %s...\n', baseName);
    
    % Filter and gridify for curvature calculation
    F = halfQuadraticFilter;
    F.alpha = 3;
    ebsdSGND = gridify(smooth(ebsd('indexed'), F));
    kappa = ebsdSGND.curvature;
    
    % HCP dislocation systems
    dS = dislocationSystem.hcp(ebsd.CS); 
    
    nu = 0.3;
    dS(dS.isEdge).u = 1 / (1 - nu);
    dS(dS.isScrew).u = 1;
    
    dSRot = ebsdSGND.orientations * dS;
    [rho, factor] = fitDislocationSystems(kappa, dSRot);
    
    % Plot Total GND Density
    figure;
    gnd_density_sum = factor * sum(abs(rho), 2, 'omitnan');
    plot(ebsdSGND, log10(gnd_density_sum), 'micronbar', 'off');
    hold on;
    plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    mtexColorMap parula;
    setColorRange([13.5, 16.5]); % Typical range for log10(m^-2), adjust if needed
    mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
    title('Total GND Density');
    
    exportScaledFigure(gcf, fullfile(outDir, [baseName, '_Total_GND.jpg']), 'target_px', 4000, 'noDateString');
    
    % Plot Edge and Screw GND Density
    figure;
    newMtexFigure('layout', [1, 2]);
    
    GND_E = factor * sum(abs(rho(:, dS.isEdge)), 2, 'omitnan');
    GND_S = factor * sum(abs(rho(:, dS.isScrew)), 2, 'omitnan');
    
    nextAxis(1, 1);
    plot(ebsdSGND, log10(GND_E), 'micronbar', 'off');
    hold on;
    plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    mtexColorMap parula;
    title('Total Edge GND Density');
    
    nextAxis(1, 2);
    plot(ebsdSGND, log10(GND_S), 'micronbar', 'off');
    hold on;
    plot(grains.boundary, 'linewidth', 1.5, 'linecolor', 'w');
    mtexColorMap parula;
    title('Total Screw GND Density');
    
    setColorRange([13.5, 16.5]); 
    mtexColorbar('title', 'log_{10}(GND Density m^{-2})');
    drawNow(gcm, 'figSize', 'large');
    
    exportScaledFigure(gcf, fullfile(outDir, [baseName, '_Edge_Screw_GND.jpg']), 'target_px', 4000, 'noDateString');
    
    close all;
end

fprintf('Analysis complete.\n');
