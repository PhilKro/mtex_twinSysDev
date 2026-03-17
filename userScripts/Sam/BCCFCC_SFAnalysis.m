%% initialize mtex

clear 
close all
load('CS_cS_andPlanes.mat')


%%
% Mtex directions (important when working with EBSD)
% setMTEXpref('xAxisDirection','east');
% setMTEXpref('zAxisDirection','intoPlane');

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% Load data: Nanomegas .ang

% path = 'K:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240725_Re_7dc_tkdPTP\06_TEM_Observation_20240814_PTP_7degcaxis_Re\Nanomegas\';
% file = '20240814_full_ref_res.ang';

% path = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240625_Re_c_tkdPTP\06_TEM_Observation_20240726_PTP_Insitu_Re_twin_Philipp_01\Nanomegas\';
file = 'userScripts\Sam\BCCFCC Deformed pillar scan3 rescan without refinement.ang';


% % define tilt of TEM sample
% %default
% flipud_bool = false;
% zrot = 0 * degree;
% 
% temAlphaTilt = 0 * degree;
% temBetaTilt = 0 * degree;

% settings for 20240728_1432_Full.ang
flipud_bool = true;
zrot = 0 * degree;

temAlphaTilt = 0 * degree;
temBetaTilt = 0; % 32 * degree;

% settings for 20240728_1240_Full.ang
% flipud_bool = true;
% zrot = 0 * degree;
% 
% temAlphaTilt = 0 * degree;
% temBetaTilt = 0 * degree;

% settings for 20240814_full_ref_res.ang
% flipud_bool = true;
% zrot = 90 * degree;
% 
% temAlphaTilt = -58 * degree;
% temBetaTilt = 0 * degree;


% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.10;

% crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% used
CS = {crystalSymmetry('m-3m','mineral','bcc'), crystalSymmetry('m-3m','mineral','fcc')};


%%% perform corrections for tkd tilt and misaligned axes %%%
ebsd = EBSD.load(file,CS,'convertEuler2SpatialReferenceFrame','setting 4');

%ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';

% crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% used
% CS = crystalSymmetry('6/mmm',[2.7 2.7 4.4],'x||a','mineral','Re');


%%% perform corrections for tem tilt and misaligned axes %%%
% ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,temAlphaTilt));
% ebsd = rotate(ebsd,rotation.byAxisAngle(yvector,temBetaTilt));

if (flipud_bool)
    rot = rotation.byAxisAngle(xvector,180*degree);
    ebsd = rotate(ebsd,rot);
end

% fine adjustment (get parentOrientation from TKD scan)
% ebsd = rotate(ebsd, inv(rotation.byAxisAngle(axis(parentOrientation_TEM, parentOrientation_TKD), angle(parentOrientation_TEM, parentOrientation_TKD))));

%use for 20240728_1432_Full.ang to be consistent with Re_cAxis_ptp1_insitu Specimen 1 Site 2 Map Data 2_BW158.ang
% reforientation = orientation.byEuler(184.137*degree,79.9861*degree,241*degree);
% ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);

ebsd0=ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;
%%
ipf1=ipfColorKey(ebsd('fcc').CS);
ipf2=ipfColorKey(ebsd('bcc').CS);

figure;
dir = yvector;
ipf1.inversePoleFigureDirection=dir;
colors = ipf1.orientation2color(ebsd('fcc').orientations);
plot(ebsd('fcc'),colors)
hold on
ipf2.inversePoleFigureDirection=dir;
colors = ipf2.orientation2color(ebsd('bcc').orientations);
plot(ebsd('bcc'),colors)
title('IPF Y')

figure;
dir = zvector;
ipf1.inversePoleFigureDirection=dir;
colors = ipf1.orientation2color(ebsd('fcc').orientations);
plot(ebsd('fcc'),colors)
hold on
ipf2.inversePoleFigureDirection=dir;
colors = ipf2.orientation2color(ebsd('bcc').orientations);
plot(ebsd('bcc'),colors)
title('IPF Y')
%% grain calc

% min angle used for grain reconstruction
grainThresh = 5*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
grains0 = grains;

%% SF calc fcc
% ebsd = ebsd0;
% grains = grains0;
% ebsd = ebsd('fcc');
% grains = grains('fcc');
% 
% sS = slipSystem.fcc(ebsd.CS);
% sS = sS.symmetrise;
% sSLocal = grains.meanOrientation * sS;
% 
% sigma = stressTensor.uniaxial(yvector)
% SF = abs(sSLocal.SchmidFactor(sigma));
% 
% % take the maximum along the rows
% [SFMax,active] = max(SF,[],2);
% 
% % plot the maximum Schmid factor
% figure;
% plot(grains,SFMax,'micronbar','off','linewidth',2)
% mtexColorbar location southoutside
% title('')
% figure;
% histogram(SFMax)

%% Calc all the SFs 
vectors = {xvector, yvector, zvector};
vNames = {'X', 'Y', 'Z'};
phases = {'bcc', 'fcc'};

for p = 1:length(phases)
    phaseName = phases{p};
    
    % Filter data for current phase
    ebsd_p = ebsd0(phaseName);
    grains_p = grains0(phaseName);
    
    if isempty(grains_p)
        continue; % Skip if phase doesn't exist in dataset
    end
    
    % Define Slip Systems based on phase
    if strcmp(phaseName, 'bcc')
        % Standard BCC slip systems {110}<111>
        sS = slipSystem.bcc(ebsd_p.CS);
    else
        % Standard FCC slip systems {111}<110>
        sS = slipSystem.fcc(ebsd_p.CS);
    end
    sS = sS.symmetrise;
    
    % Calculate local slip systems based on mean grain orientation
    sSLocal = grains_p.meanOrientation * sS;
    
    for v = 1:length(vectors)
        currentVec = vectors{v};
        currentVName = vNames{v};
        
        % Define Uniaxial Stress Tensor
        sigma = stressTensor.uniaxial(currentVec);
        
        % Calculate SF and find the maximum per grain
        SF = abs(sSLocal.SchmidFactor(sigma));
        [SFMax, ~] = max(SF, [], 2);
        
        % Plot 1: Spatial Map
        figMap = figure('Visible', 'on');
        plot(grains_p, SFMax, 'micronbar', 'off', 'linewidth', 1);
        mtexColorbar('title', 'Max Schmid Factor');
        colormap('viridis');
        
        titleName = sprintf('Max SF Map: %s (Dir: %s)', upper(phaseName), currentVName);
        title(titleName);
        
        % Export Map
        exportNameMap = sprintf('SF_Map_%s_%s.png', phaseName, currentVName);
        exportgraphics(figMap, exportNameMap, 'Resolution', 300);
        
        % Plot 2: Histogram with Statistics
        figHist = figure('Visible', 'on');
        h = histogram(SFMax, 20, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'w');
        grid on;
        
        % Calculate Statistics
        stats_mean = mean(SFMax);
        stats_max = max(SFMax);
        stats_min = min(SFMax);
        stats_std = std(SFMax);
        
        % Add statistics text box
        statsStr = {['Mean: ', num2str(stats_mean, '%.3f')], ...
                    ['Max: ', num2str(stats_max, '%.3f')], ...
                    ['Min: ', num2str(stats_min, '%.3f')], ...
                    ['Std: ', num2str(stats_std, '%.3f')]};
        
        annotation('textbox', [0.15, 0.65, 0.2, 0.2], 'String', statsStr, ...
                   'FitBoxToText', 'on', 'BackgroundColor', 'white');
        
        title(sprintf('SF Distribution: %s (%s)', upper(phaseName), currentVName));
        xlabel('Schmid Factor');
        ylabel('Frequency (Grains)');
        
        % Export Histogram
        exportNameHist = sprintf('SF_Hist_%s_%s.png', phaseName, currentVName);
        exportgraphics(figHist, exportNameHist, 'Resolution', 300);
        
        % Optional: Close figures to save memory if processing many files
        % close(figMap); close(figHist);
    end
end
