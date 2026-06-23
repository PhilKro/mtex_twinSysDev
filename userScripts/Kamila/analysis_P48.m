%% initialize mtex
clear 
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% filename
dir = 'userScripts/Kamila';
file = 'Sample48 Site 3 Map Data 6_Map_Data_6_CircMask_EntireMap_NoGB_BW158.ang';
% dir = 'onedriveData/Kamila';
% file = 'pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';

path = fullfile(dir, file);

% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.10;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Re
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');
ebsd = rotate(ebsd,rot);

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
indebsd = ebsd('indexed');
refebsd = indebsd(1);
reforientation = refebsd.orientations;
clear indebsd
reforientation = orientation.byEuler(101.3*degree,10.3*degree,260.7*degree);
% reforientation = orientation.byEuler(270.2*degree,79.7*degree,92.3*degree);

ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
ebsd0=ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;

%% some initial test plot %%%

ipf=ipfColorKey(ebsd.CS);

figure;
ipf.inversePoleFigureDirection=xvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors)
title('IPF X')

ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
nextAxis
plot(ebsd('indexed'),colors)
title('IPF Y')

ipf.inversePoleFigureDirection = zvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
nextAxis
plot(ebsd('indexed'),colors)
title('IPF Z')

nextAxis
cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
cibw = cibw.xyz;
ebsd.prop.ciColor = cibw;
plot(ebsd,cibw)
colormap gray


% ipf with ci greyscale
figure;
plot(ebsd,cibw,'micronbar','off')
hold on
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
title('IPF Y')

%% grain reconstruction
% min angle used for grain reconstruction
grainThresh = 8*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 40;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 2;

%%% perform reconstruction %%%
ebsd = ebsd0;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

% take only boundaries between indexed grains
gB = grains.boundary(CS.mineral,CS.mineral);

%%% test plot %%%

figure;
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(grains.meanOrientation);
% plot(ebsd,cibw,'micronbar','off')
% hold on
plot(grains,colors)%,'FaceAlpha',0.7)

% create and plot crystal shapes
cS1 = crystalShape.hex(ebsd.CS);

hold on
%plottablecS = grains.meanOrientation*cS1*0.5
plot(grains,0.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% plot(grains(~isParentGrain),0.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% plot(grains(isParentGrain),0.2*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
% % plot boundaries
plot(grains.boundary)
text(grains,grains.id)
%% fill ebsd data
ebsd_filled = fill(ebsd, grains);
figure; plot(ebsd,ebsd.orientations); nextAxis; plot(ebsd_filled, ebsd_filled.orientations); 

%% Graph-Based Twin Analysis (11-21 Twins)

% 1. Define the 11-21 twin system
tS_1121 = twinSystem.hexagonal_1121(CS);

% 2. Initialize and run the TwinNetwork engine
tN = TwinNetwork(grains, gB, tS_1121, 'twinThresh', 5*degree, 'maxGenerationLimit', 2);
tN = tN.buildTopology();
tN = tN.propagateReferences();
tN = tN.calculateDeformation();

% 3. Plotting
tN.plotQuotientGraph();
tN.plotGenerations();
tN.plotTwinVariants();
