%% %% RD33Chl
%  Import Script for EBSD Data

%% CLEAR ALL
close all
clear
clc

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {...
    'notIndexed',...
    crystalSymmetry('2/m', [5.237 9.07 14.285], [90,97.90,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Chlorite', 'color', 'light green'),...
    crystalSymmetry('2/m',  [6.5540 8.7080 7.0690],[90, 113.9, 90]*degree, 'X||a*','Y||b','Z||c','mineral', 'Titanite', 'color', 'magenta'),...
    crystalSymmetry('-3',  [5.0691 5.0691 13.9849],[90, 90, 120]*degree, 'X||a*', 'Z||c','mineral', 'Ilmenite', 'color', [0.8500 0.3250 0.0980]),...
    crystalSymmetry('4/mmm',  [4.5940 4.5940 2.9580],[90, 90, 90]*degree, 'X||a*','Y||b','Z||c','mineral', 'Rutile', 'color', 'yellow')};

%% Miller Indices

h_chl = [Miller(1,0,0,CS{2},'hkl'), Miller(0,1,0,CS{2},'hkl'),Miller(0,0,1,CS{2},'hkl')];
h_tit = [Miller(2,2,1,CS{3}), Miller(-1,-1,2,CS{3}), Miller(1,0,2,CS{3})];
h_ilm = [Miller(1,1,0,CS{4}), Miller(0,1,0,CS{4}), Miller(0,0,1,CS{4})];
h_rut = [Miller(1,1,0,CS{5}), Miller(0,1,0,CS{5}), Miller(0,0,1,CS{5})];


%% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names
fname = 'userScripts/RD33Chl Specimen 1 Site 1 Map Data 1.crc';

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','crc',...
  'convertEuler2SpatialReferenceFrame','setting 2');

ebsd = rotate(ebsd,rotation('axis',vector3d.X,'angle',180*degree));

%% Clean-up
% --- Chlorite (CS{2}) ---
ps_chl_a120 = orientation.byAxisAngle(Miller(1,0,0,CS{2}), 120*degree);
ps_chl_b120 = orientation.byAxisAngle(Miller(0,1,0,CS{2}), 120*degree);
ps_chl_c120 = orientation.byAxisAngle(Miller(0,0,1,CS{2}), 120*degree);

pseudoSym = [ps_chl_a120, ps_chl_b120, ps_chl_c120];

%% Simple grain construction and cleanup (probably want to modify)
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',[1 10]*degree);

figure;
plot(ebsd('Chlorite'), ebsd('Chlorite').orientations)
nextAxis
plot(grains)
saveas(gcf,'preCleanup.jpg');

%%
[ebsd_cleaned,grains_cleaned, numChanged] = cleanUpPseudoSym(ebsd,grains,pseudoSym, 'threshold', 1.3);
numChanged

%%
figure;
plot(ebsd_cleaned('Chlorite'), ebsd_cleaned('Chlorite').orientations)
nextAxis
plot(grains_cleaned)
saveas(gcf,'postCleanup.jpg');
%% 
ebsd = ebsd_cleaned;
grains = grains_cleaned;
%% Remove potentially spurious "grains"
grains = grains(grains.numPixel>5);

% reset the EBSD data to filtered grains
ebsd = ebsd(grains);

%% Recompute grain boundaries with the remaining EBSD data
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',[1 10]*degree,'Boundary','tight');
%% Filtering

% Filtering
F = halfQuadraticFilter;
F.threshold = 1.5*degree;
F.eps = 1e-2;
F.alpha = 0.01;

% Remove un-indexed grains
% ebsd(grains(grains.grainSize<=3))=[];
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'threshold',Min_gB_angle*degree);
ebsd = smooth(ebsd('indexed'),F,'fill',grains);

%% SMOOTHING
grains = grains.smooth;

%% PHASE MAP
figure,
plot(grains('indexed'),'figSize','large')
hold on
plot(ebsd('notIndexed'),'FaceColor','white')
plot(grains.boundary,'lineWidth',1.0)
hold off

%% PHASE MAP_02
figure,
plot(grains('Chlorite'),'FaceColor',[0.38, 0.84, 0.55],'figSize','large')
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend
hold off
saveas(gcf,'Phase_map.png');

%% Orientation maps X, Y, and Z

% Create the colourkey
oM1=ipfHSVKey(ebsd('Chlorite'));

% Plot the figure
figure;
oM1.inversePoleFigureDirection=xvector; %IPFx wrt X
plot(ebsd('Chlorite'),oM1.orientation2color(ebsd('Chlorite').orientations),'figSize','large');
mtexTitle('Chlorite x-axis')
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
hold off
saveas(gcf,'Orientation_X.png');

figure
nextAxis
oM1.inversePoleFigureDirection=yvector; %IPFx wrt Y
plot(ebsd('Chlorite'),oM1.orientation2color(ebsd('Chlorite').orientations),'figSize','large');
mtexTitle('Chlorite y-axis')
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
hold off
saveas(gcf,'Orientation_Y.png');

figure
nextAxis
oM1.inversePoleFigureDirection=zvector; %IPFx wrt Z
plot(ebsd('Chlorite'),oM1.orientation2color(ebsd('Chlorite').orientations),'figSize','large');
mtexTitle('Chlorite y-axis')
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
hold off
saveas(gcf,'Orientation_Z.png');

%% Orientation - Ilmenite 

% Create the colourkey
oM1=ipfHSVKey(ebsd('Titanite'));

% Plot the figure
figure;
oM1.inversePoleFigureDirection=yvector; %IPFx wrt X
plot(ebsd('Titanite'),oM1.orientation2color(ebsd('Titanite').orientations),'figSize','large');
mtexTitle('Titanite x-axis')
hold on
plot(ebsd('Chlorite'),'FaceColor',[0.8, 0.8, 0.8])
plot(ebsd('Ilmenite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.6);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.6);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
% scale('off')
hold off
saveas(gcf,'Orientation_Ilmenite_y.png');


%% ANOTHER WAY - Chlorite 
grod = ebsd('Chlorite').calcGROD(grains);
figure
plot(ebsd('Chlorite'),grod.angle./degree,'micronbar','on','figSize','large')
parulaColorMap
mtexColorbar('Title','GROD [degree]');
setColorRange([0 30]);
% set(gca,'Color','k')
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',2.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',2.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',2.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',2.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',2.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
hold on
alpha = grains('Chlorite').innerBoundary.misorientation.angle / (10*degree);
plot(grains('Chlorite').innerBoundary,'linewidth',1.8,'edgeAlpha',alpha,'linecolor','red');
hold off
saveas(gcf,'M2M_Chl.png');

%% Misorientation angle of grain boundaries 

figure
gB_Chl = grains.boundary('Chlorite','Chlorite');
plot(grains('Chlorite'),'translucent',0.3,'micronbar','off')
legend off
hold on
plot(gB_Chl,gB_Chl.misorientation.angle./degree,'linewidth',1.8)
hold off
mtexColorbar('title','misorientation angle')
saveas(gcf,'M2M_gB_Chl.png');

%% Misorientation at Subgrain Boundaries

mori = grains('Chlorite').innerBoundary.misorientation;

figure
plot(mori.axis,'fundamentalRegion','contourf','figSize','small')

mtexColorbar


%% ANOTHER WAY - Titanite 
grodt= ebsd('Titanite').calcGROD(grains);
figure
plot(ebsd('Titanite'),grodt.angle./degree,'micronbar','on','figSize','large')
parulaColorMap
mtexColorbar('Title','GROD [degree]');
setColorRange([0 10]);
hold on
plot(ebsd('Chlorite'),'FaceColor',[0.8, 0.8, 0.8])
plot(ebsd('Ilmenite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains.boundary,'linewidth',1.0);

alpha_tit = grains('Titanite').innerBoundary.misorientation.angle / (10*degree);
plot(grains('Titanite').innerBoundary,'linewidth',1.8,'edgeAlpha',alpha_tit,'linecolor','red');

hold off
saveas(gcf,'M2M_Ttn.png');

%% GOS (Grain Orientations Spread)
mis2mean = calcGROD(ebsd, grains); % This is for all the grains
GOS = ebsd.grainMean(mis2mean.angle);

figure,
plot(grains('Chlorite'),grains('Chlorite').GOS./degree,'figSize','large')
cb = mtexColorbar('Title','GOS  [degree]');
hold on
plot(ebsd('Ilmenite'),'FaceColor',[0.25, 0.25, 0.25])
plot(ebsd('Titanite'),'FaceColor',[0.35, 0.35, 0.35])
plot(ebsd('Rutile'),'FaceColor',[0.45, 0.45, 0.45])
plot(ebsd('notIndexed'),'FaceColor','white');
hold on
plot(grains('Chlorite','Chlorite').boundary,'linewidth',2.0);
plot(grains('Chlorite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Titanite').boundary,'linewidth',1.0);
plot(grains('Chlorite','Rutile').boundary,'linewidth',1.0);
plot(grains('Chlorite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Ilmenite').boundary,'linewidth',1.0);
plot(grains('Titanite','Titanite').boundary,'linewidth',1.0);
plot(grains('Rutile','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Titanite','notIndexed').boundary,'linewidth',1.0);
plot(grains('Rutile','notIndexed').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Titanite').boundary,'linewidth',1.0);
plot(grains('Titanite','Rutile').boundary,'linewidth',1.0);
plot(grains('Ilmenite','Rutile').boundary,'linewidth',1.0);
legend('off')
hold off
saveas(gcf,'GOS_Chl.png');

%% Rudigen 
% MDF = calcDensity(grains.boundary('Chlorite','Chlorite').misorientation,'halfwidth',2*degree);
% [v,pos] = max(MDF,'numLocal', 3); 


% [v,pos] = max(MDF,'numLocal',3); 
% [v,pos] = max(MDF) % the 5 largest local maxima
% [v,pos] = max(MDF,'numLocal',5) % the 5 largest local maxima


%% CPO (Area Weighted)

o_All = ebsd('Chlorite').orientations;
odf_All = calcODF(o_All,'halfwidth',10*degree);

% Plot
figure
plotPDF(odf_All,h_chl,'colorrange',[0 10],'hemisphere','Upper','antipodal','figSize','medium');
setColorRange([0 10]);
cb = mtexColorbar('Title','M.U.D.- Area');
mtexColorMap(WhiteJetColorMap)
set(gca,'Color',[.5 .5 .5])
saveas(gcf,'CPO_Chl_area_weighted.png');

%% CPO (One Point Per Grain) - TWO Ways

o_chl = ebsd('Chlorite').orientations;
odf_chl = calcODF(o_chl,'halfwidth',10*degree);

% FIRST WAY
odf_1ppg = calcDensity(ebsd('Chlorite').orientations,'halfwidth',10*degree);

% Number of grains
ng = length(grains('Chlorite')); % Number of grains
ngrains = "N = " + ng; % Variable for box 2

% Pole figure ODF
figure
plotPDF(odf_1ppg, h_chl, 'antipodal', 'MarkerSize', 2, 'micronBar', 'off', 'figSize', 'medium');
annotation('textbox','String',ngrains,'FitBoxToText','on','LineStyle','none');
% mtexColorMap(parula);
hold on;
setColorRange([0 20]);
mtexColorbar('Title', 'M.U.D. - 1 Point per grain');
set(gca, 'Color', [.5 .5 .5]);
hold off
saveas(gcf,'CPO_Chl_1PPG_counter.png');

% CPO by point per grains
% cs = ebsd('Chlorite').CS;
% ori = ebsd('Chlorite').orientations;

% 1 point per grain 
figure
plotPDF(grains('Chlorite').meanOrientation,[Miller(1,0,0,CS{2},'hkl'),Miller(0,1,0,CS{2},'hkl'),Miller(0,0,1,CS{2},'hkl')], 'MarkerSize',6,'MarkerFaceColor','green','figSize', 'medium')
saveas(gcf,'CPO_1PPG_map.png');

%% CPO (ANOTHER WAY - One Point Per grain)
% Variables

psi = calcKernel(grains('Chlorite').meanOrientation);
odf_1ppg = calcDensity(ebsd('Chlorite').orientations,'kernel',psi);

% Plot
figure
plotPDF(odf_1ppg,h_chl,'antipodal','MarkerSize',2,'micronBar','off','figSize','medium');
annotation('textbox','String',ngrains,'FitBoxToText','on','lineStyle','none');
setColorRange([0 10]);
cb = mtexColorbar('Title','M.U.D.- 1 Point per grain');
mtexColorMap(WhiteJetColorMap)
set(gca,'Color',[.5 .5 .5])
saveas(gcf,'CPO_Chl_1PPG_kernel.png');

%% GROD ANALYSIS

% Grain construction
ebsd_chl = ebsd ('Chlorite');
[grains_chl,ebsd_chl.grainId] = calcGrains(ebsd_chl('indexed'));
ebsd_chl(grains_chl(grains_chl.grainSize <= 2)) = [];
[grains_chl,ebsd_chl.grainId] = calcGrains(ebsd_chl('indexed'),'threshold',[1*degree, 10*degree]);

% Smooth and filter
grains_chl = smooth(grains_chl,5);
F = halfQuadraticFilter;
% ebsd_chl = smooth(ebsd_chl,F,grains_chl,'fill');

% GROD axes
grod = ebsd_chl.calcGROD(grains_chl);
axSpecimen = ebsd_chl.orientations .* grod.axis('noSymmetry');

% GROD Sample reference frame (SS)
figure
plot(axSpecimen,'contourf','fundamentalRegion','halfwidth',10*degree,'hemisphere','upper')
DataNum = numel(axSpecimen.x);
mtexTitle(sprintf('Chl %d data points', DataNum))
cb = mtexColorbar('Title','M.U.D.');
mtexColorMap(parula);
hold on 
saveas(gcf,'MOA_GROD_Sample_scale.png');

% GROD Crystal reference frame (CS)
figure
axCrystal = grod.axis;
plot(axCrystal,'contourf','fundamentalRegion','antipodal','figSize','large')
DataNum = numel(axCrystal.h);
mtexTitle(sprintf('Chl %d data points', DataNum))
cb = mtexColorbar('Title','M.U.D.');
mtexColorMap(parula);
hold on 
annotate(Miller({1,0,0}, {0,1,0}, {0,0,1},{1,0,1}, {-1,0,1},CS{2},'hkl'),'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_GROD_Crystal_scale.png');

%% SPATIAL DISTRIBUTION OF THE MISORIENTATION AXES
% Colour key
figure
colorKey = HSVDirectionKey(ebsd_chl.CS,'antipodal');
plot(colorKey,'figSize','small')
color = colorKey.direction2color(axCrystal);
alpha = min(grod.angle/degree/7.5,1);

% plot the data
figure
plot(ebsd_chl,color,'micronbar','off','faceAlpha',alpha,'figSize','large')
hold on
plot(grains.boundary,'linewidth',1.0)
hold off
saveas(gcf,'Spatial_distribution_of_MOA.png');

%% ********* CPO of Titanite *************

o_tit = ebsd('Titanite').orientations;
odf_tit = calcODF(o_tit,'halfwidth',10*degree);
odf_1ppg = calcDensity(ebsd('Titanite').orientations,'halfwidth',10*degree);
ng = length(grains('Titanite')); % Number of grains
ngrains = "N = " + ng; % Variable for box 2
figure
plotPDF(odf_1ppg, h_tit, 'antipodal', 'MarkerSize', 2, 'micronBar', 'off', 'figSize', 'medium');
annotation('textbox','String',ngrains,'FitBoxToText','on','LineStyle','none');
% mtexColorMap(parula);
hold on;
setColorRange([0 10]);
mtexColorbar('Title', 'M.U.D. - 1 Point per grain');
set(gca, 'Color', [.5 .5 .5]);
hold off
saveas(gcf,'CPO_Tit_1PPG_counter.png');

%% GROD ANALYSIS - TITANITE

% Grain construction
ebsd_tit = ebsd ('Titanite');
[grains_tit,ebsd_tit.grainId] = calcGrains(ebsd_tit('indexed'));
ebsd_tit(grains_tit(grains_tit.grainSize <= 2)) = [];
[grains_tit,ebsd_tit.grainId] = calcGrains(ebsd_tit('indexed'),'threshold',[1*degree, 10*degree]);

% Smooth and filter
grains_tit = smooth(grains_tit,5);
F = halfQuadraticFilter;
% ebsd_tit = smooth(ebsd_tit,F,grains_tit,'fill');

% GROD axes
grod = ebsd_tit.calcGROD(grains_tit);
axSpecimen = ebsd_tit.orientations .* grod.axis('noSymmetry');

% GROD Sample reference frame (SS)
figure
plot(axSpecimen,'contourf','fundamentalRegion','halfwidth',10*degree,'hemisphere','upper')
DataNum = numel(axSpecimen.x);
mtexTitle(sprintf('Tit %d data points', DataNum))
cb = mtexColorbar('Title','M.U.D.');
mtexColorMap(parula);
hold on 
saveas(gcf,'MOA_GROD_Sample_scale_tit.png');

% GROD Crystal reference frame (CS)
figure
axCrystal = grod.axis;
plot(axCrystal,'contourf','fundamentalRegion','antipodal','figSize','large')
DataNum = numel(axCrystal.h);
mtexTitle(sprintf('Tit %d data points', DataNum))
cb = mtexColorbar('Title','M.U.D.');
mtexColorMap(parula);
hold on 
annotate(Miller({1,0,0}, {0,1,0}, {0,0,1},{1,0,1}, {-1,0,1},CS{2},'hkl'),'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_GROD_Crystal_scale_tit.png');

%% MOA - Misorientation angle distribution (correlated)
mori = calcMisorientation(ebsd('Chlorite'),ebsd('Chlorite'))
odf_chl = calcODF(ebsd('Chlorite').orientations,'halfwidth',10*degree)
mdf_uncor = calcODF(mori)
mdf = calcMDF(odf_chl)

figure
CMo_All = grains.boundary('Chlorite','Chlorite').misorientation;
plotAngleDistribution(CMo_All,'dispName','Correlated misorientations','FaceColor','light green','figSize','large');
% hold on
% plotAngleDistribution(mdf)
hold on
plotAngleDistribution(odf_chl.CS,odf_chl.CS,'displayName','Random','linewidth',2,'Color','black')
legend
hold off
saveas(gcf,'MOA_Hist_map_02.png');

%% **************** Twinning - Titanite *************

% % Highlight twinning boundaries
mori = orientation.map(Miller(1,0,0,CS{3}),Miller(0,1,0,CS{3}),...
  Miller(0,0,1,CS{3})),Miller(1,0,1,CS{3}),Miller(-1,0,1,CS{3}),Miller(2,2,1,CS{3}),Miller(-1,-1,2,CS{3}),Miller(1,0,2,CS{3});

% the rotational axis
round(mori.axis)

% the rotational angle
mori.angle / degree

% consider only Magnesium to Magnesium grain boundaries
gB_Tn = grains.boundary('Titanite','Titanite');
isTwinning = angle(gB_Tn.misorientation,mori) < 7.5*degree;
% plot the grains and highlight the twinning boundaries
plot(grains('Titanite'),grains('Titanite').meanOrientation,'micronbar','off')
hold on
plot(gB_Tn(isTwinning),'edgecolor','red','linewidth',3.0)
hold off

% % Misorientation angle distribution function 
MOA_Tn = grains.boundary('Titanite','Titanite').misorientation;
odf_Tn = calcODF(ebsd('Titanite').orientations,'halfwidth',10*degree)

figure
plotAngleDistribution(MOA_Tn,'dispName','Correlated misorientations','FaceColor','light green','figSize','medium');
hold on
plotAngleDistribution(odf_Tn.CS,odf_Tn.CS,'displayName','Random','linewidth',2,'Color','black')
legend
hold off
saveas(gcf,'MOA_Hist_Titanite.png');

% % MOA of twin boundaries 
figure
plotAxisDistribution(gB_Tn.misorientation,'contour')
hold on 
annotate(Miller({1,0,0}, {0,1,0}, {0,0,1},{1,0,1}, {-1,0,1},CS{3},'hkl'),'all','labeled','BackgroundColor','white','figSize','medium');
hold off
saveas(gcf,'MOA_counter_Titanite.png');

%% Boundary axis distribution

bnd_chlchl = grains.boundary('Chlorite','Chlorite'); 

figure
plotAxisDistribution(bnd_chlchl.misorientation,'smooth', 'figSize','large')
mtexTitle('Chl-Chl boundary axis distribution');
mtexColorbar('Title','M.U.D.');
mtexColorMap(parula);
hold on 
annotate(Miller({1,0,0}, {0,1,0}, {0,0,1},{1,0,1}, {-1,0,1}, CS{2},'hkl'),'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_Chl_Chl.png');

%%  **********************************************************************************
%                              Grain Boundary Propoerties
% *************************************************************************************

% A) %%%%%%%%%%% MOA on Grain Boundaries %%%%%%%%%%% 
figure
plot(grains('Chlorite'),'translucent',.5,'micronbar','off')
legend off
hold on
% this reorders the boundary segement a a connected graph which results in
% a smoother plot
gB_Chl = gB_Chl.reorder;
% ipfKey = PatalaColorKey(gB_Chl);F
oM1=ipfHSVKey(ebsd('Chlorite'));
plot(gB_Chl,'linewidth',1.5)
% on my computer setting the renderer to painters gives a much more
% pleasent result
%set(gcf,'Renderer','painters')
hold on
plot(gB_Chl,oM1.orientation2color(gB_Chl.misorientation),'linewidth',1.5)
hold off
saveas(gcf,'MOA_on_grain_boundaries_map_05.png');

% C) %%%%%%%%%%% Grain Boundaries with different angles %%%%%%%%%%% 

% Compute mean misorientation of each Chlorite boundary
angle = mean(gB_Chl.misorientation.angle,2) / degree;  % mean for each boundary

% Categorize boundaries
id = zeros(size(angle));
id(angle < 10) = 1;            
id(angle >= 10 & angle < 20) = 2;
id(angle >= 20 & angle <= 50) = 3;  
id(angle >= 50 & angle <= 100) = 4;  
id(angle > 100) = 5;

% Plot
figure
gB_all = grains.boundary('Chlorite');
plot(gB_all,'linecolor','k','linewidth',1)  % all boundaries in black
hold on
plot(gB_Chl(id==1),'linecolor','yellow','linewidth',1.5,'DisplayName','<10°')
plot(gB_Chl(id==2),'linecolor','green','linewidth',1.5,'DisplayName','10°–20°')
plot(gB_Chl(id==3),'linecolor','red','linewidth',1.5,'DisplayName','20°–50°')
plot(gB_Chl(id==4),'linecolor','blue','linewidth',1.5,'DisplayName','50°–100°')
plot(gB_Chl(id==5),'linecolor',[1.00 0.55 0.00],'linewidth',1.5,'DisplayName','>100°')
hold off
legend
saveas(gcf,'Grain_boundaries_with_different_angles_map_05.png');

%% MOA of specific Chl-Chl boundaries 
gB_all = grains.boundary('Chlorite','Chlorite');  % Chl-Chl boundaries only
mis = gB_all.misorientation;

% misorientation angle (mean of each segment)
% angle = mean(mis.angle,2) / degree;
angle = mean(gB_Chl.misorientation.angle,2) / degree;  % mean for each boundary

% Create four ID groups
id = zeros(size(angle));
id(angle < 20) = 1;
id(angle >= 20 & angle < 50) = 2;
id(angle >= 50 & angle <= 100) = 3;
id(angle >= 100 & angle <= 150) = 4;
id(angle > 150) = 5;

% % --- 2. Split boundary sets ---
bnd1 = gB_all(id==1);
bnd2 = gB_all(id==2);
bnd3 = gB_all(id==3);
bnd4 = gB_all(id==4);
bnd5 = gB_all(id==5);

% % --- 3. Direction / hkl annotations for axis-distribution plots ---
annotHKL = Miller({1,0,0}, {0,1,0}, {0,0,1}, ...
                  {1,0,1}, {-1,0,1}, ...
                  {1,1,0}, {-1,1,0}, CS{2}, 'hkl');

% % Group 1: <20°
figure
plotAxisDistribution(bnd1.misorientation,'smooth','figSize','large');
mtexTitle('MOA — <20°');
% mtexColorMap(parula);
hold on
annotate(annotHKL,'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_less_than_20.png');

% Group 2: 20-50°
figure
plotAxisDistribution(bnd2.misorientation,'smooth','figSize','large');
mtexTitle('MOA — 20°–50°');
% mtexColorMap(parula);
hold on
annotate(annotHKL,'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_20_50.png');

% % Group 3: 50–100°
figure
plotAxisDistribution(bnd3.misorientation,'smooth','figSize','large');
mtexTitle('MOA — 50°–100°');
% mtexColorMap(parula);
hold on
annotate(annotHKL,'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_50_100.png');

% % Group 4: 100-150°
figure
plotAxisDistribution(bnd4.misorientation,'smooth','figSize','large');
mtexTitle('MOA — 100-150°');
% mtexColorMap(parula);
hold on
annotate(annotHKL,'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_100_150.png');

% % Group 5: >150°
figure
plotAxisDistribution(bnd5.misorientation,'smooth','figSize','large');
mtexTitle('MOA — >150°');
% mtexColorMap(parula);
hold on
annotate(annotHKL,'all','labeled','BackgroundColor','white');
hold off
saveas(gcf,'MOA_greater_than-150.png');

%% 1) Grain size analysis 

% Chlorite 
chl_grain_sizes = sqrt(area(grains('Chlorite'))*4/pi); 
avg_grain_size = mean(chl_grain_sizes)  % average grain size                                                           
chl_max = max(chl_grain_sizes) % Maximum grain 

% --- Plot histogram ---
chl_mean = mean(chl_grain_sizes);
chl_median = median(chl_grain_sizes);
chl_mode = mode(chl_grain_sizes);
figure
histogram(chl_grain_sizes);
xlabel('Grain size (µm)');
ylabel('Frequency (number of grains)');
title('Chlorite grain size');
stats_text = sprintf('Mean = %.2f µm\nMedian = %.2f µm\nMode = %.2f µm', chl_mean, chl_median, chl_mode);
x_limits = xlim;  % get x-axis limits
y_limits = ylim;  % get y-axis limits
x_pos = x_limits(2)*1.0;  % 70% of x-axis for padding from the right
y_pos = y_limits(2)*1.5;  % 90% of y-axis to be near the top
text(x_pos, y_pos, stats_text, 'FontSize', 12, 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
saveas(gcf,'Grain_size_vs_frequency_Chl.png');

%% Titanite
tit_grain_sizes = sqrt(area(grains('Titanite'))*4/pi); 
tit_avg_grain = mean(tit_grain_sizes)  % average grain size                                
tit_max = max(tit_grain_sizes) % Maximum grain 

figure
histogram(tit_grain_sizes);
xlabel('Grain size (µm)');
ylabel('Frequency (number of grains)');
title('Titanite grain size');
saveas(gcf,'Grain_size_vs_frequency_Titn.png');

%% 2) Aspect ratio
% Chlorite
chl_asp = grains('Chlorite').aspectRatio;
chl_mean_asp = mean(chl_asp)

% Titanite
tit_asp = grains('Titanite').aspectRatio;
tit_mean_asp = mean(tit_asp)

%% 3) PARIS FACTOR
% Chlorite
chl_paris = grains('Chlorite').paris;
chl_paris_avg = mean(chl_paris)

% Titanite
tit_paris = grains('Titanite').paris;
tit_paris_avg = mean(tit_paris)

%% 4) Shape Factor 
% Chlorite
chl_sfactor = shapeFactor(grains('Chlorite'));

% Titanite
tit_sfactor_avg = mean(shapeFactor(grains('Titanite')))

%% Avg. GOS - Chlorite and Titanire 
% Chlorite
GOS_chl = ebsd.grainMean(mis2mean.angle, grains ('Chlorite'));
GOS_deg_chl = GOS_chl*57.2958;
Mean_GOS_chl = mean(GOS_deg_chl)

% Titanite
GOS_tit = ebsd.grainMean(mis2mean.angle, grains ('Titanite'));
GOS_deg_tit = GOS_tit*57.2958;
Mean_GOS_tit = mean(GOS_deg_tit)

%% SPO -
large_grains = chl_grain_sizes > 30;
small_grains = chl_grain_sizes < 30;

numBin = 100;

figure
subplot (1,3,1)
weights = large_grains .* (grains('Chlorite').aspectRatio-1);
histogram(grains('Chlorite').longAxis,numBin, 'weights', weights,'FaceColor', 'red')
title('> 30 micron')

subplot (1,3,2)
weights = small_grains .* (grains('Chlorite').aspectRatio-1);
histogram(grains('Chlorite').longAxis,numBin, 'weights', weights,'FaceColor', 'green')
title('< 30 micron')

subplot (1,3,3)
weights = grains('Chlorite').aspectRatio-1;
histogram(grains('Chlorite').longAxis,numBin, 'weights', weights,'FaceColor', 'blue')
title('All grains')
saveas(gcf,'SPO.png');

%% INTRAGRAIN MISORIENTATION ANALYSIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% %%%% %%%% %%%% %%%% INTRAGRAIN MISORIENTATION ANALYSIS %%%% %%%% %%%% %%%% %%%%

%% Selecting grains by mouse

% % remove all not indexed pixels
% ebsd = ebsd('Chlorite');
%
% % reconstruct grains
% % [grains, ebsd.grainId] = calcGrains(ebsd,'angle',10*degree);
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',[1 10]*degree,'Boundary','tight');
%
% % smooth them
% grains = smooth(grains,5);

% plot the orientation data of the Chlorite phase
% plot(ebsd('Chlorite'),ebsd('Chlorite').orientations)
oM1.inversePoleFigureDirection=yvector; %IPFx wrt Y
plot(ebsd('Chlorite'),oM1.orientation2color(ebsd('Chlorite').orientations));


% plot the grain boundary on top of it
hold on
plot(grains.boundary,'lineWidth',0.7)
hold off

% Selecting grains by mouse
selectInteractive(grains,'lineColor','blue','lineWidth',1.6)

% STEP 1 (Select grains with your wish)
   % figure
   % global indSelected;
   % grains(indSelected)
   % hold on
   % plot(grains(indSelected).boundary,'lineWidth',2,'lineColor','gold')
   % hold off

%% Misorientation Profile for selected grain

[x,y] = ginput(2);
lineSec = [x(1:1),y(1,1);x(2,1),y(2,1)]; %enter coordinates
line(lineSec(:,1),lineSec(:,2),'linewidth',1.0,'Color','white');

hold off
global indSelected;
grain_selected = grains(indSelected);
ebsd_line = spatialProfile(ebsd(grain_selected),lineSec);

%% Angles of misorientation CHECK!
figure,
plot(sqrt(power((ebsd_line.x - ebsd_line.x(1)),2)+power((ebsd_line.y - ebsd_line.y(1)),2)),angle(ebsd_line(1).orientations,ebsd_line.orientations)./degree,'LineWidth', 2,'DisplayName', 'to reference orientation')
hold on
% plot(0.5*(ebsd_line.y(1:end-1)+ebsd_line.y(2:end)),angle(ebsd_line(1:end-1).orientations,ebsd_line(2:end).orientations)./degree,'LineWidth', 3,'DisplayName', 'orientation gradient')
hold off
xlabel('Distance in micrometer'); ylabel('Misorientation angle in degree')
legend('show','Location','northwest','fontSize', 12)

%% Smoothing 
% denoise orientation data
F = halfQuadraticFilter;
F.threshold = 1.5*degree;
F.eps = 1e-2;
F.alpha = 0.01;

%% SMOOTHING 
grain_selected=grain_selected.smooth;

%% Figures for single grain analysis

% M2M Map
% figure,
% plot(ebsd(grain_selected),ebsd(grain_selected).mis2mean.angle./degree)
% hold on
% parulaColorMap
% setColorRange([0 10]);
% plot(grain_selected.boundary,'linewidth',2)
% hold off


grodx = ebsd(grain_selected).calcGROD(grains);
figure
plot(ebsd(grain_selected),grodx.angle./degree,'micronbar','on')
parulaColorMap
setColorRange([0 20]);
hold on
plot(grain_selected.boundary,'linewidth',1.0)
hold off


%% CPO - One point per grain 
figure
% plotPDF(grain_selected.meanOrientation,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs), Miller(1,0,1,cs)], 'MarkerSize',8,'MarkerFaceColor','light green','figSize','medium')

plotPDF(grain_selected.meanOrientation,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs), Miller(1,0,1,cs)], 'MarkerSize',8,'MarkerFaceColor',[ 0.75 0.75 0],'figSize','medium')

% plotPDF(grain_selected.meanOrientation,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs), Miller(1,0,1,cs)], 'MarkerSize',6,'MarkerFaceColor','orange','figSize','medium')

% plotPDF(grain_selected.meanOrientation,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs), Miller(1,0,1,cs)], 'MarkerSize',6,'MarkerFaceColor','red','figSize','medium')

% plotPDF(grain_selected.meanOrientation,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs), Miller(1,0,1,cs)], 'MarkerSize',6,'MarkerFaceColor','dark red','figSize','medium')

%% Single grain PDF (To know the slip system)
figure
cs = ebsd(grain_selected).CS;
ori = ebsd(grain_selected).orientations;
plotPDF(ori,[Miller(1,0,0,cs),Miller(0,1,0,cs),Miller(0,0,1,cs)], 'MarkerSize',3,'property',ebsd(grain_selected).mis2mean.angle./degree,'antipodal')
set(findall(gcf,'type','text'),'visible','on');
mP.micronBar.visible = 'on'  ;
mtexColorbar
setColorRange([0 20]);
% mtexColorMap LaboTeX


% Angles of misorientation CHECK!
figure,
plot(sqrt(power((ebsd_line.x - ebsd_line.x(1)),2)+power((ebsd_line.y - ebsd_line.y(1)),2)),angle(ebsd_line(1).orientations,ebsd_line.orientations)./degree,'LineWidth', 2,'DisplayName', 'to reference orientation')
hold on
% plot(0.5*(ebsd_line.y(1:end-1)+ebsd_line.y(2:end)),angle(ebsd_line(1:end-1).orientations,ebsd_line(2:end).orientations)./degree,'LineWidth', 3,'DisplayName', 'orientation gradient')
hold off
xlabel('Distance in micrometer'); ylabel('Misorientation angle in degree')
legend('show','Location','northwest','fontSize', 12)

%% GROD
% GROD Calculations
ebsd_chl = ebsd(grain_selected);
[grains_chl,ebsd_chl.grainId] = calcGrains(ebsd_chl('indexed'));
ebsd_chl(grains_chl(grains_chl.grainSize <= 2)) = [];
[grains_chl,ebsd_chl.grainId] = calcGrains(ebsd_chl('indexed'),'threshold',[1*degree, 10*degree]);
grains_chl = smooth(grains_chl,5);
F = halfQuadraticFilter;

grod = ebsd_chl.calcGROD(grains_chl);
axSpecimen = ebsd_chl.orientations .* grod.axis('noSymmetry');

% GROD Sample reference frame (SS)
figure(5)
plot(axSpecimen,'contourf','fundamentalRegion','halfwidth',10*degree,'hemisphere','upper')%,'colorrange', [0 4])%
mtexColorMap(parula);
DataNum = numel(axSpecimen.x);
mtexTitle(sprintf('chl %d data points', DataNum))
mtexColorbar('MUD')


% GROD Crystal reference frame (CS)
figure(6)
axCrystal = grod.axis;
plot(axCrystal,'contourf','fundamentalRegion','antipodal','figSize','large')%,'colorrange', [0 4])%,'contourf'
mtexColorMap(parula);
DataNum = numel(axCrystal.h);
mtexTitle(sprintf('chl %d data points', DataNum))
mtexColorbar('MUD')

%% CPO of selcted grains

% Area Weighted
o_Sel = ebsd(grain_selected);ebsd(grain_selected).orientations;
odf_Sel = calcODF(o_Sel,'halfwidth',10*degree);
ng_Sel = length(grain_selected); %number of grains


calPF = calcPoleFigure(odf_Sel,h_chl);
maxPF = ['Max ',num2str(round(max(calPF),2)*100/100)];
minPF = ['Min ',num2str(round(min(calPF),2)*100/100)];

% Plot
figure
plotPDF(odf_Sel,h_chl,'colorrange',[0 22],'hemisphere','Upper','antipodal');
colorbar
mtexColorbar('title','m.u.d. - Area');

text([vector3d.X,vector3d.Y,vector3d.Z],{'Z','X','Y'},'BackgroundColor','w','tag','axesLabels','hemisphere','Upper');
% writing axes not by original scan definition
dim = [0.46 0.48 0.5 0.5]; % Box dimentions on top-right corner

set(findall(gcf,'type','text'),'visible','on');


n_grains = "N = " + ng_Sel; % Variable for box 2
dim2 = [0.02 0.48 0.5 0.5]; % Top-left corner
annotation('textbox',dim2,'String',n_grains,'FitBoxToText','on','fontsize',18); % Box with grains number

%% CPO of one point per grain

grains_Sel=grains('Chlorite');
odf_1ppg_Sel = calcODF(grains_Sel.mecpxrientation,'halfwidth',10*degree);
Sel = grains_Sel.mecpxrientation;
ng_Sel = length(grain_selected); %number of grains


calPF = calcPoleFigure(odf_Sel,h_chl);
maxPF = ['Max ',num2str(round(max(calPF),2)*100/100)];
minPF = ['Min ',num2str(round(min(calPF),2)*100/100)];

% Points per grain
figure
plotPDF(Sel,h_chl,'antipodal','MarkerSize',4,'micronBar','off')
text([vector3d.X,vector3d.Y,vector3d.Z],{'Z','X','Y'},'BackgroundColor','w','tag','axesLabels','hemisphere','Upper');
% writing axes not by original scan definition
dim = [0.46 0.48 0.5 0.5]; % Box dimentions on top-right corner

set(findall(gcf,'type','text'),'visible','on');
% saveas(gcf, ['PDF_PPG' fname(end-4)], 'png')


figure % Countor
plotPDF(odf_1ppg_Sel,h_chl,'direction','colorrange', [0 10],'hemisphere','Upper','antipodal');
colorbar
text([vector3d.X,vector3d.Y,vector3d.Z],{'Z','X','Y'},'BackgroundColor','w','tag','axesLabels','hemisphere','Upper');
% writing axes not by original scan definition
dim = [0.46 0.48 0.5 0.5]; % Box dimentions on top-right corner

set(findall(gcf,'type','text'),'visible','on');

n_grains = "N = " + ng_Sel; % Variable for box 2
dim2 = [0.02 0.48 0.5 0.5]; % Top-left corner
annotation('textbox',dim2,'String',n_grains,'FitBoxToText','on','fontsize',18); % Box with grains number

mtexColorbar('title','m.u.d. - Point per grain');

%% IPF for selected grains

gB_sel = grain_selected.boundary('Chlorite','Chlorite');
ind1 = angle(gB_sel.misorientation)<10*degree;
lowangleBoundary = gB_sel(ind1);

gBsub_sel=lowangleBoundary;
gBinn_sel=grain_selected.innerBoundary('Chlorite','Chlorite');

% M2M Map
figure,
plot(ebsd(grain_selected),ebsd(grain_selected).mis2mean.angle./degree,'micronBar','on')
% mtexTitle('M2M Map')
mtexColorMap WhiteJet
mtexColorbar
set(gca,'clim',[0 10])

hold on
plot(gBinn_sel.innerBoundary,'edgecolor','y','linewidth',1.0)

hold on
plot(gBsub_sel,'edgecolor','r','linewidth',1.0)

hold on
plot(grain_selected.boundary,'linewidth',2,'micronBar','off')
hold off

% IPF of Subgrain Boundaries
figure,
plotAxisDistribution(gBsub_sel.misorientation,'contourf')
title('Subgrain boundaries','fontsize',14)
mtexColorbar

% IPF of Inner Boundaries
figure,
plotAxisDistribution(gBinn_sel.misorientation,'contourf')
title('Inner boundaries','fontsize',14)
mtexColorbar

%% Subset polygon

oM = ipdfHSVOrientationMapping(ebsd('Chlorite').CS.Laue);
oM.inversePoleFigureDirection = yvector;

% Plot EBSD
figure;
plot(ebsd('Chlorite'),oM.orientation2color(ebsd('Chlorite').orientations));

% Select Polygon
poly = selectPolygon;
ebsd = ebsd(inpolygon(ebsd,poly));

%% Subset polygon for Titanite

oM = ipdfHSVOrientationMapping(ebsd('Titanite').CS.Laue);
oM.inversePoleFigureDirection = yvector;

% Plot EBSD
figure;
plot(ebsd('Titanite'),oM.orientation2color(ebsd('Titanite').orientations));

% Select Polygon
poly = selectPolygon;
ebsd = ebsd(inpolygon(ebsd,poly));

%% ****************** END *********************************