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
% %% mdg import
% mdgdir = 'C:\Users\phkr\Documents\MTEX\data\Rhenium_PushToPull\20241218_Re_7dc_tkdPTP_3\';
% mdgfile = 'dp_Re_4_ppt.ang';
% mdgpath = append(mdgdir, mdgfile);
% 
% CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');
% ebsd = EBSD.load(mdgpath,CS,'convertEuler2SpatialReferenceFrame','setting 2');
% ebsd0=ebsd;
% ebsd=ebsd.gridify;
% ebsd0=ebsd0.gridify;
%% Load data from .ang file

% filename
dir = 'data\Ti_PTP\Ti_cAxis_PTP_1\';
% file = 'Insitu Specimen 1 Site 9 Map Data 8_BW158.ang';
% file = 'Insitu Specimen 1 Site 8 Map Data 7_BW123.ang';
% file = 'Insitu Specimen 1 Site 7 Map Data 6_BW158.ang';
file = 'Insitu Specimen 1 Site 6 Map Data 5_BW123.ang';
% file = 'Insitu Specimen 1 Site 5 Map Data 4_BW123.ang';
% file = 'Insitu Specimen 1 Site 4 Map Data 3_BW123.ang';
% file = 'Insitu Specimen 1 Site 3 Map Data 2_BW123.ang';
% file = 'Insitu Specimen 1 Site 2 Map Data 1_BW123.ang';
% file = 'Insitu Specimen 1 Site 2 Map Data 1_BW263.ang';

path = fullfile(dir, file);

% define tilt, only for TKD
tkd_tilt = -11 * degree;
%20240730 & 20240624: tkd_tilt = -15 * degree;
%20241007: tkd_tilt = -15 * degree;

% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.12;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Ti
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

%Re
% CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');
%CS = crystalSymmetry.load("Rhenium.cif");

%nice to see bending and warping
%CS = crystalSymmetry('622',[2.761 2.761 4.458],'x||a','mineral','Re'); %



%Mg
%CS = crystalSymmetry('6/mmm',[1.570796370506287 1.570796370506287 2.094395102393195],'x||a','mineral','Magnesium');
%CS = crystalSymmetry.load("Mg-Magnesium.cif");



% perform corrections for tkd tilt and misaligned axes
ebsd = EBSD.load(path,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

% reference orientaiton for project2FundamentalRegion see https://github.com/mtex-toolbox/mtex/discussions/2311
indebsd = ebsd('indexed');
refebsd = indebsd(1);
reforientation = refebsd.orientations;
clear indebsd
reforientation = reforientation.project2FundamentalRegion();
reforientation = orientation.byEuler(184.137*degree,79.9861*degree,180.076*degree);



ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
% I think this is the better command
%project2FundamentalRegion(ebsd.orientations,refebsd.orientations)

ebsd0=ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;

%restrict to region [xmin ymin xmax-xmin ymax-ymin]
%region = [0.9 2.9 2.2-0.9 3.8-2.9];
%condition = inpolygon(ebsd,region);
%ebsd = ebsd(condition);
%ebsd0=ebsd;
% %% Load data from .h5oina file
% 
% % crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% % used
% dir = 'K:\EBSD-Analysis-2\Philipp EBSD analysis 2\0003_Re_cAxis_2\Experiment\Re_c_Axis_PTP_2\h5oina\';
% fname = 'Re_c_Axis_PTP_2 Specimen 2 Site 3 Map Data 6.h5oina';
% path = append(dir, fname);
% 
% 
% % crystal symmetry
% %Ti
% %CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');
% 
% %Re
% CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');
% %CS = crystalSymmetry.load("Rhenium.cif");
% 
% %Mg
% %CS = crystalSymmetry('6/mmm',[1.570796370506287 1.570796370506287 2.094395102393195],'x||a','mineral','Magnesium');
% %CS = crystalSymmetry.load("Mg-Magnesium.cif");
% 
% %if using a predifined CS
% %ebsd = loadEBSD_h5oina(fname,CS,'convertEuler2SpatialReferenceFrame');
% 
% %if loading the CS from the h5oina file
% ebsd = loadEBSD_h5oina(path,'convertEuler2SpatialReferenceFrame');
% %CS = ebsd.CS;
% 
% %perform corrections for misaligned axes
% rot = rotation.byAxisAngle(xvector,180*degree);
% ebsd = rotate(ebsd,rot,'keepXY');
% 
% ebsd.orientations = project2FundamentalRegion(ebsd.orientations);
% 
% ebsd0=ebsd;
% ebsd=ebsd.gridify;
% ebsd0=ebsd0.gridify;
% 
% %% load from .ctf file
% 
% % file location
% dir = 'C:\Users\phkr\Documents\MTEX\data\Rhenium_PushToPull\20240624_Test1\';
% file = 'Re_cAxis_ptp1_insitu Specimen 1 Site 1 Map Data 1.ctf';
% path = append(dir, file);
% 
% 
% % define tilt, only for TKD
% tkd_tilt = -11 * degree;
% 
% %imports a ebsd with the correct reference frame in the CS
% %CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');
% ctf_ebsd = loadEBSD_ctf(path,'convertEuler2SpatialReferenceFrame');
% %ctf_ebsd = rotate(ctf_ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
% ctf_ebsd = rotate(ctf_ebsd,rotation.byAxisAngle(xvector,180*degree),'keepXY');
% ctf_CS = ctf_ebsd.CS;
% CS = ctf_CS;
% ebsd = ctf_ebsd;
% 
% ebsd.orientations = project2FundamentalRegion(ebsd.orientations);
% 
% ebsd0=ebsd;
% ebsd=ebsd.gridify;
% ebsd0=ebsd0.gridify;
% 
% %restrict to region [xmin ymin xmax-xmin ymax-ymin]
% %region = [0.9 2.9 2.2-0.9 3.8-2.9];
% %condition = inpolygon(ebsd,region);
% %ebsd = ebsd(condition);
% %ebsd0=ebsd;
%% Load data: Nanomegas .ang

% % path = 'K:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240725_Re_7dc_tkdPTP\06_TEM_Observation_20240814_PTP_7degcaxis_Re\Nanomegas\';
% % file = '20240814_full_ref_res.ang';
% 
% path = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240625_Re_c_tkdPTP\06_TEM_Observation_20240726_PTP_Insitu_Re_twin_Philipp_01\Nanomegas\';
% file = '20240728_1432_Full.ang';
% 
% file = append(path, file);
% 
% % % define tilt of TEM sample
% % %default
% % flipud_bool = false;
% % zrot = 0 * degree;
% % 
% % temAlphaTilt = 0 * degree;
% % temBetaTilt = 0 * degree;
% 
% % settings for 20240728_1432_Full.ang
% flipud_bool = true;
% zrot = 0 * degree;
% 
% temAlphaTilt = 0 * degree;
% temBetaTilt = 0; % 32 * degree;
% 
% % settings for 20240728_1240_Full.ang
% % flipud_bool = true;
% % zrot = 0 * degree;
% % 
% % temAlphaTilt = 0 * degree;
% % temBetaTilt = 0 * degree;
% 
% % settings for 20240814_full_ref_res.ang
% % flipud_bool = true;
% % zrot = 90 * degree;
% % 
% % temAlphaTilt = -58 * degree;
% % temBetaTilt = 0 * degree;
% 
% 
% % define threshold for CI, 0.15 is usually safe
% thresh_ci = 0.12;
% 
% % crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% % used
% CS = crystalSymmetry('6/mmm',[2.7 2.7 4.4],'x||a','mineral','Re');
% 
% 
% %%% perform corrections for tkd tilt and misaligned axes %%%
% ebsd = EBSD.load(file,CS,'convertEuler2SpatialReferenceFrame','setting 4');
% 
% %ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
% ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
% 
% % crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% % used
% CS = crystalSymmetry('6/mmm',[2.7 2.7 4.4],'x||a','mineral','Re');
% 
% 
% %%% perform corrections for tem tilt and misaligned axes %%%
% % ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,temAlphaTilt));
% % ebsd = rotate(ebsd,rotation.byAxisAngle(yvector,temBetaTilt));
% 
% if (flipud_bool)
%     rot = rotation.byAxisAngle(xvector,180*degree);
%     ebsd = rotate(ebsd,rot);
% end
% 
% % fine adjustment (get parentOrientation from TKD scan)
% % ebsd = rotate(ebsd, inv(rotation.byAxisAngle(axis(parentOrientation_TEM, parentOrientation_TKD), angle(parentOrientation_TEM, parentOrientation_TKD))));
% 
% %use for 20240728_1432_Full.ang to be consistent with Re_cAxis_ptp1_insitu Specimen 1 Site 2 Map Data 2_BW158.ang
% reforientation = orientation.byEuler(184.137*degree,79.9861*degree,241*degree);
% ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);
% 
% ebsd0=ebsd;
% ebsd=ebsd.gridify;
% ebsd0=ebsd0.gridify;

%% some initial test plot %%%

% ebsd = rotate(rotate(ebsd0,rotation.byAxisAngle(zvector,180*degree)),rotation.byAxisAngle(yvector,180*degree));
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

saveFigure(append(file,'ebsd_ipf_XYZ.png'))

% ipf with ci greyscale
figure;
plot(ebsd,cibw,'micronbar','off')
hold on
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
title('IPF Y')
saveFigure(append(file,'ebsd_Yipf_IQ.png'))
%% cut to subset
% % Plot the EBSD map
% ipf = ipfColorKey(ebsd.CS);
% figure;
% ipf.inversePoleFigureDirection = xvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% plot(ebsd('indexed'), colors);
% title('Select a rectangular region (2 points)');
% 
% selection_confirmed = false;
% 
% while ~selection_confirmed
%     % Use ginput to get two points (diagonal corners)
%     [x, y] = ginput(2);
% 
%     % Compute the rectangle bounds
%     xmin = min(x);
%     xmax = max(x);
%     ymin = min(y);
%     ymax = max(y);
% 
%     % Define the rectangle
%     dx = xmax - xmin;
%     dy = ymax - ymin;
% 
%     % Select EBSD data within the defined rectangle
%     ind = inpolygon(ebsd, [xmin, ymin, dx, dy]);
%     selected_ebsd = ebsd(ind);
% 
%     % Plot the selected EBSD region
%     figure;
%     plot(selected_ebsd, ipf.orientation2color(selected_ebsd.orientations));
%     title('Selected EBSD Region');
% 
%     % Ask for confirmation
%     choice = questdlg('Is the selection fine?', 'Selection Confirmation', ...
%                       'Yes', 'No, reselect', 'Yes');
% 
%     if strcmp(choice, 'Yes')
%         ebsd = selected_ebsd;
%         selection_confirmed = true;
%         ebsd0=ebsd;
%         ebsd=ebsd.gridify;
%         ebsd0=ebsd0.gridify;
%         cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
%         cibw = cibw.xyz;
%     end
% end
%% perform grain reconstruction

% min angle used for grain reconstruction
grainThresh = 5*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;

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

saveFigure(append(file,'Grains_crystalShapes1.png'))

figure;
plot(grains,colors)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)
saveFigure(append(file,'grain_ID.png'))
%% Set Parent grain

% 1 determined as largest grain
 isParentGrain = grains.numPixel == max(grains.numPixel);
 parentOrientation = grains.meanOrientation(isParentGrain);
 
% 2 determine by id
isParentGrain = [10];
parentOrientation = mean(grains.meanOrientation(isParentGrain));

% 3 give a custon orientation
%parentOrientation = orientation.byEuler(183.2* degree, 92.5* degree, 143.4* degree, CS);

% quiver(grains(isParentGrain), xvector)
% quiver(grains(isParentGrain), yvector)
% %quiver(grains(isParentGrain), zvector)

% show all grains that are parent grain
parentErrorInDegrees = 20;
areParentGrains = angle(grains.meanOrientation, parentOrientation)/degree < parentErrorInDegrees;

%%% test plot %%%
figure;
plot(grains,areParentGrains)
c = mtexColorbar;
% Set the tick marks to be at 0 and 1
c.Ticks = [0, 1];
% Label the tick marks
c.TickLabels = {'Not Parent Grain', 'Is Parent Grain'};

%% plot misorientations ( see here: https://mtex-toolbox.github.io/TwinningBoundaries.html )

figure;

% plot ebsd data and misorientation of boundaries
plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
legend off
hold on
% plot shapes
%plot(grains,0.3*cS1,'facecolor','r')
% plot misorientation profile
plot(gB,gB.misorientation.angle./degree,'linewidth',5)
hold off
mtexColorbar('title','misorientation angle')
colormap jet
caxis([0 90])

saveFigure(append(file,'misorientation_profile.png'))


%%% plot histogram %%%
binEdges = 0:1:90;

figure;
histogram(gB.misorientation.angle./degree, 'BinEdges', binEdges, 'Normalization', 'probability', 'FaceColor', 'r','EdgeColor', 'w','FaceAlpha', 1)
xlabel('misorientation angle (degree)')
%axis tight
pbaspect([1.618 1 1])  % Set aspect ratio (x:y:z = 1:1:1)
exportgraphics(gcf, append(file,'histogramMisori.png'), 'Resolution', 300);

%% nice twin boundary plot
tS = twinningSystem.hexagonal_1012_1012(ebsd.CS);
ipf=ipfColorKey(ebsd.CS);
gB = grains.boundary('Ti','Ti')
twinninggB = grains.boundary(angle(grains.boundary.misorientation, tS.parentTwinMisorientation)<5*degree)
 
cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
cibw = cibw.xyz;
ebsd.prop.ciColor = cibw;

f1 = figure;
plot(ebsd,cibw,'micronbar','off')
hold on
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'),colors,'FaceAlpha',0.5)
plot(gB, 'linecolor',[0.2,0.2,0.2],'linewidth',1.5)
plot(twinninggB, 'linecolor', 'b','linewidth',1.5)
exportgraphics(f1, append(file,'ebsd_Yipf_IQ_gB.png'), 'Resolution', 300);
%% plot tilt (0) and twist (90) boundaries
% smoothGrains = smooth(grains);
% smoothgB = smoothGrains.boundary(CS.mineral,CS.mineral);
% oriGB = ebsd('id',smoothgB.ebsdId).orientations;
% axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');
% 
% figure;
% plot(smoothGrains,'FaceAlpha',0.2)
% hold on
% plot(smoothgB,angle(smoothgB.direction,axS)./degree,'linewidth',2)
% mtexColorMap hsv
% mtexColorbar
% hold off
%% plot selected crystal shapes
% figure;
% ipf=ipfColorKey(ebsd.CS);
% ipf.inversePoleFigureDirection = zvector;
% colors = ipf.orientation2color(ebsd.orientations);
% nextAxis
% plot(ebsd,colors)
% ax = gca;
% 
% cSfigure = figure;
% cS1 = crystalShape.hex(ebsd.CS);
% g_collector = [];
% while true
%     axes(ax); 
%     p = ginput(1);
%     gn = findByLocation(ebsd,p);
%     g = ebsd(gn);
%     g_collector = [g_collector;g];
%     % create and plot crystal shapes
%     figure(cSfigure)
%     plot(g.orientations * cS1, 'FaceColor',colors(gn,:), 'FaceAlpha',0.5)
%     hold on
%     disp(g.orientations.Euler)
%     if length(g_collector)>1
%         basalPlane = Miller(0,0,0,1,CS);
%         g1 = g_collector(length(g_collector)-1);
%         g2 = g_collector(length(g_collector));
%         disp("Angle (°) between basal planes between last and second to last orientation selected:")
%         disp(angle(g1.orientations*basalPlane,g2.orientations*basalPlane,'noSymmetry')/degree)
%     end
% end
%% plot selected pole figures of twinning relationships
possibleTwins = [twinningSystem.hexagonal_1012_1012(ebsd.CS), twinningSystem.hexagonal_110K(1,ebsd.CS), twinningSystem.hexagonal_1122_1124(ebsd.CS), twinningSystem.hexagonal_1121_0001(ebsd.CS)]

selectFig = figure;
plot(grains.boundary); hold on; plot(ebsd,ebsd.orientations);
hold on
title('Select multiple points and press Enter when done')
[x_selected, y_selected] = ginput;

% --- FIND NEAREST INDICES ---
n = length(x_selected);
indices = zeros(n, 1);
id = zeros(n, 1);
% indices_up1 = zeros(n, 1);

for i = 1:n
    distances = sqrt((ebsd.x - x_selected(i)).^2 + (ebsd.y - y_selected(i)).^2);
    [~, indices(i)] = min(distances, [], 'all');
    id(i) = ebsd.id(indices(i));
    text(ebsd(indices(i)), id(i));
    % Find up1 indices
    % indices_up1(i) = find(ebsdGIdT == id(i));
end

figure;
for i = 1:length(possibleTwins)
    plotPDF(ebsd(indices).orientations,ipf.orientation2color(ebsd(indices).orientations),[round(possibleTwins(i).k1),round(possibleTwins(i).eta1)],'MarkerEdgeColor','k')
    if i~=length(possibleTwins)
        nextAxis
    end
end
% 
% misori = inv(ebsd(indices(row)).orientations)*ebsd(indices(column)).orientations
% 
% for twin in possibleTwins
%     var = twin.variantDetermination2([all misoris in the table], deviationInDegrees)
%     if var
%         print twin.char in the cell
%     end
% end
% show the table
%% Grain relative orientation deviation
%cutoffForGrodAxis = 1; %DOESNT work yet

ebsd.prop.grod = ebsd.calcGROD(grains)

figure;
plot(grains.boundary); hold on
plot(ebsd, ebsd.prop.grod.angle()/degree)
colormap magma
title("Grain relative orientation deviation")
mtexColorbar('title',"Deviation in °")
saveFigure(append(file,'GROD_angle.png'))

ipf=ipfColorKey(ebsd.CS);

figure;
plot(grains.boundary); hold on
plot(ebsd,ipf.Miller2Color(ebsd.prop.grod.axis))
colormap magma
title("Grain relative orientation deviation axis. Miller to ipf color")
saveFigure(append(file,'GROD_axis.png'))


%% in-grain misorientation axis density (inverse pole figure)
% nonParentEBSD = ebsd(ebsd.grainId ~= 17 & ebsd.grainId ~= 0);
% plot(nonParentEBSD)
% IGMA = Miller.empty;
% for j = unique(nonParentEBSD.grainId)'
%     IGMA = [IGMA; axis(inv(nonParentEBSD(nonParentEBSD.grainId ==j).orientations)* grainMeanOris(j),CS)];
% end
% 
% miller_DF = calcDensity(IGMA);
%% Grain relative orientation deviation pole figure
% grainID = 10
% grainEBSD = ebsd(ebsd.grainId == grainID);
% axis(inv(ebsdInGrains(ebsdInGrains.grainId ==x).orientations)* grainMeanOris(x),CS); 
%% Misorientation
% % spot1 = 63868
% % spot2 = 61576
% 
% % spot1 = 27298
% % spot2 = 22725
% 
% % spot1 = 45845
% % spot2 = 53471
% 
% spot1 = 45845
% spot2 = 41020
% round(axis(inv(ebsd(spot1).orientations)*ebsd(spot2).orientations,CS))
% angle(ebsd(spot1).orientations,ebsd(spot2).orientations)/degree
%% Kernel average misorientation
%a bit of smoothing
F = halfQuadraticFilter;
F.alpha = 0.2;
ebsdS = smooth(ebsd,F,'fill',grains);
%ebsdS = ebsd;
%KAM Map
figure;
neighborOrder = 1; %how many neighbours are included in the calculation (higher i.e. more robust but a bit blurred
plot(ebsdS,ebsdS.KAM('threshold',2.5*degree,'order',neighborOrder) ./ degree,'micronbar','off')
mtexColorbar
title("Kernel Average Misorientation")
mtexColorMap parula %LaboTeX
saveFigure(append(file,'KAM.png'))
%% Geometrically necessary dislocations 1: Curvature, https://mtex-toolbox.github.io/GND.html
% 
% F = halfQuadraticFilter;
% ebsdSGND = smooth(ebsd('indexed'),F,'fill',grains);
% 
% compute the curvature tensor
% kappa = ebsd.curvature
% 
% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
% kappa(2,3)
% 
% newMtexFigure('nrows',3,'ncols',3);
% 
% cycle through all components of the tensor
% for i = 1:3
%   for j = 1:3
% 
%     nextAxis(i,j)
%     plot(ebsd,kappa{i,j},'micronBar','off')
%     hold on; plot(grains.boundary,'linewidth',2); hold off
% 
%   end
% end
% 
% unify the color rage  - you may also use setColoRange equal
% setColorRange equal
% setColorRange([-0.2,0.2])
% drawNow(gcm,'figSize','large')
% mtexColorbar
% saveFigure(append(file,'_Curvature.png'))

%% Geometrically necessary dislocations 2: GND calculation
% FILTER THE DATA FIRST IF WANTED
% F = halfQuadraticFilter;
% ebsdSGND = smooth(ebsd,F,'fill',grains);
% %ebsdSGND = ebsd;
% 
% %dislocationSystem
% dS = dislocationSystem.hcp(CS)
% %dSSymm = symmetrise(dS,'unique')
% [norm(dS(1).b), norm(dS(12).b),1/3 * norm(ebsd.CS.aAxis)]
% [norm(dS(13).b), norm(dS(end).b),1/3 * sqrt(norm(ebsd.CS.aAxis)^2+norm(ebsd.CS.cAxis)^2)] %dS(end) is a screw dislocation idk what the burger there is
% 
% %Assumption about energy of screw dislocations based on poisson ratio
% nu = 0.3;
% % energy of the screw dislocations arbitrarily set to 1
% dS(dS.isEdge).u = 1/(1 - 0.3);
% % energy of the screw dislocations
% dS(dS.isScrew).u = 1;
% 
% dSRot = ebsdSGND.orientations * dS;
% 
% [rho,factor] = fitDislocationSystems(kappa,dSRot);
% % the restored dislocation density tensors
% alpha = sum(dSRot.tensor .* rho,2);
% 
% % we have to set the unit manually since it is not stored in rho
% alpha.opt.unit = '1/um';
%% Geometrically necessary dislocations 3: Complete Curvature Tensor plotting
% % restore the full curvature tensor
% kappa = alpha.curvature;
% newMtexFigure('nrows',3,'ncols',3);
% 
% % cycle through all components of the tensor
% for i = 1:3
%   for j = 1:3
% 
%     nextAxis(i,j)
%     plot(ebsd,kappa{i,j},'micronBar','off')
%     hold on; plot(grains.boundary,'linewidth',2); hold off
% 
%   end
% end
% 
% setColorRange([-0.2,0.2])
% drawNow(gcm,'figSize','large');
% h = mtexColorbar
% default_position = h.Position;
% % Set height as a percentage of the default (e.g., 70% of the default height)
% percentage = 0.89;  % 70% of default height
% new_height = default_position(4) * percentage;
% % Update the colorbar position with the new height
% h.Position = [default_position(1), default_position(2), default_position(3), new_height];
% % from https://github.com/mtex-toolbox/mtex/discussions/2258
% saveFigure(append(file,'_Curvature_Full.png'))
% 


%% Geometrically necessary dislocations 4: Plot total GND density
% %set limits for logarithmic scale in m^-2
% limits = [14.5 16.8];
% 
% figure;
% % rho .* dSRot.u SHOULD BE AN ENERGY J/m^3 NOT A Dislocation density https://github.com/mtex-toolbox/mtex/discussions/921 
% ebsd.prop.totalDislocationEnergy = factor*sum(abs(rho .* dSRot.u),2);
% ebsd.prop.totalGNDdensity = sum(abs(rho),2);
% 
% plot(ebsd, log10(ebsd.prop.totalGNDdensity), 'micronbar', 'off','logarithmic');
% title('Total GND dentsiy')
% mtexColorMap parula
% saveFigure(append(file,'total_GND_density.png'))
%% Geometrically necessary dislocations 5: Plot GND density for the seperate dislocation systems
% %calculate the screw and edge dislocation density
% GND_sS = abs(rho .* dSRot.u);
% 
% %edge GNDs
% GND_E = GND_sS;
% GND_E(dSRot.isScrew) = nan;
% GND_E = factor*nansum(GND_E,2);
% 
% %screw GNDs
% GND_S = GND_sS;
% GND_S(dSRot.isEdge) = nan;
% %GND_S = factor*nansum(GND_S,2);
% 
% %basal GNDs
% basalEdge = dS.b==Miller(1,1,-2,0,CS)&dS.l==Miller(1,-1,0,0,CS);
% basalEdge = basalEdge';
% basalEdge = repmat(basalEdge, length(ebsd), 1);
% basalScrew = dS.isScrew & dS.b==Miller(1,1,-2,0,CS);
% basalScrew = basalScrew';
% basalScrew = repmat(basalScrew, length(ebsd), 1);
% 
% GND_basal = GND_sS;
% GND_basal(~(basalEdge | basalScrew)) = nan;
% GND_basal = factor*nansum(GND_basal,2);
% 
% %prismatic GNDs
% prismaticEdge = dS.b==Miller(1,1,-2,0,CS)&dS.l==Miller(0,0,0,1,CS);
% prismaticEdge = prismaticEdge';
% prismaticEdge = repmat(prismaticEdge, length(ebsd), 1);
% 
% GND_prismatic = GND_sS;
% GND_prismatic(~(prismaticEdge | basalScrew)) = nan;
% GND_prismatic = factor*nansum(GND_prismatic,2);
% 
% %pyramidal 1 and 2 GNDs
% GND_pyramidal = GND_sS;
% GND_pyramidal((basalEdge | basalScrew | prismaticEdge)) = nan;
% GND_pyramidal = factor*nansum(GND_pyramidal,2);
% 
% 
% plot(ebsd, log10(ebsd.prop.totalGNDdensity), 'micronbar', 'off','logarithmic');
% title('Total GND dentsiy')
% caxis(limits);
% %nextAxis
% % plot(ebsd, log10(GND_S), 'micronbar', 'off','logarithmic');
% % title('Screw GND dentsiy')
% % caxis(limits);
% % nextAxis
% % plot(ebsd, log10(GND_E), 'micronbar', 'off','logarithmic');
% % title('Edge GND dentsiy')
% % caxis(limits);
% nextAxis
% plot(ebsd, log10(GND_basal), 'micronbar', 'off','logarithmic');
% title('Basal GND dentsiy')
% caxis(limits);
% 
% nextAxis
% plot(ebsd, log10(GND_prismatic), 'micronbar', 'off','logarithmic');
% title('Prismatic GND dentsiy')
% caxis(limits);
% 
% nextAxis
% plot(ebsd, log10(GND_pyramidal), 'micronbar', 'off','logarithmic');
% title('Pyramidal 1 and 2 GND dentsiy')
% caxis(limits);
% 
% 
% 
% % set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% % set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
% h = mtexColorbar('title', append('log_{10}', 'rho_{GND}', '\', 'm^{-2}'));
% default_position = h.Position;
% % Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)
% percentage = 0.89;  % 70% of default height
% new_height = default_position(4) * percentage;
% % Update the colorbar position with the new height
% h.Position = [default_position(1), default_position(2), default_position(3), new_height];
% mtexColorMap parula
% % saveFigure(append(file,'seperated_GND_density.png'))
%% get angle between c-Axis of dominant orientataion to the y direction (axis of tension)
odf = calcDensity(ebsd.orientations);
figure;
plotPDF(odf, Miller({0,0,0,1}, CS),'grid','grid_res',15*degree)
nextAxis
plotPDF(odf, Miller({1,-1,0,0}, CS),'grid','grid_res',15*degree)
nextAxis
plotPDF(odf, Miller({1,1,-2,0}, CS),'grid','grid_res',15*degree)
[~,maxOri] = max(odf);
angle_caxis_yDir = angle (yvector, maxOri * Miller({0,0,0,1}, CS))/degree
saveFigure(append(file,'ODF.png'))

%% Quaternion kmeans clustering with initial centers
% cannot define my own distance function easily quatAngle
% 
% % Extract orientations as quaternions
% quats = quaternion(ebsd.orientations);
% quats = quats.reshape([],1);
% n = quats.norm;
% %check if unit quat or NaN
% assert(all(n(~isnan(n))-1.0<1e-6));
% quats = [quats.a;quats.b;quats.c;quats.d]';
% 
% % Define expected orientations and convert to quaternions
% expectedQuats = quaternion(DTsWant); % Modify as needed
% n = expectedQuats.norm;
% assert(all(n(~isnan(n))-1.0<1e-6));
% expectedQuats = [expectedQuats.a;expectedQuats.b;expectedQuats.c;expectedQuats.d]';
% 
% quatAngle = @(q1, q2) acos(2 * dot(q1, q2)^2 - 1);
% 
% % Perform k-means clustering using quaternion distance
% [idx, clusterCenters] = kmeans(quats, length(expectedQuats), ...
%     'Start', expectedQuats, ...
%     'Distance', quatAngle);
% 
% % Convert cluster centers back to orientations
% clusterOris = orientation(clusterCenters);
% 
% % Plot results
% plot(ebsd, idx) % Color grains by cluster
% hold on
% plot(clusterOris, 'markerSize', 20, 'linewidth', 2) % Show cluster centers
% hold off
% 

%% Cluster orientations and plot them in Sigma section pole figures
%Definitions
tS = tS_1121_Re;
clusterNum = 0; % amount of clusters, 0 if unconstrained
clusterHalfWidth = 10*degree; %size of a cluster


odf = calcDensity(grains.meanOrientation);
%odf = calcDensity(ebsd.orientations);

%ODF plot
figure;
plot(odf,'sections',6,'silent','sigma','grid','grid_res',15*degree);
saveFigure(append(file,'orientationDensity_Sigma.png'));


%plotPDF(odf,directionODF,'antipodal','grid'); %'logarithmic'
%plotPDF(ebsd.orientations, Miller({0,0,0,1},CS),'grid')
%plotPDF(grains.meanOrientation, Miller({0,0,0,1},CS),'grid', 'MarkerColor','k', 'MarkerSize',2)
%sigma Sections

%---------
%Calculate the expected double twin orientations
%mark out the expected orientations for the 1121 double twins starting from
%the matrix orientation
%11-21 -> 11-21
tS_1121_Re_variants = tS.symmetrise;
PT = tS_1121_Re_variants.parentTwinMisorientation;

uniqueDT = unique(tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation);
doubleTwinVars = tS_1121_Re_variants(1).parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation;

%uniqueDT 1 and 3 are the DT that is close to another primary twin and matrix respectively
dontWant = [uniqueDT(1),uniqueDT(3)]
DTsWant = doubleTwinVars(~any(angle(dontWant,doubleTwinVars)<1e-5,1));
%Doest cat easily
for i = 1:length(PT)
    DTsWant = [DTsWant, PT(i)];
end
unDTsWant= unique(DTsWant)
expectedVarCol = parulaColorMap(length(unDTsWant))
%--------

%0001 pole figure
pf = figure;
plotPDF(grains.meanOrientation, Miller({0,0,0,1},CS),'grid', 'MarkerColor','k', 'MarkerSize',4);
hold on
%Sigma sections
sf = figure;
oS = sigmaSections(odf.CS,odf.SS);
plotSection(grains.meanOrientation,oS,'grid', 'MarkerColor','k', 'MarkerSize',4) 
hold on



%[cId,center] = calcCluster(grains.meanOrientation,'numCluster',6);
%cluster the grains according to their orientation
if clusterNum == 0
    [cId,center] = calcCluster(grains.meanOrientation,'halfwidth',clusterHalfWidth);
else
    [cId,center] = calcCluster(grains.meanOrientation,'halfwidth',clusterHalfWidth,'numCluster',clusterNum);
end
colmap = hsvColorMap(length(center));
%plotSection(center,oS,'grid', 'MarkerColor','k','Marker','x', 'MarkerSize',6)
plotSection(parentOrientation,oS,'grid', 'MarkerColor','k','Marker','x', 'MarkerSize',10)

figure(pf.Number); %plots fig 2 in parallel
plotPDF(parentOrientation, Miller({0,0,0,1},CS),'grid', 'MarkerColor','k','Marker','x', 'MarkerSize',10);
figure(sf.Number);

hold on
%plot the expected primary twin and double locations
for i = 1:length(unDTsWant)
    hold on    
    plotSection(parentOrientation.*DTsWant(angle(DTsWant, unDTsWant(i))<1e-5),oS,'grid', 'Marker','diamond', 'MarkerColor',expectedVarCol(i,:), 'MarkerSize',10) 
    figure(pf.Number); %plots fig 2 in parallel
    hold on
    plotPDF(parentOrientation.*DTsWant(angle(DTsWant, unDTsWant(i))<1e-5), Miller({0,0,0,1},CS),'grid', 'Marker','diamond', 'MarkerColor',expectedVarCol(i,:), 'MarkerSize',10);
    figure(sf.Number);

end
hold on
%plot the clusters
for i = 1:length(center)
    plotSection(grains.meanOrientation(cId==i),oS,'grid', 'MarkerColor',colmap(i,:), 'MarkerSize',4)
    figure(pf.Number); %plots fig 2 in parallel
    hold on
    plotPDF(grains.meanOrientation(cId==i), Miller({0,0,0,1},CS),'grid', 'MarkerColor',colmap(i,:), 'MarkerSize',4);
    figure(sf.Number);
end
figure(sf.Number);
saveFigure(append(file,'Cluster_Sigma.png'));
figure(pf.Number); 
saveFigure(append(file,'Cluster_0001_pole.png'));

figure;
% plot the grains according to their color in the sigma section
plot(grains,'faceAlpha',0)
hold on
for i = 1:length(center)
    plot(grains(cId==i), 'FaceColor',colmap(i,:));
    hold on
end
saveFigure(append(file,'Cluster_Map.png'));
%plot(ebsd,colors)
%mtexColorbar

% %plot twinning variants in ODF
% hold on
% plotVariants(tS_1.parentTwinMisorientation,parentOrientation,char(tS_1));
% hold on
% plotVariants(tS_2.parentTwinMisorientation,parentOrientation,char(tS_2));
% hold on
% plotPDF(parentOrientation, directionODF,'label','SX-orientation')
% %saveFigure(append(file,'texture2.png')

% Define symbols and labels
legendLabels = {'ParentOrientation', 'Theoretical twin variants, primary twins and two types of secondary twins(!!not all possibilities!!)', 'Orientations', [ 'Number of calculated Clusters: ',num2str(length(center))] };
colors = {'k', expectedVarCol(1,:) , colmap(1,:), colmap(1,:)};  % Colors for markers
markerType = {'x', 'd', 'o','o'};  % Different marker styles

% Create a new figure for the legend
figure;
hold on
axis([1 2 0 1]); % Fix x-axis starting at 1, set reasonable y-axis range
axis off; % Hide axes for a clean legend-only plot

xPos = 1; % Fixed X position for markers
yPos = linspace(0.8, 0.2, length(legendLabels)); % Y positions for each 

% Loop through each legend entry
for i = 1:length(legendLabels)
    % Plot marker for legend
    plot(xPos, yPos(i), markerType{i}, 'Color', colors{i}, 'MarkerSize', 10, 'LineWidth', 2);
    
    % Add text next to the marker
    text(xPos + 0.1, yPos(i), legendLabels{i}, 'FontSize', 12, 'Interpreter', 'none');
end
saveFigure(append(file,'Cluster_Legend.png'));
hold off;

% plot traces of the rotation Axis
%plot(parentOrientation*tS_1121_Re_variants(5).rotAxis,'plane','linecolor','r','linewidth',2,'add2all')
%plot((parentOrientation*tS_1121_Re_variants(5).parentTwinMisorientation)*tS_1121_Re_variants(4).rotAxis,'plane','linecolor','b','linewidth',2,'add2all')
%% cluster by using misorientation to parentOrientation

% angular threshold for considering results (miso of grain to expected
% orientations of twins


%calculate the primary twin and secondary twin variants
tS = tS_1121_Re;
tS_1121_Re_variants = tS.symmetrise;
PT = tS_1121_Re_variants.parentTwinMisorientation;


uniqueDT = unique(PT*PT);
%uniqueDT 1 and 3 are the DT that is close to another primary twin and matrix respectively
uDTs_dontWant = [uniqueDT(1),uniqueDT(3)];
uDTs_Want = [uniqueDT(2),uniqueDT(4)];

% Color setup
parentColor = [0.5 0.5 0.5];
baseColors = lines(6);  % Distinct primary colors for T1–T6
nonCozoneHighlight = [1, 0.2, 1];


% Initialize arrays
%allGuessOrientations = orientation.empty;
allLabels = strings(0,1);
twinGenerationClass = [];
LabelColor = [];

% Parent (matrix) orientation
allGuessOrientations = orientation.byEuler(0,0,0,CS,CS);
allLabels(end+1) = "P";
twinGenerationClass(end+1) = 0;
LabelColor(end+1, :) = parentColor;

% ---- Primary twins: T1 through T6 ----
for i = 1:length(tS_1121_Re_variants)
    % Generate or import the orientation for twin i
    o = tS_1121_Re_variants(i).parentTwinMisorientation;
    allGuessOrientations(end+1) = o;
    twinGenerationClass(end+1) = 1;
    allLabels(end+1) = "T" + i;
    LabelColor(end+1, :) = baseColors(i, :);
end

% ---- Secondary Twins T11–T66 ----
for i = 1:length(PT)
    for j = 1:length(PT)
        mo = PT(i) * PT(j);
        [angleToWant, matchIdx] = min(angle(uDTs_Want, mo));

        if angleToWant < 1e-5
            allGuessOrientations(end+1) = mo;
            allLabels(end+1) = "T" + i + j;
            twinGenerationClass(end+1) = 1 + matchIdx;

            if matchIdx == 2
                % Darker shade of Ti
                darkColor = 0.5 * baseColors(i, :);
                LabelColor(end+1, :) = darkColor;
            else
                % Bright red
                LabelColor(end+1, :) = nonCozoneHighlight;
            end
        end
    end
end



% ---- Create final table ----
guessTable = table(allGuessOrientations(:), allLabels(:), twinGenerationClass(:), LabelColor(:,:), ...
                   'VariableNames', {'Misorientation', 'Label', 'GenerationClass', 'LabelColor'});

% ---- assign to the table ----
A = angle_outer(grains.meanOrientation, parentOrientation* guessTable.Misorientation);
[minAngle, idx] = min(A, [], 2);
matchedLabels = guessTable.Label(idx);
matchedTwinGeneration = guessTable.GenerationClass(idx);
matchedColors = guessTable.LabelColor(idx, :);

resultTable = table(grains.id(:), grains.meanOrientation(:), matchedLabels, matchedColors , matchedTwinGeneration, minAngle, ...
                    'VariableNames', {'GrainId','GrainOrientation', 'Label', 'LabelColors', 'GenerationClass', 'Angle'});

% ---- FILTER: exclude poor matches (angle > 8.5°) ----
validMask = resultTable.Angle < angleThreshold;
filteredTable = resultTable(validMask, :);


answer = input('Do you want to save the results to a MAT file? (y/n): ', 's');
if lower(answer) == 'y'
    save(append(file,'TwinClassification.mat'), 'ebsd', 'grains', 'resultTable', 'guessTable', 'file', 'parentOrientation');
end

%% ---- plot Results of misorientation clustering ----
angleThreshold = 17 * degree;




CS = grains.CS;
% figure; plotPDF(resultTable.GrainOrientation, resultTable.GenerationClass, Miller(0,0,0,1,CS),'lower')
% figure; plotPDF(resultTable.GrainOrientation, resultTable.LabelColors, Miller(0,0,0,1,CS),'lower')
% exportgraphics(gcf, append(file,'grouping_map.png'), 'Resolution', 300);
% figure; plot(grains(resultTable.GrainId), resultTable.GenerationClass)
figure; plot(grains(resultTable.GrainId), resultTable.LabelColors); hold on
plot(grains.boundary, 'lineWidth',1.5);
exportgraphics(gcf, append(file,'unfiltered_grouping_PDF.png'), 'Resolution', 300);

% ---- plot filtered Results ----
% figure with all results (not filtered for the angle threshold
% figure; plotPDF(resultTable.GrainOrientation,Miller(0,0,0,1,CS), 'MarkerColor', 'k', 'MarkerSize', 1.5);hold on; 
% exportgraphics(gcf, append(file,'unfiltered_grouping_PDF.png'), 'Resolution', 300);
% plotPDF(filteredTable.GrainOrientation, filteredTable.GenerationClass, Miller(0,0,0,1,CS))

% figure with filtered results
figure; plotPDF(resultTable.GrainOrientation,Miller(0,0,0,1,CS), 'MarkerColor', 'k', 'MarkerSize', 5); hold on;
plotPDF(filteredTable.GrainOrientation, filteredTable.LabelColors, Miller(0,0,0,1,CS), 'MarkerSize', 7)
% doesn't work -> annotate(filteredTable.GrainOrientation, 'label',filteredTable.Label)
setMTEXpref('annotationStyle', {'marker'  'o'  'MarkerEdgeColor'  'w'  'MarkerFaceColor'  'k' 'MarkerSize' 1 'hitTest'  'off'})
[uniqueLabels, ia] = unique(filteredTable.Label, 'first');
for i = 1: length(ia)
    annotate(filteredTable.GrainOrientation(ia(i)), 'label',filteredTable.Label(ia(i)));
end
exportgraphics(gcf, append(file,'filtered_grouping_PDF.png'), 'Resolution', 300);

% figure; plot(grains(filteredTable.GrainId), filteredTable.GenerationClass)
figure; plot(grains(filteredTable.GrainId), filteredTable.LabelColors); hold on;
plot(grains.boundary, 'lineWidth',1.5)
exportgraphics(gcf, append(file,'filtered_grouping_map.png'), 'Resolution', 300);

% figure with the non-cozonal twins labled
onlyDTTable = filteredTable(filteredTable.GenerationClass == 2,:);
figure; plot(grains(filteredTable.GrainId), filteredTable.LabelColors); hold on
text(grains(onlyDTTable.GrainId), onlyDTTable.Label)
exportgraphics(gcf, append(file,'filtered_grouping_map_labelled.png'), 'Resolution', 300);

% pole figure with crystal shapes
figure; plotPDF(parentOrientation*(guessTable.Misorientation),Miller({0,0,0,1},ebsd.CS),'MarkerColor','w','MarkerSize', 0.1,'grid');
hold on
plotPDF(grains.meanOrientation, Miller({0,0,0,1},ebsd.CS),'MarkerColor','k','MarkerSize',3)
hold on
plotPDF(filteredTable.GrainOrientation, filteredTable.LabelColors, Miller(0,0,0,1,CS),'MarkerSize', 3)

% plotCs = (filteredTable.GrainOrientation)* cS_Rhenium * 0.2;
plotCs = 0.2*rotate(cS_Rhenium,filteredTable.GrainOrientation);
for i = 1:length(filteredTable.GrainOrientation)
    plot(filteredTable.GrainOrientation(i), plotCs(i),'FaceColor',filteredTable.LabelColors(i,:), 'add2all')
end

plottedCs = (parentOrientation*(guessTable.Misorientation)) * cS_Rhenium * 0.4;
for i = 1:length(guessTable.Misorientation)
    hold on
    plot(parentOrientation*guessTable.Misorientation(i),plottedCs(i),'FaceColor',guessTable.LabelColor(i,:),'FaceAlpha',1,'add2all')
end
exportgraphics(gcf, append(file,'filtered_grouping_cS.png'), 'Resolution', 300);


%% select crystal shapes interactively
% !!!!! CHANGE   %scaling = sqrt(grains.area); TO scaling = 1; in edit
% grains/plot line 111 crystalShape case

figSelect = figure;
plot(grains(resultTable.GrainId), resultTable.LabelColors);
%hold on; plot(grains(resultTable.GrainId), cS_Rhenium);

while true
    global indSelected;
    indSelected = [];

    fprintf('Select a grain (click in the figure, or close the figure to stop)...\n');
    selectInteractive(grains, 'lineColor', 'gold');

    % Wait for a selection
    while isempty(indSelected)
        if ~isvalid(figSelect)
            disp('Figure closed. Exiting selection loop.');
            break;
        end
        pause(0.1);
    end

    if isempty(indSelected) || ~isvalid(figSelect)
        break;  % Exit loop if no selection or figure closed
    end

    selectedGrain = grains(indSelected);
    fig1 = figure;
    plot(selectedGrain, 'lineColor', 'w', 'faceColor', 'w', 'micronbar', 'off'); hold on;

    grainColor = resultTable.LabelColors(resultTable.GrainId == indSelected, :);
    plot(selectedGrain, cS_Rhenium, 'faceColor', grainColor, 'faceAlpha', 0.9, 'lineWidth', 3);
    
    exportgraphics(fig1, append(file, 'cS_Grain', num2str(indSelected), '.png'), 'Resolution', 300);
    close(fig1);
end

close(figSelect);

%% Primary twin system detection
%Definitions
errorInDegrees = 10; %allowed misorientation from the theoretical twin orientation
%
Systems_112K = arrayfun(@(x_int) twinningSystem.hexagonal_211K(x_int, CS), [1,2]);
Systems_110K = arrayfun(@(x_int) twinningSystem.hexagonal_110K(x_int, CS), [1,2]);

twinningSystems = [Systems_112K, Systems_110K];

%parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
%parentOrientation=grains(9).meanOrientation;
misOrientation = inv(grains.meanOrientation).*parentOrientation;

twinningsystemdeviation_logical = angle(misOrientation,twinningSystems.parentTwinMisorientation())/degree < errorInDegrees;

figure;
plot(ebsd, cibw)
hold on
plot(grains,'FaceAlpha',0)
hold on

colmap = hsvColorMap(length(twinningSystems));

for i = 1:length(twinningSystems)
    plot(grains(twinningsystemdeviation_logical(:,i)),'faceColor','r')%, 'faceColor', colmap(i,:))
    hold on
end
legend off
mtexColorMap(colmap)
setColorRange([1 length(twinningSystems)]);
tickDistance = (length(twinningSystems)-1)/length(twinningSystems);
mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems))
myColorbar = mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems),'TickLabels', cellfun(@(x) strjoin(x), num2cell(string(twinningSystems.n.hkil),2)), 'TickLength', 0);
saveFigure(append(file,'TwinningSystems.png'))
hold off
% 
% % misorientation from the ideal twin systems (within the given
% % errorInDegrees)
% % Compute the mean orientation for each grain
% meanOrientations = grains(twinningsystemdeviation_logical(:,2)).meanOrientation;
% 
% % Initialize array for misorientations
% ebsdMisorientations = zeros(size(ebsd.orientations));
% tsMisorientations = twinningSystems.parentTwinMisorientation;
% 
% figure;
% plot(grains,'FaceAlpha',0)
% hold on
% 
% % Compute misorientations
% for j = 1:length(tsMisorientations)
%     grainIdInTwinOrientation = grains(twinningsystemdeviation_logical(:,j)).id;
%     for i = 1:length(grainIdInTwinOrientation)
%         plot(ebsd(ebsd.grainId==grainIdInTwinOrientation(i)), angle(inv(ebsd(ebsd.grainId==grainIdInTwinOrientation(i)).orientations).*parentOrientation, tsMisorientations(j))/degree)
%         hold on
%     end
% end
% legend off
% % Plot misorientation as a heatmap
% mtexColorbar('title','Misorientation from expected Twin Orientation in (degrees)');
% caxis([1 10]); % Adjust color axis as needed
% xlabel('X (\mum)');
% ylabel('Y (\mum)');
% saveFigure(append(file,'misorientationOffGivenTwinningSystem.png'))
%% 11-21 retwin variants
% %11-21 -> 11-21
% tS = tS_1121_Re;
% % initial grains, of the untwinned microstrucutre
% %parentGrains = grains(areParentGrains);
% 
% % twin system which is to be analysed
% tS = twinningSystem.hexagonal_1121_0001(CS);
% 
% %allowed misorientation from the theoretical twin orientation in degrees
% errorInDegrees = 10;
% 
% tS_1121_Re_variants = tS_1121_Re.symmetrise;
% %angle(tS_1121_Re.parentTwinMisorientation)/degree
% %angle(tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation)/degree
% retwinning_Variants_misorientation = unique(tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation);
% 
% disp(append('There are ', string(length(retwinning_Variants_misorientation)),' geometrically unique double twin types. Of which ', num2str(any(angle(retwinning_Variants_misorientation) < 1e-5)), ' is/are the back-transformation to the parent orientation.'))
% retwinning_Variants_misorientation
% 
% misOrientation = inv(grains.meanOrientation).*parentOrientation;
% 
% colmapsize = prod(size(retwinning_Variants_misorientation));
% colmap = viridisColorMap(colmapsize);
% 
% figure;
% plot(grains,'FaceAlpha',0)
% hold on
% %twinningsystemdeviation_logical = angle(misOrientation,twinningSystems.parentTwinMisorientation())/degree < errorInDegrees;
% 
% for i = 1:length(retwinning_Variants_misorientation)
%     twinningsystemdeviation_logical= angle(misOrientation, retwinning_Variants_misorientation(i))/degree < errorInDegrees; %, 'noSymmetry'
%     plot(grains(twinningsystemdeviation_logical), 'faceColor', colmap(i,:))
%     text(grains(twinningsystemdeviation_logical),append(char(i)))
%     pause(1)
%     % find(squeeze(twinningsystemdeviation_logical(i,j,:)))
%     % plot(grains(squeeze(twinningsystemdeviation_logical(i,j,:))), 'faceColor', colmap(i,j,:))
%     hold on
% end
% legend off
% mtexColorMap(colmap)
% setColorRange([1 colmapsize]);
% myColorbar = mtexColorbar('title',append('Double twins of the ', char(tS),'-system'));
%% Twin system detection by GB misorientation
% errorInDegrees = 10; %allowed misorientation from the theoretical twin orientation
% %
% Systems_112K = arrayfun(@(x_int) twinningSystem.hexagonal_211K(x_int, CS), 1:2);
% %Systems_110K = arrayfun(@(x_int) twinningSystem.hexagonal_110K(x_int, CS), 2:2);
% twinningSystems = [Systems_112K];% Systems_110K];
% 
% twinningsystemdeviation = angle(gB.misorientation,twinningSystems.parentTwinMisorientation())/degree;
% twinningsystemdeviation_logical = twinningsystemdeviation < errorInDegrees;
% 
% figure;
% plot(ebsd, cibw)
% hold on
% plot(grains,'FaceAlpha',0.5)
% hold on
% 
% colmap = hsvColorMap(length(twinningSystems));
% 
% for i = 1:length(twinningSystems)
%     plot(gB(twinningsystemdeviation_logical(:,i)), 'linecolor',colmap(i,:),'linewidth',2.5) %'faceColor',
%     hold on
% end
% legend off
% mtexColorMap(colmap)
% setColorRange([1 length(twinningSystems)]);
% tickDistance = (length(twinningSystems)-1)/length(twinningSystems);
% %mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems))
% myColorbar = mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems),'TickLabels', cellfun(@(x) strjoin(x), num2cell(string(twinningSystems.n.hkil),2)), 'TickLength', 0);
% title('Twinning System')
% saveFigure(append(file,'TwinningSystems_alongGBs.png'))
% hold off
% 
% % misorientation from the ideal twin systems (within the given
% % errorInDegrees)
% % Compute the mean orientation for each grain
% deviationValues = twinningsystemdeviation_logical .* twinningsystemdeviation;
% 
% reducedDeviationValues = zeros(size(deviationValues, 1), 1);                    % Initialize the output array
% for i = 1:size(deviationValues, 1)
%     idx = find(deviationValues(i, :), 1);                  % Find the first nonzero element in the row
%     if ~isempty(idx)
%         reducedDeviationValues(i) = deviationValues(i, idx);                   % Assign the first nonzero value to B
%     else
%         reducedDeviationValues(i) = NaN;
%     end
% end
% 
% figure;
% % plot ebsd data and misorientation of boundaries
% plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
% legend off
% hold on
% % plot shapes
% plot(grains,0.3*cS1,'facecolor','r')
% % plot misorientation profile
% plot(gB,reducedDeviationValues,'linewidth',5)
% hold off
% mtexColorbar('title','misorientation angle (°)')
% title('Misorientation from the Twinning System misorientation')
% colormap jet
% caxis([0 10])
% 
% saveFigure(append(file,'TwinningSystems_alongGBs_misorientation.png'))
%% 11-22 and 11-21 twinning on ipf map
% errorInDegrees = 10; %allowed misorientation from the theoretical twin orientation
% %
% Systems_112K = arrayfun(@(x_int) twinningSystem.hexagonal_211K(x_int, CS), 1:2);
% twinningSystems = [Systems_112K];
% 
% twinningsystemdeviation = angle(gB.misorientation,twinningSystems.parentTwinMisorientation())/degree;
% twinningsystemdeviation_logical = twinningsystemdeviation < errorInDegrees;
% 
% ipf=ipfColorKey(ebsd.CS);
% 
% figure;
% plot(ebsd, cibw)
% hold on
% ipf.inversePoleFigureDirection=yvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% plot(ebsd('indexed'),colors,'FaceAlpha',0.8)
% hold on
% 
% colmap = gray(3);
% colmap = colmap(1:2,:)
% 
% for i = 1:length(twinningSystems)
%     plot(gB(twinningsystemdeviation_logical(:,i)), 'linecolor',colmap(i,:),'linewidth',2.5) %'faceColor',
%     hold on
% end
% legend off
% mtexColorMap(colmap)
% setColorRange([1 length(twinningSystems)]);
% tickDistance = (length(twinningSystems)-1)/length(twinningSystems);
% %mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems))
% myColorbar = mtexColorbar('title','Twin systems','Ticks', 1+tickDistance/2:tickDistance:length(twinningSystems),'TickLabels', cellfun(@(x) strjoin(x), num2cell(string(twinningSystems.n.hkil),2)), 'TickLength', 0);
% title('Twinning System')
% saveFigure(append(file,'11-22and11-21_twinningSystems.png'))
% hold off
% 
% % misorientation from the ideal twin systems (within the given
% % errorInDegrees)
% % Compute the mean orientation for each grain
% deviationValues = twinningsystemdeviation_logical .* twinningsystemdeviation;
% 
% reducedDeviationValues = zeros(size(deviationValues, 1), 1);                    % Initialize the output array
% for i = 1:size(deviationValues, 1)
%     idx = find(deviationValues(i, :), 1);                  % Find the first nonzero element in the row
%     if ~isempty(idx)
%         reducedDeviationValues(i) = deviationValues(i, idx);                   % Assign the first nonzero value to B
%     else
%         reducedDeviationValues(i) = NaN;
%     end
% end
% 
% figure;
% % plot ebsd data and misorientation of boundaries
% plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
% legend off
% hold on
% % plot shapes
% plot(grains,0.3*cS1,'facecolor','r')
% % plot misorientation profile
% plot(gB,reducedDeviationValues,'linewidth',5)
% hold off
% mtexColorbar('title','misorientation angle (°)')
% title('Misorientation from the Twinning System misorientation')
% colormap jet
% caxis([0 10])
% 
% saveFigure(append(file,'11-22and11-21_twinningSystems_misorientation.png'))
%% Retwinning variant map along GBs

% initial grains, of the untwinned microstrucutre
parentGrains = grains(areParentGrains);

% twin system which is to be analysed
tS = twinningSystem.hexagonal_1121_0001(CS);
tSvariants = tS.symmetrise;

% iteration depth, i_depth=1 only primary twins are detected
i_depth = 2;

%allowed misorientation from the theoretical twin orientation in degrees
errorInDegrees = 8;

%grainIDs
grainIds = grains.id;

%parent grain logical
parentGrainLogical = false([length(grains),1]);
parentGrainLogical(parentGrains.id) = true;

% iterator
currentGrainsLogical = parentGrainLogical;

% record analysed grains to avoid recursion
prevGrainsLogical = parentGrainLogical;

% prepare blank grain plot
varBoundaries = figure;
% plot(grains,grains.meanOrientation)
% hold on
% plot boundaries
plot(ebsd, cibw,'micronbar','off')
hold on
plot(grains.boundary)

varDevFromTwin = figure;
% plot(grains,grains.meanOrientation)
% hold on
% plot boundaries
plot(ebsd, cibw,'micronbar','off')
hold on
plot(grains.boundary)

varHabitPlanes = figure;
plot(ebsd, cibw,'micronbar','off')
hold on
plot(grains.boundary)

colmap = hsvColorMap(tS.numberOfVariants);

% % prepare crystal directions plotting
% % dotted lines are positive z (into plane with the right plotting convention)
% % and solid lines are negative z
neg_posMarker = {'-',':'};
% % [1,1,-2,0], [1,-1,0,0] and [0,0,0,1] directions get plotted
% grainDirs =  arrayfun(@(x) {x, rotate(Miller(1,1,-2,0,CS,'UTVW'), x.meanOrientation), rotate(Miller(1,-1,0,0,CS,'UTVW'), x.meanOrientation),rotate(Miller(0,0,0,1,CS,'UTVW'), x.meanOrientation)},grains(currentGrainsLogical),'UniformOutput',false);
% grainDirs = vertcat(grainDirs{:});
% 
% % plot crystal directions of parents
% for row = 1:size(grainDirs,1)
%     thisgrain = grainDirs{row,1};
%     vec_a1 = grainDirs{row,2};
%     vec_a2 = grainDirs{row,3};
%     vec_c = grainDirs{row,4};
%     %both 'linecolor' and 'color' do not work, you have to edit the
%     %@grain2d/quiver.m method:
%     % copy line:  h = optiondraw(quiver(xy(:,1),xy(:,2),dir.x,dir.y),varargin{:});
%     % 8 lines up and comment out varargin = delete_option(varargin,'color',1);
%     quiver(thisgrain, vec_a1,'color','r','LineStyle',neg_posMarker{ispos(vec_a1.z)+1},'noScaling')%'AutoScaleFactor',0.5)
%     quiver(thisgrain, vec_a2,'color','g','LineStyle',neg_posMarker{ispos(vec_a2.z)+1},'noScaling')%'AutoScaleFactor',0.5)
%     quiver(thisgrain, vec_c,'color','b','LineStyle',neg_posMarker{ispos(vec_c.z)+1},'noScaling')%'AutoScaleFactor',0.5)
% end
% 
% % plot crystal shapes of parents
% cS1 = crystalShape.hex(ebsd.CS);
% hold on
% hGrains = plot(grains(currentGrainsLogical),0.5*cS1,'faceAlpha', 0.1);
% hText = text(grains(currentGrainsLogical), 'parent', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');


for i = 1:i_depth

    currentgBSegments = gB(gB.I_FG(:,currentGrainsLogical));
    %currentgBSegments = gB(grains(currentGrainsLogical));
    currentGrainIds = grainIds(currentGrainsLogical);

    % CHUGG OUT PREV GRAINS
    possibleChildGrainIds = grainIds(~(parentGrainLogical | prevGrainsLogical));

    % Create a grid of all combinations
    [X, Y] = ndgrid(currentGrainIds,possibleChildGrainIds);

    % Reshape into a two-column array with all combinations of (x, y)
    currentToChildCombinations = [X(:), Y(:)];
    
    %selectByGrainId alway yields the gB as parent -> twin (X to Y) so no
    %inv misorientation
    currentgBSegments = gB.selectByGrainId(currentToChildCombinations);
    %currentgBSegments = gB.selectByGrainId([107,77;93,77]);
    previousgBSegments = gB(grains(prevGrainsLogical));
    if isempty(currentgBSegments)
        continue
    end
    %parentgBSegments_Twinning = parentgBSegments(angle(parentgBSegments.misorientation,tS.parentTwinMisorientation())/degree < errorInDegrees);
    [gBtwinVariants] = tS.variantDetermination2((currentgBSegments.misorientation), errorInDegrees);
    [angleFromPerfectTwin] = tS.anglefromTwin(currentgBSegments.misorientation);
    
    twinninggBsegments = currentgBSegments(logical(gBtwinVariants));
    

    % excludes the 0 for not a twinvariant 
    gBactualTwinVariants = gBtwinVariants(logical(gBtwinVariants));
    angleFromPerfectTwin = angleFromPerfectTwin(logical(gBtwinVariants));
    
    %twinninggBsegments with a unique parent twin combination to quiver the
    % habit planes
    [~, uniqueIdx] = unique(twinninggBsegments.grainId, 'rows', 'stable');
    uniqueTwinninggBsegments = twinninggBsegments(uniqueIdx);

    % plot gb segments with color depending on the variant
    figure(varBoundaries.Number);
    plot(twinninggBsegments,colmap(gBactualTwinVariants,:),'linewidth',2.5,'HandleVisibility','off','linewidth',4)
    
    %plot the misorientation from perfect twin
    figure(varDevFromTwin.Number);
    plot(twinninggBsegments,mapToParula(angleFromPerfectTwin,errorInDegrees),'linewidth',2.5,'HandleVisibility','off','linewidth',4)
    
    % quiver expected habit planes with color depending on the variant
    figure(varHabitPlanes.Number);
    k1Array = tSvariants(gBactualTwinVariants(uniqueIdx)).k1;

    for uni = 1:length(uniqueIdx)
        thisTwinninggBsegments = uniqueTwinninggBsegments(uni);
        parent = thisTwinninggBsegments.grainId(1);

        twin = thisTwinninggBsegments.grainId(2);
        k1inSpecimenCoordinates = grains(parent).meanOrientation*k1Array(uni);
        k1inSpecimenCoordinates.antipodal=1;
        quiver(grains(twin),cross(k1inSpecimenCoordinates,zvector),'autoScale','off','color',colmap(gBactualTwinVariants(uniqueIdx(uni)),:),'linewidth',4)
        hold on
    end


    % if i==1
    %     plot(twinninggBsegments,'r','linewidth',2.5,'HandleVisibility','off')
    % elseif i==2
    %     plot(twinninggBsegments,'b','linewidth',2.5,'HandleVisibility','off')
    % else
    %     plot(twinninggBsegments,'g','linewidth',2.5,'HandleVisibility','off')
    % end
    %plot(twinninggBsegments,repmat([1 0.5*i-0.5 0],[length(twinninggBsegments),1]),'linewidth',2.5,'HandleVisibility','off')
    % Flatten grainIds of twinningGbsegments into a single column and get unique values
    allAnalysedGrains = unique(twinninggBsegments.grainId(:));

    % Find unique elements in allAnalysed that are not current grains
    twinnedGrains = setdiff(allAnalysedGrains, grains(currentGrainsLogical).id);
    

    
    % % plot crystal directions of current grains
    % grainDirs =  arrayfun(@(x) {x, rotate(Miller(1,1,-2,0,CS,'UTVW'), x.meanOrientation), rotate(Miller(1,-1,0,0,CS,'UTVW'), x.meanOrientation),rotate(Miller(0,0,0,1,CS,'UTVW'), x.meanOrientation)},grains(currentGrainsLogical),'UniformOutput',false);
    % grainDirs = vertcat(grainDirs{:});
    % for row = 1:size(grainDirs,1)
    %     thisgrain = grainDirs{row,1};
    %     vec_a1 = grainDirs{row,2};
    %     vec_a2 = grainDirs{row,3};
    %     vec_c = grainDirs{row,4};
    %     %both 'linecolor' and 'color' do not work
    %     quiver(thisgrain, vec_a1,'color','r','LineStyle',neg_posMarker{ispos(vec_a1.z)+1},'noScaling')%'AutoScaleFactor',1.3)
    %     quiver(thisgrain, vec_a2,'color','g','LineStyle',neg_posMarker{ispos(vec_a2.z)+1},'noScaling')%'AutoScaleFactor',1.3)
    %     quiver(thisgrain, vec_c,'color','b','LineStyle',neg_posMarker{ispos(vec_c.z)+1},'noScaling')%'AutoScaleFactor',1.3)
    % end

    % %quiver the habit planes
    % k1Array = tSvariants(variantIndex).k1;
    % k1inSpecimenCoordinates = parentOrientation*k1Array;
    % k1inSpecimenCoordinates.antipodal=1;
    % 
    % for i = 1:length(grains(isTwinning))
    %     quiver(twinGrains(i),cross(k1inSpecimenCoordinates(i),zvector),'autoScale','off','color',grainColors(i,:))
    %     hold on
    % end

    % track already analysed grains

    
    currentGrainsLogical = false([length(grains),1]);
    currentGrainsLogical(twinnedGrains) = true;
    prevGrainsLogical = prevGrainsLogical | currentGrainsLogical;

    %plot crystal shapes
    hGrains = plot(grains(currentGrainsLogical),0.7*cS1, 'faceAlpha', 0.1)%, 'linewidth',2.5);
    hText = text(grains(currentGrainsLogical), i, 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 90);
    figure(varBoundaries.Number);
    hGrains = plot(grains(currentGrainsLogical),0.7*cS1, 'faceAlpha', 0.1)%, 'linewidth',2.5);
    hText = text(grains(currentGrainsLogical), i, 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 90);
    
    
end


figure(varBoundaries.Number);
mtexColorMap(colmap)
setColorRange([1 tS.numberOfVariants]);
myColorbar = mtexColorbar('title',append('Variants of the ', char(tS),'-system'), 'TickLength', 0);
tickDistance = (length(tS.symmetrise)-1)/length(tS.symmetrise);
myColorbar.Ticks = 1+tickDistance/2:tickDistance:length(tS.symmetrise);
myColorbar.TickLabels = arrayfun(@(x)(convertCharsToStrings(char(x))), tS.symmetrise);
saveFigure(append(file,'RetwinningVariants_alongGBs.png'))

figure(varDevFromTwin.Number);
colormap(parula);
colorbar;
caxis([0 errorInDegrees]);
exportgraphics(figure(varDevFromTwin.Number), append(file,'deviationFromPerfectTwinningrelation.png'), 'Resolution', 300);

figure(varHabitPlanes.Number);
mtexColorMap(colmap)
setColorRange([1 tS.numberOfVariants]);
myColorbar = mtexColorbar('title',append('K1 planes of the variants of the ', char(tS),'-system'), 'TickLength', 0);
tickDistance = (length(tS.symmetrise)-1)/length(tS.symmetrise);
myColorbar.Ticks = 1+tickDistance/2:tickDistance:length(tS.symmetrise);
myColorbar.TickLabels = arrayfun(@(x)(convertCharsToStrings(char(x))), tS.symmetrise);
legend off
saveFigure(append(file,'HabitPlanesOfRetwinningVariants_alongGBs.png'))

%% common Zone axis of multiple twinning
% parentId = 107;
% primaryTwinId = 77;
% secondaryTwinId = 99;
% 
% % twin system which is to be analysed
% tS = twinningSystem.hexagonal_1121_0001(CS);
% 
% %allowed misorientation from the theoretical twin orientation in degrees
% errorInDegrees = 8;
% 
% %stuff
% tSsym = tS.symmetrise();
% parent = grains(parentId);
% t1 = grains(primaryTwinId);
% t2 = grains(secondaryTwinId);
% 
% %variant determination
% twinVar1 = tS.variantDetermination2(inv(t1.meanOrientation) * parent.meanOrientation, errorInDegrees);
% twinVar2 = tS.variantDetermination2(inv(t2.meanOrientation) * t1.meanOrientation, errorInDegrees);
% 
% twinVar1 = tSsym(twinVar1);
% twinVar2 = tSsym(twinVar2);
% 
% %the twin rotation axis as specimen coordinate vector3d
% primaryRotAxis = rotate(twinVar1.rotAxis, parent.meanOrientation);
% secondaryRotAxis = rotate(twinVar2.rotAxis, t1.meanOrientation);
% 
% %checking if the two vector3ds are aligned
% if angle(primaryRotAxis, secondaryRotAxis)/degree<5 || angle(primaryRotAxis,-secondaryRotAxis)/degree<5
%     disp("Primary and secondary Twinning have the same zone axis")
% else
%     disp("Primary and secondary Twinning DO NOT share the same zone axis")
% end
%% Variant Map
%Definitions
tS = twinningSystem.hexagonal_1121_0001(CS);

%tS = twinningSystem.hexagonal_1122_1124(CS);

errorInDegrees = 10;

%parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
misOrientation = inv(grains.meanOrientation).*parentOrientation;



colmap = hsvColorMap(tS.numberOfVariants);

[variantIndex] = tS.variantDetermination2(misOrientation,errorInDegrees);

figure;
f = plot(ebsd,cibw)%,'micronbar','off')
legend off

% BROKEN
%sth like this
hold on
plot(grains(variantIndex~=0), variantIndex(variantIndex~=0))
mtexColorbar
colormap viridis
%

%plot(grains(isParentGrain), cS1*0.2, 'facecolor', 'gray', 'DisplayName', 'Parent grain')
figure;
plot(ebsd,cibw)
hold on
plot(grains,colors,'FaceAlpha',0.3)
%plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
hold on
plot(grains.boundary)
isTwinning = variantIndex~=0;
twinGrains = grains(isTwinning);
grainColors = hsvColorMap(tS.numberOfVariants);

filteredIndex = variantIndex(isTwinning);

warning('off', 'all');
for i = 1:length(grains(isTwinning))
    plot(twinGrains(i), cS1,'facecolor', grainColors(filteredIndex(i),:),'HandleVisibility','off')
    hold on
end
warning('on', 'all');

mtexColorMap(colmap)
setColorRange([1 tS.numberOfVariants]);
myColorbar = mtexColorbar('title',append('Variants of the ', char(tS),'-system'), 'TickLength', 0);
tickDistance = (length(tSvariants)-1)/length(tSvariants);
myColorbar.Ticks = 1+tickDistance/2:tickDistance:length(tSvariants);
myColorbar.TickLabels = arrayfun(@(x)(convertCharsToStrings(char(x))), tSvariants);
%saveFigure(char(append(file,tS.filename_char(),'Variants.png')))
exportgraphics(gcf, char(append(file,tS.filename_char(),'Variants.png')), 'Resolution', 300)

figure;
plot(ebsd,cibw)
hold on
plot(grains,colors,'FaceAlpha',0.3)
%plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
hold on
plot(grains.boundary)
k1Array = tSvariants(filteredIndex).k1;
k1inSpecimenCoordinates = parentOrientation*k1Array;
k1inSpecimenCoordinates.antipodal=1;

for i = 1:length(grains(isTwinning))
    quiver(twinGrains(i),cross(k1inSpecimenCoordinates(i),zvector),'autoScale','off','color',grainColors(filteredIndex(i),:),'LineWidth',4)
    hold on
end

mtexColorMap(colmap)
setColorRange([1 tS.numberOfVariants]);
myColorbar = mtexColorbar('title',append('Variants of the ', char(tS),'-system'));
tickDistance = (length(tSvariants)-1)/length(tSvariants);
myColorbar.Ticks = 1+tickDistance/2:tickDistance:length(tSvariants);
myColorbar.TickLabels = arrayfun(@(x)(convertCharsToStrings(char(x))), tSvariants);
%saveFigure(char(append(file,tS.filename_char(),'Variant_Traces.png')))
exportgraphics(gcf, char(append(file,tS.filename_char(),'Variant_Traces.png')), 'Resolution', 300)
%% Schmid Factors of different variants based on a global stress state (e.g. uniaxial tension)

%parentOrientation - orientation (CS → xyz)
%parentOrientation;
% parentOrientation = maxOri;
%stress state in sampleCoordinates
errorInDegrees = 10;
stressTensor = stressTensor.uniaxial(yvector);
%twin System definition
tS = twinningSystem.hexagonal_1121_0001(CS);

% twinSystemName = '11-22';
% tS = twinningSystem.hexagonal_1122_1124(CS);

r = inv(parentOrientation)*yvector;

matrixStressTensor = rotate(stressTensor, inv(parentOrientation));

misOrientation = inv(grains.meanOrientation).*parentOrientation;

[isTwinning, variantIndex,tSvariants] = tS.variantDetermination(misOrientation,errorInDegrees);
twinGrains = grains(isTwinning);
disp(['possible Schmid factors'])
schmidFactors = tSvariants.SchmidFactor(matrixStressTensor)

fig = figure;
plot(grains.boundary)
hold on
posSave = fig.Position;
for i = 1:length(grains(isTwinning))
    plot(twinGrains(i), schmidFactors(variantIndex(i)), 'label', num2str(variantIndex(i)))
    text(twinGrains(i), num2str(round(schmidFactors(variantIndex(i)),2)),'uniformOutput',false)
    hold on
end
fig.Position = posSave;
setColorRange([0 0.5])
mtexColorbar('title',append('Schmid factor of the primary ', tS.char(),'-twins'))
saveFigure(append(file,tS.char(),['SchmidFactorOfPrimaryTwins.png']))

%% SAME FOR THE RETWIN
% %isPrimaryTwin = [25,12];
% isPrimaryTwin = 25;
% primaryTwinOrientation = mean(grains.meanOrientation(isPrimaryTwin));
% 
% errorInDegrees = 10;
% 
% %parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
% misOrientation = inv(grains.meanOrientation).*primaryTwinOrientation;
% 
% twinStressTensor = rotate(stressTensor, inv(primaryTwinOrientation));
% 
% [isTwinning, variantIndex,tSvariants] = tS.variantDetermination(misOrientation,errorInDegrees);
% twinGrains = grains(isTwinning);
% 
% schmidFactors = tSvariants.SchmidFactor(twinStressTensor);
% slipSys = slipSystem.basal(CS);
% slipSys = slipSys.symmetrise;
% slipSys.SchmidFactor(twinStressTensor)
% 
% fig = figure;
% plot(grains.boundary)
% posSave = fig.Position;
% for i = 1:length(grains(isTwinning))
%     plot(twinGrains(i), schmidFactors(variantIndex(i)), 'label', num2str(variantIndex(i)))
%     text(twinGrains(i), num2str(round(schmidFactors(variantIndex(i)),2)),'uniformOutput',false)
%     hold on
% end
% %FOR FIGURE RESIZE SET A STOPPOINT IN THE NEXT LINE
% fig.Position = posSave;
% setColorRange([0 0.5])
% mtexColorbar('title',append('Schmid factor for the retwinning of the ', tS.char(),''))
% saveFigure(append(file,tS.char(),['SchmidFactorOfRetwins.png']))
% %% FUNCTIONS
% function result = cell2ndlogical(C)
%     % Get the size of the cell array (n-dimensional)
%     cellSize = size(C);
% 
%     % Determine the length of the logical arrays in the cells
%     % Assuming all cells contain logical arrays of the same length
%     firstLogicalArray = C{1};  % Get the first logical array
%     arrayLength = length(firstLogicalArray);  % Determine its length
% 
%     % Convert the n-dimensional cell array into a 2D matrix
%     % Use cell2mat to concatenate all cells into a big 2D matrix
%     concatenated_matrix = cell2mat(C);
% 
%     % Reshape the concatenated matrix into an (arrayLength, cellSize(1), cellSize(2), ...)
%     reshapeDims = [arrayLength, cellSize];  % The target dimensions
%     reshaped_matrix = reshape(concatenated_matrix, reshapeDims);
% 
%     % Permute the reshaped matrix so that it has the form (cellSize(1), cellSize(2), ..., arrayLength)
%     numDims = length(cellSize);  % Get the number of dimensions of the original cell array
%     permuteOrder = [2:(numDims+1), 1];  % Permutation order to bring the logical array length to the last dimension
%     result = permute(reshaped_matrix, permuteOrder);
% end
% 
% function angleValue = conditionalAngle(x)
%     if x.grainId == 0
%         angleValue = NaN;
%     else
%         angleValue = angle(grainMeanOris(x.grainId), x.orientation);
%     end
% end
%% double twinning ori relationship
% % https://www.sciencedirect.com/science/article/pii/S0921509314006261?via%3Dihub#f0030
% %10-11 -> 10-12
% doubleTwin_Mg_1011_1012 = (tS_1011_Mg_variants.parentTwinMisorientation*tS_1012_Mg_variants.parentTwinMisorientation);
% angle(unique(doubleTwin_Mg_1011_1012))/degree
% round(axis(unique(doubleTwin_Mg_1011_1012)))
% 
% %11-21 -> 11-21
% tS_1121_Re_variants = tS_1121_Re.symmetrise;
% angle(tS_1121_Re.parentTwinMisorientation)/degree
% angle(tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation)/degree
% unique(tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation)
% 
% doubleTwinVars = tS_1121_Re_variants.parentTwinMisorientation*tS_1121_Re_variants.parentTwinMisorientation
% cozoneDoubletwins = doubleTwinVars(abs(angle(doubleTwinVars)/degree - 68.8243) <10e-5); %for Re with 68 degs
% % %g_collector from %% plot selected crystal shapes
% % 
% % lastEBSD = g_collector(end);
% % twinOrientaitons = lastEBSD.orientations *tS_1121_Re_variants.parentTwinMisorientation;
% % 
% % figure(cSfigure)
% % plot(twinOrientaitons * cS1, 'FaceColor','k', 'FaceAlpha',0.2)
% % hold on
%%
% % invariant plane of the co-zonal double twin, in which direction do I
% % rotate around the zone axis? Crocker double twinning
% %x = Miller(-1,-1,2,6,CS_Rhenium,'UVTW')
% x = Miller (1,1,-2,1,CS_Rhenium,'UVTW')
% % x parallel to the twinning direction of the primary twin, phi has been
% % calculated according to crocker
% 
% inv1_1 = rotate(x, rotation.byAxisAngle(Miller(-1,1,0,0,CS_Rhenium,'UVTW'), -1.012526951));
% % inv1_1.dispStyle= 'hkil'
% round(inv1_1)
% inv1_2 = rotate(x, rotation.byAxisAngle(Miller(1,-1,0,0,CS_Rhenium,'UVTW'), -1.012526951));
% % inv1_2.dispStyle= 'hkil'
% round(inv1_2)
% 
% inv2_1 = rotate(x, rotation.byAxisAngle(Miller(-1,1,0,0,CS_Rhenium,'UVTW'), 0.005258291422));
% % inv2_1.dispStyle= 'hkil'
% round(inv2_1)
% % inv1_2.dispStyle= 'hkil'
% inv2_2 = rotate(x, rotation.byAxisAngle(Miller(1,-1,0,0,CS_Rhenium,'UVTW'), 0.005258291422));
% % inv2_2.dispStyle= 'hkil'
% round(inv2_2)
%%
% x = Miller(1,0,-1,1,CS_Magnesium,'hkil')
% inv1_1 = rotate(x, rotation.byAxisAngle(Miller(-1,2,-1,0,CS_Magnesium,'UVTW'), 0.5222354181))
% round(inv1_1)
% 
% %correct same result as in crocker double twinning 
% inv1_2 = rotate(x, rotation.byAxisAngle(Miller(1,-2,1,0,CS_Magnesium,'UVTW'), 0.5222354181))
% round(inv1_2)
% 
% inv2_1 = rotate(x, rotation.byAxisAngle(Miller(-1,2,-1,0,CS_Magnesium,'UVTW'), -0.9798112))
% round(inv2_1)
% 
% %correct same result as in crocker double twinning 
% inv2_2 = rotate(x, rotation.byAxisAngle(Miller(1,-2,1,0,CS_Magnesium,'UVTW'), -0.9798112))
% round(inv2_2)
% 
% 
% % % inv1_1.dispStyle= 'hkil'
% % % 
% % inv1_2 = rotate(x, rotation.byAxisAngle(Miller(1,-1,0,0,CS_Rhenium,'UVTW'), -1.012526951))
% % inv1_2.dispStyle= 'hkil'
% % round(inv1_2)
% % inv2_1 = rotate(x, rotation.byAxisAngle(Miller(-1,1,0,0,CS_Rhenium,'UVTW'), 0.005258291422))
% % inv2_1.dispStyle= 'hkil'
% % round(inv2_1)
% % inv1_2 = rotate(x, rotation.byAxisAngle(Miller(1,-1,0,0,CS_Rhenium,'UVTW'), -1.012526951))
% % inv1_2.dispStyle= 'hkil'
% % inv2_2 = rotate(x, rotation.byAxisAngle(Miller(1,-1,0,0,CS_Rhenium,'UVTW'), 0.005258291422))
% % inv2_2.dispStyle= 'hkil'
% % round(inv2_2)
