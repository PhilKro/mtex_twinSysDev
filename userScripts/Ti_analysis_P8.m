%% initialize mtex

clear 
close all

% Mtex directions (important when working with EBSD)
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');


%% Load data

% filename
fname = '2024-04-30_Kamila TKD Ti Pillar8 Site 1 Map Data 31_BW158.ang';

% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is safe
thresh_ci = 0.12;

% crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% used
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');




%%% perform corrections for tkd tilt and misaligned axes %%%
ebsd = EBSD.load(fname,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');


ebsd0=ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;

%%% some initial test plot %%%
ipf=ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection=yvector;
colors = ipf.orientation2color(ebsd.orientations);

figure;
plot(ebsd,colors)
saveFigure([fname,'ebsd_ipf_Y.png'])


%% perform grain reconstruction

% min angle used for grain reconstruction
grainThresh = 15*degree;
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during reconstruction
alpha = 3;


%%% perform reconstruction %%%
ebsd = ebsd0;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',grainThresh,'alpha',alpha);
ebsd(grains(grains.grainSize<=minSize)) = 'notIndexed';
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',grainThresh,'alpha',alpha);


%%% test plot %%%

figure;
plot(grains,grains.meanOrientation)

% create and plot crystal shapes
cS1 = crystalShape.hex(ebsd.CS);

hold on
plot(grains,0.5*cS1,'facecolor','r')

% plot boundaries
plot(grains.boundary)

saveFigure([fname,'Grains_crystalShapes1.png'])

figure;
plot(grains,grains.meanOrientation)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)
saveFigure([fname 'grain_ID.png'])



%% twinning


plotColor = {'r','b','c','Yellow'};

K1 = Miller({1,1,-2,2},{1,1,-2,1},{1,0,-1,2},{1,0,-1,1},CS);
K1s = Miller({1,1,-2,2},{1,1,-2,1},{1,0,-1,2},{1,0,-1,1},CS);

L = Miller({1,-1,0,0},{1,-1,0,0},{1,-2,1,0},{1,-2,1,0},CS);
Ls = Miller({-1,1,0,0},{-1,1,0,0},{-1,2,-1,0},{-1,2,-1,0},CS);

clear twinSys
for i=1:length(K1)
    twinSys{i} = orientation.map(K1(i),K1s(i),L(i),Ls(i));
end

v = variants(twinSys{1});
t2t = twinSys{1} .* inv(v(:));

for i=1:length(t2t)
    twinSys{end+1} = t2t(i);
    plotColor{end+1} = 'Orange';
end

%% calc and plot twinBounds


% take only boundaries between Ti grains
gB = grains.boundary('Ti','Ti');


twinMatch = zeros(length(twinSys),length(gB));
for i=1:length(twinSys)
    twinMatch(i,:) = angle(gB.misorientation,twinSys{i});
end

[twinBoundFit,twinBoundSys] = min(twinMatch,[],1);

isTwinning = twinBoundFit < 10*degree;



figure;
plot(ebsd,ebsd.ci)
mtexColorMap(gray)
hold on
plot(gB,'linecolor','k','LineWidth',1);
plot(grains,0.3*cS1,'facecolor','r')

for i=1:length(twinSys)
    hold on
    % plot(gB(isTwinning & twinBoundSys == i),'lineColor',plotColor{i},'LineWidth',3,'DisplayName',char(K1(i)),'latex');
    plot(gB(isTwinning & twinBoundSys == i),'lineColor',plotColor{i},'LineWidth',3);
end


saveFigure([fname,'full_twinBound.png'])















%% plot misorientations ( see here: https://mtex-toolbox.github.io/TwinningBoundaries.html )

figure;

% take only boundaries between Ti grains
gB = grains.boundary('Ti','Ti');

% plot ebsd data and misorientation of boundaries
plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
legend off
hold on
% plot shapes
plot(grains,0.3*cS1,'facecolor','r')
% plot misorientation profile
plot(gB,gB.misorientation.angle./degree,'linewidth',5)
hold off
mtexColorbar('title','misorientation angle')
colormap jet
caxis([0 90])

saveFigure('misorientation_profile.png')


%%% plot histogram %%%
figure;
histogram(gB.misorientation.angle./degree,40)
xlabel('misorientation angle (degree)')
saveFigure('F_histogramMisori.png')

%% check for twinning orientation relationship (https://mtex-toolbox.github.io/TwinningBoundaries.html)

% misorientations that might be twin boundaries
ind = gB.misorientation.angle > 64*degree & gB.misorientation.angle < 65*degree;
mori = gB.misorientation(ind);

% quick plot of rotation axis in IPF
% figure;
% plotAxisDistribution(mori,'contourf')

% determine the mean of the cluster
mori_mean = mean(mori,'robust');

% determine the closest special orientation relation ship
[a,b,c,d]=round2Miller(mori_mean)

twinning3 = orientation.map(round(a),round(b),...
  round(c),round(d));


%% 


%11-2-2
K1 = Miller(1,1,-2,2,CS);
K1s = Miller(1,1,-2,2,CS);

% L = Miller(1,-1,0,0,CS,'uvw');
% Ls = Miller(1,0,-1,0,CS,'uvw');
L = Miller(1,-1,0,0,CS);
Ls = Miller(-1,1,0,0,CS);


%11-21
K1_2 = Miller(1,1,-2,1,CS);
K1s_2 = Miller(1,1,-2,1,CS);

L_2 = Miller(1,-1,0,0,CS,'uvw');
Ls_2 = Miller(-1,1,0,0,CS,'uvw');

%10-12
K1_3 = Miller(1,0,-1,2,CS);
K1s_3 = Miller(1,0,-1,2,CS);

L_3 = Miller(1,-2,1,0,CS,'uvw');
Ls_3 = Miller(-1,2,-1,0,CS,'uvw');


% K2 = Miller(-1,-1,2,4,CS);
e1 = Miller(1,1,-2,-3,CS,'uvw');
% e2 = Miller(2,2,-4,3,CS,'uvw');
% 
% L2 = Miller(-1,0,1,0,CS,'uvw');
% 
% L3 = Miller(1,-1,0,0,CS);
% 
% % % 10-12
% % K1 = Miller(1,0,-1,2,CS);
% % K2 = Miller(-1,0,1,-2,CS);
% % e1 = Miller(1,0,-1,-1,CS,'uvw');
% % e2 = Miller(1,0,-1,1,CS,'uvw');
% % 
% % lambda = Miller(1,-2,1,0,CS,'uvw');
% % lambda2 = Miller(-1,2,-1,0,CS,'uvw');





% define twinning orientation relationship (11-22) twin
twinning = orientation.map(K1,K1s,L,Ls);

% define twinning orientation relationship (10-11) twin
twinning2 = orientation.map(K1_2,K1s_2,L_2,Ls_2);

twinning3 = orientation.map(K1_3,K1s_3,L_3,Ls_3);

% % define twinning orientation relationship (11-22) twin
% twinning2 = orientation.map(Miller(1,0,-1,2,CS),Miller(1,0,-1,-2,CS),...
%   Miller(1,0,-1,-1,CS,'uvw'),Miller(1,0,-1,1,CS,'uvw'));

% check twinning axis and angle

% axis
rot = round(twinning.axis)
rot2 = round(twinning2.axis)
rot3 = round(twinning3.axis) 

% the rotational angle
twinning.angle / degree
twinning2.angle / degree
twinning3.angle / degree


%% plot polefigures
ipf = ipfColorKey(CS);
ipf.inversePoleFigureDirection = xvector;
colors = ipf.orientation2color(grains.meanOrientation);

% ids=[106,105,93,91];
ids=[93,105];

inds = id2ind(grains,ids);

figure;
% plotPDF(grains.meanOrientation,colors,[K1,e1,lambda2])
plotPDF(grains(inds).meanOrientation,colors(inds,:),[L K1 K1_2],'markerSize',15,'grid')
hold on
plot(grains(inds(2)).meanOrientation*L.symmetrise,'plane')







%%

% take only boundaries between Ti grains
gB = grains.boundary('Ti','Ti');

twinMatch = angle(gB.misorientation,twinning);
isTwinning = twinMatch < 5*degree;
twinBound = gB(isTwinning);

figure;
% newMtexFigure('layout',[1,2])
% % plot ebsd data and misorientation of boundaries
% plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
% legend off
% hold on
% % plot shapes
% plot(grains,0.3*cS1,'facecolor','r')
% % plot misorientation profile
% plot(twinBound,'linecolor','k','linewidth',2)
% hold off
% 
% nextAxis

% plot ebsd data and misorientation of boundaries
plot(ebsd,ebsd.orientations,'faceAlpha',0.3,'micronbar','off')
legend off
hold on
% plot shapes
plot(grains,0.3*cS1,'facecolor','r')
% plot misorientation profile
plot(gB,twinMatch/degree,'linewidth',5)
hold off
% mtexColorbar('title','fit to 10-12 system')
mtexColorbar('title','fit to 11-22 system')
colormap jet
caxis([0 45])








