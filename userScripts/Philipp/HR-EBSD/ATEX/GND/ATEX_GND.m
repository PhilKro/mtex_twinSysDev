%% initialize mtex

clear 
close all
load('CS_cS_andPlanes.mat')

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;
%% Load data from .ang file

% filename
ang_dir = 'data\Ti_PTP\Ti_cAxis_PTP_1\';
ang_file = 'Insitu Specimen 1 Site 2 Map Data 1_BW263.ang';
path = fullfile(ang_dir, ang_file);

% define tilt, only for TKD
tkd_tilt = -10 * degree;

% define threshold for CI, 0.15 is usually safe
thresh_ci = 0.12;

% crystal symmetry, make sure to use correct space group and c/a ratio

%Ti
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

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
%% import the ATEX data into the ebsd properties

ATEX_data_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0007_Ti_cAxis_PTP_1\ATEX';
ATEX_data_file = 'Insitu Specimen 1 Site 2 Map Data 1_SphericallyIndexed_EXPORT_PIXELS.txt';

ebsd = read_atex_data_complete(fullfile(ATEX_data_dir,ATEX_data_file), ebsd);

%% gridify the data and calcGrains
ebsd0 = ebsd;
ebsd=ebsd.gridify;
ebsd0=ebsd0.gridify;

% min angle used for grain reconstruction
grainThresh = 5*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
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

% ipf with ci greyscale
figure;
plot(ebsd,ebsd.prop.x_NumberOfIterationsOfTheIC_GNAlgorithm__average_5_361e1_S_Dev_,'micronbar','off')
title('IC GN Iterations')
mtexColorbar

%% classical GND calculation
% ----
% Geometrically necessary dislocations 1: Curvature, https://mtex-toolbox.github.io/GND.html
% ----

F = halfQuadraticFilter;
ebsdSGND = smooth(ebsd,F,'fill',grains);

% compute the curvature tensor
kappa = ebsd.curvature;

% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
% kappa(2,3)

newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsd,kappa{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end

% unify the color rage  - you may also use setColoRange equal
setColorRange equal
setColorRange([-0.2,0.2])
drawNow(gcm,'figSize','large')
mtexColorbar
% saveFigure(append(file,'_Curvature.png'))
% ----
% Geometrically necessary dislocations 2: GND calculation
% ----

%dislocationSystem
dS = dislocationSystem.hcp(CS)
%dSSymm = symmetrise(dS,'unique')
[norm(dS(1).b), norm(dS(12).b),1/3 * norm(ebsd.CS.aAxis)]
[norm(dS(13).b), norm(dS(end).b),1/3 * sqrt(norm(ebsd.CS.aAxis)^2+norm(ebsd.CS.cAxis)^2)] %dS(end) is a screw dislocation idk what the burger there is

%Assumption about energy of screw dislocations based on poisson ratio
nu = 0.3;
% energy of the screw dislocations arbitrarily set to 1
dS(dS.isEdge).u = 1/(1 - 0.3);
% energy of the screw dislocations
dS(dS.isScrew).u = 1;

dSRot = ebsdSGND.orientations * dS;

[rho,factor] = fitDislocationSystems(kappa,dSRot);
% the restored dislocation density tensors
alpha = sum(dSRot.tensor .* rho,2);

% we have to set the unit manually since it is not stored in rho
alpha.opt.unit = '1/um';

% ----
% Geometrically necessary dislocations 3: Complete Curvature Tensor plotting
% ----

% restore the full curvature tensor
kappa = alpha.curvature;
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsd,kappa{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end

setColorRange([-0.2,0.2])
drawNow(gcm,'figSize','large');
h = mtexColorbar
default_position = h.Position;
% Set height as a percentage of the default (e.g., 70% of the default height)
percentage = 0.89;  % 70% of default height
new_height = default_position(4) * percentage;
% Update the colorbar position with the new height
h.Position = [default_position(1), default_position(2), default_position(3), new_height];
% from https://github.com/mtex-toolbox/mtex/discussions/2258
%saveFigure(append(file,'_Curvature_Full.png'))

% ----
% Geometrically necessary dislocations 4: Plot total GND density
% ----

%set limits for logarithmic scale in m^-2
% limits = [14.5 16.8];

figure;
% rho .* dSRot.u SHOULD BE AN ENERGY J/m^3 NOT A Dislocation density https://github.com/mtex-toolbox/mtex/discussions/921 
ebsd.prop.totalDislocationEnergy = factor*sum(abs(rho .* dSRot.u),2);
ebsd.prop.totalGNDdensity = sum(abs(rho),2);

plot(ebsd, ebsd.prop.totalGNDdensity, 'micronbar', 'off');
title('Total GND dentsiy')
mtexColorMap parula
mtexColorbar

figure;
plot(ebsd, ebsd.prop.totalDislocationEnergy);
title('Total GND dentsiy')
mtexColorMap parula
mtexColorbar

%% WBV calculation from spherically indexed data
% calculate WBV with integral method
windowSize = 5; % some sort of smoothening parameter
W = ebsd.weightedBurgersVec('integral',windowSize);

% filter out data with low "dislocation amplitude"
% W(W.norm<0.01) = NaN;

Wcr = inv(ebsd.orientations) .* W;



ipf = ipfColorKey(ebsd.CS);
colors = ipf.Miller2Color(Wcr);


signW = sign(dot(W,yvector));


% transparancy scaling
transpRange = [0.001 0.05];
transpField = W.norm;
transpField(transpField<transpRange(1)) = transpRange(1);
transpField(transpField>transpRange(2)) = transpRange(2);
transpField = rescale(transpField);

% take the c+a burgers vector as longest b in hcp to estimate the
% geometrically necessary disc density as in:
% The weighted Burgers vector: a new quantity for constraining dislocation densities and types using electron backscatter diffraction on 2D sections through crystalline materials - WHEELER - 2009 - Journal of Microscopy - Wiley Online Library

ca_b = Miller(0,-1,1,'UVW',CS_Rhenium); % = 1/3*[1,1,-2,3] but not really bc of the factor 3 hidden somewhere
norm(ca_b) % should be roughly 5.2437 Anström
gndDensity = W.norm/((10^-4)*ca_b.norm); %W.norm in um^-1 and ca_b.norm in Angström

% ----
% plotting
% ----
figure; %WBV intensity
% plot(ebsd,signW.*W.norm)
plot(ebsd,W.norm, 'micronbar','off')
caxis([0.001 0.1])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Intensity (\mu m^{-1})'; % Add units to the colorbar label

figure; %GND density in um-2
% plot(ebsd,signW.*W.norm)
plot(ebsd,gndDensity, 'micronbar','off')
caxis([0.001/((10^-4)*ca_b.norm) 0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('GND density estimation')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (\mu m^{-2})'; % Add units to the colorbar label

figure; %GND density in m-2
% plot(ebsd,signW.*W.norm)
plot(ebsd,(10^12)*gndDensity, 'micronbar','off')
caxis([(10^12)*0.001/((10^-4)*ca_b.norm) (10^12)*0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('GND density estimation')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (m^{-2})'; % Add units to the colorbar label
exportScaledFigure(gcf, [fname,'_GNDDensity_m.png'])
%% ATEX data WBV calculation from 

% from edit EBSD.weightedBurgersVec
% % the incomplete Nye tensor
% alpha = dislocationDensity(kappa);
% 
% % the weighted Burgers vector is simply its last column
% W = vector3d(alpha.M(1,3,:,:),alpha.M(2,3,:,:),alpha.M(3,3,:,:));
% W = reshape(W,size(ebsd));
% 

W_Atex = vector3d(ebsd.prop.x_a13ComponentOfTheNye_sTensor_curvaturesAlone_InSampleFrame_CS*10^-6,ebsd.prop.x_a23ComponentOfTheNye_sTensor_curvaturesAlone_InSampleFrame_CS*10^-6,ebsd.prop.x_a33ComponentOfTheNye_sTensor_curvaturesAlone_InSampleFrame_CS*10^-6);


Wcr_ATEX = inv(ebsd.orientations) .* W_Atex;



ipf = ipfColorKey(ebsd.CS);
colors = ipf.Miller2Color(Wcr_ATEX);


signW = sign(dot(W_Atex,yvector));


% transparancy scaling
transpRange = [0.001 0.05];
transpField = W_Atex.norm;
transpField(transpField<transpRange(1)) = transpRange(1);
transpField(transpField>transpRange(2)) = transpRange(2);
transpField = rescale(transpField);

% take the c+a burgers vector as longest b in hcp to estimate the
% geometrically necessary disc density as in:
% The weighted Burgers vector: a new quantity for constraining dislocation densities and types using electron backscatter diffraction on 2D sections through crystalline materials - WHEELER - 2009 - Journal of Microscopy - Wiley Online Library

ca_b = Miller(0,-1,1,'UVW',ebsd.CS); % = 1/3*[1,1,-2,3] but not really bc of the factor 3 hidden somewhere
fprintf('Norm of the ca burgers vector should be 5.53 A (Ti) or 5.24 A (Re). The calculated value is:\n %d \n',round(norm(ca_b),4))
gndDensity = W_Atex.norm/((10^-4)*ca_b.norm); %W.norm in um^-1 and ca_b.norm in Angström
%% ATEX data WBV plotting
figure; %WBV intensity
plot(ebsd,W_Atex.norm, 'micronbar','off')
%caxis([0.001 0.1])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Intensity (\mu m^{-1})'; % Add units to the colorbar label

figure; %GND density in um-2
plot(ebsd,gndDensity, 'micronbar','off')
%caxis([0.001/((10^-4)*ca_b.norm) 0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (\mu m^{-2})'; % Add units to the colorbar label

figure; %GND density in m-2
% plot(ebsd,signW.*W.norm)
plot(ebsd,(10^12)*gndDensity, 'micronbar','off')
%caxis([(10^12)*0.001/((10^-4)*ca_b.norm) (10^12)*0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (m^{-2})'; % Add units to the colorbar label

%% comparison ATEX and spherical indexing WBV

warning('The scale for the ATEX WBV is shifted because of the wierd "bug" with ATEX')
range = [0.001 0.1]

figure; %WBV intensity
plot(ebsd,W_Atex.norm, 'micronbar','off')
caxis(range.*10^-3)
set(gca,'ColorScale','log')
colormap viridis
title('ATEX WBV Intensity [\mu m^{-1}]')

nextAxis

% plot(ebsd,signW.*W.norm)
plot(ebsd,W.norm, 'micronbar','off')
caxis(range)
set(gca,'ColorScale','log')
colormap viridis
title('Spherical indexing WBV Intensity [\mu m^{-1}]')
cb2 = mtexColorbar; % Get the colorbar handle
%cb2.Label.String = 'Intensity (\mu m^{-1})'; % Add units to the colorbar label

% --- Add annotation at bottom of figure ---
annotation('textbox', [0 0 1 0.05], ...  % [x y width height] in normalized figure units
    'String', 'The scale for the ATEX WBV is shifted because of the weird "bug" with ATEX.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 10, ...
    'Interpreter', 'tex');  % allows superscripts etc.

exportScaledFigure(gcf, [ang_file,'_GNDDensity_m.png'],'target_px', 3000)