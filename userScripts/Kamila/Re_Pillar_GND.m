%% Geometrically necessary dislocations 1: Curvature, https://mtex-toolbox.github.io/GND.html
% 
clc
clear 
close all

% Mtex plotting directions
plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;


%% Load data

% read data folder (extract .ang files)
% dirPath = 'onedriveData/Philipp/Ti_PTP/Ti_cAxis_PTP_2'; % The folder containing your .ang files
% filePattern = fullfile(dirPath, '*.ang');
% files = dir(filePattern);
% fname = files(5).name;

dirPath = "userScripts\Kamila\";
fname = "pillar29 Sample 23 Site 2 Map Data 7_Map_Data_7_BW158.ang";

file = fullfile(dirPath, fname);

% define tilt, only for TKD
tkd_tilt = -20 * degree;

% define threshold for CI, 0.15 is safe
thresh_ci = 0.1;

% crystal symmetry, make sure to use correct space group and c/a ratio
CS = crystalSymmetry('6/mmm',[2.76 2.76 4.45],'x||a','mineral','Re');

%%% perform corrections for tkd tilt and misaligned axes %%%
ebsd = EBSD.load(file,CS,'convertEuler2SpatialReferenceFrame','setting 2');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
ebsd(ebsd.ci < thresh_ci) = 'notIndexed';
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');
ebsd0 = ebsd.gridify;


%% cut to subset
% % Plot the EBSD map
ipf = ipfColorKey(ebsd.CS);
figure;
ipf.inversePoleFigureDirection = xvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd('indexed'), colors);
title('Select a rectangular region (2 points)');

selection_confirmed = false;

while ~selection_confirmed
    % Use ginput to get two points (diagonal corners)
    [x, y] = ginput(2);

    % Compute the rectangle bounds
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);

    % Define the rectangle
    dx = xmax - xmin;
    dy = ymax - ymin;

    % Select EBSD data within the defined rectangle
    ind = inpolygon(ebsd, [xmin, ymin, dx, dy]);
    selected_ebsd = ebsd(ind);

    % Plot the selected EBSD region
    figure;
    plot(selected_ebsd, ipf.orientation2color(selected_ebsd.orientations));
    title('Selected EBSD Region');

    % Ask for confirmation
    choice = questdlg('Is the selection fine?', 'Selection Confirmation', ...
                      'Yes', 'No, reselect', 'Yes');

    if strcmp(choice, 'Yes')
        ebsd = selected_ebsd;
        selection_confirmed = true;
        ebsd0=ebsd;
        ebsd=ebsd.gridify;
        ebsd0=ebsd0.gridify;
        cibw = (reshape(vector3d([1,1,1]).*((ebsd.ci-min(ebsd.ci,[],'all'))./(max(ebsd.ci,[],'all')-min(ebsd.ci,[],'all'))),[],1));
        cibw = cibw.xyz;
    end
end

%% grain calc
reforientation = orientation.byEuler(211.5*degree,75.4*degree,203.7*degree, ebsd.CS);
%reforientation = grains(34).meanOrientation;
ebsd.orientations = ebsd.orientations.project2FundamentalRegion(reforientation);

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree,'alpha',3,'minPixel',3);

grod = ebsd.calcGROD(grains);

ebsd.prop.ci_rgb = repmat(mat2gray(ebsd.prop.ci), 1, 1, 3);


%% Geometrically necessary dislocations 1: Smoothening and initial curvature calc


% smoothening
F = halfQuadraticFilter;
F.alpha=3; % amount of smoothening (not sure about defintition)


ebsdSGND = gridify(smooth(ebsd('indexed'),F));

% compute the curvature tensor
kappa = ebsdSGND.curvature;

% 
% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
% kappa(2,3)
% 
figure;
newMtexFigure('nrows',3,'ncols',3);
% 
% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsdSGND,kappa{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end
% 
% unify the color rage  - you may also use setColoRange equal
setColorRange equal
setColorRange([-0.2,0.2])
drawNow(gcm,'figSize','large')
mtexColorbar
saveFigure(append(file,'_Curvature.png'))

%% Geometrically necessary dislocations 2: GND calculation


% define all hcp slip systems
sSBasal = symmetrise(slipSystem.basal(CS),'antipodal');
sSPrism = symmetrise(slipSystem.prismaticA(CS),'antipodal');
sSPyraCA = symmetrise(slipSystem.pyramidalCA(CS),'antipodal');
sSPyraCA2 = symmetrise(slipSystem.pyramidal2CA(CS),'antipodal');
sSPyraA = symmetrise(slipSystem.pyramidalA(CS),'antipodal');

% combine systems
sSAll = [sSBasal;sSPrism;sSPyraCA;sSPyraCA2;sSPyraA];
% sSAll = [sSPyraA];

% create dislocation systems
dS = dislocationSystem(sSAll);



%dSSymm = symmetrise(dS,'unique')
% [norm(dS(1).b), norm(dS(12).b),1/3 * norm(ebsdSGND.CS.aAxis)]
% [norm(dS(13).b), norm(dS(end).b),1/3 * sqrt(norm(ebsdSGND.CS.aAxis)^2+norm(ebsdSGND.CS.cAxis)^2)] %dS(end) is a screw dislocation idk what the burger there is


% Check if <c+a> edge and screw dislocations have the same Burgers vector norm
ca_screw_idx = dS.isScrew & angle(dS.b,sSPyraCA2(1).b)<10^-5;
ca_edge_idx = dS.isEdge & angle(dS.b,sSPyraCA2(1).b)<10^-5;
if any(ca_screw_idx) && any(ca_edge_idx)
    ca_screw_norm = max(dS(ca_screw_idx).b.norm);
    ca_edge_norm = max(dS(ca_edge_idx).b.norm);
    if abs(ca_screw_norm - ca_edge_norm) > 1e-4
        warning('The norm of the <c+a> edge dislocations (%g) does not match the norm of the <c+a> screw dislocations (%g).', ca_edge_norm, ca_screw_norm);
    end
end

% Assumption about energy of dislocations based on b and poisson ratio
nu = 0.3;
% energy of the screw dislocations, scales with b^2
dS(dS.isScrew).u = dS(dS.isScrew).b.norm.^2;

% energy of edge dislocations, scales with b^2/(1-nu)
dS(dS.isEdge).u = dS(dS.isEdge).b.norm.^2/(1 - nu);

dSRot = ebsdSGND.orientations * dS;
%% fit dislocation systems
[rho,factor] = fitDislocationSystems(kappa,dSRot);

% the restored dislocation density tensors
alpha = sum(dSRot.tensor .* rho,2);

% we have to set the unit manually since it is not stored in rho
alpha.opt.unit = '1/um';


%% Geometrically necessary dislocations 3: Complete Curvature Tensor plotting

% restore the full curvature tensor
kappa = alpha.curvature;
figure
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsdSGND,kappa{i,j},'micronBar','off')
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
saveFigure(append(file,'_Curvature_Full.png'))



%% Geometrically necessary dislocations 4: Plot total GND density
%set limits for logarithmic scale in m^-2
limits = [10^(14.5) 10^(16.5)];

figure;
% rho .* dSRot.u SHOULD BE AN ENERGY J/m^3 NOT A Dislocation density https://github.com/mtex-toolbox/mtex/discussions/921 
ebsdSGND.prop.totalDislocationEnergy = factor*sum(abs(rho .* dSRot.u),2);
ebsdSGND.prop.totalGNDdensity = factor*sum(abs(rho),2);

plot(ebsdSGND, ebsdSGND.prop.totalGNDdensity, 'micronbar', 'off');
set(gca,'ColorScale','log')
title('Total GND dentsiy')
caxis(limits)
mtexColorbar
mtexColorMap parula
saveFigure(append(file,'_total_GND_density.png'))




%% Geometrically necessary dislocations 5: Plot GND density for the seperate dislocation systems

GND_sS = rho;

n = size(GND_sS,2);

%%% screw and edge, not used currently
%edge GNDs
GND_E = GND_sS;
GND_E(dSRot.isScrew) = nan;
GND_E = factor*nansum(GND_E,2);

%screw GNDs
GND_S = GND_sS;
GND_S(dSRot.isEdge) = nan;
GND_S = factor*nansum(abs(GND_S),2);
%%%

% define a cell array to store the name of the dislocation system for later
% plotting
DislNames = cell(n,1);

% <a> screw dislocations
aScrew0 = dS.isScrew & dS.b==Miller(1,1,-2,0,CS,'UVTW');
DislNames(aScrew0) = {'<a> Screw'};

aScrew = aScrew0';
aScrew = repmat(aScrew, length(ebsdSGND), 1);

GND_aScrew = GND_sS;
GND_aScrew(~(aScrew)) = NaN;
GND_aScrew = factor*nansum(abs(GND_aScrew),2);


%basal edge GNDs
basalEdge0 = dS.b==Miller(1,1,-2,0,CS,'UVTW') & dS.l==Miller(1,-1,0,0,CS,'UVTW');
DislNames(basalEdge0) = {'Basal Edge'};

basalEdge = basalEdge0';
basalEdge = repmat(basalEdge, length(ebsdSGND), 1);

GND_basal = GND_sS;
GND_basal(~(basalEdge)) = nan;
GND_basal = factor*nansum(abs(GND_basal),2);


%prismatic edge GNDs
prismaticEdge0 = dS.b==Miller(1,1,-2,0,CS,'UVTW')&dS.l==Miller(0,0,0,1,CS,'UVTW');
DislNames(prismaticEdge0) = {'Prism Edge'};

prismaticEdge = prismaticEdge0';
prismaticEdge = repmat(prismaticEdge, length(ebsdSGND), 1);

GND_prismatic = GND_sS;
GND_prismatic(~(prismaticEdge)) = nan;
GND_prismatic = factor*nansum(abs(GND_prismatic),2);


%pyr1a edge GNDs 
Pyr1aEdge0 = dS.b==Miller(1,1,-2,0,CS)&dS.l==Miller(1,-1,0,2,CS,'UVTW');
DislNames(Pyr1aEdge0) = {'Pyr1a Edge'};

Pyr1aEdge = Pyr1aEdge0';
Pyr1aEdge = repmat(Pyr1aEdge, length(ebsdSGND), 1);

GND_pyr1a = GND_sS;
GND_pyr1a(~(Pyr1aEdge)) = nan;
GND_pyr1a = factor*nansum(abs(GND_pyr1a),2);


%<c+a> screw
caScrew0 = dS.isScrew & dS.b==Miller(-1,-1,2,-3,CS,'UVTW');
DislNames(caScrew0) = {'<c+a> screw'};

caScrew = caScrew0';
caScrew = repmat(caScrew, length(ebsdSGND), 1);

GND_caScrew = GND_sS;
GND_caScrew(~(caScrew)) = nan;
GND_caScrew = factor*nansum(abs(GND_caScrew),2);


%pyramidal2 c+a edge
Pyr2caEdge0 = dS.b==Miller(-1,-1,2,-3,CS,'UVTW') & dS.l==Miller(1,-1,0,0,CS,'UVTW');
DislNames(Pyr2caEdge0) = {'Pyr2 c+a Edge'};

Pyr2caEdge = Pyr2caEdge0';
Pyr2caEdge = repmat(Pyr2caEdge, length(ebsdSGND), 1);

GND_pyr2ca = GND_sS;
GND_pyr2ca(~(Pyr2caEdge)) = nan;
GND_pyr2ca = factor*nansum(abs(GND_pyr2ca),2);


%pyramidal c+a edge
PyrcaEdge0 = dS.isEdge & dS.b==Miller(-1,-1,2,-3,CS,'UVTW') & ~Pyr2caEdge0;
DislNames(PyrcaEdge0) = {'Pyr1 c+a Edge'};

PyrcaEdge = PyrcaEdge0';
PyrcaEdge = repmat(PyrcaEdge, length(ebsdSGND), 1);

GND_pyrca = GND_sS;
GND_pyrca(~(PyrcaEdge)) = nan;
GND_pyrca = factor*nansum(abs(GND_pyrca),2);



%% % plotting, per slip system and edge/screw

limits = [1e13 1e16];

gndAlpha = 0.9;

fig1a = figure;
newMtexFigure('layout',[2,2]);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_aScrew, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('<a> Screw GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_basal, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('Basal Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_prismatic, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('Prism Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_pyr1a, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('Pyr <a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
h = mtexColorbar('title', append('rho_{GND}', '\', 'm^{-2}'));
mtexColorMap parula
exportScaledFigure(fig1a, append(file, '_seperated_GND_density1_part1.jpg'), 'noDateString')


fig1b = figure;
newMtexFigure('layout',[2,2]);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb, 'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_caScrew, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('<c+a> Screw GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_pyrca, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('Pyr1 <c+a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsd, ebsd.prop.ci_rgb,'micronbar', 'off'); 
axis image off
hold on
plot(ebsdSGND, GND_pyr2ca, 'micronbar', 'off', 'FaceAlpha',gndAlpha);
title('Pyr2 <c+a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
h = mtexColorbar('title', append('rho_{GND}', '\', 'm^{-2}'));
mtexColorMap parula
exportScaledFigure(fig1b, append(file, '_seperated_GND_density1_part2.jpg'), 'noDateString')


%% plot individual dislocation systems

% re-sort dislocation densities, for nicer plotting
sorting = [find(basalEdge0);
           find(prismaticEdge0);
           find(Pyr1aEdge0);
           find(aScrew0);
           find(PyrcaEdge0);
           find(Pyr2caEdge0);
           find(caScrew0)];

DislNamesSort = DislNames(sorting);
GND_sS_sorted = GND_sS(:,sorting);
dS_sorted = dS(sorting);

% define limits, this time from negative to positive to also show the sign
% of b... 
% !!! So no log plotting for now!!!!
limits = [-5e14 5e14];

n=size(GND_sS_sorted,2);
max_per_fig = 4;
fig_idx = 1;

% loop over all individual systems
for i=1:n
    if mod(i-1, max_per_fig) == 0
        if i > 1
            mtexColorMap jet
            mtexColorbar
            exportScaledFigure(gcf, append(file, '_seperated_GND_density2_part', num2str(fig_idx), '.jpg'), 'noDateString')
            fig_idx = fig_idx + 1;
        end
        figure;
        newMtexFigure('layout',[2,2]);
    end

    nextAxis
    % plot individualn system GND
    plot(ebsdSGND, factor*GND_sS_sorted(:,i), 'micronbar', 'off');
    
    title(DislNamesSort{i})
    caxis(limits);

    % plot burgers vector and slip plane trace(s) (not 100% this is correct, but probably
    % it is
    hold on
    % plotting is different for edge and screw. if edge, there is only 1
    % possible slip plane (cross(b,l))
    % if screw, we plot all the possible slip planes of the corresponding
    % slip systems
    if dS_sorted(i).isEdge
        quiver(grains,grains.meanOrientation*dS_sorted(i).l,'color','w','antipodal','project2plane')
        slipPlanes = grains.meanOrientation*cross(dS_sorted(i).b,dS_sorted(i).l);
        quiver(grains,cross(slipPlanes,zvector),'antipodal','color','r');
    else
        sS_screw = sSAll(vector3d(sSAll.b)==vector3d(dS_sorted(i).b));
        slipPlanes = trace(grains.meanOrientation*sS_screw);
        for j=1:size(slipPlanes,2)
            quiver(grains,slipPlanes(:,j),'antipodal','color','r');
        end
    end

    quiver(grains,grains.meanOrientation*dS_sorted(i).b,'color','k','project2plane')
end

mtexColorMap jet
mtexColorbar

% set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
% h = mtexColorbar('title', append('log_{10}', 'rho_{GND}', '\', 'm^{-2}'));
% default_position = h.Position;
% Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)

% percentage = 0.89;  % 70% of default height
% new_height = default_position(4) * percentage;
% % Update the colorbar position with the new height
% h.Position = [default_position(1), default_position(2), default_position(3), new_height];

% mtexColorMap parula
% saveFigure(append(fname,'fully_seperated_GND_density.png'))
exportScaledFigure(gcf, append(file, '_seperated_GND_density2_part', num2str(fig_idx), '.jpg'), 'noDateString')




%% plot combination of GNDs per burgers vector!

% re-sort dislocation densities
sorting = [find(basalEdge0);
           find(prismaticEdge0);
           find(Pyr1aEdge0);
           find(aScrew0);
           find(PyrcaEdge0);
           find(Pyr2caEdge0);
           find(caScrew0)];

DislNamesSort = DislNames(sorting);
GND_sS_sorted = GND_sS(:,sorting);
dS_sorted = dS(sorting);

limits = [-5e14 5e14];

% find all unique burgers vectors
b_unique = round(unique(dS_sorted.b.round,'noSymmetry','antipodal'));
n=size(b_unique,1);

[~,I]=sort(abs(b_unique.z));
b_unique = b_unique(I);

max_per_fig = 4;
fig_idx = 1;

for i=1:n
    if mod(i-1, max_per_fig) == 0
        if i > 1
            mtexColorMap jet
            mtexColorbar
            exportScaledFigure(gcf, append(file, '_seperated_GND_density3_part', num2str(fig_idx), '.jpg'), 'noDateString')
            fig_idx = fig_idx + 1;
        end
        figure;
        newMtexFigure('layout',[2,2]);
    end

    nextAxis
    % combine all activity from this b (not abs values!)
    GND_b = sum(GND_sS_sorted(:,eq(vector3d(dS_sorted.b.round),vector3d(b_unique(i)),'antipodal') ),2);
    plot(ebsdSGND, factor*GND_b, 'micronbar', 'off');
    
    
    % title(DislNamesSort{i})
    caxis(limits);

    % plot burgers and line vector
    hold on
    % if dS_sorted(i).isEdge
    %     quiver(grains,grains.meanOrientation*dS_sorted(i).l,'color','w','antipodal')
    %     slipPlanes = grains.meanOrientation*cross(dS_sorted(i).b,dS_sorted(i).l);
    %     quiver(grains,cross(slipPlanes,zvector),'antipodal','color','r');
    % else
    %     sS_screw = sSAll(vector3d(sSAll.b)==vector3d(dS_sorted(i).b));
    %     slipPlanes = trace(grains.meanOrientation*sS_screw);
    %     for j=1:size(slipPlanes,2)
    %         quiver(grains,slipPlanes(:,j),'antipodal','color','r');
    %     end
    % end

    quiver(grains,grains.meanOrientation*b_unique(i),'color','k','project2plane')
    % quiver(grains,grains.meanOrientation*b_unique(i),'color','k')
    mtexTitle(char(b_unique(i),'LaTeX'))
end

mtexColorMap jet
mtexColorbar

% set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
% h = mtexColorbar('title', append('log_{10}', 'rho_{GND}', '\', 'm^{-2}'));
% default_position = h.Position;
% Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)

% percentage = 0.89;  % 70% of default height
% new_height = default_position(4) * percentage;
% % Update the colorbar position with the new height
% h.Position = [default_position(1), default_position(2), default_position(3), new_height];

% mtexColorMap parula
% saveFigure(append(fname,'b_seperated_GND_density.png'))
exportScaledFigure(gcf, append(file, '_seperated_GND_density3_part', num2str(fig_idx), '.jpg'), 'noDateString')



%%

% % % %% twins check
% % % 
% % % tS1 = symmetrise(slipSystem.twinC2(CS),'antipodal');
% % % tS2 = symmetrise(slipSystem.twinC1(CS),'antipodal');
% % % tS3 = symmetrise(slipSystem.twinT1(CS),'antipodal');
% % % tS4 = symmetrise(slipSystem.twinT2(CS),'antipodal');
% % % 
% % % tS1Local = grains.meanOrientation * tS1;
% % % tS2Local = grains.meanOrientation * tS2;
% % % tS3Local = grains.meanOrientation * tS3;
% % % tS4Local = grains.meanOrientation * tS4;
% % % 
% % % 
% % % how2plot=plottingConvention;
% % % how2plot.north = yvector;
% % % 
% % % figure;
% % % plot(ebsd0,ebsd0.prop.ci,how2plot)
% % % hold on
% % % for i=1:length(tS1)
% % %     quiver(grains,tS1Local(:,i).trace,'color','k')
% % % end
% % % mtexTitle(tS1(1).n.char)
% % % 
% % % nextAxis
% % % 
% % % plot(ebsd0,ebsd0.prop.ci,how2plot)
% % % hold on
% % % for i=1:length(tS2)
% % %     quiver(grains,tS2Local(:,i).trace,'color','k')
% % % end
% % % mtexTitle(tS2(1).n.char)
% % % 
% % % nextAxis
% % % 
% % % plot(ebsd0,ebsd0.prop.ci,how2plot)
% % % hold on
% % % for i=1:length(tS3)
% % %     quiver(grains,tS3Local(:,i).trace,'color','k')
% % % end
% % % mtexTitle(tS3(1).n.char)
% % % 
% % % nextAxis
% % % 
% % % plot(ebsd0,ebsd0.prop.ci,how2plot)
% % % hold on
% % % for i=1:length(tS4)
% % %     quiver(grains,tS4Local(:,i).trace,'color','k')
% % % end
% % % mtexTitle(tS4(1).n.char)

%% Save GND Data
disp('Saving GND data to mat file...');
save(append(file, '_GND_Data.mat'), 'rho', 'factor', 'alpha', 'ebsdSGND', 'kappa', 'dS', 'GND_sS', 'b_unique', 'CS', 'ebsd0', 'ebsd', 'grains', 'grod', '-v7.3');
disp('Done!');