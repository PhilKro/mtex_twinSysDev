%% Geometrically necessary dislocations 1: Curvature, https://mtex-toolbox.github.io/GND.html
% 
if ~isempty(evalin('base', 'whos')) || ~isempty(get(0, 'Children'))
    choice = questdlg('Variables or figures detected. Run script from the top?', ...
        'Run Confirmation', ...
        'Yes', 'No (Abort)', 'Yes');
    
    if strcmp(choice, 'No (Abort)') || isempty(choice)
        error('Script execution stopped by user.');
    end
    % If 'Yes' is chosen, the script simply continues to the next line
end

clear 
close all

% % Mtex directions (important when working with EBSD)
% setMTEXpref('xAxisDirection','east');
% setMTEXpref('zAxisDirection','intoPlane');

load('CS_cS_andPlanes.mat')

% Mtex directions (important when working with EBSD)
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');


%% Load data

% fname
% fname = '2024-04-30_Kamila TKD Ti Pillar 55 Site 3 Map Data 41 new_BW158c.ang';
dir = 'C:\Users\phkr\OneDrive - empa.ch\mtex-6.2.beta.3\data\Rhenium_PushToPull\20241007_Re_cAxis_PTP2\';
fname = 'Re_c_Axis_PTP_2 Specimen 2 Site 3 Map Data 6_BW263.ang';

% define tilt, only for TKD
tkd_tilt = -11 * degree;

% define threshold for CI, 0.15 is safe
thresh_ci = 0.15;

% crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% used
CS = crystalSymmetry('6/mmm',[2.76 2.76 4.45],'x||a','mineral','Re');

%%% perform corrections for tkd tilt and misaligned axes %%%
ebsd = EBSD.load(fullfile(dir, fname),CS,'convertEuler2SpatialReferenceFrame','setting 2');
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
ebsd0 = ebsd;

% min angle used for grain reconstruction
grainThresh = 6*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;

%%% perform reconstruction %%%

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',grainThresh,'alpha',alpha);
ebsd(grains(grains.grainSize<=minSize)) = 'notIndexed';
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

gBs = grains.boundary(CS.mineral,CS.mineral);


%% restrict to CI of 0.3 (done after grain calculation, such that gBs look nice)

restrictiveCI = 0.3;
ebsd(ebsd.ci < restrictiveCI) = 'notIndexed';
ebsd0 = ebsd.gridify;
ebsd = ebsd0;

%% Geometrically necessary dislocations 1: Smoothening and initial curvature calc


% smoothening
F = halfQuadraticFilter;
F.alpha=1.5; % amount of smoothening (not sure about defintition)


ebsdSGND = gridify(smooth(ebsd('indexed'),F));

% compute the curvature tensor
kappa = ebsdSGND.curvature;

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
saveFigure(append(fname,'_Curvature.png'))

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
% sSAll = [sSBasal;sSPrism;sSPyraCA;sSPyraA];

% create dislocation systems
dS = dislocationSystem(sSAll);

%dSSymm = symmetrise(dS,'unique')
% [norm(dS(1).b), norm(dS(12).b),1/3 * norm(ebsdSGND.CS.aAxis)]
% [norm(dS(13).b), norm(dS(end).b),1/3 * sqrt(norm(ebsdSGND.CS.aAxis)^2+norm(ebsdSGND.CS.cAxis)^2)] %dS(end) is a screw dislocation idk what the burger there is

%Assumption about energy of dislocations based on b and poisson ratio
nu = 0.3;
% energy of the screw dislocations, scales with b^2
dS(dS.isScrew).u = dS(dS.isScrew).b.norm.^2;

% energy of edge dislocations, scales with b^2/(1-nu)
dS(dS.isEdge).u = dS(dS.isEdge).b.norm.^2/(1 - nu);

% rotate dislocation systems
dSRot = ebsdSGND.orientations * dS;

% fit dislocation systems
% [rho,factor] = fitDislocationSystems(kappa,dSRot,'coneprog',1e-6);
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
saveFigure(append(fname,'_Curvature_Full.png'))



%% Geometrically necessary dislocations 4: Plot total GND density
%set limits for logarithmic scale in m^-2
limits = [10^(14.5) 10^(17)];

figure;
% rho .* dSRot.u SHOULD BE AN ENERGY J/m^3 NOT A Dislocation density https://github.com/mtex-toolbox/mtex/discussions/921 
ebsdSGND.prop.totalDislocationEnergy = factor*sum(abs(rho .* dSRot.u),2);
ebsdSGND.prop.totalGNDdensity = factor*sum(abs(rho),2);

plot(ebsdSGND, ebsdSGND.prop.totalGNDdensity, 'micronbar', 'off');
set(gca,'ColorScale','log')
title('Total GND dentsiy')
caxis(limits)
mtexColorbar
mtexColorMap jet
saveFigure(append(fname,'total_GND_density.png'))





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
% PyrcaEdge0 = dS.isEdge & dS.b==Miller(-1,-1,2,-3,CS,'UVTW') & ~Pyr2caEdge0;
PyrcaEdge0 = dS.isEdge & dS.b==Miller(-1,-1,2,-3,CS,'UVTW');
DislNames(PyrcaEdge0) = {'Pyr1 c+a Edge'};

PyrcaEdge = PyrcaEdge0';
PyrcaEdge = repmat(PyrcaEdge, length(ebsdSGND), 1);

GND_pyrca = GND_sS;
GND_pyrca(~(PyrcaEdge)) = nan;
GND_pyrca = factor*nansum(abs(GND_pyrca),2);



%%% plotting, per slip system and edge/screw

limits = [10^14 10^16.5];
width = 5;
grainsSmooth = smooth(grains);
bound = grainsSmooth.boundary;

figure
newMtexFigure('layout',[2,3])

nextAxis
plot(ebsdSGND, GND_aScrew, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('<a> Screw GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_basal, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('Basal Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_prismatic, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('Prism Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_pyr1a, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('Pyr <a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_caScrew, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('<c+a> Screw GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_pyrca, 'micronbar', 'off');
hold on
plot(bound,'lineWidth',width,'linecolor','w')
title('Pyr1 <c+a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

nextAxis
plot(ebsdSGND, GND_pyr2ca, 'micronbar', 'off');
title('Pyr2 <c+a> Edge GND dentsiy')
set(gca,'ColorScale','log')
caxis(limits);

mtexColorMap jet
h = mtexColorbar('title', append('rho_{GND}', '\', 'm^{-2}'));

% set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
% default_position = h.Position;
% % % Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)
% percentage = 0.89;  % 70% of default height
% new_height = default_position(4) * percentage;
% % Update the colorbar position with the new height
% h.Position = [default_position(1), default_position(2), default_position(3), new_height];



saveFigure(append(fname,'seperated_GND_density2.png'))




%% plot individual dislocation systems

% re-sort dislocation densities, for nicer plotting
sorting = [find(basalEdge0);
           find(prismaticEdge0);
           find(Pyr1aEdge0);
           find(aScrew0);
           find(PyrcaEdge0);
           % find(Pyr2caEdge0);
           find(caScrew0)];

DislNamesSort = DislNames(sorting);
GND_sS_sorted = GND_sS(:,sorting);
dS_sorted = dS(sorting);

% plotting params
limits = [10^13.5 10^16];
width = 3;

figure;
newMtexFigure('layout',[4,9])
n=size(GND_sS_sorted,2);
% loop over all individual systems
for i=1:n

    nextAxis
    % plot individualn system GND
    plot(ebsdSGND, factor*abs(GND_sS_sorted(:,i)), 'micronbar', 'off');
    hold on
    plot(bound,'lineWidth',width,'linecolor','w')
    
    burgerString = ['[' num2str(dS_sorted(i).b.round.UVTW) ']'];
    title([DislNamesSort{i}, burgerString])
    set(gca,'ColorScale','log')
    caxis(limits);

    % plot burgers vector and slip plane trace(s) (not 100% sure this is correct, but probably
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


saveFigure(append(fname,'fully_seperated_GND_density2.png'))



%% plot combination of GNDs per burgers vector!

% % re-sort dislocation densities
% sorting = [find(basalEdge0);
%            find(prismaticEdge0);
%            find(Pyr1aEdge0);
%            find(aScrew0);
%            find(PyrcaEdge0);
%            % find(Pyr2caEdge0);
%            find(caScrew0)];
% 
% DislNamesSort = DislNames(sorting);
% GND_sS_sorted = GND_sS(:,sorting);
% dS_sorted = dS(sorting);
% 
% limits = [-5e14 5e14];
% 
% figure;
% % newMtexFigure('layout',[7,6])
% 
% % find all unique burgers vectors
% b_unique = round(unique(dS_sorted.b.round,'noSymmetry','antipodal'));
% n=size(b_unique,1);
% 
% [~,I]=sort(abs(b_unique.z));
% b_unique = b_unique(I);
% 
% 
% for i=1:n
%     nextAxis
%     % combine all activity from this b (not abs values!)
%     GND_b = sum(GND_sS_sorted(:,eq(vector3d(dS_sorted.b.round),vector3d(b_unique(i)),'antipodal') ),2);
%     plot(ebsdSGND, factor*GND_b, 'micronbar', 'off');
% 
% 
%     % title(DislNamesSort{i})
%     caxis(limits);
% 
%     % plot burgers and line vector
%     hold on
%     % if dS_sorted(i).isEdge
%     %     quiver(grains,grains.meanOrientation*dS_sorted(i).l,'color','w','antipodal')
%     %     slipPlanes = grains.meanOrientation*cross(dS_sorted(i).b,dS_sorted(i).l);
%     %     quiver(grains,cross(slipPlanes,zvector),'antipodal','color','r');
%     % else
%     %     sS_screw = sSAll(vector3d(sSAll.b)==vector3d(dS_sorted(i).b));
%     %     slipPlanes = trace(grains.meanOrientation*sS_screw);
%     %     for j=1:size(slipPlanes,2)
%     %         quiver(grains,slipPlanes(:,j),'antipodal','color','r');
%     %     end
%     % end
% 
%     quiver(grains,grains.meanOrientation*b_unique(i),'color','k','project2plane')
%     % quiver(grains,grains.meanOrientation*b_unique(i),'color','k')
%     mtexTitle(char(b_unique(i),'LaTeX'))
% end
% 
% mtexColorMap jet
% mtexColorbar
% 
% % set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% % set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
% % h = mtexColorbar('title', append('log_{10}', 'rho_{GND}', '\', 'm^{-2}'));
% % default_position = h.Position;
% % Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)
% 
% % percentage = 0.89;  % 70% of default height
% % new_height = default_position(4) * percentage;
% % % Update the colorbar position with the new height
% % h.Position = [default_position(1), default_position(2), default_position(3), new_height];
% 
% % mtexColorMap parula
% 
% saveFigure(append(fname,'b_seperated_GND_density.png'))

%% Plot Grain Boundaries and Colored {11-21} Traces with Formatted Legend

% 1. Setup Canvas
figH = figure;
plot(grains.boundary, 'linewidth', 1.5, 'lineColor', [0.2 0.2 0.2]); 
hold on;

% 2. Define the plane family (using 'antipodal' to get the 6 unique trace lines)
h_fam = Miller(1,1,-2,1, grains.CS).symmetrise('antipodal', 'unique');

% Define a distinct color for each variant
variantColors = lines(length(h_fam)); 

% 3. Loop through variants
for i = 1:length(h_fam)
    
    % Current variant
    h_variant = h_fam(i);
    
    % Calculate Trace: Cross product of Plane Normal (n) and Surface Normal (z)
    n_spec = grains.meanOrientation * h_variant;
    traceDir = cross(n_spec, zvector);
    traceDir = normalize(traceDir);
    
    % Create the formatted label string: [ h k i l ]
    % We use num2str on the rounded hkil values
    labelStr = ['[' num2str(round(h_variant.hkil)) ']'];
    
    % Plot the trace (Positive direction)
    % We set 'DisplayName' here so the legend picks it up automatically
    quiver(grains, traceDir, 'color', variantColors(i,:), ...
        'linewidth', 2, 'autoScale', 'off', ...
        'DisplayName', labelStr);
    
    % Plot the trace (Negative direction) - HandleVisibility off to avoid duplicate legend entries
    quiver(grains, -traceDir, 'color', variantColors(i,:), ...
        'linewidth', 2, 'autoScale', 'off', ...
        'HandleVisibility', 'off');
end

% 4. Create Legend
% The legend will display the colored line next to the formatted text
lgd = legend('show', 'Location', 'bestoutside');
title(lgd, 'Trace Variant');
title('Traces of \{11-21\} planes');

% Final adjust
axis equal tight;
hold off;

%% twinvariant det

cS1 = crystalShape.hex(grains.CS);

% first get the parent grain
% 1 determined as largest grain
 isParentGrain = grains.grainSize == max(grains.grainSize);
 parentOrientation = grains.meanOrientation(isParentGrain);
 
% 2 determine by id
%isParentGrain = [93];
%parentOrientation = mean(grains.meanOrientation(isParentGrain));

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

% then figure out the twinning
tS = twinningSystem.hexagonal_1121_0001(CS);
tSvariants = tS.symmetrise;
%tS = twinningSystem.hexagonal_1122_1124(CS);

errorInDegrees = 10;

%parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
misOrientation = inv(grains.meanOrientation).*parentOrientation;



colmap = hsvColorMap(tS.numberOfVariants);

[variantIndex] = tS.variantDetermination2(misOrientation,errorInDegrees);

%plot(grains(isParentGrain), cS1*0.2, 'facecolor', 'gray', 'DisplayName', 'Parent grain')
figure;
plot(ebsd,cibw)
hold on
plot(grains,'FaceAlpha',0.3)
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

%% Continuation: Plotting Specific GND Densities

% --- Definitions and Settings ---
% Plotting limits (Density in m^-2, adjust if necessary based on your data range)
plotLimits = [5*10^14 3*10^16]; 
cmapHotInverted = flipud(hot);

% Ensure factor is applied for m^-2 conversion if not already in the previous variables
% (Based on the previous script, 'factor' appears to handle the conversion to m^-2 
% or the variables like GND_basal are already scaled. We use the variables directly.)

% --- 1. Total GND Density (Viridis) ---
figure;

plot(ebsdSGND, ebsdSGND.prop.totalGNDdensity, 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
title('Total GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_Total_GND_Viridis.png']);


% --- 2. Basal GND Density (Hot) ---
% Using the previously calculated 'GND_basal'
figure;

plot(ebsdSGND, abs(GND_basal), 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
title('Basal GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_Basal_GND_Hot.png']);


% --- 3. Prismatic <a> GND Density (Hot) ---
% Using the previously calculated 'GND_prismatic'
figure;

plot(ebsdSGND, abs(GND_prismatic), 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
title('Prismatic <a> GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_PrismA_GND_Hot.png']);


% --- 4. Pyramidal I <a> GND Density (Hot) ---
% Using the previously calculated 'GND_pyr1a'
figure;

plot(ebsdSGND, abs(GND_pyr1a), 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
title('Pyramidal I <a> GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_Pyr1A_GND_Hot.png']);


% --- 5. Total <c+a> GND Density (Hot) ---
% Summing available c+a components (Screw + Pyr1 Edge). 
% Note: Pyr2 Edge was commented out in the fitting section of the original script.
% If it was included, add 'abs(GND_pyr2ca)' to the sum below.

% Calculate total c+a (ignoring NaNs)
GND_ca_total = nansum([abs(GND_caScrew), abs(GND_pyrca), abs(GND_pyr2ca)], 2); 

% Handle cases where all were NaN (optional, depending on preference)
% GND_ca_total(all(isnan([GND_caScrew, GND_pyrca]), 2)) = NaN;

figure;

plot(ebsdSGND, GND_ca_total, 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
title('Total <c+a> GND Density (Pyr1_c+a_Edge + Pyr2_c+a_Edge + Screw)');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_TotalCA_GND_Hot.png']);


% --- 6. Specific Basal Edge Dislocation b=[-2110] (Hot) ---
% Target Burgers vector: [-2 1 1 0]
target_b = Miller(-2, 1, 1, 0, CS, 'UVTW');

% Identify indices in dS that match this Burgers vector AND are Edge type
% We use angle check for robustness against rounding errors
isTargetB = angle(dS.b, target_b, 'noSymmetry') < 0.1*degree;
isTargetSystem = isTargetB & dS.isEdge;

% Check if we found any systems
if ~any(isTargetSystem)
    warning(['No basal edge system with b=[', num2str(target_b.UVTW), '] found in the dislocation list.']);
    GND_specific = zeros(length(ebsdSGND), 1);
else
    % Sum rho for the specific systems (handling potential antipodal variants if unique wasn't used)
    GND_specific = factor * nansum(abs(rho(:, isTargetSystem)), 2);
end

figure;

plot(ebsdSGND, GND_specific, 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(cmapHotInverted);
caxis(plotLimits);
% Create a title with the latex string for the burgers vector
title(['Basal Edge GND Density, b=', char(target_b)]);
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_SpecificBasalEdge_GND_Hot.png']);


% --- 7. Screw <a> GND Density (Inverted Hot) ---
% Using the previously calculated 'GND_aScrew'
figure;
plot(ebsdSGND, abs(GND_aScrew), 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(flipud(hot)); % Inverted Hot
caxis(plotLimits);
title('Screw <a> GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_ScrewA_GND_HotInv.png']);


% --- 8. Prismatic <a> GND Density (Inverted Hot) ---
% Using the previously calculated 'GND_prismatic'
figure;
plot(ebsdSGND, abs(GND_prismatic), 'micronbar', 'off'); hold on;
plot(grains.boundary);
set(gca, 'ColorScale', 'log');
colormap(flipud(hot)); % Inverted Hot
caxis(plotLimits);
title('Prismatic <a> GND Density');
cb = mtexColorbar;
cb.Label.String = 'Density (m^{-2})';
exportScaledFigure(gcf, [fname, '_PrismA_GND_HotInv.png']);