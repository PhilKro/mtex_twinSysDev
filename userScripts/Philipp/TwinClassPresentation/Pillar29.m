%% Twin Analysis of Rhenium EBSD Data
% This script performs a twin analysis on a Rhenium EBSD dataset. It
% identifies a parent grain, calculates primary and secondary twin
% orientations based on the {11-21} twinning system, and visualizes the
% parent, primary, and secondary twins in the grain map.

% Clear workspace and close figures
clear all
close all

% Set MTEX preferences for axis directions
plotx2east
plotzIntoPlane

% Create output directory for figures
outDir = fullfile(pwd, 'Output_Pillar29');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% 1. Load and Pre-process Data
disp('Loading and pre-processing EBSD data...');

% Define file name and crystal symmetry for Rhenium
% Please ensure the .ang file is in the MATLAB path.
fname = 'onedriveData/Kamila/pillar29 Specimen 3 Site 1 Map Data 3_BW158.ang';
CS = crystalSymmetry('6/mmm',[2.761 2.761 4.458],'x||a','mineral','Re');

% Load EBSD data
ebsd = EBSD.load(fname, CS, 'convertEuler2SpatialReferenceFrame','setting 2');

% Clean data based on a confidence index (CI) threshold
ebsd(ebsd.ci < 0.15) = 'notIndexed';

tkd_tilt = -20 * degree;
% Perform corrections for tkd tilt and misaligned axes
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
rot = rotation.byAxisAngle(zvector,180*degree);
ebsd = rotate(ebsd,rot,'keepXY');

%% 1.5 Crop EBSD Data to Rectangular Subset
disp('Cropping EBSD data to region of interest...');
% Define a rectangular region of interest: [xmin, ymin, width, height]
% Adjust these coordinates to match your desired subset area
roi = [3.3, 5, 1.3, 1.6]; 
ebsd = ebsd(inpolygon(ebsd, roi));

% Gridify the data for faster processing and visualization
% (Note: subsetting removes the grid property, so we gridify afterwards)
ebsd = ebsd.gridify;

% some initial test plot
ipf=ipfColorKey(ebsd.CS);
ipf.inversePoleFigureDirection=yvector;
colors = ipf.orientation2color(ebsd.orientations);

figure;
plot(ebsd,colors)
exportScaledFigure(gcf, fullfile(outDir, 'ipf_map.jpg'));
%% 2. Grain Reconstruction
disp('Reconstructing grains...');

% Define parameters for grain reconstruction
grainThresh = 5*degree; 
minPixel = 5; % minimum number of pixels per grain
alpha = 1; % alpha parameter for filling during reconstruction

% Calculate grains from the indexed EBSD data
[grains, ebsd.grainId] = calcGrains(ebsd, 'angle', grainThresh, 'minPixel', minPixel, 'alpha', alpha);

% Smooth grain boundaries
grains = smooth(grains, 5);

% Plot reconstructed grains
figure;
plot(grains, grains.meanOrientation);
title('Reconstructed Grains in Rhenium Sample');
exportScaledFigure(gcf, fullfile(outDir, 'grain_map.jpg'));
%% 3. Identify Parent Grain and Define Twin Systems
disp('Identifying parent grain and defining twin systems...');

% Identify the largest grain as the parent grain for the analysis
[~, parentId] = max(grains.numPixel);
parentGrain = grains(parentId);
parentGrains = grains(angle(parentGrain.meanOrientation, grains.meanOrientation) < 8*degree);

% Create a twinSystem object for the {11-21} twinning in Rhenium
tS_1121 = twinSystem.hexagonal_1121_0001(CS);

% Symmetrise to get all unique twin variants
tS_sym = tS_1121.symmetrise;

% Generate the first layer of twins from the parent grain
grain_sym_twin = parentGrain * tS_sym;

% Generate the second layer of twins from the first layer
analysis_twin = grain_sym_twin * tS_sym;

%% 4. Twin Identification and Analysis
disp('Identifying primary and secondary twins...');

% Get orientations of the primary (first layer) twins
% implement sth like this
% primary_twins = unique(analysis_twin.parent);
% primary_twin_oris = primary_twins.orientation;
% this is the shortcut
primary_twin_oris = grain_sym_twin.orientation;

% Get orientations of the secondary (second layer) twins
secondary_twin_oris = analysis_twin.orientation;

% Define tolerance for orientation matching
tolerance = 8 * degree;

% Initialize logical arrays to classify grains
isPrimary = false(length(grains), 1);
isSecondary = false(length(grains), 1);
primaryVariantId = zeros(length(grains), 1);
secondaryVariantId = zeros(length(grains), 2);

% --- Identify Primary Twins ---
% Iterate through all grains to find primary twins
for i = 1:length(grains)
    % Skip the parent grain
    if grains(i).id == parentGrain.id, continue; end
    
    % Find the best matching primary twin variant for the current grain
    [minAngle, bestFit] = min(angle(grains(i).meanOrientation, primary_twin_oris));
    
    % If the match is within tolerance, classify as a primary twin
    if minAngle < tolerance
        isPrimary(i) = true;
        primaryVariantId(i) = grain_sym_twin(bestFit).variantId;
    end
end

% --- Identify Secondary Twins ---
% Iterate through all grains to find secondary twins among the remaining ones
for i = 1:length(grains)
    % Skip parent and primary grains
    if any(grains(i).id == parentGrains.id) || isPrimary(i), continue; end
    
    % Check if the grain orientation matches any secondary twin orientation
    [minAngle, bestFit] = min(angle(grains(i).meanOrientation, secondary_twin_oris));
    
    % If the match is within tolerance, classify as a secondary twin
    if minAngle < tolerance
        isSecondary(i) = true;
        secondaryVariantId(i,:) = [analysis_twin(bestFit).parent.variantId, analysis_twin(bestFit).variantId];
    end
end

%% 5. Visualization
disp('Visualizing the twin analysis results...');

% Create a new figure for the final plot
figure;
% Plot the base grain map with some transparency
plot(grains, grains.meanOrientation, 'facealpha', 0.3)
hold on


% Define colormap for primary twin variants
colors = lines(length(primary_twin_oris));

% Plot the parent grain in a grey tone
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
% text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

% Plot the identified primary twins with variant-specific colors
primaryGrains = grains(isPrimary);
variantIdsOfPrimaryGrains = primaryVariantId(isPrimary);
uniqueVariants = unique(variantIdsOfPrimaryGrains);

for v_idx = 1:length(uniqueVariants)
    variant = uniqueVariants(v_idx);
    if variant == 0, continue; end % Skip if variant ID is 0
    
    % Select grains belonging to the current variant
    grainsOfVariant = primaryGrains(variantIdsOfPrimaryGrains == variant);
    
    if ~isempty(grainsOfVariant)
        hold on
        plot(grainsOfVariant, 'facecolor', colors(variant,:));
        % text(grainsOfVariant, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
    end
end

% Plot the identified secondary twins in pink
secondaryGrains = grains(isSecondary);
if ~isempty(secondaryGrains)
    hold on
    plot(secondaryGrains, 'facecolor', 'm'); % 'm' for magenta/pink
    secIds = secondaryVariantId(isSecondary, :);
end

hold on
title('Primary and Secondary Twin Analysis of Rhenium');
exportScaledFigure(gcf, fullfile(outDir, 'twinAnalysis_woText.jpg'));
legend off % Hide legend for clarity

text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
for v_idx = 1:length(uniqueVariants)
    variant = uniqueVariants(v_idx);
    if variant == 0, continue; end % Skip if variant ID is 0
    
    % Select grains belonging to the current variant
    grainsOfVariant = primaryGrains(variantIdsOfPrimaryGrains == variant);
    if ~isempty(grainsOfVariant)
        text(grainsOfVariant, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
    end
end
if ~isempty(secondaryGrains)
    labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    text(secondaryGrains, labels, 'HorizontalAlignment', 'center');
end
exportScaledFigure(gcf, fullfile(outDir, 'twinAnalysis_withText.jpg'));
%% 6. Quiver K1 Plane Traces
disp('Plotting traces of K1 planes for primary and secondary twins...');

% Create a new figure for the traces
figure;
plot(grains, grains.meanOrientation, 'facealpha', 0.3)
hold on
plot(grains.boundary, 'lineWidth', 1.5);
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);

% --- Plot Primary Twin K1 Traces ---
primaryGrains = grains(isPrimary);
varIds_P = primaryVariantId(isPrimary);

for i = 1:length(primaryGrains)
    varID = varIds_P(i);
    if varID == 0, continue; end
    
    % The k1 plane in crystal coordinates
    k1_crystal = tS_sym(varID).k1;
    
    % Transform to specimen coordinates using parent orientation
    k1_specimen = parentGrain.meanOrientation * k1_crystal;
    k1_specimen.antipodal = 1;
    
    % Calculate trace direction (cross product with zvector)
    trace_dir = cross(k1_specimen, zvector);
    
    % Plot the trace inside the primary twin grain
    hold on;
    quiver(primaryGrains(i), trace_dir, 'autoScale', 'off', 'color', colors(varID,:), 'LineWidth', 2);
end

% --- Plot Secondary Twin K1 Traces ---
secondaryGrains = grains(isSecondary);
varIds_S = secondaryVariantId(isSecondary, :);

for i = 1:length(secondaryGrains)
    p_varID = varIds_S(i, 1);
    s_varID = varIds_S(i, 2);
    if p_varID == 0 || s_varID == 0, continue; end
    
    % The theoretical orientation of the primary twin (parent of this secondary)
    primaryOri = parentGrain.meanOrientation * tS_sym(p_varID).parentTwinMisorientation;
    
    % The k1 plane of the secondary twin in crystal coordinates
    k1_crystal = tS_sym(s_varID).k1;
    
    % Transform to specimen coordinates using the primary twin orientation
    k1_specimen = primaryOri * k1_crystal;
    k1_specimen.antipodal = 1;
    
    % Calculate trace direction
    trace_dir = cross(k1_specimen, zvector);
    
    % Plot the trace inside the secondary twin grain
    quiver(secondaryGrains(i), trace_dir, 'autoScale', 'off', 'color', 'm', 'LineWidth', 2);
end

title('K1 Plane Traces of Identified Twins');
exportScaledFigure(gcf, fullfile(outDir, 'twinAnalysisTraces_woText.jpg'));
text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
exportScaledFigure(gcf, fullfile(outDir, 'twinAnalysisTraces_withText.jpg'));
%% 7. GND calc

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
exportScaledFigure(gcf, fullfile(outDir, 'Curvature.jpg'));

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
exportScaledFigure(gcf, fullfile(outDir, 'Curvature_Full.jpg'));



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
exportScaledFigure(gcf, fullfile(outDir, 'totalGND_density.jpg'));





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



exportScaledFigure(gcf, fullfile(outDir, '7_Separated_GND_density.jpg'));


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

exportScaledFigure(gcf, fullfile(outDir, '7_Fully_Separated_GND_density.jpg'));