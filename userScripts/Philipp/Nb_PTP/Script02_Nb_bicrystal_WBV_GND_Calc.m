%% initialize mtex

clear 
close all
load('CS_cS_andPlanes.mat')

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector; 
%% USE THE IMPORTER_Nb_Data_TKD_DIC_SE_alignment.m TO CREATE THE ebsd object
% THE FITTED EBSD DATA IS NOT GOOD FOR GND (LIKELY ALSO WBV), ARTIFACTS
% FROM INTERPOLATION PROBABLY


% % Increment 1
% loadOptions = struct( ...
%     matName = 'insitu Specimen 1 0 Map Data 1_BW263_20251127_092547_aligned_newEBSD_withDIC.mat',...
%     saveName = 'insitu Specimen 1 0 Map Data 1_BW263')

% % Increment 2
% loadOptions = struct( ...
%     matName = 'insitu Specimen 1 1 Map Data 2_BW263_20251127_091912_aligned_newEBSD_withDIC.mat',...
%     saveName = 'insitu Specimen 1 1 Map Data 2_BW263')

% % Increment 3
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 2 Map Data 3_BW263_20251126_221700_aligned_newEBSD_withDIC.mat', ...
%     saveName = 'insitu Specimen 1 2 Map Data 3_BW263')

% % Increment 4
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 3 Map Data 4_BW263_20251127_154633_aligned_newEBSD_withDIC.mat', ...
%     saveName = 'insitu Specimen 1 3 Map Data 4_BW263')

% Increment 5
loadOptions = struct(...
    matName = 'insitu Specimen 1 4 Map Data 5_BW263__aligned_newEBSD_withDIC.mat', ...
    saveName = 'insitu Specimen 1 4 Map Data 5_BW263')

% % Increment 6
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 5 Map Data 6_BW263_20251201_105153_aligned_newEBSD_withDIC.mat', ...
%     saveName = 'insitu Specimen 1 5 Map Data 6_BW263')

fileVars = whos('-file', loadOptions.matName);
varNames = {fileVars.name};

if ~ismember('originalEBSD', varNames)
    warning('The specified MAT file does not contain ''originalEBSD''. Please run the importer script first. Without the original EBSD data, the WBV and GND calculations do not work properly.');
end

load(loadOptions.matName)


if ~isfield(newEBSD.prop,'originalId') || ~isfield(originalEBSD.prop,'originalId')
  warning('Missing ''originalId'' field in EBSD data properties. GND and WBV alignment will not work properly.');
  choice = questdlg('Continue anyway?', 'Missing originalId', 'Yes', 'No', 'Yes');
  if strcmp(choice, 'No')
    return
  end
end

if isstruct(newEBSD.prop) && isfield(newEBSD.prop,'idNearest')
    %remove my self created field bc it breaks the 
  newEBSD.prop = rmfield(newEBSD.prop,'idNearest');
end



newEBSD.prop.preAlignmentIndices = originalEBSD.id2ind(newEBSD.prop.oldId);

%%
ebsd = newEBSD;
ebsd0 = ebsd;
CS = ebsd.CS;
clear newEBSD

originalEBSD = originalEBSD.gridify;
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

figure;
plot(ebsd,ebsd.prop.ci)
title('Spherical Indexing CI')
colormap gray

nextAxis
plot(ebsd,ebsd.prop.bc)
title('Band Contrast')

nextAxis
plot(ebsd,ebsd.prop.bs)
title('Band Slope')

nextAxis
plot(ebsd,ebsd.prop.error)
title('Error')

nextAxis
plot(ebsd,ebsd.prop.mad)
title('Mean Angular Deviation')

nextAxis
plot(ebsd,ebsd.prop.bands)
title('Bands')

% exportScaledFigure(gcf, [loadOptions.matName, '_qualityMaps','.png'],'nodatestring')
%% cut to subset
% Plot the EBSD map
ipf = ipfColorKey(ebsd.CS);
figure;
ipf.inversePoleFigureDirection = xvector;
colors = ipf.orientation2color(ebsd('indexed').orientations);
plot(ebsd, ebsd.prop.SE_raw);
hold on
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
    plot(selected_ebsd, selected_ebsd.prop.SE_raw);
    hold on
    plot(selected_ebsd('indexed'), colors);
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
%% perform grain reconstruction
ipf = ipfColorKey(ebsd.CS);
% min angle used for grain reconstruction
grainThresh = 7*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 10;
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

%exportScaledFigure(gcf, [loadOptions.matName, '_Grains_crystalShapes1','.png'],'nodatestring')

figure;
plot(grains,colors)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)
exportScaledFigure(gcf, [loadOptions.matName, '_grain_ID','.png'],'nodatestring')
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

exportScaledFigure(gcf, [loadOptions.matName, '_misorientation_profile','.png'],'nodatestring')


%%% plot histogram %%%
binEdges = 0:1:90;

figure;
histogram(gB.misorientation.angle./degree, 'BinEdges', binEdges, 'Normalization', 'probability', 'FaceColor', 'r','EdgeColor', 'w','FaceAlpha', 1)
xlabel('misorientation angle (degree)')
%axis tight
pbaspect([1.618 1 1])  % Set aspect ratio (x:y:z = 1:1:1)
exportScaledFigure(gcf, [loadOptions.matName,'histogramMisori.png'],'nodatestring');
%% plot tilt (90° angle between gb normal and rotation axis) and twist (0° angle between gb normal and rotation axis) boundaries
smoothGrains = smooth(grains, 1000);
smoothgB = smoothGrains.boundary(CS.mineral,CS.mineral);
oriGB = ebsd('id',smoothgB.ebsdId).orientations;
axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');

figure;
plot(smoothGrains,'FaceAlpha',0.2)
hold on
plot(smoothgB,angle(cross(smoothgB.direction, zvector),axS)./degree,'linewidth',2)
title('Angle rotation axis to GB normal (in plane)')
mtexColorMap parula
mtexColorbar
exportScaledFigure(gcf, [loadOptions.matName, '_Angle_MisorientationAxis_GBnormal.png'],'nodatestring')

figure;
nextAxis
plot(smoothGrains,'FaceAlpha',0.2)
hold on
plot(smoothgB,angle(zvector,axS,'antipodal')./degree,'linewidth',2)
title('Angle rotation axis to out of plane direction (zvector)')
mtexColorbar
exportScaledFigure(gcf, [loadOptions.matName, '_Angle_MisorientationAxis_zvector.png'],'nodatestring')

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
%% Grain relative orientation deviation
%cutoffForGrodAxis = 1; %DOESNT work yet

ebsd.prop.grod = ebsd.calcGROD(grains)

figure;
plot(grains.boundary); hold on
plot(ebsd, ebsd.prop.grod.angle()/degree)
colormap magma
title("Grain relative orientation deviation")
mtexColorbar('title',"Deviation in °")
exportScaledFigure(gcf, [loadOptions.matName, '_GROD_angle','.png'],'nodatestring')

ipf=ipfColorKey(ebsd.CS);

figure;
plot(grains.boundary); hold on
plot(ebsd,ipf.Miller2Color(ebsd.prop.grod.axis))
colormap magma
title("Grain relative orientation deviation axis. Miller to ipf color")
exportScaledFigure(gcf, [loadOptions.matName, '_GROD_axis','.png'],'nodatestring')
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
exportScaledFigure(gcf, [loadOptions.matName, '_KAM','.png'],'nodatestring')
%% Slip system traces

sS_bbc = slipSystem.bcc(ebsd.CS);

%line [211], burger = [111], plane=(110) edge GNDs
slip_110 = sS_bbc(1).symmetrise('antipodal');

%line [101], burger = [111], plane=(211) edge GNDs
slip_211 = sS_bbc(2).symmetrise('antipodal');

%line [451], burger = [111], plane=(213) edge GNDs
slip_321 = sS_bbc(3).symmetrise('antipodal');

% chose slip System to analyse
chosensS = slip_110;

% rotate symmetrically equivalent slip systems into specimen coordinates of
% each grain
sSLocal = grains.meanOrientation * chosensS;

% compute Schmid factor according Z axis
sigma = stressTensor.uniaxial(vector3d.Y);
SF = abs(sSLocal.SchmidFactor(sigma));

% sort the slips system according to their Schmid factor from largest to
% smallest
[SF, reorder] = sort(SF,2,"descend");
reorder = sub2ind(size(reorder),repmat((1:size(reorder,1))',1,size(reorder,2)),reorder);
sSLocal = sSLocal(reorder);
% Finally, we plot the trace of the slip planes together with the slip direction in the EBSD map.

order = (1:3);
for i = 1:length(order)
    figure;
    ipf.inversePoleFigureDirection = zvector;
    plot(ebsd, ebsd.prop.Strains_MaximumShearStrain_RC)
    caxis([0,0.05])

    hold on
    plot(grains.boundary,'lineColor','r','linewidth',3)
    hold off

    hold on
    % visualize the trace of the slip plane
    quiver(grains,sSLocal(:,order(i)).trace,'color','red','linewidth',5)

    % and the slip direction
    quiver(grains,sSLocal(:,order(i)).b,'color','w','linewidth',5)
    hold on
    text(grains, num2str(SF(:,order(i))))
    title(sprintf('%s slip system with %d-th highest SF', chosensS.char{1}, order(i)));
    exportScaledFigure(gcf, sprintf('%s_%s_Rank_%d.png', loadOptions.matName, chosensS.char{1}, order(i)), 'extraheight','nodatestring')
end
%% --------- Grain calc originalEBSD for WBV and GND ----------
% WBV and GND calculations have to be done on the originalEBSD object and then can be fused into
% the newEBSD (ebsd) object based on ebsd.prop.oldID

[originalGrains,originalEBSD.grainId,originalEBSD.mis2mean] = calcGrains(originalEBSD,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

% take only boundaries between indexed grains
originalGB = originalGrains.boundary(CS.mineral,CS.mineral);

% if ~all(all(originalEBSD.id == originalEBSD.prop.originalId))
%     warning('ebsd.id has changed during grain reconstruction. I hoped this wouldnt happen I dont know if reassigning the ebsd.prop.originalId to ebsd.id breaks the grain to ebsd referencing.... But if you continue like this it wont work! Good luck!')
% end

%% --------- WBV calculation ---------
%% plot {110} plane traces

% for each ebsd data point, calculate the 110 trace
B = Miller(1,1,0,originalEBSD.CS,'hkl');

Bvec = originalEBSD.orientations .* B;%B.symmetrise('unique');
Btrace = cross(zvector,Bvec);

% define spacing for plotting (in data points)
spacing = 5;

% scaling factor for size of lines
scaling = 0.5;

ipf = ipfColorKey(CS);
ipf.inversePoleFigureDirection = xvector;
colors1 = ipf.orientation2color(originalEBSD.orientations);

figure;
plot(originalEBSD,colors1);
colormap gray
hold on

[m,n] = size(originalEBSD);
ebsdPlot = originalEBSD(1:spacing:m,1:spacing:n);
ebsdPlot = ebsdPlot(:);
BtracePlot = Btrace(1:spacing:end,1:spacing:end);
BtracePlot = BtracePlot(:);

quiver(ebsdPlot,BtracePlot,'antipodal','color','k','autoScaleFactor',scaling,'linewidth',2);

%exportScaledFigure(gcf, [loadOptions.matName, '_BasalPlaneTraces','.png'],'nodatestring')

%% plot [111] direction traces
% for each ebsd data point, calculate the basal trace
C = Miller(1,1,1,originalEBSD.CS,'uvw');

Cvec = originalEBSD.orientations .* C;


% define spacing for plotting (in data points)
spacing = 7;

% scaling factor for size of lines
scaling = 0.5;

figure;
plot(originalEBSD,originalEBSD.ci);
colormap gray
hold on

[m,n] = size(originalEBSD);
ebsdPlot = originalEBSD(1:spacing:m,1:spacing:n);
ebsdPlot = ebsdPlot(:);

CPlot = Cvec(1:spacing:end,1:spacing:end);
CPlot = CPlot(:);

quiver(ebsdPlot,CPlot,'antipodal','color','r','autoScaleFactor',scaling,'linewidth',2);
%% weighted burgers vector calculations
% calculate WBV with integral method
windowSize = 5; % some sort of smoothening parameter
W = originalEBSD.weightedBurgersVec('integral',windowSize);

% filter out data with low "dislocation amplitude"
% W(W.norm<0.01) = NaN;

Wcr = inv(originalEBSD.orientations) .* W;



ipf = ipfColorKey(originalEBSD.CS);
colors = ipf.Miller2Color(Wcr);


signW = sign(dot(W,yvector));


% transparancy scaling
transpRange = [0.0001 0.005];
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

%% setup CI thresholding
ci_tresh = 0.01;

if exist("W0",'var')
    W = W0;
    Wcr = Wcr0;
    transpField = transpField0;
    gndDensity = gndDensity0;
end

W0 = W;
Wcr0 = Wcr;
transpField0 = transpField;
gndDensity0 = gndDensity;

W(originalEBSD.prop.ci <ci_tresh)=NaN;
Wcr(originalEBSD.prop.ci <ci_tresh)=NaN;
transpField(originalEBSD.prop.ci <ci_tresh)=NaN;
gndDensity(originalEBSD.prop.ci <ci_tresh)=NaN;

%% WBV plotting
figure; %WBV intensity
% plot(originalEBSD,signW.*W.norm)
plot(originalEBSD,W.norm, 'micronbar','off')
caxis([0.001 0.1])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Intensity (\mu m^{-1})'; % Add units to the colorbar label

figure; %GND density in um-2
% plot(originalEBSD,signW.*W.norm)
plot(originalEBSD,gndDensity, 'micronbar','off')
caxis([0.001/((10^-4)*ca_b.norm) 0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (\mu m^{-2})'; % Add units to the colorbar label

figure; %GND density in m-2
% plot(originalEBSD,signW.*W.norm)
plot(originalEBSD,(10^12)*gndDensity, 'micronbar','off')
caxis([(10^12)*0.001/((10^-4)*ca_b.norm) (10^12)*0.1/((10^-4)*ca_b.norm)])
set(gca,'ColorScale','log')
colormap viridis
title('Intensity')
cb = mtexColorbar; % Get the colorbar handle
cb.Label.String = 'Density (m^{-2})'; % Add units to the colorbar label
exportScaledFigure(gcf, [loadOptions.matName,'_GNDDensity_m.png'],'nodatestring')
%% check distance to real b

% define a and ca burgers vectors
b = Miller(1,1,1,'uvw',CS);

[b_sym, a_l,~] = b.symmetrise('unique','antipodal');

b_sym.dispStyle = 'UVTW';


% calculate local b
b_sym_local = originalEBSD.orientations * b_sym;


% calculate angle between everything
b_angle = angle(W(:),b_sym_local,'antipodal');


min_b_angle = min(b_angle,[],2);



figure;

plot(originalGrains.boundary, 'micronbar','off'); hold on;
pl=plot(originalEBSD,min_b_angle/degree);
pl.AlphaDataMapping = "scaled";
pl.AlphaData = transpField;
caxis([0 45])
mtexColorbar
title(['angle to closest to ' char(round(b(1)))])

% exportScaledFigure(gcf, [loadOptions.matName, '_seperated_GND_density','.png'],'nodatestring')
exportScaledFigure(gcf, [loadOptions.matName,'_angleToB.png'],'nodatestring','extraheight')

%% Import WBV into newEBSD (ebsd) object based on originalId in both ebsd and originalEBSD
% originalEBSD.prop.WBV = W;
% originalEBSD.prop.WBV_cr = Wcr;
originalEBSD.prop.WBV_gndDens = gndDensity;


sW = size(W);
W_withNanRow = [W,repmat(vector3d.nan,sW(1),1)];
Wcr_withNanRow = [Wcr,repmat(Miller.nan(CS),sW(1),1)];
gndDensity_withNanRow = [gndDensity,NaN(sW(1),1)];
lastNanEntryId = length(W_withNanRow);

% extract the indices for a given originalId which point to the same data
% in ebsd and originalEBSD. i.e. originalEBSD(originalEBSD.prop.originalId)
% ==  ebsd(ebsd.prop.originalId) {CAN BE MULTIPLE}
listInds = originalEBSD.prop2ind('originalId',ebsd.prop.originalId);

% Where the listinds are 0 default to the nans which we have written into
% the last datapoint of the lists
listInds(listInds==0) = lastNanEntryId;

aligned_W = W_withNanRow(listInds);
aligned_Wcr = Wcr_withNanRow(listInds);
aligned_Wcr_GNDdens = gndDensity_withNanRow(listInds);


figure;
plot(ebsd, aligned_Wcr_GNDdens)
%% ----- WBV Done
%% ----- Classical GND

%% Geometrically necessary dislocations 1: Curvature, https://mtex-toolbox.github.io/GND.html

F = halfQuadraticFilter;
ebsdSGND = smooth(originalEBSD('indexed'),F,'fill',originalGrains);

% compute the curvature tensor
kappa = originalEBSD.curvature;

% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
% kappa(2,3)

newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(originalEBSD,kappa{i,j},'micronBar','off')
    hold on; plot(originalGrains.boundary,'linewidth',2); hold off

  end
end

% unify the color rage  - you may also use setColoRange equal
setColorRange equal
setColorRange([-0.2,0.2])
drawNow(gcm,'figSize','large')
mtexColorbar
exportScaledFigure(gcf, [loadOptions.matName, '_Curvature','.png'],'nodatestring')

%% Geometrically necessary dislocations 2: GND calculation
% FILTER THE DATA FIRST IF WANTED
F = halfQuadraticFilter;
ebsdSGND = smooth(originalEBSD,F,'fill',originalGrains);
%ebsdSGND = originalEBSD;

%dislocationSystem
dS = dislocationSystem.bcc(originalEBSD.CS);
%dSSymm = symmetrise(dS,'unique')
[norm(dS(1).b), norm(dS(12).b),1/3 * norm(originalEBSD.CS.aAxis)]
[norm(dS(13).b), norm(dS(end).b),1/3 * sqrt(norm(originalEBSD.CS.aAxis)^2+norm(originalEBSD.CS.cAxis)^2)] %dS(end) is a screw dislocation idk what the burger there is

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

%% Geometrically necessary dislocations 3: Complete Curvature Tensor plotting
% restore the full curvature tensor
kappa = alpha.curvature;
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(originalEBSD,kappa{i,j},'micronBar','off')
    hold on; plot(originalGrains.boundary,'linewidth',2); hold off

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
exportScaledFigure(gcf, [loadOptions.matName, '_Curvature_Full','.png'],'nodatestring')
% 


%% Geometrically necessary dislocations 4: Plot total GND density
%set limits for logarithmic scale in m^-2
limits = [14.5 16.8];

figure;
% rho .* dSRot.u SHOULD BE AN ENERGY J/m^3 NOT A Dislocation density https://github.com/mtex-toolbox/mtex/discussions/921 
originalEBSD.prop.totalDislocationEnergy = factor*sum(abs(rho .* dSRot.u),2);
originalEBSD.prop.totalGNDdensity = sum(abs(rho),2);

plot(originalEBSD, log10(originalEBSD.prop.totalGNDdensity), 'micronbar', 'off','logarithmic');
title('Total GND dentsiy')
mtexColorMap parula
exportScaledFigure(gcf, [loadOptions.matName, '_total_GND_density','.png'],'nodatestring')
%% Geometrically necessary dislocations 5: Plot GND density for the seperate dislocation systems
%calculate the screw and edge dislocation density
GND_sS = abs(rho .* dSRot.u);

%line [211], burger = [111], plane=(110) edge GNDs
Edge_211 = dS.b==Miller(1,1,1,CS)&dS.l==Miller(2,1,1,CS);
Edge_211 = Edge_211';
Edge_211 = repmat(Edge_211, length(originalEBSD), 1);

GND_211 = GND_sS;
GND_211(~(Edge_211 | dSRot.isScrew)) = nan;
GND_211 = factor*nansum(GND_211,2);

%line [101], burger = [111], plane=(211) edge GNDs
Edge_101 = dS.b==Miller(1,1,-2,0,CS)&dS.l==Miller(0,0,0,1,CS);
Edge_101 = Edge_101';
Edge_101 = repmat(Edge_101, length(originalEBSD), 1);

GND_101 = GND_sS;
GND_101(~(Edge_101 | dSRot.isScrew)) = nan;
GND_101 = factor*nansum(GND_101,2);

%line [451], burger = [111], plane=(213) edge GNDs
GND_541 = GND_sS;
GND_541((Edge_211 | dSRot.isScrew | Edge_101)) = nan;
GND_541 = factor*nansum(GND_541,2);


plot(originalEBSD, log10(originalEBSD.prop.totalGNDdensity), 'micronbar', 'off','logarithmic');
title('Total GND density')
caxis(limits);

nextAxis
plot(originalEBSD, log10(GND_211), 'micronbar', 'off','logarithmic');
title('{110}<111> GND density')
caxis(limits);

nextAxis
plot(originalEBSD, log10(GND_101), 'micronbar', 'off','logarithmic');
title('{211}<111> GND density')
caxis(limits);

nextAxis
plot(originalEBSD, log10(GND_541), 'micronbar', 'off','logarithmic');
title('{213}<111> GND dentsiy')
caxis(limits);



% set(gcf, 'Position', get(0, 'Screensize'))  % Maximize the figure window
% set(gca, 'Position', get(gca, 'OuterPosition'))  % Adjust axes to fit everything
h = mtexColorbar('title', append('log_{10}', 'rho_{GND}', '\', 'm^{-2}'));
default_position = h.Position;
% Set height of the colorbar as a percentage of the default (e.g., 70% of the default height)
percentage = 0.89;  % 70% of default height
new_height = default_position(4) * percentage;
% Update the colorbar position with the new height
h.Position = [default_position(1), default_position(2), default_position(3), new_height];
mtexColorMap parula
exportScaledFigure(gcf, [loadOptions.matName, '_seperated_GND_density','.png'],'nodatestring')
%% Geometrically necessary dislocations 5: Plot GND density edge and screw
%edge GNDs
GND_E = GND_sS;
GND_E(dSRot.isScrew) = nan;
GND_E = factor*nansum(GND_E,2);

%screw GNDs
GND_S = GND_sS;
GND_S(dSRot.isEdge) = nan;
GND_S = factor*nansum(GND_S,2);

plot(originalEBSD, log10(originalEBSD.prop.totalGNDdensity), 'micronbar', 'off','logarithmic');
title('Total GND density')
caxis(limits);
nextAxis
plot(originalEBSD, log10(GND_S), 'micronbar', 'off','logarithmic');
title('Screw GND dentsiy')
caxis(limits);
nextAxis
plot(originalEBSD, log10(GND_E), 'micronbar', 'off','logarithmic');
title('Edge GND dentsiy')
caxis(limits);

%% Import GNDs into newEBSD (ebsd) object based on originalId in both ebsd and originalEBSD 

% It is easy to plot the doubles with the same method as before
% sth like GND_sS is not that

sW = size(GND_E);
GND_E_withNanRow = [GND_E,NaN(sW(1),1)];
lastNanEntryId = length(GND_E_withNanRow);

% extract the indices for a given originalId which point to the same data
% in ebsd and originalEBSD. i.e. originalEBSD(originalEBSD.prop.originalId)
% ==  ebsd(ebsd.prop.originalId) {CAN BE MULTIPLE}
listInds = originalEBSD.prop2ind('originalId',ebsd.prop.originalId);

% Where the listinds are 0 default to the nans which we have written into
% the last datapoint of the lists
listInds(listInds==0) = lastNanEntryId;

aligned_GND_E = GND_E_withNanRow(listInds);

figure;
plot(ebsd, aligned_GND_E)


figure;
% Step 1: Convert SE_raw to grayscale RGB
SE_gray = mat2gray(ebsd.prop.SE_raw); % Normalize to [0,1]
SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image

% Step 2: Create a new MTEX figure and plot SE_raw as background
plot(ebsd, SE_rgb); 
axis image off
hold on

% Step 3: Plot Fyx with jet colormap and transparency
plot(ebsd, log10(aligned_GND_E), 'micronbar', 'off','logarithmic', 'FaceAlpha', 0.7);
colormap(viridis);
caxis(limits);
mtexColorbar('title', 'GND_E');
exportScaledFigure(gcf,[loadOptions.matName,'AlignedSE_GND_E.png'],'extraheight','nodatestring')
%% --- GND done

%% --- Save GND and WBV in .mat file
save([loadOptions.saveName, '_Alinged_EBSD_DIC_GND.mat'],'rho', 'factor', 'dSRot', 'dS', 'ebsd', 'originalEBSD', 'W', 'Wcr')

%% select grains
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
