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
%     matName = 'insitu Specimen 1 0 Map Data 1_BW263_Alinged_EBSD_DIC_GND.mat',...
%     saveName = 'insitu Specimen 1 0 Map Data 1_BW263');

% % Increment 2
% loadOptions = struct( ...
%     matName = 'insitu Specimen 1 1 Map Data 2_BW263_Alinged_EBSD_DIC_GND.mat',...
%     saveName = 'insitu Specimen 1 1 Map Data 2_BW263');

% % Increment 3
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 2 Map Data 3_BW263_20251126_221700_aligned_newEBSD_withDIC.mat', ...
%     saveName = 'insitu Specimen 1 1 Map Data 2_BW263')

% % Increment 4
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 3 Map Data 4_BW263_Alinged_EBSD_DIC_GND.mat', ...
%     saveName = 'insitu Specimen 1 3 Map Data 4_BW263',...
%     dateString = datestr(now,'yyyymmdd_HHMMSS'))

% Increment 5
loadOptions = struct(...
    matName = 'insitu Specimen 1 4 Map Data 5_BW263_Alinged_EBSD_DIC_GND.mat', ...
    saveName = 'insitu Specimen 1 4 Map Data 5_BW263',...
    dateString = datestr(now,'yyyymmdd_HHMMSS'))

% % Increment 6
% loadOptions = struct(...
%     matName = 'insitu Specimen 1 4 Map Data 5_BW263_Alinged_EBSD_DIC_GND.mat', ...
%     saveName = 'insitu Specimen 1 4 Map Data 5_BW263',...
%     dateString = datestr(now,'yyyymmdd_HHMMSS'))

fileVars = whos('-file', loadOptions.matName);
varNames = {fileVars.name};

if ~ismember('originalEBSD', varNames)
    warning('The specified MAT file does not contain ''originalEBSD''. Please run the importer script first. Without the original EBSD data, the WBV and GND calculations do not work properly.');
end

load(loadOptions.matName)
CS = ebsd.CS;
%% calcGrains for original and aligned ebsd

ipf = ipfColorKey(ebsd.CS);
% min angle used for grain reconstruction
grainThresh = 5*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
% min grainsize, in pixels
minSize = 1000;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
ebsd(grains(grains.numPixel<=minSize)) = 'notIndexed';
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

% take only boundaries between indexed grains
gB = grains.boundary(CS.mineral,CS.mineral);

[originalGrains,originalEBSD.grainId,originalEBSD.mis2mean] = calcGrains(originalEBSD,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
originalEBSD(originalGrains(originalGrains.numPixel<=minSize)) = 'notIndexed';
[originalGrains,originalEBSD.grainId,originalEBSD.mis2mean] = calcGrains(originalEBSD,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

% take only boundaries between indexed grains
originalGB = originalGrains.boundary(CS.mineral,CS.mineral);
%% Seperate GNDs

% Geometrically necessary dislocations 5: Plot GND density for the seperate dislocation systems
%calculate the screw and edge dislocation density
GND_sS = abs(rho .* dSRot.u);

%line [211], burger = [111], plane=(110) edge GNDs
Edge_110 = dS.b==Miller(1,1,1,CS)&dS.l==Miller(2,1,1,CS);
Edge_110 = Edge_110';
Edge_110 = repmat(Edge_110, length(originalEBSD), 1);

GND_110 = GND_sS;
GND_110(~(Edge_110 | dSRot.isScrew)) = nan;
GND_110 = factor*nansum(GND_110,2);

%line [101], burger = [111], plane=(211) edge GNDs
Edge_211 = dS.b==Miller(1,1,-2,0,CS)&dS.l==Miller(0,0,0,1,CS);
Edge_211 = Edge_211';
Edge_211 = repmat(Edge_211, length(originalEBSD), 1);

GND_211 = GND_sS;
GND_211(~(Edge_211 | dSRot.isScrew)) = nan;
GND_211 = factor*nansum(GND_211,2);

%line [451], burger = [111], plane=(213) edge GNDs
GND_123 = GND_sS;
GND_123((Edge_110 | dSRot.isScrew | Edge_211)) = nan;
GND_123 = factor*nansum(GND_123,2);


%edge GNDs
GND_E = GND_sS;
GND_E(dSRot.isScrew) = nan;
GND_E = factor*nansum(GND_E,2);

%screw GNDs
GND_S = GND_sS;
GND_S(dSRot.isEdge) = nan;
GND_S = factor*nansum(GND_S,2);

%% plot GNDS
plotGNDQuantities = {'GND_E','GND_S','GND_110','GND_211','GND_123'};

% limits = [14.5 16.8];
limits = [15 16.4];

for i = 1:length(plotGNDQuantities)
    quant = plotGNDQuantities{i};
    figure;
    % Step 1: Convert SE_raw to grayscale RGB
    SE_gray = mat2gray(ebsd.prop.SE_raw); % Normalize to [0,1]
    SE_rgb = repmat(SE_gray, 1, 1, 3);       % Convert to RGB image
    
    % Step 2: Create a new MTEX figure and plot SE_raw as background
    plot(ebsd, SE_rgb); 
    axis image off
    hold on
    
    % maxStrainTest = ebsd.prop.Strains_MaximumShearStrain_RC;
    % maxStrainTest(ebsd.prop.Others_CorrelationValue>0.9) = 0;
    
    % Step 3: Plot Fyx with jet colormap and transparency
    plot(ebsd, log10(transformDataOriginal2Aligned(eval(quant),originalEBSD,ebsd)), 'FaceAlpha', 0.7);
    colormap(viridis);
    caxis(limits);
    mtexColorbar('title', ['GND density', quant]);
    exportScaledFigure(gcf,[loadOptions.saveName,loadOptions.dateString,quant,'GNDs.png'],'extraheight','nodatestring')
end
%% total seperation of GNDs

for i = 1:14 %size(GND_sS,2)
    current_b = dS(i).b;
    current_l = dS(i).l;
    current_GND = GND_sS(:,i);
    
    figure;
    plot(originalEBSD,log10(current_GND))
    mtexColorbar
    title(['GNDs b:', char(round(current_b)), 'and l:', char(round(current_l))])
end

%% GND densities for given slip systems

dS_Edge = dS(dS.isEdge);
equiv_sS = slipSystem(dS_Edge.b,cross(dS_Edge.b,dS_Edge.l)); %slipSystem(b,n)GND_E = GND_sS;
GND_Edge = GND_sS(:,dS.isEdge);

%have to be the same size
assert(all(size(equiv_sS)==size(dS_Edge)) && size(GND_Edge,2)==size(dS_Edge,1))

% rotate slip systems that correspond to the dislocation systems into specimen coordinates of
% each grain
sSLocal = originalGrains.meanOrientation * equiv_sS;

% compute Schmid factor according Z axis
sigma = stressTensor.uniaxial(vector3d.Y);
SF = abs(sSLocal.SchmidFactor(sigma));

% sort the slips systems/GNDs/dS_Edges according to the Schmid factor on the slip system from largest to
% smallest
[~, reorder] = sort(SF,2,"descend");
reorder = sub2ind(size(reorder),repmat((1:size(reorder,1))',1,size(reorder,2)),reorder);
%sSLocal = sSLocal(reorder);

order = (1:5);
indices = zeros(length(order),1);
for j = 1:length(order)
    [~, indices(j)] = find(ismember(reorder,order(j)));
end


for i = 1:length(order)
    figure;
    plot(originalEBSD, log10(factor*GND_Edge(:,indices(i))))
    mtexColorbar
    caxis(limits)
    hold on
    plot(originalGrains.boundary,'lineColor','r','linewidth',3)

    % visualize the trace of the slip plane
    quiver(originalGrains,sSLocal(:,indices(i)).trace,'color','red','linewidth',5)

    % and the slip direction
    quiver(originalGrains,sSLocal(:,indices(i)).b,'color','w','linewidth',5)
    hold on
    text(originalGrains, num2str(SF(:,indices(i))))
    title(sprintf('%s slip system with %d-th highest SF', chosensS.char{i}, order(i)));
    exportScaledFigure(gcf, sprintf('%s_%s_Rank_%d.png', loadOptions.matName, chosensS.char{i}, order(i)), 'extraheight','nodatestring')
end

%% Slip system traces

sS_bbc = slipSystem.bcc(ebsd.CS);

%line [211], burger = [111], plane=(110) edge GNDs
slip_110 = sS_bbc(1).symmetrise('antipodal');

%line [101], burger = [111], plane=(211) edge GNDs
slip_211 = sS_bbc(2).symmetrise('antipodal');

%line [451], burger = [111], plane=(213) edge GNDs
slip_321 = sS_bbc(3).symmetrise('antipodal');

% chose slip System to analyse
chosensS = slip_321;

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
sSLocalSorted = sSLocal(reorder);
% Finally, we plot the trace of the slip planes together with the slip direction in the EBSD map.

order = (1:5);
for i = 1:length(order)
    figure;
    ipf.inversePoleFigureDirection = zvector;
    plot(ebsd, ebsd.prop.Strains_MaximumShearStrain_RC)
    caxis([0,0.1])

    hold on
    plot(grains.boundary,'lineColor','r','linewidth',3)
    hold off

    hold on
    % visualize the trace of the slip plane
    quiver(grains,sSLocalSorted(:,order(i)).trace,'color','red','linewidth',5)

    % and the slip direction
    quiver(grains,sSLocalSorted(:,order(i)).b,'color','w','linewidth',5)
    hold on
    text(grains, num2str(SF(:,order(i))))
    title(sprintf('%s slip system with %d-th highest SF', chosensS.char{1}, order(i)));
    exportScaledFigure(gcf, sprintf('%s_%s_Rank_%d.png', loadOptions.matName, chosensS.char{1}, order(i)), 'extraheight','nodatestring')
end

figure;
plot(grains)
hold on
% plot boundaries
plot(grains.boundary)
text(grains,grains.id)
%% GrainB

% active slip [grainId1, activeSlip_index_in_sSLocal1; grainId2, activeSlip_index_in_sSLocal2]


smoothGrains = smooth(grains, 1000);
smoothgB = smoothGrains.boundary(CS.mineral,CS.mineral);
oriGB = ebsd('id',smoothgB.ebsdId).orientations;
axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');


average_gB_Dir = mean(smoothGrains.boundary.direction);
average_gB_Normal = cross(average_gB_Dir,zvector);

figure;
plot(smoothGrains,'FaceAlpha',0.2)
hold on
quiver(grains, average_gB_Dir, 'r');
quiver(grains, average_gB_Normal, 'b');
title('avg. gB-normal in blue and gB-dir in red')
exportScaledFigure(gcf, [loadOptions.matName, 'avg_gB_normalAndDir.png'],'nodatestring')

%% m' calc of slip transfer, angle between planes

% % active slip [grainId1, activeSlip_index_in_sSLocal1; grainId2, activeSlip_index_in_sSLocal2]
% 
% % insitu Specimen 1 4 Map Data 5_BW263_Alinged_EBSD_DIC_GND
% activeSlip = [1, 1; 2, 3];
% % activeSlip = [1, 3; 2, 4];
% 
% mprime = @(sSLocal_1, sSLocal_2) cos(angle(sSLocal_1.b,sSLocal_2.b))*cos(angle(sSLocal_1.n,sSLocal_2.n));
% 
% % https://www.sciencedirect.com/science/article/pii/S1359028614000205
% lrb_crit = @(sSLocal_1, sSLocal_2, gBnormal) cos(angle(sSLocal_1.b,sSLocal_2.b)) * cos(angle(cross(gBnormal,sSLocal_1.n),cross(gBnormal,sSLocal_2.n)));
% 
% sSLocal_grain1 = sSLocal(activeSlip(1,1),:);
% sSLocal_grain2 = sSLocal(activeSlip(2,1),:);
% 
% grain1 = grains(activeSlip(1,1));
% grain2 = grains(activeSlip(2,1));
% 
% transmittingSs = [sSLocal_grain1(activeSlip(1,2)),sSLocal_grain2(activeSlip(2,2))];
% correspondingGrains = [grain1,grain2];
% 
% figure;
% ipf.inversePoleFigureDirection = zvector;
% plot(ebsd, ebsd.prop.Strains_MaximumShearStrain_RC)
% caxis([0,0.05])
% 
% hold on
% plot(grains.boundary,'lineColor','r','linewidth',3)
% hold off
% 
% hold on
% % visualize the trace of the slip plane
% quiver(correspondingGrains(1),transmittingSs(1).trace,'color','red','linewidth',5)
% quiver(correspondingGrains(2),transmittingSs(2).trace,'color','red','linewidth',5)
% 
% % and the slip direction
% quiver(correspondingGrains(1),transmittingSs(1).b,'color','w','linewidth',5)
% quiver(correspondingGrains(2),transmittingSs(2).b,'color','w','linewidth',5)
% 
% hold on
% 
% txt = sprintf('LRB-Criterion: %s \n mPrime: %s', num2str(round(lrb_crit(transmittingSs(1),transmittingSs(2),average_gB_Normal),2)) ,num2str(round(mprime(transmittingSs(1),transmittingSs(2)),2)));
% 
% text(grains, txt)
% 
% exportScaledFigure(gcf, sprintf('%s_%s_Rank_%d.png', loadOptions.matName, chosensS.char{1}, order(i)), 'extraheight','nodatestring')

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
sSLocalSorted = sSLocal(reorder);

% m' calc of slip transfer, angle between planes
% Define active slip systems as input: [grainId1, ssIndex1; grainId2, ssIndex2; ...]
activeSlips = [1, 1; 2, 2];

% Define criteria functions
mprime = @(sSLocal_1, sSLocal_2) abs(cos(angle(sSLocal_1.b,sSLocal_2.b))*cos(angle(sSLocal_1.n,sSLocal_2.n)));
lrb_crit = @(sSLocal_1, sSLocal_2, gBnormal) abs(cos(angle(sSLocal_1.b,sSLocal_2.b)) * cos(angle(cross(gBnormal,sSLocal_1.n),cross(gBnormal,sSLocal_2.n))));

% Loop through all slip system pairs
for pairIdx = 1:size(activeSlips,1)-1
    grain1Id = activeSlips(pairIdx,1);
    ss1Idx = activeSlips(pairIdx,2);
    grain2Id = activeSlips(pairIdx+1,1);
    ss2Idx = activeSlips(pairIdx+1,2);
    
    % Get slip systems and grains
    sSLocal_1 = sSLocalSorted(grain1Id, ss1Idx);
    sSLocal_2 = sSLocalSorted(grain2Id, ss2Idx);
    grain1 = grains(grain1Id);
    grain2 = grains(grain2Id);
    SF1 = SF(grain1Id, ss1Idx);
    SF2 = SF(grain2Id, ss2Idx);

    % Calculate criteria
    m_prime_val = mprime(sSLocal_1, sSLocal_2);
    lrb_val = lrb_crit(sSLocal_1, sSLocal_2, average_gB_Normal);
    
    % Plot
    figure;
    ipf.inversePoleFigureDirection = zvector;
    plot(ebsd, ebsd.prop.Strains_MaximumShearStrain_RC)
    caxis([0,0.1])
    hold on
    plot(grains.boundary,'lineColor','r','linewidth',3)
    
    % Visualize slip traces and directions
    quiver(grain1, sSLocal_1.trace, 'color', 'red', 'linewidth', 5)
    quiver(grain2, sSLocal_2.trace, 'color', 'red', 'linewidth', 5)
    quiver(grain1, sSLocal_1.b, 'color', 'w', 'linewidth', 5)
    quiver(grain2, sSLocal_2.b, 'color', 'w', 'linewidth', 5)
    
    txt = sprintf('Grain %d→%d \n LRB: %.2f | m'': %.2f \n SF grain %d: %.2f (%dth rank) \n SF grain %d: %.2f (%dth rank)', grain1Id, grain2Id, lrb_val, m_prime_val, grain1Id, SF1, ss1Idx, grain2Id, SF2,ss2Idx);
    text(grains, txt)
    
    title(txt)
    exportScaledFigure(gcf, sprintf('%s_SlipTransfer_G%d_rank%d_G%d_rank%d.png', loadOptions.saveName, grain1Id,ss1Idx, grain2Id, ss2Idx), 'extraheight','nodatestring')
end
