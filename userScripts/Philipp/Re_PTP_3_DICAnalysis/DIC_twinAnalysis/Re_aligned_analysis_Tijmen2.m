%% initialize mtex

clear
close all

% Mtex directions (important when working with EBSD)
plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

% inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3_Subset19_Step6_relative2First_firstOrderSF';
% outputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3_Subset19_Step6_relative2First_firstOrderSF\MTEX_Analysis_output';

inputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3';
outputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3\MTEX_output';
%% Load data

% tkd_tilt = -15*degree;
% thresh_ci = 0.1;
% 
% 
% % old files sum of differential
% % fnames = {fullfile(inputPath,'inc_0_newEBSD__withStrain_05-08-2025_10-27'),...
% %     fullfile(inputPath,'inc_4_newEBSD__withStrain_05-08-2025_23-58'),...
% %     fullfile(inputPath,'inc_6_newEBSD__withStrain_05-12-2025_15-07'),...
% %     fullfile(inputPath,'inc_7_newEBSD__withStrain_05-12-2025_15-16'),...
% %     fullfile(inputPath,'inc_8_newEBSD__withStrain_05-12-2025_15-23')};
% 
% % garbage strain files relative to first
% % fnames = {fullfile(inputPath,'inc_0_newEBSD__withStrain_04-14-2026_20-45'),...
% %     fullfile(inputPath,'inc_4_newEBSD__withStrain_04-14-2026_19-52'),...
% %     fullfile(inputPath,'inc_6_newEBSD__withStrain_04-14-2026_19-40'),...
% %     fullfile(inputPath,'inc_7_newEBSD__withStrain_04-14-2026_20-41'),...
% %     fullfile(inputPath,'inc_8_newEBSD__withStrain_04-14-2026_20-37')};
% 
% % relative to first subset 19 step 6
% fnames = {fullfile(inputPath,'inc_4_newEBSD__withStrain_04-15-2026_16-48'),...
%     fullfile(inputPath,'inc_6_newEBSD__withStrain_04-15-2026_17-01'),...
%     fullfile(inputPath,'inc_7_newEBSD__withStrain_04-15-2026_17-06')};
% 
% 
% %load('imported_data_converted_mtex6.1_ZnCS.mat','area');
% 
% N=length(fnames);
% for i=1:N
%     data = load(fnames{i}); % Load all variables to ensure the new file contains them
%     ebsd0{i} = data.newEBSD;
% 
%     names = fieldnames(ebsd0{i}.prop.DIC);
%     for j=1:length(names)
%         ebsd0{i}.prop.(names{j}) = ebsd0{i}.prop.DIC.(names{j});
%     end
%     ebsd0{i}.prop = rmfield(ebsd0{i}.prop, 'DIC');
% 
%     % fix mtex version difference
%     CSList = {'notIndexed',crystalSymmetry('6/mmm',[2.7 2.7 4.4],'x||a','mineral','Re')};
%     ebsd{i} = EBSD(vector3d(ebsd0{i}.x,ebsd0{i}.y,zeros(size(ebsd0{i}.x))),...
%         ebsd0{i}.orientations,...
%         ebsd0{i}.phaseId,...
%         CSList,...
%         ebsd0{i}.prop);
% 
% 
%     ebsd{i} = rotate(ebsd{i},rotation.byAxisAngle(xvector,-70*degree + tkd_tilt),'keepXY');
%     ebsd{i}(ebsd{i}.prop.ci < thresh_ci) = 'notIndexed';
%     rot = rotation.byAxisAngle(zvector,180*degree);
% 
%     ebsd{i} = rotate(ebsd{i},rot,'keepXY');
% 
%     %%% select area or use area data
%     figure;
%     plot(ebsd{i}, ebsd{i}.orientations);
%     title(['Select area for dataset ', num2str(i)]);
%     [~,area{i},~] = selectInteractive(ebsd{i});
% 
%     % areaTemp = [min(area{i}(:,1)),...
%     %             min(area{i}(:,2)),...
%     %             max(area{i}(:,1))-min(area{i}(:,1)),...
%     %             max(area{i}(:,2))-min(area{i}(:,2))];
%     ebsd{i} = ebsd{i}(inpolygon(ebsd{i},area{i}));
%     ebsd{i}(ebsd{i}.orientations.isnan) = 'notIndexed';
%     % ebsd{i}(ebsd{i}.prop.specimenlocation) = 'notIndexed';
% 
%     % num pixel, adjust for small stepsize due to alignment
%     numPix = (4*1.5)^2;
% 
%     [grains{i}, ebsd{i}.grainId] = calcGrains(ebsd{i},'alpha',3,'angle',15*degree,'minPixel',numPix);
%     % ebsd{i}(grains{i}(grains{i}.numPixel < numPix)) = 'notIndexed';
%     % [grains{i}, ebsd{i}.grainId] = calcGrains(ebsd{i},'alpha',3,'angle',15*degree);
% 
%     ebsd{i} = ebsd{i}.gridify;
% 
%     % Apply interactive mask 
%     if isfield(ebsd{i}.prop, 'ci')
%         ci_matrix = ebsd{i}.prop.ci;
%     else
%         ci_matrix = ebsd{i}.ci;
%     end
%     ci_matrix(isnan(ci_matrix)) = 0;
%     % [mask_2d, ~] = interactive_specimen_mask(ci_matrix);
%     % ebsd{i}.prop.specimenlocation = mask_2d(:);
% 
%     [~,maxId] = max(grains{1}.numPixel);
%     ori = grains{1}(maxId).meanOrientation;
% 
%     grains{i}.meanOrientation = grains{i}.meanOrientation.project2FundamentalRegion(ori);
%     ebsd{i}.orientations = ebsd{i}.orientations.project2FundamentalRegion(ori);
% 
%     % Save all fully processed relevant variables to skip this loop in the future
%     data.newEBSD = ebsd{i};
%     data.area = area{i};
%     data.grains = grains{i};
%     data.ori = ori;
%     save([fnames{i}, '_MTEX6.1.mat'], '-struct', 'data');
% end
% 
% % Save all selected areas into a separate appropriately named .mat file as well
% save(fullfile(inputPath, 'selected_areas_all.mat'), 'area');

%% --- QUICK LOAD SECTION ---
% old files sum of differential
fnames = {fullfile(inputPath,'inc_0_newEBSD__withStrain_05-08-2025_10-27'),...
    fullfile(inputPath,'inc_4_newEBSD__withStrain_05-08-2025_23-58'),...
    fullfile(inputPath,'inc_6_newEBSD__withStrain_05-12-2025_15-07'),...
    fullfile(inputPath,'inc_7_newEBSD__withStrain_05-12-2025_15-16'),...
    fullfile(inputPath,'inc_8_newEBSD__withStrain_05-12-2025_15-23')};
% garbage strain files relative to first
% fnames = {fullfile(inputPath,'inc_0_newEBSD__withStrain_04-14-2026_20-45'),...
%     fullfile(inputPath,'inc_4_newEBSD__withStrain_04-14-2026_19-52'),...
%     fullfile(inputPath,'inc_6_newEBSD__withStrain_04-14-2026_19-40'),...
%     fullfile(inputPath,'inc_7_newEBSD__withStrain_04-14-2026_20-41'),...
%     fullfile(inputPath,'inc_8_newEBSD__withStrain_04-14-2026_20-37')};
% fnames = {fullfile(inputPath,'inc_4_newEBSD__withStrain_04-15-2026_16-48'),...
%     fullfile(inputPath,'inc_6_newEBSD__withStrain_04-15-2026_17-01'),...
%     fullfile(inputPath,'inc_7_newEBSD__withStrain_04-15-2026_17-06')};
N=length(fnames);

for i = 1:N
    loaded_data = load([fnames{i}, '_MTEX6.1.mat']);
    % loaded_data = load(fnames{i});
    ebsd{i}   = loaded_data.newEBSD;
    reforientation = orientation.byEuler(0*degree,0*degree,0*degree, ebsd{i}.CS);
    ebsd{i}.orientations = ebsd{i}.orientations.project2FundamentalRegion(reforientation);
    
    
    ebsd{i}(~ebsd{i}.prop.specimenlocation) = 'notIndexed';
    % ----
    % possible grain calc
    % num pixel, adjust for small stepsize due to alignment
    numPix = (4*1.5)^2;
    [grains{i}, ebsd{i}.grainId] = calcGrains(ebsd{i},'alpha',3,'angle',15*degree,'minPixel',numPix);
    % ----
    % grains{i} = loaded_data.grains;
    % ----
    
    area{i}   = loaded_data.area;
end
% ---------------------------------

CS = ebsd{1}.CS;
clear data ebsd0


% load('imported_data_correct.mat')


%% plot ipf

lineWidth = 2;


% f=newMtexFigure('layout',[1 N]);
for i=1:N
    figure;
    % nextAxis
    % plot(ebsd{i},ebsd{i}.orientations,'ipfd',yvector)
    plot(ebsd{i},ebsd{i}.orientations,'ipfd',yvector);
    hold on
    plot(grains{i}.boundary,'LineWidth',lineWidth,'lineColor','k');
    
    saveFigure(fullfile(outputPath,['inc_',num2str(i),'_ipf.png']))
end

% saveFigure('IPFs_aligned.png')

for i=1:N
    figure;
    % nextAxis
    % plot(ebsd{i},ebsd{i}.orientations,'ipfd',yvector)
    plot(ebsd{i},ebsd{i}.ci);
    colormap gray
    hold on
    % plot(grains{i}.boundary,'LineWidth',lineWidth,'lineColor','k');
    
    saveFigure(fullfile(outputPath,['inc_',num2str(i),'_ci.png']))
end



%% recalc strains
filterSize = 3;
%%% recalculate strains
for i=1:N
    data.X = ebsd{i}.x;
    data.Y = ebsd{i}.y;
    data.U = ebsd{i}.prop.u * 0.005;
    data.V = ebsd{i}.prop.v * 0.005;

    [x,y,straindata]=calcStrainsFromRawDisp(data,filterSize);

    ebsd{i}.prop.exy2 = straindata.Exy;
    ebsd{i}.prop.exx2 = straindata.Exx;
    ebsd{i}.prop.eyy2 = straindata.Eyy;


end
clear data

%% plot overlap Exy

lineWidth = 2;
lineColor = 'w';
strainthresh = 1;

% figure;
% newMtexFigure('layout',[1 N-1])
for i=2:N
    figure;
    plotData = ebsd{i}.prop.eyy2;

    plotData(ebsd{i}.prop.geff>strainthresh) = NaN;
    plotData(ebsd{i}.prop.geff==0) = NaN;
    plot(ebsd{i},plotData)
    hold on
    % plot gray for NaNs
    ebsdNaN = ebsd{i}(isnan(plotData));
    colors = repmat([0.5 0.5 0.5],[length(ebsdNaN),1]);
    plot(ebsdNaN,colors)
    hold on
    % plot(grains{i}.boundary,'LineWidth',lineWidth,'lineColor',lineColor);
    % caxis([-0.3 0.3])
    caxis([-0.05 0.3])
    mtexColorMap("viridis")
    mtexColorbar
    saveFigure(fullfile(outputPath,['inc_',num2str(i),'_eyy2_alignedNewFilt_Filt_noBound',num2str(filterSize),'.png']))
end


%% plot eyy

lineWidth = 2;
lineColor = 'w';
strainthresh = 1;

% figure;
% newMtexFigure('layout',[1 N-1])
for i=2:N
    figure;
    plotData = ebsd{i}.prop.eyy2;

    plotData(ebsd{i}.prop.geff>strainthresh) = NaN;
    plotData(ebsd{i}.prop.geff==0) = NaN;
    plot(ebsd{i},plotData)
    hold on
    % plot gray for NaNs
    ebsdNaN = ebsd{i}(isnan(plotData));
    colors = repmat([0.5 0.5 0.5],[length(ebsdNaN),1]);
    plot(ebsdNaN,colors)
    hold on
    plot(grains{i}.boundary,'LineWidth',lineWidth,'lineColor',lineColor);
    % caxis([-0.3 0.3])
    caxis([-0.05 0.4])
    mtexColorMap("viridis")
    mtexColorbar
    saveFigure(fullfile(outputPath,['inc_',num2str(i),'_eyy2_alignedNewFilt_Filt',num2str(filterSize),'.png']))
end
