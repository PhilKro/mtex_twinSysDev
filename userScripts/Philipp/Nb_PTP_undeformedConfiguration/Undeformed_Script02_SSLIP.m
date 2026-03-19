%% Initiallize Mtex

clear
close all

plotx2east
plotzIntoPlane

plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

%% define a crystal orientation and slip systems
dir = fullfile('onedriveData','Philipp','Nb_PTP', 'Processed_Output');
matname = 'insitu Specimen 1 0 Map Data 1_BW263_20260319_144702Processing_subset=25_step=8_undeformedEBSD_withAllDIC';
load(fullfile(dir, [matname,'.mat']))

%% calc grains
if exist('grains', 'var')
    disp('Variable "grains" already exists in the workspace. Skipping grain calculation.');
else    
    thresh_ci = 0.13;
    newEBSD(newEBSD.ci < thresh_ci) = 'notIndexed';
    
    % min angle used for grain reconstruction
    grainThresh = 5*degree; %BRANDON CRITERIA, significance 15°/sqrt(Σ); Σ being the CSL number https://link.springer.com/article/10.1007/s10853-006-0665-8
    % min grainsize, in pixels
    minSize = 2000;
    % aplha parameters, defines how much to fill during recons
    % truction
    alpha = 1;
    
    [grains,newEBSD.grainId,newEBSD.mis2mean] = calcGrains(newEBSD,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
    i = 0;
    while length(grains.id)>2
        if i > 10
            error('MTEX doesnt calculate less than 3 grains. I guess you have more?')
        else
            newEBSD(grains(grains.numPixel<=minSize)) = 'notIndexed';
            [grains,newEBSD.grainId,newEBSD.mis2mean] = calcGrains(newEBSD,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
            i = i+1;
        end
    end
    
    % plot grains with grain ID numbers
    figure;
    plot(grains,grains.meanOrientation)
    hold on
    text(grains,grains.id)
    hold off
    
    save(fullfile(dir,[matname,'_withGrains.mat']), 'newEBSD', 'originalEBSD', 'grains')
    fprintf('Saved as: %s\n', [matname,'_withGrains.mat'])
end
%% plot some things to show the data

% plot effective strain map with grain boundary overlay to check alignment
figure;
plot(newEBSD,newEBSD.prop.B0007_Displacements_x_displacement)
colormap viridis
set(gca,'colorscale')

hold on
plot(grains.boundary,'linewidth',2,'linecolor','r')
text(grains,grains.id,'color','w')
hold off
% exportScaledFigure(gcf, )
%% extract the prefixes of the DIC steps
all_fields = fields(newEBSD.prop);

% 1. Filter: Find all field names containing the unique displacement keyword.
% The 'Displacements_' part is unique to the fields of interest.
displacement_fields_logical = contains(all_fields, '_Displacements_');
displacement_fields = all_fields(displacement_fields_logical);

% 2. Split: Extract the prefix by splitting the string at the first underscore '_'.
% The prefix is the part before the first '_'.
prefixes_with_duplicates = cellfun(@(x) x(1:find(x=='_', 1)-1), ...
                                   displacement_fields, ...
                                   'UniformOutput', false);

% 3. Unique & Sort: Get a unique, sorted list of prefixes.
prefix_list = unique(prefixes_with_duplicates);
clear all_fields prefixes_with_duplicates prefixes_with_duplicates
%% choose a grain and define slip systems
% grainId = 1;

% define BCC slip systems
sS = slipSystem.bcc(newEBSD.CS);

% %select only (110)
sS = sS(:);

sS = symmetrise(sS,'antipodal');

% % get local slip systems
% sSLocal = grains(grainId).meanOrientation * sS;
% 
% % extract ebsd data (which includes the disp data) for the choosen grain
% ebsd = newEBSD(grains(grainId));
% ebsd = ebsd.gridify;


%% slip ID
clear IDoptions
close all

% CHOOSE THE DIC INCREMENT FROM WHICH TO TAKE THE DISPLACEMENT DATA
% increment = 8;
increments = [7];
grains_ids = [2];
% loop through DIC increments
for i_i = 1:length(increments)
    increment = increments(i_i);
    % loop through grains
    for g_i = 1:length(grains_ids)
        clear IDoptions
        close all
        grainId = grains_ids(g_i);
        
        % get local slip systems
        sSLocal = grains(grainId).meanOrientation * sS;
        
        % constrain slip systems to those with an SF>0.2
        threshold_SF = 0.2;
        stressState = stressTensor.uniaxial(yvector);
        SF_crit = abs(sSLocal.SchmidFactor(stressState))>threshold_SF;
        if threshold_SF>0
            warning('Slip systems constrained to those with a SF higher than: %f', threshold_SF)
        end
        
        sSLocal = sSLocal(SF_crit);
        sS = sS(SF_crit);

        % extract ebsd data (which includes the disp data) for the choosen grain
        ebsd = newEBSD(grains(grainId));
        ebsd = ebsd.gridify;

        
        % define the displacement field which should be used
        U = ebsd.prop.([prefix_list{increment},'_Displacements_y_displacement']) * newEBSD.dPos; % displacement * pixel size (dPos) in um to get proper displacement field
        V = ebsd.prop.([prefix_list{increment},'_Displacements_y_displacement']) * newEBSD.dPos; % displacement * pixel size (dPos) in um to get proper displacement field
        
        % Batches for identification (name of batch and slip system numbers to use
        % for ID). Can be used to easily try some things. Not currently capable of
        % running multiple batches in a loop.
        
        % idName = name used to same figures (if IDoptions.saveFig == 1)
        % idName = name used to same figures
        idName{1} = ['Nb_bicrystal_DICincrement',prefix_list{increment}]; NoSs_batch{1} = 1:length(sS);
        idName{2} = ['Nb_bicrystal_sS_2_1_DICincrement0',prefix_list{increment}]; NoSs_batch{2} = [1 7 12 17 23 29 33 35 44 48]; % the systems which are actually active...
        
        % choose which "case/batch" to run
        Ss_batch = [1];
        
        %%%%%%% INPUT OPTIONS %%%%%%%
        %%% for more information, see the SSLIP.m file, in which most options are
        %%% explained a bit better
        
        %%%
        % set SSLIP method
        % 1: constrained and minimized slip ID (As used predominantly in the SSLIP paper. 
        % Contraint: || H^exp - H^their || < H_thresh. Minimzation of sum of absolute value of slip activities ) 
        
        % 2: constrained slip ID (only solve problem based on the || H^exp - H^their || < H_thresh  constraint, no minimization)
        
        
        % 3: single slip system ID (check for each pixel if a SINGLE system fits the activity. Used to "initialize" the SSLIP id om Figure 11 of the paper)
        % Recommended as trial for uncertain/complex situations.
        IDoptions.IDMethod = 1;
        %%%
        
        % filtering and coarse graining options (performed on displacement fields)
        IDoptions.filterSize = 4;
        IDoptions.coarsegrain = 4;
        
        % define threshold under which the "residual" should be to stop solving
        % (note that this is an "absolute" value, so in case of low strains (e.g.
        % 0.02), the default value (0.01) is rather high.
        % IDoptions.threshResidual = 0.005;
        IDoptions.threshResidual = 0.02;
        
        % set minimum strain to perform SSLIP at a datapoint (i.e. skip
        % pixels with low strain)
        IDoptions.minEeff = 0.01;
        
        % Positive constraint: choose to constrain the slip amplitudes to be positive or not. This
        % only works well if the slip systems are "configured" to have a positive
        % amplitude under a certain load (which is normally assured in the main
        % script, assuming e.g. uniaxial tension).
        % How to "reconfigure" the slip system under complex loads is T.B.D.
        IDoptions.posConstr = 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%% PLOTTING OPTIONS %%%%%%%
        
        % plot SSLIP results?
        IDoptions.plotSSLIP = 1;
        
        % max strain to plot, do not specify if no max is required (only required
        % in some cases where something weird happens)
        % IDoptions.maxE = 0.2;
        
        % layout of plotting (leave empty if unknown) [rows columns]
        IDoptions.layout = [3 4];
        
        % plot slip traces and slip directions over the activity fields
        IDoptions.plotTraces = 1;
        
        % provide stress state for SF calculation
        IDoptions.stress = stressState;
        
        % give this a name for saving
        IDoptions.casename = [matname,num2str(grainId),'_sS_batch',num2str(idName{Ss_batch}),'DIC_inc',num2str(increment)];
        
        % specify wether or not to plot defGrad fields and Eff shear strain fields,
        % before performing SSLIP
        IDoptions.plotDefGrad = 1;
        
        % save the figures?
        IDoptions.saveFig = 1;
        IDoptions.plotname = IDoptions.casename;
        
        % extra saved figure with different extension? (only used when saveFig==1)
        % IDoptions.saveExt = '.eps';
        
        % change "size" of figure, e.g. to enlarge font etc (only used for
        % individual plots)
        IDoptions.sizeAdjust = 0.5;
        
        % plot residual field after ID
        IDoptions.plotResidual = 1;
        
        % use logscale for plotting
        IDoptions.logscale = 1;
        IDoptions.logmin = 0.01;
        %IDoptions.maxE = 0.1;
        
        
        % assign colormap
        IDoptions.cmap = viridis(512);
        
        % put in a comment, will be used in filenames when saving figures
        IDoptions.comment = '';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% END OF INPUT %%%%%%%%%
        
        
        
        
        %%% set some variables
        
        %%% choose which analysis to do, based on batch number (could be used to try
        % a lot of things in a loop)
        IDoptions.NoSs = NoSs_batch{Ss_batch};
        
        
        %%% perform SSLIP analysis
        [ebsdID,optOut] = SSLIP(ebsd,U,V,sSLocal,IDoptions);
        
        %%% save the results in a matfile
        save(fullfile(dir,optOut.plotname),'ebsdID','sSLocal','optOut');
        
        SF = abs(SchmidFactor(sSLocal,IDoptions.stress));
        SF = round(SF*100)/100;
        
        is_in_batch = false(1, length(sS));
        is_in_batch(NoSs_batch{Ss_batch}) = true;
        
        T = table((1:length(sS))', arrayfun(@(x) x.b.char, sS, 'UniformOutput', false), arrayfun(@(x) x.n.char, sS, 'UniformOutput', false), SF', is_in_batch',...
                  'VariableNames', {'SlipSystem_ID', 'Burgers_Vector_b', 'Plane_Normal_n', 'Schmid_Factor', 'Contained_in_batch'});
        
        % 3. Write the Table to a CSV file.
        filename = fullfile(dir,['SSLIP_CGR_' num2str(IDoptions.coarsegrain),'_Filt_' num2str(IDoptions.filterSize), '_',matname,num2str(grainId),'_sS_batch',num2str(idName{Ss_batch}),'_SlipSystem_Data.csv']);        
        writetable(T, filename);
    end
end
% plot one activity field (system 1), just to demonstrate how data is
% structured
% figure;
% plot(ebsdID,ebsdID.prop.slipIDcor(find(optOut.NoSs==1),:))
% mtexColorbar




% %%% potentially, for replotting:
% plotSSLIP(ebsdID.prop.slipIDcor,ebsdID.prop.residualEeff,ebsdID,sSLocal,optOut)
