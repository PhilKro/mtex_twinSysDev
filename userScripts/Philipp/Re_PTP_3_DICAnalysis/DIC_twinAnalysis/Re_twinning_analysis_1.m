
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% twinning analysis Re
%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% load('imported_data_correct.mat')
close all

plotx2east
plotzIntoPlane


inc = 5;
CS = grains{inc}.CS;

comment = '_filt3';
% outputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\MTEX_Analysis_output';
outputPath = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Re_PTP\DIC_Aligned_Data\data_Re_7deg_PTP3\MTEX_output';
%% might be obsolete but run if the entire ebsd map appears as indexed
ebsd{inc}(ebsd{inc}.orientations.isnan) = 'notIndexed';
% ebsd{inc}(~ebsd{inc}.prop.specimenlocation) = 'notIndexed';
reforientation = orientation.byEuler(0*degree,0*degree,0*degree, ebsd{inc}.CS);
ebsd{inc}.orientations = ebsd{inc}.orientations.project2FundamentalRegion(reforientation);


numPix = (4*1.5)^2;
[grains{inc}, ebsd{inc}.grainId] = calcGrains(ebsd{inc},'alpha',3,'angle',15*degree,'minPixel',numPix);
%% plot grains
figure;
ipf=ipfColorKey(ebsd{inc}.CS);
ipf.inversePoleFigureDirection = yvector;
colors = ipf.orientation2color(grains{inc}.meanOrientation);
plot(grains{inc},colors)
% create and plot crystal shapes
cS1 = crystalShape.hex(ebsd{inc}.CS);
hold on
plot(grains{inc},1.5*cS1,'facecolor','b','faceAlpha',1,'linewidth',2.5)
plot(grains{inc}.boundary)


%% 3. Identify Parent Grain and Define Twin Systems
disp('Identifying parent grain and defining twin systems...');

% Identify the largest grain as the parent grain for the analysis
[~, parentId] = max(grains{inc}.numPixel);
parentGrain = grains{inc}(parentId);
parentGrains = grains{inc}(angle(parentGrain.meanOrientation, grains{inc}.meanOrientation) < 8*degree);

% Create a twinSystem object for the {11-21} twinning in Rhenium
tS_1121 = twinSystem.hexagonal_1121_0001(CS);

% Symmetrise to get all unique twin variants
tS_sym = tS_1121.symmetrise;

% Generate the first layer of twins from the parent grain
primaryTwin = parentGrain * tS_sym;

% Generate the second layer of twins from the first layer
secondTwin = primaryTwin * tS_sym;

%% 4. Twin Identification and Analysis
disp('Identifying primary and secondary twins...');

% Get orientations of the primary (first layer) twins
% implement sth like this
% primary_twins = unique(analysis_twin.parent);
% primary_twin_oris = primary_twins.orientation;
% this is the shortcut
primary_twin_oris = primaryTwin.orientation;

% Get orientations of the secondary (second layer) twins
secondary_twin_oris = secondTwin.orientation;

% Define tolerance for orientation matching
tolerance = 8 * degree;

% Initialize logical arrays to classify grains
isPrimary = false(length(grains{inc}), 1);
isSecondary = false(length(grains{inc}), 1);
primaryVariantId = zeros(length(grains{inc}), 1);
secondaryVariantId = zeros(length(grains{inc}), 2);

% --- Identify Primary Twins ---
% Iterate through all grains to find primary twins
for i = 1:length(grains{inc})
    % Skip the parent grain
    if grains{inc}(i).id == parentGrain.id, continue; end
    
    % Find the best matching primary twin variant for the current grain
    [minAngle, bestFitPrimary] = min(angle(grains{inc}(i).meanOrientation, primary_twin_oris));
    
    % If the match is within tolerance, classify as a primary twin
    if minAngle < tolerance
        isPrimary(i) = true;
        primaryVariantId(i) = primaryTwin(bestFitPrimary).variantId;
    end
end

% --- Identify Secondary Twins ---
% Iterate through all grains to find secondary twins among the remaining ones
secondaryVariantIdTotal = [];
for i = 1:length(grains{inc})
    % Skip parent and primary grains
    if any(grains{inc}(i).id == parentGrains.id) || isPrimary(i), continue; end
    
    % Check if the grain orientation matches any secondary twin orientation
    [minAngle, bestFitPrimary] = min(angle(grains{inc}(i).meanOrientation, secondary_twin_oris));
    
    % If the match is within tolerance, classify as a secondary twin
    if minAngle < tolerance
        isSecondary(i) = true;
        secondaryVariantId(i,:) = [secondTwin(bestFitPrimary).parent.variantId, secondTwin(bestFitPrimary).variantId];
        
        secondaryVariantIdTotal(end+1) = bestFitPrimary;
        
    end
end

%% 5. Visualization
disp('Visualizing the twin analysis results...');

% Create a new figure for the final plot
figure;
% Plot the base grain map with some transparency
plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.3)
title('Primary and Secondary Twin Analysis of Rhenium');
hold on


% Define colormap for primary twin variants
colors = lines(length(primary_twin_oris));

% Plot the parent grain in a grey tone
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

% Plot the identified primary twins with variant-specific colors
primaryGrains = grains{inc}(isPrimary);
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
        text(grainsOfVariant, ['PV' num2str(variant)], 'HorizontalAlignment', 'center');
    end
end

% Plot the identified secondary twins in pink
secondaryGrains = grains{inc}(isSecondary);
if ~isempty(secondaryGrains)
    hold on
    plot(secondaryGrains, 'facecolor', 'm'); % 'm' for magenta/pink
    secIds = secondaryVariantId(isSecondary, :);
    labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    text(secondaryGrains, labels, 'HorizontalAlignment', 'center');
end

hold off

legend off % Hide legend for clarity

saveFigure(fullfile(outputPath,['inc_',num2str(inc),'_twinId.png']))


%% plot k1 traces

% calculate k1 traces of primary twins
pT_k1 = primaryTwin.parent.meanOrientation .* primaryTwin.k1;
pT_k1_vars = pT_k1(variantIdsOfPrimaryGrains);

pT_k1_traces = cross(pT_k1_vars,zvector);
pT_k1_traces.antipodal = true;


% calculate k1 traces of secondary twins
sT_k1 = secondTwin.parent.orientation .* secondTwin.k1;
sT_k1_vars = sT_k1(secondaryVariantIdTotal);

sT_k1_traces = cross(sT_k1_vars,zvector);
sT_k1_traces.antipodal = true;



figure;
% Plot the base grain map with some transparency
plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.3)
% title('Primary and Secondary Twin Analysis of Rhenium');
hold on


% Define colormap for primary twin variants
colors = lines(length(primary_twin_oris));

% Plot the parent grain in a grey tone
plot(parentGrains, 'facecolor', [0.8 0.8 0.8]);
% text(parentGrains, 'P', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

% Plot the identified primary twins with variant-specific colors
primaryGrains = grains{inc}(isPrimary);
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

% plot traces of primary twins
hold on
quiver(primaryGrains,pT_k1_traces,'color','k')


% Plot the identified secondary twins in pink
secondaryGrains = grains{inc}(isSecondary);
if ~isempty(secondaryGrains)
    hold on
    plot(secondaryGrains, 'facecolor', 'm'); % 'm' for magenta/pink
    secIds = secondaryVariantId(isSecondary, :);
    % labels = arrayfun(@(x) sprintf('T%d,%d', secIds(x,1), secIds(x,2)), 1:size(secIds,1), 'UniformOutput', false);
    % text(secondaryGrains, labels, 'HorizontalAlignment', 'center');
end

hold on
quiver(secondaryGrains,sT_k1_traces,'color','k')

hold off

legend off % Hide legend for clarity

saveFigure(fullfile(outputPath,['inc_',num2str(inc),'_twin_k1.png']))


%% try to plot some deformation stuff


H = primaryTwin.displacementGradient;

Hvars = H(variantIdsOfPrimaryGrains);
F = Hvars + tensor.eye(size(Hvars));
F = deformationGradientTensor(F);

E = GreenStrain(F);


Hxx = squeeze(Hvars.M(1,1,:));
Hxy = squeeze(Hvars.M(1,2,:));
Hyx = squeeze(Hvars.M(2,1,:));
Hyy = squeeze(Hvars.M(2,2,:));

Exx = squeeze(E.M(1,1,:));
Exy = squeeze(E.M(1,2,:));
Eyy = squeeze(E.M(2,2,:));

cl = [-0.15 0.15];

figure;
newMtexFigure('layout',[1,3])
% Plot the base grain map with some transparency
plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
hold on
plot(primaryGrains,Exx);
caxis(cl)
title('E_{xx}')

nextAxis
plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
hold on
plot(primaryGrains,Exy);
caxis(cl)
title('E_{xy}')

% 
% nextAxis
% plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
% hold on
% plot(primaryGrains,Eyx);
% caxis(cl)

nextAxis
plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
hold on
plot(primaryGrains,Eyy);
caxis(cl)
title('E_{yy}')

colormap viridis
mtexColorbar

hold on

saveFigure(fullfile(outputPath,['inc_',num2str(inc),'_strain_twins.png']))


%%%% try variant identification from strains

%% primary variants considered
% pVars = [1 2 3];
% pVars = [1,2,4];
pVars = 1:6;
N = numel(pVars);

F_id = H(pVars) + tensor.eye(size(H(pVars)));
F_id = deformationGradientTensor(F_id);
% consider Green strain tensor to disregard rigid rotations in the
% comparison
E_id = F_id.GreenStrain;

E_id2d = E_id.M(1:2,1:2,:);

try
    E_exp = zeros(2,2,numel(ebsd{inc}.prop.exx2));
    E_exp(1,1,:) = ebsd{inc}.prop.exx2(:);
    E_exp(2,1,:) = ebsd{inc}.prop.exy2(:);
    E_exp(1,2,:) = ebsd{inc}.prop.exy2(:);
    E_exp(2,2,:) = ebsd{inc}.prop.eyy2(:);
    
    DefDiff = zeros(2,2,numel(ebsd{inc}.prop.exx2),N);
    NormDiff = zeros(numel(ebsd{inc}.prop.exx2),N);
catch
    warning('Using exx/exy/eyy instead of exx2/exy2/eyy2')
    E_exp = zeros(2,2,numel(ebsd{inc}.prop.exx));
    E_exp(1,1,:) = ebsd{inc}.prop.exx(:);
    E_exp(2,1,:) = ebsd{inc}.prop.exy(:);
    E_exp(1,2,:) = ebsd{inc}.prop.exy(:);
    E_exp(2,2,:) = ebsd{inc}.prop.eyy(:);
    
    DefDiff = zeros(2,2,numel(ebsd{inc}.prop.exx),N);
    NormDiff = zeros(numel(ebsd{inc}.prop.exx),N);
end

% set NaNs on bad strain data based on effective shear strain
try
    nanField = ebsd{inc}.prop.geff > 2 | ebsd{inc}.prop.geff == 0;
catch ME
    if any(strcmp(ME.identifier, {'MATLAB:noSuchMethodOrField', 'MATLAB:nonExistentField'}))
        warning('prop.geff not found. Calculating on the fly.');
        ebsd{inc}.prop.geff = sqrt(((ebsd{inc}.prop.Fxx - ebsd{inc}.prop.Fyy)/2).^2 + ((ebsd{inc}.prop.Fxy - ebsd{inc}.prop.Fyx)/2).^2);
        nanField = ebsd{inc}.prop.geff > 2 | ebsd{inc}.prop.geff == 0;
    else
        rethrow(ME);
    end
end

for i=1:N
    DefDiff(:,:,:,i) = E_exp - E_id2d(:,:,i);
    NormDiff(:,i) = vecnorm(reshape(squeeze(DefDiff(:,:,:,i)), [], size(squeeze(DefDiff(:,:,:,i)),3)));
    NormDiff(nanField(:),i) = NaN;
end


%%% identify primary twins based on difference
% norm threshold
threshP = 0.2;

% identify best fitting variant
[fitPrimary,bestFitPrimary] = min(NormDiff,[],2);

% create color array based on bestFit and threshold
plotColor = ones(length(bestFitPrimary),3);
plotColor(fitPrimary<threshP,:) = colors(bestFitPrimary(fitPrimary<threshP),:);


%%% plot the differences to pTwins
cl = [0 1];
figure;
newMtexFigure('layout',[1,N])
% Plot the base grain map with some transparency
% plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
pT_k1_trace = cross(pT_k1,zvector);
pT_k1_trace.antipodal = true;
for i=1:N
    nextAxis
    plot(ebsd{inc},NormDiff(:,i));
    caxis(cl)
    title(['Difference to pV',num2str(pVars(i))]);
    hold on
    quiver(parentGrain, pT_k1_trace(pVars(i)), 'color','w')
end
colormap viridis
mtexColorbar

saveFigure(fullfile(outputPath,['inc_',num2str(inc),'_normRed',comment,'.png']))


figure;
plot(ebsd{inc},plotColor)
colormap(colors)
caxis([1 N])
mtexColorbar
saveFigure(fullfile(outputPath,['inc_',num2str(inc),'_bestFitPrimaryRed',comment,'.png']))




%% secondary variants considered
% 
% % sVars = [2:6,7:12,14:18];
% sVars = 7:18;
% sVars([2,9])=[];
% N = numel(sVars);
% 
% % calc primary twin def
% parentTwins = secondTwin(sVars).parent;
% Hp = tensor.zeros(N,1);
% for i=1:length(parentTwins); Hp(i) = parentTwins(i).displacementGradient; end
% 
% % calc secondary twin def
% Hs = secondTwin.displacementGradient;
% Hs = Hs(sVars);
% 
% % calculate separate def grad tensors before multiplying for proper solid
% % mechanics shit
% F_idP = deformationGradientTensor( Hp + tensor.eye(size(Hp)) );
% F_idS = deformationGradientTensor( Hs + tensor.eye(size(Hs)) );
% 
% F_id = deformationGradientTensor (mtimes(F_idS,F_idP) );
% 
% E_id = F_id.GreenStrain;
% 
% E_id2d = E_id.M(1:2,1:2,:);
% 
% E_exp = zeros(2,2,numel(ebsd{inc}.prop.exx2));
% E_exp(1,1,:) = ebsd{inc}.prop.exx2(:);
% E_exp(2,1,:) = ebsd{inc}.prop.exy2(:);
% E_exp(1,2,:) = ebsd{inc}.prop.exy2(:);
% E_exp(2,2,:) = ebsd{inc}.prop.eyy2(:);
% 
% DefDiff = zeros(2,2,numel(ebsd{inc}.prop.exx2),N);
% NormDiff = zeros(numel(ebsd{inc}.prop.exx2),N);
% 
% % set NaNs on bad strain data
% nanField = ebsd{inc}.prop.geff > 2 | ebsd{inc}.prop.geff==0;
% 
% for i=1:N
%     DefDiff(:,:,:,i) = E_exp - E_id2d(:,:,i);
%     NormDiff(:,i) = vecnorm(reshape(squeeze(DefDiff(:,:,:,i)), [], size(squeeze(DefDiff(:,:,:,i)),3)));
%     NormDiff(nanField(:),i) = NaN;
% end
% 
% 
% 
% % cl = [0 1];
% % figure;
% % newMtexFigure('layout',[3,6])
% % % Plot the base grain map with some transparency
% % % plot(grains{inc}, grains{inc}.meanOrientation, 'facealpha', 0.1)
% % for i=1:N
% %     nextAxis
% %     plot(ebsd{inc},NormDiff(:,i));
% %     caxis(cl)
% %     title(['Difference to sV',num2str(sVars(i))]);
% % 
% % end
% % colormap viridis
% % mtexColorbar
% % 
% % saveFigure(['inc_',num2str(inc),'_normSecond_FpFs',comment,'.png'])
% 
% 
% 
% 
% %%% identify secondary twins based on difference, only if not primary
% %%% found!
% 
% % norm threshold
% threshS = 0.2;
% 
% % identify best fitting variant
% [fitSecond,bestFitSecond] = min(NormDiff,[],2);
% 
% % create color array based on bestFit and threshold
% plotColorS = plotColor;
% plotColorS(fitSecond<threshS & fitPrimary>threshP,:) = repmat([1 0 1],[sum(fitSecond<threshS & fitPrimary>threshP),1]);
% 
% figure;
% plot(ebsd{inc},plotColorS)
% colormap(colors)
% caxis([1 6])
% mtexColorbar
% saveFigure(['inc_',num2str(inc),'_bestFitSecond_Reduced',comment,'.png'])

%% Tijmens kitchen

% v1 = 3
% v2 = 6
% point = 404756
% Ev1 = E_id2d(:,:,v1)
% Ev2 = E_id2d(:,:,v2)
% Ev1 = Ev1(:)
% Ev2 = Ev2(:)
% Ee = E_exp(:,:,point)
% Ee = Ee(:)
% % still including strain magnitude
% x1 = Ev1\Ee
% x2 = Ev2\Ee
% % this is a shape parameter
% norm(Ee-x1*Ev1)
% norm(Ee-x2*Ev2)

%% 0. Explicitly Define Primary Variants
disp('Re-initializing Primary Twin Tensors...');
pVars = 1:6;
% pVars = [1,2,4];
N = numel(pVars);

% Recalculate E_id2d specifically for the 6 primary variants 
% so it doesn't accidentally use the 10 secondary variants
F_id = H(pVars) + tensor.eye(size(H(pVars)));
F_id = deformationGradientTensor(F_id);
E_id = F_id.GreenStrain;
E_id2d = E_id.M(1:2,1:2,:);

disp('Checking similarity of 2D projections of the theoretical twin shear tensors...');
E_check2d = reshape(E_id.M(1:2, 1:2, :), 4, []);
E_check2d_norm = E_check2d ./ vecnorm(E_check2d, 2, 1);

simMatrix = E_check2d_norm' * E_check2d_norm;
simThresh = 0.95;
[simRow, simCol] = find(abs(triu(simMatrix, 1)) > simThresh);

if ~isempty(simRow)
    warning('Very similar theoretical twin shear tensors detected (similarity > %.2f)! This might cause variant misidentification.', simThresh);
    for k = 1:length(simRow)
        fprintf(' -> Variant %d and Variant %d: similarity = %.3f\n', ...
            pVars(simRow(k)), pVars(simCol(k)), simMatrix(simRow(k), simCol(k)));
    end
else
    disp('All theoretical twin shear tensors are sufficiently distinct in shape.');
end

% Make sure E_exp is defined (from earlier in your script)
E_exp = zeros(2,2,numel(ebsd{inc}.prop.exx2));
E_exp(1,1,:) = ebsd{inc}.prop.exx2(:);
E_exp(2,1,:) = ebsd{inc}.prop.exy2(:);
E_exp(1,2,:) = ebsd{inc}.prop.exy2(:);
E_exp(2,2,:) = ebsd{inc}.prop.eyy2(:);

%% 1. Map-wide Shape and Magnitude Analysis
disp('Calculating shape mismatch and magnitude across the entire map...');

% Define cutoff criteria
magCutoff = 0.27; 
shapeCutoff = 0.5; 

% Flatten the experimental strain tensors
% E_exp is size [2, 2, nPoints], reshape it to [4, nPoints]
try
    nPoints = numel(ebsd{inc}.prop.exx2);
catch
    nPoints = numel(ebsd{inc}.prop.exx);
end
Ee_flat = reshape(E_exp, 4, nPoints); 

% Flatten the theoretical twin strain tensors
Ev_flat = reshape(E_id2d, 4, N);

% Initialize arrays to store the results
shapeMismatch = zeros(N, nPoints);
magnitudeX = zeros(N, nPoints);

% Calculate x and the shape mismatch for every point and every variant
for v = 1:N
    Ev = Ev_flat(:, v); 
    
    % Calculate the scalar projection factor 'x' 
    x = (Ev' * Ee_flat) ./ (Ev' * Ev);
    magnitudeX(v, :) = x;
    
    % Calculate the projected 'best fit' tensor
    Proj = Ev * x; 
    
    % Calculate the residual mismatch vector
    Mismatch = Ee_flat - Proj;
    
    % The shape parameter is the norm of this mismatch vector
    shapeMismatch(v, :) = vecnorm(Mismatch, 2, 1);
end


% Apply nanField to results before finding the minimum
shapeMismatch(:, nanField) = NaN;
magnitudeX(:, nanField) = NaN;

% Find the variant with the lowest shape mismatch for each point
[minMismatch, bestVariantIdx] = min(shapeMismatch, [], 1);

% Extract the magnitude 'x' associated with that specific best-fit variant
linearIndices = sub2ind(size(magnitudeX), bestVariantIdx, 1:nPoints);
bestMagnitudeX = magnitudeX(linearIndices);

% Create a logical mask of valid points based on criteria
% Twinning is unidirectional, so we only want positive x values
isValid = (bestMagnitudeX > magCutoff) & (minMismatch < shapeCutoff);

% Get indices of valid pixels and remove NaNs
validIndices = find(isValid);
nanMask = isnan(minMismatch); 
validIndices = validIndices(~nanMask(validIndices));

% Ensure the 'colors' colormap has enough rows for all N variants
if size(colors, 1) < N
    extra_colors = lines(N);
    colors = [colors; extra_colors(size(colors, 1) + 1 : N, :)];
end


% 2. Plot Best Fit Variants
disp('Plotting Best Fit Variants ...');

% Initialize plotting colors as pure white for the background
plotColorNoBackdrop = ones(nPoints, 3); 

% Assign the pure variant colors to the valid pixels
plotColorNoBackdrop(validIndices, :) = colors(bestVariantIdx(validIndices), :);

% Plotting the results
figure;
plot(ebsd{inc}, plotColorNoBackdrop) 
hold on 
% overlay outside boundaries
plot(grains{inc}.boundary('indexed','notIndexed'), 'lineColor', 'k','linewidth',4)

title('Best Fit Variant (No Backdrop)')

% Apply the colormap for the colorbar
colormap(colors)
caxis([1 N])
mtexColorbar
saveFigure(fullfile(outputPath,['inc_', num2str(inc), '_shape_magnitude_nobackdrop', comment, '.png']))

%% 3. Plot the "scalar projection factor 'x'". 1 means shear fits prefectly (in mag and dir)

figure;
plot(ebsd{inc}, bestMagnitudeX)
colormap turbo
setColorRange([0.2,1.35])
mtexColorbar
title('Scalar Projection Factor of Best Fit Variant')
saveFigure(fullfile(outputPath,['inc_', num2str(inc), '_shape_magnitude_projectionFactor', comment, '.png']))
return
%% --- GARBAGE SECTION ---

%% 3. Create Greyscale Backdrop from Max(x) and Plot
% alphaVariant = 0.6; % Define your alpha here (0 = fully transparent, 1 = fully solid)
% 
% disp('Plotting Best Fit Variants with Greyscale Backdrop...');
% 
% % Find the maximum magnitude across ALL variants for each point
% maxX_all = max(magnitudeX, [], 1);
% 
% % Handle any potential NaNs in the data to avoid breaking the normalization
% maxX_all(isnan(maxX_all)) = 0;
% 
% % Normalize the max(x) values to a 0 to 1 range for greyscale plotting.
% % Using the 80th percentile prevents a few hot pixels from washing out the contrast.
% c_max = prctile(maxX_all, 85);
% if c_max == 0; c_max = 1; end % Prevent division by zero
% normX = maxX_all / c_max;
% 
% % Clamp values to ensure they stay strictly within [0, 1]
% normX(normX > 1) = 1;
% normX(normX < 0) = 0;
% 
% % Invert so that low strain is white and high strain is dark 
% greyValues = normX(:); %1 - normX(:); 
% 
% % Convert the 1D array into an Nx3 RGB triplet matrix
% bg_color_rgb = repmat(greyValues, 1, 3);
% 
% % Initialize plotting colors with the greyscale backdrop
% plotColorShape = bg_color_rgb; 
% 
% % *** ALPHA BLENDING CALCULATIONS ***
% % Get the pure variant colors for the valid pixels
% baseVariantColors = colors(bestVariantIdx(validIndices), :);
% 
% % Get the existing background colors for those exact pixels
% baseBackgroundColors = bg_color_rgb(validIndices, :);
% 
% % Blend them using the alpha formula: C_out = alpha*C_front + (1-alpha)*C_back
% blendedColors = (alphaVariant * baseVariantColors) + ((1 - alphaVariant) * baseBackgroundColors);
% 
% % Apply the new blended colors to the final array
% plotColorShape(validIndices, :) = blendedColors;
% % ***********************************
% 
% % Plotting the results with the backdrop
% figure;
% plot(ebsd{inc}, plotColorShape) % No 'FaceAlpha' flag needed here
% title('Best Fit Variant on Max(x) Greyscale Backdrop')
% 
% % Apply the colormap for the colorbar
% colormap(colors)
% caxis([1 N])
% mtexColorbar
% saveFigure(['inc_', num2str(inc), '_shape_magnitude_backdrop', comment, '.png'])
% 
% %% 4. Specific Variant Shear Magnitude Analysis
% disp('Plotting shear magnitude for a specific variant...');
% 
% % --- User Settings for this block ---
% targetVariant = 4;         % Which variant ID to isolate (1 to N)
% plotDirectedShear = false;  % True = shear in twin direction, False = absolute max 2D shear
% % ------------------------------------
% 
% % Calculate expected 2D max engineering shear for the target theoretical variant
% % Using Mohr's circle for strain: gamma_max = 2 * sqrt(((exx - eyy)/2)^2 + exy^2)
% E_target_theo = E_id2d(:,:,targetVariant);
% expected_gamma_2d = 2 * sqrt(((E_target_theo(1,1) - E_target_theo(2,2))/2)^2 + E_target_theo(1,2)^2);
% 
% % Calculate experimental 2D max engineering shear for ALL points
% exx_exp = squeeze(E_exp(1,1,:));
% eyy_exp = squeeze(E_exp(2,2,:));
% exy_exp = squeeze(E_exp(1,2,:));
% gamma_max_exp = 2 * sqrt(((exx_exp - eyy_exp)/2).^2 + exy_exp.^2);
% 
% % Create the mask for pixels that were identified as the target variant
% targetMask = isValid & (bestVariantIdx == targetVariant);
% targetIndices = find(targetMask);
% 
% % Initialize the plot data with NaNs (so non-target grains remain white/empty)
% plotShearData = NaN(nPoints, 1);
% 
% if plotDirectedShear
%     % Plot the shear magnitude *projected* onto the twinning mode
%     % magnitudeX contains the scale factor 'x'. x * expected_gamma_2d gives the directed shear.
%     plotShearData(targetIndices) = magnitudeX(targetVariant, targetIndices)' .* expected_gamma_2d;
%     plotTitle = ['Directed 2D Shear Magnitude for Variant ', num2str(targetVariant)];
% else
%     % Plot the absolute max 2D shear of the experimental tensor
%     plotShearData(targetIndices) = gamma_max_exp(targetIndices);
%     plotTitle = ['Absolute Max 2D Shear for Variant ', num2str(targetVariant)];
% end
% 
% % Plotting
% figure;
% plot(ebsd{inc}, plotShearData);
% title(plotTitle);
% 
% % Customize colormap and colorbar for shear magnitude
% colormap viridis; % Good sequential colormap for magnitudes
% cb = mtexColorbar;
% setColorRange([expected_gamma_2d*0.8, expected_gamma_2d*1.2]);
% ylabel(cb, 'Shear Strain \gamma');
% 
% % Add the annotation text at the bottom of the figure
% annotationText = sprintf('Expected 2D Twinning Shear (\\gamma): %.4f', expected_gamma_2d);
% annotation('textbox', [0.1, 0.01, 0.8, 0.05], 'String', annotationText, ...
%            'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
%            'FontSize', 12, 'FontWeight', 'bold');
% 
% % saveFigure(['inc_', num2str(inc), '_shear_mag_variant_', num2str(targetVariant), comment, '.png'])