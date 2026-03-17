%% Stress Tensors
map_s11 = ATEX_stress.Stresss11insampleframe;
map_s22 = ATEX_stress.Stresss22insampleframe;
map_s33 = ATEX_stress.Stresss33insampleframe;
map_s12 = ATEX_stress.Stresss12insampleframe;
map_s23 = ATEX_stress.Stresss23insampleframe;
map_s31 = ATEX_stress.Stresss13insampleframe;


stressMatrices = cat(4, cat(3, map_s11, map_s12, map_s31),cat(3, map_s12, map_s22, map_s23),cat(3,map_s31, map_s23, map_s33));

% Put this in in the static methods of the stressTensor class
    % function sigma = zeros(varargin)
    %     r = get_option(varargin,'rank',2);
    %     varargin = delete_option(varargin,'rank',1);
    %     [cs,varargin] = getClass(varargin,'symmetry');
    %     d = [repmat(3,1,r),varargin{:},1];
    %     sigma = stressTensor(zeros(d),'rank',r);
    %     if ~isempty(cs), sigma.CS = cs; end
    % end
stressTensors = stressTensor.zeros(size(map_s11));
%if stressMatrices are NaN it throws a warning
warning('off')
for i = 1:size(stressTensors, 1)
    for j = 1:size(stressTensors, 2)
        stressTensors(i,j) = stressTensor(squeeze(stressMatrices(i,j,1:end,1:end)));
    end
end
warning('on')
%% stress tensor rotation
%param
grainOfinterestID = 24;
% Crosscourt: 
% The rotation axes are:-
% 1 – Horizontal, to the right.
% 2 – Vertical, up the map.
% 3 – Normal, out of the surface.

%MTEX
% 1 – ('xAxisDirection') ,('east') to the right.
% 2 – Vertical, down the map.
% 3 – ('zAxisDirection') ,('intoPlane') Normal, into the surface
Coordiffrence_CC_MTEXrot = rotation.byAxisAngle(xvector,180*degree);


%parent grain is determined as largest grain
isParentGrain = grains.grainSize == max(grains.grainSize);
parentOrientation = grains.meanOrientation(isParentGrain);

%


stressTensors0 = stressTensors;


%CRITICAL!!!
%parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
%CRITICAL!!!
MTEXStressTensors = rotate(stressTensors,Coordiffrence_CC_MTEXrot);
%matrixStressTensors = rotate(MTEXStressTensors, rotation(parentOrientation));
matrixStressTensors = inv(parentOrientation) .* MTEXStressTensors;

%determine k1 plane of twin
tS = twinningSystem.hexagonal_1121_0001(CS)
tS_sym = tS.symmetrise;
errorInDegrees = 10;




%define parent grain on your own
%parentOrientation = orientation.byEuler(183.2* degree, 92.5* degree, 143.4* degree, CS);

%parent child misorientations https://mtex-toolbox.github.io/MisorientationTheory.html
misOrientation = inv(grains.meanOrientation(grainOfinterestID)).*parentOrientation;

%[variantIndex] = tS.variantDetermination2(misOrientation,errorInDegrees);
variantIndex = 1;

variant = tS_sym(variantIndex)
 
%CRITICAL!!!
%%IS THIS CORRECT?????
%habitPlaneOri = orientation.byMiller(variant.k1,variant.eta1, CS, CS);


% ROTATION INTO THE ETA1 and ROTAXIS DIRECTIONS
%habitPlaneOri = orientation.map(variant.eta1, inv(parentOrientation)*xvector, variant.rotAxis, inv(parentOrientation)*yvector)
%mapping the global coordinate representation of the eta1 direction onto the x vector,
%and mapping the global coordinate representation of the k1 direction onto the y vector
%habitPlaneRot = rotation.map(parentOrientation*variant.eta1, xvector,parentOrientation*variant.rotAxis,yvector)

%does this change anything?
%mapping rot to e.g. y and k1 to z
habitPlaneRot = rotation.map(parentOrientation*variant.rotAxis, yvector,parentOrientation*variant.k1,zvector)
%habitPlaneRot = rotation.map(parentOrientation*variant.eta1, xvector,parentOrientation*variant.k1,zvector)


% this is a orientation
%testToSee = orientation.map(variant.eta1, xvector,variant.rotAxis,yvector)
% the mapping creates a orientation Re -> xyz, thats why the inverse has to
% be applied to the StressTensors, thus also the inverse of the rotation
% (habitPlaneRot) should be applied

%twinningPlaneStressTensors = rotate(matrixStressTensors, rotation(habitPlaneOri));
twinningPlaneStressTensors = inv(habitPlaneRot) .* MTEXStressTensors;
%% coordinate plotting
% all the:
% Warning: Possibly applying an orientation to an object in specimen coordinates!
% are created by the plotting of the cS1
% (plot(grains,0.1*cS1,'facecolor','r'))
figure;
plot(ebsd,ebsd.ci, 'figureTitle','MTEX coordinates')
hold on
quiver(grains(isParentGrain), xvector)
quiver(grains(isParentGrain), yvector)
quiver(grains(isParentGrain), zvector)
hold off

figure;
plot(ebsd,ebsd.ci, 'figureTitle','ATEX coordinates')
hold on
quiver(grains(isParentGrain), rotate(xvector, Coordiffrence_CC_MTEXrot))
quiver(grains(isParentGrain), rotate(yvector, Coordiffrence_CC_MTEXrot))
quiver(grains(isParentGrain), rotate(zvector, Coordiffrence_CC_MTEXrot))
% quiver(grains(isParentGrain), inv(Coordiffrence_CC_MTEXrot)*xvector)
% quiver(grains(isParentGrain), inv(Coordiffrence_CC_MTEXrot)*yvector)
% quiver(grains(isParentGrain), inv(Coordiffrence_CC_MTEXrot)*zvector)
hold off

figure;
plot(ebsd,ebsd.ci, 'figureTitle','ParentGrain/matrix coordinates')
hold on


plot(grains,0.1*cS1,'facecolor','r')

% parentX = rotate(xvector, rotation(parentOrientation));
% parentY = rotate(yvector, rotation(parentOrientation));
% parentZ = rotate(zvector, rotation(parentOrientation));
parentX = (parentOrientation)*CS.X;
parentY = (parentOrientation)*CS.Y;
parentZ = (parentOrientation)*CS.Z;
quiver(grains(isParentGrain), parentX)
quiver(grains(isParentGrain), parentY)
quiver(grains(isParentGrain), parentZ)
hold off

figure
eta1_uvtw = round(variant.eta1);
k1_hkil = round(variant.k1);
plot(ebsd,ebsd.ci, 'figureTitle',['x || [', char(num2str(eta1_uvtw.UVTW)), '] ; z || (', num2str(k1_hkil.hkil), ') (direction and plane of the ParentGrain/matrix)'])
hold on
% shearDir = rotate(parentX, rotation(habitPlaneOri));
% rotAxisDir = rotate(parentY, rotation(habitPlaneOri));
% shearPnormal = rotate(parentZ, rotation(habitPlaneOri));
shearDir = (parentOrientation)*variant.eta1;
% the rotation of the xvector with the habitPlaneRot should be equivalent
% to the multiplication of parentOrientation*variant.eta1
%equivalent for rotAxis and k1, y and z, respectively 
assert(abs(dot(normalize(shearDir),normalize(inv(habitPlaneRot)*xvector)))-1<1e-15)

rotAxisDir = (parentOrientation)*variant.rotAxis;
assert(abs(dot(normalize(rotAxisDir),normalize(inv(habitPlaneRot)*yvector)))-1<1e-15)

shearPnormal = (parentOrientation)*variant.k1;
assert(abs(dot(normalize(shearPnormal),normalize(inv(habitPlaneRot)*zvector))-1)<1e-15)

quiver(grains(isParentGrain), shearDir)
quiver(grains(isParentGrain), rotAxisDir)
quiver(grains(isParentGrain), shearPnormal)

plot(grains,0.1*cS1,'facecolor','r')

% nextAxis
% plot(parentOrientation*cS1,'faceColor','blue', 'faceAlpha',0)
% hold on
% plotTwinningSystem(cS1,variant)
hold off
%% matrix Stress Plots %%%
% param
stressLim = 5000;

%

matrixStressMatrices = matrixStressTensors.M;

matrix_s11 = squeeze(matrixStressMatrices(1,1,:,:));
matrix_s12 = squeeze(matrixStressMatrices(1,2,:,:));
matrix_s31 = squeeze(matrixStressMatrices(3,1,:,:));
matrix_s22 = squeeze(matrixStressMatrices(2,2,:,:));
matrix_s23 = squeeze(matrixStressMatrices(2,3,:,:));
matrix_s33 = squeeze(matrixStressMatrices(3,3,:,:));

mtexFigure('layout',[3 3], 'FigureTitle', 'Stresses in matrix coordinate system');

plot(ebsd,matrix_s11.*(abs(matrix_s11)<stressLim))
title('s11')

nextAxis(1,2)
plot(ebsd,matrix_s12.*(abs(matrix_s12)<stressLim))
title('s12')

nextAxis(1,3)
plot(ebsd,matrix_s31.*(abs(matrix_s31)<stressLim))
title('s13')

nextAxis(2,2)
plot(ebsd,matrix_s22.*(abs(matrix_s22)<stressLim))
title('s22')

nextAxis(2,3)
plot(ebsd,matrix_s23.*(abs(matrix_s23)<stressLim))
title('s23')

nextAxis(3,3)
plot(ebsd,matrix_s33.*(abs(matrix_s33)<stressLim))
title('s33')

setColorRange([-stressLim,stressLim],'current')
mtexColorbar('title', 'Stress / GPa')
colormap jet


saveFigure(append(file,['_StressesInMatrixCoordinates.png']))

%% twinning Plane Stress plots
% param
stressLim = 8000;

%
twinningPlaneStressMatrices = twinningPlaneStressTensors.M;

tP_s11 = squeeze(twinningPlaneStressMatrices(1,1,:,:));
tP_s12 = squeeze(twinningPlaneStressMatrices(1,2,:,:));
tP_s31 = squeeze(twinningPlaneStressMatrices(3,1,:,:));
tP_s22 = squeeze(twinningPlaneStressMatrices(2,2,:,:));
tP_s23 = squeeze(twinningPlaneStressMatrices(2,3,:,:));
tP_s33 = squeeze(twinningPlaneStressMatrices(3,3,:,:));

uvtw = round(variant.eta1);
rotAxis = round(variant.rotAxis);
fig = mtexFigure('layout',[3 3], 'FigureTitle', append('Stresses in the coordinate system: 1 || [', num2str((uvtw.UVTW)), '] ; 3 || (', num2str(k1_hkil.hkil), ') (direction and plane of the matrix)'));
title(append('Stresses in the coordinate system: 1 || [', num2str(uvtw.UVTW), '] ; 3 || (', num2str((k1_hkil.hkil)), ') (direction and plane of the matrix)'));

nextAxis(1,1)
plot(ebsd,tP_s11.*(abs(tP_s11)<stressLim))
title('s11')


nextAxis(1,2)
plot(ebsd,tP_s12.*(abs(tP_s12)<stressLim))
title('s12')

nextAxis(1,3)
plot(ebsd,tP_s31.*(abs(tP_s31)<stressLim))
title('s13')

nextAxis(2,2)
plot(ebsd,tP_s22.*(abs(tP_s22)<stressLim))
title('s22')

nextAxis(2,3)
plot(ebsd,tP_s23.*(abs(tP_s23)<stressLim))
title('s23')

nextAxis(3,3)
plot(ebsd,tP_s33.*(abs(tP_s33)<stressLim))
title('s33')
setColorRange([-stressLim,stressLim],'current')
mtexColorbar('title', append('Stress / GPa in the coordinate system:',newline,' 1 || [', num2str(uvtw.UVTW), '] ; 3 || (', num2str((k1_hkil.hkil)), ') (direction and plane of the matrix)'));
colormap jet
hold off
saveFigure(append(file,['_StressesInTwinPlaneCoordinates.png']))
%% with alpha
fig = mtexFigure('layout',[3 3], 'FigureTitle', ...
    append('Stresses in the coordinate system: 1 || [', num2str((uvtw.UVTW)), ...
    '] ; 2 || (', num2str(rotAxis.UVTW), ') (direction and plane of the matrix)'));

% List of stress maps and titles
stressMaps = {tP_s11, [], [], tP_s12, tP_s22, [], tP_s31, tP_s23, tP_s33};
titles = {'s11','','','s12','s22','','s13','s23','s33'};

for k = 1:9
    if isempty(stressMaps{k})
        mtexColorbar off
        colormap(gca, jet)
        continue
    end
    [r, c] = ind2sub([3 3], k);
    nextAxis(r, c)
    
    % Background: grayscale cibw image
    plot(ebsd, cibw, 'micronbar','off')
    hold on

    % Generate stress map with transparency
    stressMap = stressMaps{k};
    mask = abs(stressMap) < stressLim;
    
    % Overlay: transparent colored patch of stress values
    stressData = stressMap;
    stressData(~mask) = NaN; % suppress out-of-range

    plot(ebsd, stressData, 'FaceAlpha', 0.25)
    mtexColorbar off
    % Axis title
    title(titles{k})
    colormap(gca, jet)
    setColorRange([-stressLim, stressLim], 'current')
end

% Apply colormap only to the final stress layer
mtexColorbar('title', append('Stress / MPa in the coordinate system:',newline,...
    ' 1 || [', num2str(uvtw.UVTW), '] ; 2 || (', num2str((rotAxis.UVTW)), ') (direction and plane of the matrix)'));
%% line profile from stress
% Step 1: Plot the s12 stress component
fig = figure;
stressMap = tP_s31;
plot(ebsd, stressMap, 'FaceAlpha', 1);
colormap(jet)
setColorRange([-stressLim, stressLim])
mtexColorbar('title','s13 / MPa')
title('Draw a line to extract s13 profile');

%extract along line
h = drawline('Color', 'k', 'LineWidth', 1.5);
pos = h.Position;  % 2x2 array: [x1 y1; x2 y2]

% Extract the coordinates
xClicks = pos(:,1);
yClicks = pos(:,2);

% Step 3: Interpolate stress along the line
nPoints = 200;
xLine = linspace(xClicks(1), xClicks(2), nPoints);
yLine = linspace(yClicks(1), yClicks(2), nPoints);

% Prepare interpolation grid
F = scatteredInterpolant(ebsd.x(:), ebsd.y(:), stressMap(:), 'linear', 'none');
stressLine = F(xLine, yLine);

% Compute distance along the line
dist = sqrt((xLine - xLine(1)).^2 + (yLine - yLine(1)).^2);

% Step 4: Plot the profile
figure;
plot(dist, stressLine, 'o', ...
    'MarkerSize', 6, ...          % size of the dots
    'MarkerEdgeColor', 'b', ...   % border color
    'MarkerFaceColor', 'b')       % fill color
  %'b-', 'LineWidth', 1.5)
xlabel('Distance along line [μm]')
ylabel('s13 stress [MPa]')
title('Line profile of s13 stress')
grid on
