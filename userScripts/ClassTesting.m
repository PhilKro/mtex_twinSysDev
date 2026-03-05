cs = crystalSymmetry('6/mmm', 'mineral', 'Re');
tS_1121 = twinSystem.hexagonal_1121_0001(cs);
tS_sym = tS_1121.symmetrise;
tS_sym(2)
slippy_the_SlipSystem = slipSystem.basal(cs);

parentOri_1 = orientation.byEuler(1,1,1,cs);
parentOri_2 = orientation.byEuler(2,1,2,cs);

% Create dummy EBSD (2x2 grid)
[x,y] = meshgrid(0:1,0:1);
pos = vector3d(x(:),y(:),0);
oris = [parentOri_1, parentOri_1, parentOri_2, parentOri_2];
dummyEBSD = EBSD(pos, oris, 2*ones(numel(pos),1),cs, struct('ci', ones(numel(pos),1)));
dummyEBSD = dummyEBSD.gridify;

dummyGrains = calcGrains(dummyEBSD,'angle',10*degree);


%% testing multiple orientations, one tS

orientatedTwin = oris*tS_1121
orientatedTwin.parent
orientatedTwin.orientaiton
orientatedTwin.parentTwinMisorientation

%% testing multiple grains, one tS

grainedTwin = mtimes(dummyGrains,tS_1121)
grainedTwin.parent
grainedTwin.orientation
grainedTwin.parentTwinMisorientation

%% testing single orientation, multiple tS
ori = oris(1);
orientated_sym_Twin = ori*tS_sym
orientated_sym_Twin.parent
orientated_sym_Twin.orientation
orientated_sym_Twin.parentTwinMisorientation

%% testing single grain, multiple tS
grain = dummyGrains(1);
grained_sym_Twin = grain*tS_sym
grained_sym_Twin.parent
grained_sym_Twin.orientation
grained_sym_Twin.parentTwinMisorientation

%% testing multiple grain, multiple tS

multi_grained_sym_Twin = dummyGrains*tS_sym
multi_grained_sym_Twin.parent
multi_grained_sym_Twin.orientation
multi_grained_sym_Twin.parentTwinMisorientation

%% testing multiple ori, multiple tS

multi_orientated_sym_Twin = dummyGrains*tS_sym
multi_orientated_sym_Twin.parent
multi_orientated_sym_Twin.orientation
multi_orientated_sym_Twin.parentTwinMisorientation

%% testing nested twins

%TODO



%% next prompt going down the rabbit hole (or maybe blocking this by introducing a "isMultiplied" property that stops the user from going here)
% okay what happens if i symmetrise a activated twin system. that should be no problem right. the parent of the symmetrised twins will be the same and makes physical sense (should cause no problems, please think of problems that could arise??). right now these properties get chugged out during symmetrization. however if this twinSystem is the parent of another, how would that be handled, either some recursive stuff (or better forbid it, what do you think)
