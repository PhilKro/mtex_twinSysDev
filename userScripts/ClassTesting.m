clear all
close
%%

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
orientatedTwin.orientation
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

%% testing multiple grains, multiple tS

multi_orientated_sym_Twin = dummyGrains*tS_sym
multi_orientated_sym_Twin.parent
multi_orientated_sym_Twin.orientation
multi_orientated_sym_Twin.parentTwinMisorientation

%% testing single tS, multiple tS

twinned_Twins = tS_1121*tS_sym
% if I pass in a 1x1tS*nx1tS_sym the line 260 (parentsExpanded =
% parentObj(idx);) creates a nx1 parentExpnaded twin object which only
% contains 1x1 properties (idx in that case is [1;1;1;...]. maybe sth wrong in subasng or sth like this? are we gonna breack
% any other functions changing that?

twinned_Twins.parent
twinned_Twins.orientation

%% testing multiple tS, single tS
multitwinned_Twin = tS_sym*tS_1121
% multitwinned_Twin.parent
% multitwinned_Twin.orientation

% here it fales bc idx = [1;2;3;4;...] but sth wrong in the size of the
% PARENT (FIGURE OUT) mnb


%% testing multiple tS, multiple tS
multitwinned_Twins = tS_sym*tS_sym
multitwinned_Twins.parent
multitwinned_Twins.orientation

%% next prompt going down the rabbit hole (or maybe blocking this by introducing a "isMultiplied" property that stops the user from going here)
% okay what happens if i symmetrise a activated twin system. that should be no problem right. the parent of the symmetrised twins will be the same and makes physical sense (should cause no problems, please think of problems that could arise??). right now these properties get chugged out during symmetrization. however if this twinSystem is the parent of another, how would that be handled, either some recursive stuff (or better forbid it, what do you think)
% ja i dont think thats it the length of the object is determined correctly and it is created correctly as a 6x1 object array, but it is only populated with 1x1 object arrays. doesnt that mean sth is wrong with cat or sth or subsref but numel(tS.eta1) reliably returns the number of tS elements. somehow implementing a numel method breaks other logic bc it is called befor the subsref (sth like this i made a comment in the subsref method). please think again

