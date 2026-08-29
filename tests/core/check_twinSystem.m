function check_twinSystem
% compatibility checks for twinSystem on MTEX 7

csH = crystalSymmetry('6/mmm',[3.21 3.21 5.21], ...
  'mineral','synthetic magnesium');
csC = crystalSymmetry('m-3m',[3.6 3.6 3.6], ...
  'mineral','synthetic cubic');
csT = crystalSymmetry('4/mmm',[3 3 4], ...
  'mineral','synthetic tetragonal');

hcp = {twinSystem.hexagonal_1012(csH), ...
  twinSystem.hexagonal_1121(csH), ...
  twinSystem.hexagonal_2241(csH), ...
  twinSystem.hexagonal_1011(csH), ...
  twinSystem.hexagonal_1011_TypeI(csH), ...
  twinSystem.hexagonal_2021(csH), ...
  twinSystem.hexagonal_1013(csH), ...
  twinSystem.hexagonal_1013_TypeI(csH), ...
  twinSystem.hexagonal_1122(csH), ...
  twinSystem.hexagonal_1124(csH)};
systems = [hcp,{twinSystem.fcc_111(csC), ...
  twinSystem.bcc_112(csC),twinSystem.bct_101(csT)}];

for k = 1:numel(systems)
  tS = systems{k};
  assert(isa(tS,'twinSystem') && isscalar(tS))
  assert(isa(tS.k1,'Miller') && isa(tS.eta1,'Miller'))
  assert(isfinite(tS.shearMagnitude) && tS.shearMagnitude > 0)
end

tS = hcp{1};
[variants,id] = symmetrise(tS);
assert(length(variants) > 1 && all(id == 1))
assert(isequal(reshape(variants.variantId,1,[]),1:length(variants)))

mori = parentTwinMisorientation(variants);
assert(isa(mori,'orientation') && length(mori) == length(variants))
assert(all(isImproper(mori)))

r = Miller(0,0,0,1,csH,'UVTW');
sf = SchmidFactor(tS,r);
assert(isfinite(sf) && isscalar(sf))

F = calcDeformationSequence(variants,[1 2]);
assert(isequal(size(F),[3 3]) && all(isfinite(F(:))))

ori = orientation.byEuler(20*degree,30*degree,10*degree,csH);
rotated = ori * variants;
assert(isa(rotated,'twinSystem') && length(rotated) == length(variants))
assert(~isa(rotated.k1,'Miller') && ~isa(rotated.eta1,'Miller'))

k2 = Miller(1,1,0,csC,'hkl');
eta1 = Miller(1,1,1,csC,'uvw');
typeII = twinSystem.byK2eta1(k2,eta1,2);
assert(typeII.twinType == 2 && ~isImproper(parentTwinMisorientation(typeII)))

disp('check_twinSystem: passed')

end
