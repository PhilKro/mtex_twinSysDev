function plotShape(tS)
%WHATZEFAK Summary of this function goes here
%   Detailed explanation goes here
% go with CS.lattice == 6  -> hexagonal
if strcmp(tS.CS.char(),'"6/mmm"')
    cS = crystalShape.hex(tS.CS);
elseif strcmp(tS.CS.char(),'"m-3m"')
    cS = crystalShape.cube(tS.CS);
else
    disp('Do not have a crystal shape for this symmetry, please define yourself in @twinningSystem/plotShape'); 
    return
end

plot(cS, 'faceAlpha', 0)
hold on 
plot(cS,tS.k1,'faceColor','blue','label','K1')



d = tS.eta1;
%d = normalize(d) * min(m);
%h = [h,arrow3d(d,varargin{:})];
plot(cS,tS.eta1,'uvw','faceColor','blue','label','K1')

end
