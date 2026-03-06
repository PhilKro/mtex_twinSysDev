function tS = subSet(tS,ind)
% subindex vector3d

tS.k1 = tS.k1(ind);
tS.k2 = tS.k2(ind);
tS.eta1 = tS.eta1(ind);
tS.eta2 = tS.eta2(ind);
tS.rotAxis = tS.rotAxis(ind);
tS.CRSS = tS.CRSS(ind);
tS.twinType = tS.twinType(ind);
if ~isempty(tS.parent)
    tS.parent = tS.parent(ind);
else
    tS.parent = [];
end

end