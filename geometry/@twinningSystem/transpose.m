function tS=transpose(tS)
% transpose list of slipSystem

tS.b = tS.b.';
tS.n = tS.n.';
tS.eta1 = tS.eta1.';
tS.eta2 = tS.eta2.';
tS.k1 = tS.k1.';
tS.k2 = tS.k2.';
tS.rotAxis = tS.rotAxis.';

tS.CRSS = tS.CRSS.';