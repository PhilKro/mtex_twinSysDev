function tS = rotate_outer(tS,rot)
% rotate twin systems by every rotation

fields = {'rotAxis','k1','eta1','k2','eta2'};
for k = 1:numel(fields)
  tS.(fields{k}) = rotate_outer(tS.(fields{k}),rot);
end
tS.CRSS = repmat(tS.CRSS(:).',length(rot),1);
tS.twinType = repmat(tS.twinType(:).',length(rot),1);
tS.variantId = repmat(tS.variantId(:).',length(rot),1);

end
