function l = isMirrorPlane(m)
% check if a Miller index corresponds to a mirror plane
%
% Syntax
%   l = isMirrorPlane(m)
%
% Input
%   m - @Miller
%
% Output
%   l - logical

% get symmetry operations
symOps = rotation(m.CS);

% generate reflection operators for the planes
mir = reflection(m);

% check for each plane if it corresponds to a symmetry operation
l = false(size(m));
for i = 1:length(symOps)
  l = l | angle(mir,symOps(i)) < 1e-4;
end