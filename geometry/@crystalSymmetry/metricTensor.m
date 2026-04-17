function [G_dir] = metricTensor(cs)
% METRICTENSOR Computes the direct metric tensor for a crystalSymmetry.
%
% Syntax
%   G_dir = cs.metricTensor()
%
% Output
%   G_dir - 3x3 direct metric tensor (Gram matrix).
M = double(cs.axes);
G_dir = M' * M;
end

