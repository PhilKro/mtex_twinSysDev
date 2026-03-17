function out = mapEBSDprop(W, preId)
% mapEBSDprop  Remap an EBSD property W using index map preId.
%
% out(i) = W(preId(i)) when preId(i) is valid
% out(i) = nan-of-correct-type when preId(i) is NaN

valid = ~isnan(preId);
n = numel(preId);

% ---- allocate nan-array of correct class ----

if isnumeric(W) || islogical(W)
    out = nan(n,1, class(W));

elseif isa(W, 'Miller') || isa(W, 'Miller2d')
    CS = W.CS;
    out = repmat(Miller.nan(CS), n, 1);

elseif ismethod(W, 'nan')
    out = W.nan(n,1);   % rotation, vector3d, slipsystem, etc.

else
    error('Unsupported data type for mapping.');
end

% ---- assign only valid positions ----
out(valid) = W(preId(valid));

end

