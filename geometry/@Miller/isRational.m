function bool = isRational(m, maxRationalIndex)
    % check if the Miller index is rational
    % maxRationalIndex is the maximum allowed value for the indices to be considered rational. Defined for 3 index notation (4 index notation will be converted to 3 index notation before checking)
    if nargin < 2
        maxRationalIndex = 5;
    end
    maxRationalIndex = maxRationalIndex*2;
    [~, err] = m.round('maxHKL', maxRationalIndex,'round3Index');
    bool = err < 1e-4;
end