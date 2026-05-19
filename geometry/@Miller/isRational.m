function bool = isRational(m, maxRationalIndex)
    % check if the Miller index is rational
    % maxRationalIndex is the maximum allowed value for the indices to be considered rational. Defined for 3 index notation (4 index notation will be converted to 3 index notation before checking)
    if nargin < 2
        maxRationalIndex = 3;
    end
    format = m.dispStyle;
    switch format
        case {'hkl', 'uvw'}
            intMiller = round(m, 'maxHKL', maxRationalIndex);
            bool = angle(m, intMiller) < 1e-4;
        case 'hkil'
            mtemp = m;
            mtemp.dispStyle = 'hkl';
            intMiller = round(mtemp, 'maxHKL', maxRationalIndex);
            bool = angle(mtemp, intMiller) < 1e-4;
        case 'UVTW'
            mtemp = m;
            mtemp.dispStyle = 'uvw';
            intMiller = round(mtemp, 'maxHKL', maxRationalIndex);
            bool = angle(mtemp, intMiller) < 1e-4;
    end
end