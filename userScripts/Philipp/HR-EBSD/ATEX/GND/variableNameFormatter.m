function cleanedVarName = variableNameFormatter(varName)
    % Replace common Greek symbols and subscripts with ASCII equivalents
    replacements = {
        'σ', 's';
        'ε', 'e';
        'ω', 'w';
        'α', 'a';
        'β', 'b';
        'γ', 'g';
        'µ', 'u';
        '¹', '1';
        '²', '2';
        '₁', '1';
        '₂', '2';
        '₃', '3';
        '⁻', '-';
        '–', '-';
    };

    for i = 1:size(replacements,1)
        varName = strrep(varName, replacements{i,1}, replacements{i,2});
    end

    % Use MATLAB's native sanitization
    cleanedVarName = matlab.lang.makeValidName(varName);
end
