function ebsd = read_atex_data_complete(filename, ebsd)
    % Reads an ATEX .txt data file and assigns all columns as properties to an EBSD object.
    %
    % INPUTS:
    %   filename - string, path to the ATEX data file
    %   ebsd - EBSD object to which properties will be added
    %
    % OUTPUT:
    %   ebsd - EBSD object with added properties from the ATEX file
    
    % Check if the ebsd object is not ungridified
    if isa(ebsd, "EBSDsquare") || isa(ebsd, "EBSDhex")
        error("Please provide the ungridified ebsd object!");
    end
    
    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Initialize
    col_names = {};
    X_pixels = 0;
    Y_pixels = 0;
    
    % Parse header
    while true
        line = fgetl(fid);
        if ~ischar(line)
            error('Unexpected end of file before DATA BEGIN.');
        end
        
        if contains(line, '# X pixels:')
            X_pixels = sscanf(line, '# X pixels: %d');
        elseif contains(line, '# Y pixels:')
            Y_pixels = sscanf(line, '# Y pixels: %d');
        elseif startsWith(line, '# Column')
            % Extract column index and name
            col_idx = sscanf(line, '# Column%d:', 1);
            col_name = extractAfter(line, sprintf('# Column%d: ', col_idx));
            col_name = strtrim(col_name);
            col_names{col_idx} = variableNameFormatter(col_name);
        elseif startsWith(line, 'DATA BEGIN')
            break;
        end
    end

    % Read numeric data (whitespace-separated)
    data = fscanf(fid, '%f', [length(col_names), inf]);
    fclose(fid);
    data = data'; % rows = datapoints, columns = variables
    
    % Sanity check
    nPoints = size(data, 1);
    if ebsd.length ~= nPoints
        error('EBSD object length (%d) does not match number of data points (%d).', ebsd.length, nPoints);
    end
    
    % Assign each column to ebsd.prop
    for i = 1:length(col_names)
        col_name = col_names{i};
        ebsd.prop.(col_name) = data(:, i);
    end
    
    % ============================================
    % (1) UNIT CORRECTION: m^-2 → µm^-1 for Nye tensor components
    % ============================================
    propNames = fieldnames(ebsd.prop);

    % Find the µm⁻¹ and m⁻² variants
    field_um = find(contains(propNames, "NormOfTheNye", "IgnoreCase", true) & ...
                    contains(propNames, "CurvaturesOnly", "IgnoreCase", true) & ...
                    contains(propNames, "um_1", "IgnoreCase", true), 1);
    field_m = find(contains(propNames, "NormOfTheNye", "IgnoreCase", true) & ...
                   contains(propNames, "CurvaturesOnly", "IgnoreCase", true) & ...
                   contains(propNames, "m_2", "IgnoreCase", true), 1);

    if ~isempty(field_um) && ~isempty(field_m)
        f1 = propNames{field_um};
        f2 = propNames{field_m};

        v1 = ebsd.prop.(f1);
        v2 = ebsd.prop.(f2);

        % find first non-zero index in the µm⁻¹ variable
        idx_nonzero = find(v1 ~= 0 & ~isnan(v1), 1);
        if ~isempty(idx_nonzero)
            ratio = v1(idx_nonzero) / v2(idx_nonzero);
            fprintf('Found conversion ratio between m⁻² and µm⁻¹: %.3e\n', ratio);

            % Scale all Nye tensor components (curvatures only) in sample frame
            nyeFields = propNames(contains(propNames, "Component", "IgnoreCase", true) & ...
                                  contains(propNames, "Nye", "IgnoreCase", true) &...
                                  contains(propNames, "InSampleFr", "IgnoreCase", true));
            for i = 1:numel(nyeFields)
                ebsd.prop.(nyeFields{i}) = ebsd.prop.(nyeFields{i}) .* ratio;
            end

            fprintf('Scaled %d Nye tensor component fields by conversion factor.\n', numel(nyeFields));
            if numel(nyeFields) ~= 8
                warning('8 instead of %d Nye tensor component fields should have been rescaled.\n', numel(nyeFields))
            end
        else
            warning('No non-zero entries found in %s. No scaling applied.', f1);
        end
    else
        warning('µm⁻¹ and/or m⁻² variants of NormOfTheNye (curvatures only) not found. No unit conversion applied.');
    end
    fprintf('Assigned %d properties to EBSD object.\n', length(col_names));
    
    % ============================================
    % (2) FILTERING based on IC-GN iteration count
    % ============================================
    iterationFields = fieldnames(ebsd.prop);
    idxField = find(contains(iterationFields, "NumberOfIterationsOfThe"), 1);

    if isempty(idxField)
        warning("The number of iterations of the IC-GN algorithm weren't found, no filtering for <=0 or >200 was applied.");
    else
        iterationField = iterationFields{idxField};
        iterValues = ebsd.prop.(iterationField);

        % Identify invalid iteration entries
        invalidIdx = iterValues <= 0 | iterValues > 200;

        % Set all properties to NaN where invalid
        propNames = fieldnames(ebsd.prop);
        for i = 1:numel(propNames)
            vals = ebsd.prop.(propNames{i});
            if isnumeric(vals) && numel(vals) == nPoints
                vals(invalidIdx) = NaN;
                ebsd.prop.(propNames{i}) = vals;
            end
        end

        fprintf("Filtered %d ATEX data points based on %s (IC-GN iteration count <=0 or >200).\n", sum(invalidIdx), iterationField);
    end 
end
