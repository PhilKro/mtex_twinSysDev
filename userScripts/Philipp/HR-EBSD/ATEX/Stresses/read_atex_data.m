function data_struct = read_atex_data(filename, requested_columns)
    % Reads an ATEX .txt data file and extracts specified columns into 2D matrices.
    %
    % INPUTS:
    %   filename - string, path to the data file
    %   requested_columns - cell array of strings, names of the columns to extract
    %
    % OUTPUT:
    %   data_struct - struct with field names corresponding to requested columns,
    %                 containing 2D arrays of extracted data
    
    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read header lines and extract metadata
    col_names = {};
    X_pixels = 0;
    Y_pixels = 0;
    
    while true
        line = fgetl(fid);
        if contains(line, '# X pixels:')
            X_pixels = sscanf(line, '# X pixels: %d');
        elseif contains(line, '# Y pixels:')
            Y_pixels = sscanf(line, '# Y pixels: %d');
        elseif startsWith(line, '# Column')
            % Extract column name
            col_idx = sscanf(line, '# Column%d:', 1);
            col_name = extractAfter(line, sprintf('# Column%d: ', col_idx));
            col_name = strtrim(col_name);
            col_names{col_idx} = col_name;
        elseif startsWith(line, 'DATA BEGIN')
            break;
        end
    end
    
    % Function to clean and format variable names
    %format_name = @(s) regexprep(s, {'σ', 'ε', 'ω', '[{}()°/-]', '\s+', '₁', '₂', '₃'}, {'s', 'e', 'w', '', '', '1', '2', '3'});
    
    % Find indices of requested columns
    col_indices = zeros(size(requested_columns));
    formatted_names = cell(size(requested_columns));
    for i = 1:length(requested_columns)
        match = find(contains(col_names, requested_columns{i}), 1);
        if isempty(match)
            error('Column not found: %s', requested_columns{i});
        end
        col_indices(i) = match;
        formatted_names{i} = variableNameFormatter(requested_columns{i});
    end
    
    % Read the data
    data = fscanf(fid, '%f', [length(col_names), inf]);
    fclose(fid);
    data = data'; % Transpose so that rows correspond to entries
    
    % Reshape selected columns into 2D matrices and store in struct
    data_struct = struct();
    for i = 1:length(requested_columns)
        column_data = data(:, col_indices(i));
        data_struct.(formatted_names{i}) = transpose(reshape(column_data, [X_pixels, Y_pixels]));
    end
end
