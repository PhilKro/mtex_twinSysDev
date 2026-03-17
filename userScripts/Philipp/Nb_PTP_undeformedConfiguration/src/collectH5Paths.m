function allPaths = collectH5Paths(infoStruct, currentPath)
    if nargin < 2
        currentPath = '';
    end

    allPaths = {};

    % Normalize path (remove trailing slash)
    if ~isempty(currentPath) && currentPath(end) == '/'
        currentPath(end) = [];
    end

    % Datasets
    for i = 1:numel(infoStruct.Datasets)
        ds = infoStruct.Datasets(i).Name;

        if isempty(currentPath)
            full = ['/' ds];
        else
            full = [currentPath '/' ds];
        end

        allPaths{end+1} = full; %#ok<AGROW>
    end

    % Subgroups
    for i = 1:numel(infoStruct.Groups)
        grpName = infoStruct.Groups(i).Name;  % Already like "/B00001"
        allPaths = [allPaths, collectH5Paths(infoStruct.Groups(i), grpName)]; %#ok<AGROW>
    end
end
