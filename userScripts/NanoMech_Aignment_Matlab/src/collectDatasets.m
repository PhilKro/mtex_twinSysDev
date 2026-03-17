function datasetPaths = collectDatasets(info)
% collectDatasets Recursively collect HDF5 dataset full paths
%   datasetPaths = collectDatasets(info)
%
%   Given an output of `h5info` (the root info struct or a group struct),
%   this function returns a cell array of dataset full paths using '/'
%   separators which are suitable for use with `h5read`.

datasetPaths = {};
if isempty(info)
    return
end

% If an array of info/group structs is passed, gather results from each
if numel(info) > 1
    for i = 1:numel(info)
        sub = collectDatasets(info(i));
        if ~isempty(sub)
            datasetPaths = [datasetPaths; sub(:)]; %#ok<AGROW>
        end
    end
    return
end

% Add datasets in this group (if any)
if isfield(info,'Datasets') && ~isempty(info.Datasets)
    for d = 1:numel(info.Datasets)
        datasetPaths{end+1,1} = [info.Name '/' info.Datasets(d).Name]; %#ok<AGROW>
    end
end

% Recurse into subgroups (if any)
if isfield(info,'Groups') && ~isempty(info.Groups)
    for g = 1:numel(info.Groups)
        sub = collectDatasets(info.Groups(g));
        if ~isempty(sub)
            datasetPaths = [datasetPaths; sub(:)]; %#ok<AGROW>
        end
    end
end

end
