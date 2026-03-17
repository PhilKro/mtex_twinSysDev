function fullPath = findDatasetPath(allPaths, requiredRelPath)
req = strsplit(requiredRelPath, '/');
reqTail = req(max(end-1,1):end);
fullPath = '';
for i = 1:numel(allPaths)
p = allPaths{i};
parts = strsplit(p, '/');
if length(parts) >= length(reqTail)
if isequal(parts(end-numel(reqTail)+1:end), reqTail)
fullPath = p;
return;
end
end
end
end