function check_twinNetwork
% Regression checks for local EBSD matching and virtual twin generations.

[network,bulkDeviation] = localDoubleTwinNetwork;

assert(bulkDeviation > 5*degree, ...
  'check_twinNetwork: the fixture bulk means accidentally match a twin sequence')

G = network.buildInitialGraph('maxGenStep',2);

assert(numnodes(G) == 3 && numedges(G) == 2, ...
  'check_twinNetwork: a double twin did not create one virtual node and two edges')
assert(nnz(G.Nodes.Type == "Physical") == 2 && ...
    nnz(G.Nodes.Type == "Virtual") == 1, ...
  'check_twinNetwork: physical and virtual node types are incorrect')
assert(all(isfinite(G.Edges.AngularDeviation)) && ...
    max(G.Edges.AngularDeviation) < 1e-6*degree, ...
  'check_twinNetwork: the exact local boundary relation was not matched')
assert(ismember('MeanEbsdOri1',G.Edges.Properties.VariableNames) && ...
    ismember('MeanEbsdOri2',G.Edges.Properties.VariableNames), ...
  'check_twinNetwork: local EBSD endpoint orientations were not retained')
assert(all(cellfun(@(s) isfinite(str2double(s)),G.Edges.Sequence)), ...
  'check_twinNetwork: the virtual chain lost its per-edge variant sequence')

T = network.extractSpanningTree(G);
T = network.resolveNodeStates(T);

assert(numnodes(T) == 3 && numedges(T) == 2, ...
  'check_twinNetwork: the virtual chain was not retained in the spanning tree')
assert(all(T.Nodes.State ~= "[]") && isequal(sort(T.Nodes.Depth),[0;1;2]), ...
  'check_twinNetwork: state resolution did not walk through the virtual node')
assert(all(isfinite(T.Edges.Variant)) && ...
    all(isfinite(T.Edges.AngularDeviation)), ...
  'check_twinNetwork: per-step variant locking failed')

Tpruned = network.pruneAndRewireOutliers(T, ...
  'outlierMode','theoretical','maxErrThreshold',15*degree);
assert(isa(Tpruned,'graph') && numnodes(Tpruned) == numnodes(T), ...
  'check_twinNetwork: the restored pruning pipeline failed on a valid network')

Q = network.reduceToQuotient(T);

assert(numnodes(Q) == 3 && numedges(Q) == 2, ...
  'check_twinNetwork: quotient reduction lost a twin generation')
assert(any(Q.Nodes.Area == 0), ...
  'check_twinNetwork: the purely virtual generation is absent from the quotient graph')
assert(all(isfinite(Q.Edges.Variant)), ...
  'check_twinNetwork: quotient edges lost their locked twin variants')

disp('check_twinNetwork: passed');

end

% -------------------------------------------------------------------------
function [network,bulkDeviation] = localDoubleTwinNetwork
% local boundary pixels match exactly while grain means differ by a gradient

cs = crystalSymmetry('6/mmm',[3.21 3.21 5.21], ...
  'mineral','synthetic magnesium');
tS = twinSystem.hexagonal_1012(cs);
variants = parentTwinMisorientation(symmetrise(tS));
improper = isImproper(variants);
variants(improper) = -variants(improper);

% choose a two-step operator distinct from both the parent and every
% primary variant
bestSeparation = 0;
pair = [0 0];
for i = 1:length(variants)
  for j = 1:length(variants)
    op = variants(i) * variants(j);
    separation = min([angle(op),min(angle(op,variants))]);
    if separation > bestSeparation
      bestSeparation = separation;
      pair = [i j];
    end
  end
end
assert(bestSeparation > 20*degree, ...
  'check_twinNetwork: no distinct double-twin fixture could be constructed')
doubleTwin = variants(pair(1)) * variants(pair(2));

% Each half has a 12 degree intragranular gradient. The two orientations
% directly adjacent to the boundary obey the double-twin relation exactly.
nRows = 8;
nCols = 8;
ori = orientation.nan(nRows,nCols,cs);
offset = [-12 -8 -4 0 0 4 8 12] * degree;
for col = 1:4
  value = orientation.byEuler(offset(col),0,0,cs);
  ori(:,col) = repmat(value,nRows,1);
end
for col = 5:8
  value = orientation.byEuler(offset(col),0,0,cs) * doubleTwin;
  ori(:,col) = repmat(value,nRows,1);
end

ebsd = EBSDsquare([],ori,2*ones(nRows,nCols),[0 1], ...
  {'notIndexed',cs},'dxy',[1 1]);
[grains,ebsd] = calcGrains(ebsd,'threshold',20*degree,'minPixel',1);

assert(length(grains) == 2, ...
  'check_twinNetwork: the synthetic gradient map did not reconstruct as two grains')

% Demonstrate that bulk-mean matching would reject this boundary.
nVariants = length(variants);
previous = repmat((1:nVariants)',nVariants,1);
next = repelem((1:nVariants)',nVariants,1);
operators = [variants(:);variants(previous).*variants(next)];
bulkDeviation = min(angle(grains(2).meanOrientation, ...
  grains(1).meanOrientation * operators));

network = TwinNetwork(grains,grains.boundary,tS,ebsd,'thresh',1*degree);

end
