%% MTEX 7 twinSystem and TwinNetwork validation
% Run this script section by section to inspect the intermediate objects.

clearvars
close all

makePlots = true;
runPruning = true;
plotSigmaSections = false;
maxGeneration = 2;
twinThreshold = 3 * degree;

%% Load magnesium EBSD data and reconstruct grains

mtexdata twins silent
ebsd = ebsd('indexed');
[grains,ebsd] = calcGrains(ebsd,'threshold',5*degree);

% Use the MTEX 7 boundary-smoothing API (including staircase removal).
grains = smoothBoundary(grains,1);
gB = grains.boundary;
cs = grains(1).CS;

fprintf('Reconstructed %d grains and %d boundary segments.\n', ...
  length(grains),length(gB));

if makePlots
  figure
  plot(grains,grains.meanOrientation)
  title('Reconstructed magnesium grains')
end

%% Inspect the extension-twin family

tS = twinSystem.hexagonal_1012(cs);
[variants,sourceId] = symmetrise(tS);
twinMisorientation = parentTwinMisorientation(variants);

% EBSD orientations use a proper representative of the twin law.
improper = isImproper(twinMisorientation);
twinMisorientation(improper) = -twinMisorientation(improper);

variantSummary = table(reshape(variants.variantId,[],1), ...
  reshape(sourceId,[],1),reshape(variants.shearMagnitude,[],1), ...
  reshape(angle(twinMisorientation)/degree,[],1), ...
  'VariableNames',{'Variant','SourceMode','Shear','MisorientationDeg'});
disp(variantSummary)

assert(all(isfinite(variantSummary.Shear)))
assert(all(~isImproper(twinMisorientation)))

%% Exercise specimen-frame mechanics

parentOrientation = grains(1).meanOrientation;
specimenVariants = parentOrientation * variants;
sigma = stressTensor.uniaxial(vector3d.Z);

schmid = SchmidFactor(specimenVariants,sigma);
work = interactionWork(specimenVariants,sigma);

mechanicsSummary = table(reshape(variants.variantId,[],1), ...
  reshape(schmid,[],1),reshape(work,[],1), ...
  'VariableNames',{'Variant','SchmidFactor','InteractionWork'});
disp(mechanicsSummary)

[~,favouredVariant] = max(work);
fprintf('Largest interaction work: variant %d.\n', ...
  variants.variantId(favouredVariant));

variantPath = [1,min(2,length(variants))];
Fcrystal = calcDeformationSequence(variants,variantPath);
Fspecimen = calcDeformationSequence(variants,variantPath,parentOrientation);

fprintf('Two-step deformation determinants: crystal %.8f, specimen %.8f.\n', ...
  det(Fcrystal),det(Fspecimen));

assert(abs(det(Fcrystal)-det(Fspecimen)) < 1e-10)

%% Calculate one theoretical Type-II mode

csCubic = crystalSymmetry('m-3m',[3.6 3.6 3.6], ...
  'mineral','synthetic fcc');
etaTypeII = Miller(1,1,1,csCubic,'uvw');
[typeIISystems,typeIIDetails] = twinSystem.calculateTheoreticalTwins( ...
  csCubic,etaTypeII,1:2,'FCC','type2',1);

fprintf('The [111] FCC Type-II search returned %d mode(s).\n', ...
  length(typeIISystems));
if ~isempty(typeIIDetails)
  twinSystem.printTheoreticalTwins(typeIIDetails)
end

assert(~isempty(typeIISystems))
assert(all(ismember(typeIISystems.twinType,[0 2])))

%% Build one- and two-generation graphs

network = TwinNetwork(grains,gB,tS,ebsd,'thresh',twinThreshold);
Gprimary = network.buildInitialGraph('maxGenStep',1);
G = network.buildInitialGraph('maxGenStep',maxGeneration);

printGraphSummary('primary graph',Gprimary)
printGraphSummary('multi-generation graph',G)

numVirtual = nnz(G.Nodes.Type == "Virtual");
fprintf('The multi-generation graph contains %d virtual node(s).\n', ...
  numVirtual);

if makePlots
  network.plotGraph(Gprimary)
  title('Primary-twin graph')

  network.plotGraph(G)
  title(sprintf('Twin graph through generation %d',maxGeneration))
end

%% Compare local-boundary and grain-mean deviations

bulkDeviation = calcBulkDeviation(network,Gprimary);
localDeviation = Gprimary.Edges.AngularDeviation;
isPhysicalEdge = isfinite(bulkDeviation);

deviationComparison = table((1:numedges(Gprimary))', ...
  localDeviation/degree,bulkDeviation/degree, ...
  (bulkDeviation-localDeviation)/degree, ...
  'VariableNames',{'Edge','LocalDeg','GrainMeanDeg','MeanMinusLocalDeg'});
deviationComparison = deviationComparison(isPhysicalEdge,:);

[~,order] = sort(abs(deviationComparison.MeanMinusLocalDeg),'descend');
disp(deviationComparison(order(1:min(10,height(deviationComparison))),:))

fprintf('Median local deviation: %.3f degrees.\n', ...
  median(deviationComparison.LocalDeg));
fprintf('Median grain-mean deviation: %.3f degrees.\n', ...
  median(deviationComparison.GrainMeanDeg));

if makePlots
  figure
  scatter(deviationComparison.GrainMeanDeg, ...
    deviationComparison.LocalDeg,30,'filled')
  hold on
  lim = max([xlim ylim]);
  plot([0 lim],[0 lim],'k--')
  hold off
  axis equal
  xlim([0 lim])
  ylim([0 lim])
  xlabel('Grain-mean deviation (degree)')
  ylabel('Local-boundary deviation (degree)')
  title('Local versus grain-mean twin matching')
end

%% Resolve twin states through the graph

T = network.extractSpanningTree(G);
T = network.resolveNodeStates(T);
printGraphSummary('resolved spanning tree',T)

virtualRows = find(T.Nodes.Type == "Virtual");
nodeRows = unique([reshape(1:min(12,numnodes(T)),[],1);virtualRows]);
disp(T.Nodes(nodeRows,{'Type','GrainId','State','Depth'}))

virtualEdges = find(any(ismember(T.Edges.EndNodes,virtualRows),2));
edgeRows = unique([reshape(1:min(12,numedges(T)),[],1);virtualEdges]);
disp(T.Edges(edgeRows, ...
  {'EndNodes','Sequence','Variant','AngularDeviation'}))

assert(all(T.Nodes.State ~= "[]"))
assert(all(isfinite(T.Edges.Variant)))

if makePlots
  network.plotGraph(T)
  title('Resolved twin spanning tree')
end

%% Exercise robust orientation centering

[~,largestIndex] = max(grains.area);
largestGrainId = grains(largestIndex).id;
grainPixels = ebsd(ebsd.grainId == largestGrainId).orientations;
grainPixels = grainPixels(1:min(30,length(grainPixels)));

if length(grainPixels) > 2
  outlier = grainPixels(1) * network.variants(1);
  centerInput = [grainPixels(:);outlier];
  robustCenter = network.computeRobustCenter(centerInput);
  ordinaryCenter = mean(centerInput);

  fprintf('Robust versus ordinary center: %.3f degrees.\n', ...
    angle(robustCenter,ordinaryCenter)/degree);
end

%% Prune and reduce to the quotient graph

if runPruning
  Tfinal = network.pruneAndRewireOutliers(T, ...
    'outlierMode','theoretical','maxErrThreshold',15*degree, ...
    'debug',false);
else
  Tfinal = T; %#ok<UNRCH>
end

Q = network.reduceToQuotient(Tfinal);
printGraphSummary('quotient graph',Q)

if makePlots
  network.plotGraph(Q,'noGrainsBG',true)
  title('Twin quotient graph')

  network.plot('twinsOnly',false,'noText',true,'lineWidth',2)
  title('Twin generations in the EBSD map')
end

if plotSigmaSections
  network.plotSigmaSections(Tfinal) %#ok<UNRCH>
end

disp('MTEX7TwinningValidation: completed')

%% Local helpers

function deviation = calcBulkDeviation(network,G)

deviation = nan(numedges(G),1);
for edgeId = 1:numedges(G)
  nodes = G.Edges.EndNodes(edgeId,:);
  if any(G.Nodes.Type(nodes) ~= "Physical"), continue; end

  sourceId = G.Nodes.GrainId(nodes(1));
  targetId = G.Nodes.GrainId(nodes(2));
  source = network.grains(network.grains.id2ind(sourceId)).meanOrientation;
  target = network.grains(network.grains.id2ind(targetId)).meanOrientation;
  variant = str2double(G.Edges.Sequence{edgeId});

  deviation(edgeId) = angle(target,source * network.variants(variant));
end

end

function printGraphSummary(label,G)

if ismember('Type',G.Nodes.Properties.VariableNames)
  physical = nnz(G.Nodes.Type == "Physical");
  virtual = nnz(G.Nodes.Type == "Virtual");
else
  physical = NaN;
  virtual = NaN;
end

fprintf('%s: %d nodes, %d edges, %g physical, %g virtual.\n', ...
  label,numnodes(G),numedges(G),physical,virtual);

end
