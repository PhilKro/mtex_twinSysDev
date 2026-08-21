clear
close all
mtexdata twins silent
CS = ebsd('indexed').CS;

% reconstruct grains
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',5*degree);

grains = smooth(grains,1);
gB = grains.boundary;

% plot the grains
plot(grains,grains.meanOrientation)

%% Define {10-12} Twin System and initiate TwinNetwork
% K1 = {10-12}, eta1 = <-1011>
cs = grains(1).CS;
tS = twinSystem.hexagonal_1012(cs);

network = TwinNetwork(grains, gB, tS, ebsd, 'thresh', 3*degree);

%% 5. Execute Twin Network Topology Pipeline
disp('    Step 1: Building Initial Topology...');
G = network.buildInitialGraph();
network.plotGraph(G);
figG = gcf;
title(sprintf('[%s] Step 1: Initial Raw Topology'), 'Interpreter', 'none');
% exportScaledFigure(figG, fullfile(outDir, [defStep.stepName, '_01_Raw_Topology.png']), 'nodatestring');
% close(figG);

disp('    Step 2: Extracting Spanning Tree...');
mstG = network.extractSpanningTree(G);
network.plotGraph(mstG);
figMST = gcf;
title(sprintf('[%s] Step 2: Minimum Spanning Tree'), 'Interpreter', 'none');
% exportScaledFigure(figMST, fullfile(outDir, [defStep.stepName, '_02_Minimum_Spanning_Tree.png']), 'nodatestring');
% close(figMST);

disp('    Step 3: Resolving Node States and Variants...');
mstG = network.resolveNodeStates(mstG);
network.plotGraph(mstG);
figRes = gcf;
title(sprintf('[%s] Step 3: MST after FZ Propagation'), 'Interpreter', 'none');
% exportScaledFigure(figRes, fullfile(outDir, [defStep.stepName, '_03_Resolved_Node_States.png']), 'nodatestring');
% close(figRes);

% disp('    Step 4: Pruning and Rewiring Outliers...');
% mstG = network.pruneAndRewireOutliers(mstG, 'debug', false);
% network.plotGraph(mstG);
% figPrune = gcf;
% title(sprintf('[%s] Step 4: MST after Pruning and Rewiring'), 'Interpreter', 'none');
% exportScaledFigure(figPrune, fullfile(outDir, [defStep.stepName, '_04_Pruned_Rewired.png']), 'nodatestring');
% close(figPrune);

disp('    Step 5: Reducing to Quotient Graph...');
Q = network.reduceToQuotient(mstG);
network.plotGraph(Q,'noGrainsBG', true);
figQ = gcf;
title(sprintf('[%s] Step 5: Final Quotient Graph'), 'Interpreter', 'none');
% exportScaledFigure(figQ, fullfile(outDir, [defStep.stepName, '_05_Quotient_Graph.png']), 'nodatestring');
% close(figQ);

disp('    Step 6: Plotting grains with twin generation assignment...');
network.plot('twinsOnly', false);
figGen = gcf;
title(sprintf('[%s] Step 6: Twin Generations & Network Map'), 'Interpreter', 'none');
% exportScaledFigure(figGen, fullfile(outDir, [defStep.stepName, '_06_Twin_Generations.png']), 'nodatestring');
% close(figGen);

disp('    Step 6.5: Plotting grains with twin generation assignment with no text...');
network.plot('twinsOnly', false, 'noText', true, 'lineWidth', 2);
figGen = gcf;
title(sprintf('[%s] Step 6: Twin Generations & Network Map'), 'Interpreter', 'none');
% exportScaledFigure(figGen, fullfile(outDir, [defStep.stepName, '_06b_Twin_Generations_Nice.png']), 'nodatestring');
% close(figGen);


% close all; % Cleanup any remaining open handles
