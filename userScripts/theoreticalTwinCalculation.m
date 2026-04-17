%define crystal symmetry
CS_Mag = crystalSymmetry('6/mmm', [3.21,3.21,5.21]);

% calculate theoretical twin systems for the lattice
theoTwins = twinSystem.calculateTheoreticalTwins(CS_Mag, 2, 3);

% sort by shear magnitudes
shearValues = [theoTwins.shear];
[~, sortIdx] = sort(shearValues);
theoTwins_sorted = theoTwins(sortIdx);

% print out results
twinSystem.printTheoreticalTwins(theoTwins_sorted(1:40))