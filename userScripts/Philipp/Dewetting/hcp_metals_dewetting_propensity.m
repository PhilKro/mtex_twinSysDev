% Data for standard hcp metals
elements_hcp = {'Be', 'Mg', 'Sc', 'Ti', 'Co', 'Zn', 'Y', 'Zr', 'Tc', 'Ru', 'Cd', 'Hf', 'Re', 'Os', 'Tl', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Lu'}';
B_hcp = [130, 45, 57, 110, 180, 70, 41, 89, 151, 220, 42, 110, 370, 395, 43, 38, 38.7, 41, 40, 44, 44, 48]';
V_hcp = [4.85, 14.00, 15.00, 10.60, 6.67, 9.16, 19.80, 14.00, 8.63, 8.17, 13.00, 13.40, 8.86, 8.42, 17.20, 19.90, 19.30, 19.00, 18.70, 18.40, 18.10, 17.80]';
dH_hcp = [-609, -601, -636, -472, -238, -350, -635, -550, -219, -152, -258, -572, -224, -98, -131, -606, -621, -621, -626, -632, -629, -626]';
MeltingPoint_hcp = [1560, 923, 1814, 1941, 1768, 693, 1799, 2128, 2430, 2607, 594, 2506, 3459, 3306, 577, 1585, 1629, 1680, 1747, 1802, 1818, 1925]';

% Data for non-hcp metals (Gold, Nickel, Molybdenum)
elements_nonhcp = {'Au', 'Ni', 'Mo'}';
B_nonhcp = [180, 180, 230]';
V_nonhcp = [10.21, 6.59, 9.38]';
% dH values calculated per oxygen atom for Au2O3, NiO, and MoO3
dH_nonhcp = [-4.3, -240.0, -248.3]'; 
MeltingPoint_nonhcp = [1337, 1728, 2896]';

% Combine arrays for the full dataset
elements = [elements_hcp; elements_nonhcp];
B = [B_hcp; B_nonhcp];
V = [V_hcp; V_nonhcp];
dH = [dH_hcp; dH_nonhcp];
MeltingPoint_K = [MeltingPoint_hcp; MeltingPoint_nonhcp];

% Calculations
n = sqrt(B ./ V);
dH_eV = dH / 96.485; % Conversion factor: 1 eV/particle is roughly 96.485 kJ/mol

% Create and display a formatted MATLAB table
CrystalStructure = [repmat({'hcp'}, length(elements_hcp), 1); {'fcc'; 'fcc'; 'bcc'}];
metal_table = table(elements, CrystalStructure, B, V, n, dH, dH_eV, MeltingPoint_K, ...
    'VariableNames', {'Element', 'Structure', 'BulkModulus', 'MolarVolume', 'n', 'Enthalpy_kJ', 'Enthalpy_eV', 'MeltingPoint_K'});
disp(metal_table);

% Define marker sizes proportional to melting point (scaled for visibility)
sz_hcp = MeltingPoint_hcp / 10;
sz_nonhcp = MeltingPoint_nonhcp / 10;

% Create the plot
figure('Name', 'dH vs n for Metals', 'Color', 'w');

% Left Y-Axis (kJ/mol)
yyaxis left

% Plot hcp metals (blue) with size mapped to melting point
scatter(n(1:length(elements_hcp)), dH_hcp, sz_hcp, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'k');
hold on;
% Plot non-hcp metals (red) with size mapped to melting point
scatter(n(length(elements_hcp)+1:end), dH_nonhcp, sz_nonhcp, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
grid on;
box on;

% Label all data points with element symbols
for i = 1:length(elements)
    text(n(i) + 0.12, dH(i), elements{i}, 'FontSize', 10, 'VerticalAlignment', 'middle');
end
hold off;

ylabel('\DeltaH per Oxygen Atom (kJ/mol)');
ylim_left = ylim; % Capture the automatically generated limits

% Right Y-Axis (eV per O atom)
yyaxis right
% Scale limits to match left axis exactly
ylim(ylim_left / 96.485);
ylabel('\DeltaH per Oxygen Atom (eV)');

% Format axes colors to black to avoid default blue/red color mismatch
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% Final plot formatting and legend
title('\DeltaH of Oxide Formation vs n (Marker Area scales with Melting Point)');
xlabel('Parameter n = (B/V)^{0.5}');
legend('hcp Metals', 'Non-hcp Metals', 'Location', 'southeast');