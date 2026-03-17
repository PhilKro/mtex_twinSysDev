% Description of a, c, m and r plane
% https://www.globalsino.com/EM/image2/3917a.GIF
load('CS_cS_andPlanes.mat')

%% Define Twinning and Slip Systems

% 11-21 Extension Twin
tS_1121_Ru = twinningSystem.hexagonal_1121_0001(CS_Ru);
tS_1121_Ru_sym = tS_1121_Ru.symmetrise;

% 11-22 Contraction Twin
tS_1122_Ru = twinningSystem.hexagonal_1122_1124(CS_Ru);
tS_1122_Ru_sym = tS_1122_Ru.symmetrise;

% 10-11 Contraction Twin
tS_1011_Ru = twinningSystem.hexagonal_1101_1013(CS_Ru);
tS_1011_Ru_sym = tS_1011_Ru.symmetrise;

% Slip Systems
basal_slip_Ru = slipSystem.basal(CS_Ru);
basal_slip_Ru_sym = basal_slip_Ru.symmetrise;

prism_slip_Ru = slipSystem.prismaticA(CS_Ru);
prism_slip_Ru_sym = prism_slip_Ru.symmetrise;

%% Define Planes
aplane = Miller(1,1,-2,0,'hkil', CS_Ru);
cplane = Miller(0,0,0,1,'hkil', CS_Ru);
mplane = Miller(1,-1,0,0,'hkil', CS_Ru);
rplane = Miller(1,-1,0,2,'hkil', CS_Ru);

%% Define Stress Tensors (Uniaxial Compression Normalized)
a_Compression_stress = -stressTensor.uniaxial(aplane) / norm(stressTensor.uniaxial(aplane));
c_Compression_stress = -stressTensor.uniaxial(cplane) / norm(stressTensor.uniaxial(cplane));
m_Compression_stress = -stressTensor.uniaxial(mplane) / norm(stressTensor.uniaxial(mplane));
r_Compression_stress = -stressTensor.uniaxial(rplane) / norm(stressTensor.uniaxial(rplane));

%% Calculate All Schmid Factors

stresses = [a_Compression_stress, c_Compression_stress, m_Compression_stress, r_Compression_stress];
stress_names = {'a-axis', 'c-axis', 'm-axis', 'r-axis'};

systems = {tS_1121_Ru_sym, tS_1122_Ru_sym, tS_1011_Ru_sym, basal_slip_Ru_sym, prism_slip_Ru_sym};
system_names = {'11-21 Twin', '11-22 Twin', '10-11 Twin', 'Basal Slip', 'Prismatic Slip'};

max_SF = zeros(length(stresses), length(systems));
all_SF_data = cell(length(stresses), length(systems)); 

for i = 1:length(stresses)
    for j = 1:length(systems)
        current_sys = systems{j};
        sf_temp = zeros(length(current_sys), 1);
        is_twin = contains(system_names{j}, 'Twin');
        
        for k = 1:length(current_sys)
            if is_twin
                sf_temp(k) = stresses(i).calcShearStress(current_sys(k).k1, current_sys(k).eta1);
            else
                sf_temp(k) = stresses(i).calcShearStress(current_sys(k).n, current_sys(k).b);
            end
        end
        
        % Filter for active (positive) variants
        if is_twin
            valid_sfs = sf_temp(sf_temp > 1e-4); 
        else
            valid_sfs = abs(sf_temp);
            valid_sfs = valid_sfs(valid_sfs > 1e-4); 
        end
        
        all_SF_data{i, j} = valid_sfs;
        
        if isempty(valid_sfs)
            max_SF(i, j) = 0;
        else
            max_SF(i, j) = max(valid_sfs);
        end
    end
end

%% Plot the Results

figure;
% Plot the background bars (transparent) to show the maximums
b = bar(max_SF, 'grouped', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Customize colors 
colors = [0.2 0.6 0.8;  % Light Blue for 11-21 Twin
          0.1 0.3 0.6;  % Dark Blue for 11-22 Twin
          0.4 0.1 0.6;  % Purple for 10-11 Twin
          0.8 0.4 0.2;  % Orange for Basal
          0.4 0.8 0.4]; % Green for Prismatic

for j = 1:length(systems)
    b(j).FaceColor = colors(j, :);
end

hold on;

% Define the maximum horizontal spread for the jitter
jitter_width = 0.12; 

% Overlay every valid variant as a scatter point with jitter
for i = 1:length(stresses)
    for j = 1:length(systems)
        % Get the X-coordinate of the center of this specific bar
        x_center = b(j).XEndPoints(i);
        
        % Get the SF data points for this group
        sf_vals = all_SF_data{i, j};
        
        % Plot the points if they exist
        if ~isempty(sf_vals)
            % Generate random horizontal offsets
            num_points = length(sf_vals);
            x_offsets = (rand(1, num_points) - 0.5) * jitter_width;
            x_jittered = x_center + x_offsets;
            
            % Plot the jittered points
            scatter(x_jittered, sf_vals, 35, colors(j, :), 'filled', 'MarkerEdgeColor', 'k');
        end
    end
end

% Format the plot
set(gca, 'XTickLabel', stress_names, 'FontSize', 12);
legend(b, system_names, 'Location', 'northeast'); 
ylabel('Schmid Factor (All Active Variants)', 'FontSize', 12);
xlabel('Compression Direction', 'FontSize', 12);
title('Distribution of Schmid Factors for Ru under Uniaxial Compression', 'FontSize', 14);
grid on;
ylim([0 0.55]);

% Add vertical dividers
xline(1.5, 'k:', 'LineWidth', 1);
xline(2.5, 'k:', 'LineWidth', 1);
xline(3.5, 'k:', 'LineWidth', 1);

hold off;