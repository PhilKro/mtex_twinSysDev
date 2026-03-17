% Define schematic time points (scaled for visual clarity) 
% Assume starting and ending at room temperature (25 C)
t = [0, 1, 2, 3, 4, 5.5, 9, 11];
T = [25, 250, 250, 500, 500, 1000, 1000, 25];

figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Plot the main heating curve
plot(t, T, 'b-', 'LineWidth', 2.5);
hold on;

% Add points to highlight the transitions between ramps and holds
scatter(t, T, 60, 'k', 'filled'); 

% Format axes
title('Schematic Heating Curve for Thermal Processing', 'FontSize', 14);
ylabel('Temperature (^\circC)', 'FontSize', 12);
xlabel('Process Timeline (Schematic Scale)', 'FontSize', 12);
set(gca, 'XTick', [], 'FontSize', 11); % Hide x-ticks since the time scale is arbitrary
grid on;
ylim([0, 1150]); % Leave room at the top for the 1000 C annotation

% --- Add Annotations ---

% 1. Initial Ramp
text(0.4, 150, 'Ramp 2^\circC/min', 'HorizontalAlignment', 'right', 'FontSize', 11);

% 2. First Hold (250 C)
text(1.5, 330, {'Hold: 250^\circC for 20 min', '(step necessary if sample was in atmosphere)'}, ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

% 3. Second Ramp
text(2.4, 390, 'Ramp 2^\circC/min', 'HorizontalAlignment', 'right', 'FontSize', 11);

% 4. Second Hold (500 C)
text(3.5, 600, {'Hold: 500^\circC for 30 min', '(flow of forming gas to get rid of oxide)'}, ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

% 5. Third Ramp
text(4.6, 770, 'Ramp to 1000^\circC', 'HorizontalAlignment', 'right', 'FontSize', 11);

% 6. Final Hold (1000 C)
text(7.25, 1080, {'Hold: 1000^\circC for \geq 12h minimum'}, ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold', 'Margin', 3);

% 7. Cool Down
text(10.2, 500, 'Cool down: 2^\circC/min', 'HorizontalAlignment', 'left', 'FontSize', 11);

hold off;