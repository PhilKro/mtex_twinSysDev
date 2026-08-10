%% define .tdms filename
close all
clear all

% % Ti PTP_1
% ebsd_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0007_Ti_cAxis_PTP_1\h5oina';
% dic_dir  = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_1\test\DIC';
% dic_vid_dir = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_1\test\DIC\DIC continual';
% filename = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_1\test\AMICS\insitu tkd Ti_caxis_1 2025.10.06-11.04.44.989049_000.tdms';
% drift_time_val = 132;

% Ti PTP_2
ebsd_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0009_Ti_cAxis_PTP_2\H5OINA';
dic_dir  = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\insitu\DIC';
dic_vid_dir = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\insitu\DIC\continual';
filename = "C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2\TiPTP2_AMICSRT Log 2026.01.20-15.09.00.559488_000.tdms";
drift_time_val = 260;

% % Nb bicrystal
% ebsd_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0008_Nb_bicrystal_PTP1\insitu\h5oina';
% dic_dir  = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Nb_Bicrystal\Experiment\DIC';
% dic_vid_dir = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Nb_Bicrystal\Experiment\DIC\videos';
% filename = "N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Nb_Bicrystal\Experiment\AMICS\insitu tkd nb bicrystal 2025.11.12-15.32.28.895618_000.tdms";
% drift_time_val = 0;

% % Re PTP 3
% ebsd_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0005_Re_7deg_PTP_3\h5oina\Originals';
% dic_dir  = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20241218_Re_7dc_tkdPTP_3\04_insitu';
% dic_vid_dir = '';
% filename = "N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20241218_Re_7dc_tkdPTP_3\Alemnis\Re_7deg_3 2025.01.20-11.13.40.410848_000.tdms";
% drift_time_val = 0;

% Re cAxis PTP 2
% ebsd_dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0003_Re_cAxis_2\Experiment\Re_c_Axis_PTP_2\h5oina';
% dic_dir  = '';
% dic_vid_dir = '';
% filename = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20240913_Re_c_tkdPTP_2\Alemnis\Re_Caxis_PTP_2 2024.10.08-08.45.55.844284_000.tdms';
% drift_time_val = 60;

%% inspect the tdms file (optional)
info = tdmsinfo(filename);
disp(info.ChannelList);
%% read .tdms file

% Read all data (modern computers handle 380MB easily in memory)
data = tdmsread(filename);

relevantData = (data(:,2));
relevantData = relevantData{1,1};
%% analyse the data

disp('Column Names in your table:');
disp(relevantData.Properties.VariableNames);

t_s = relevantData.("Time [s]");
disp_m = relevantData.("Displacement [m]");
load_N = relevantData.("Load [N]");

% 2. Convert Units and Zero the Data
% Meters -> Microns (um), Newtons -> milliNewtons (mN)
% We subtract the first value (index 1) to make the curve start at (0,0)
disp_um = (disp_m - disp_m(1)) * 1e6;
load_mN = (load_N - load_N(1)) * 1e3;

% 3. Plot Load vs Displacement
figure;
plot(disp_um, load_mN, 'LineWidth', 1.5);
grid on;
xlabel('Displacement (\mu m)');
ylabel('Load (mN)');
title('Microcompression: Load-Displacement Curve');

% Optional: Reverse X-axis if the compression is going "negative"
% set(gca, 'XDir', 'reverse');
%% plot load over time

% 3. Plot Load vs Displacement
figure;
plot(t_s, load_mN, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Load (mN)');
title('Microcompression: Load-Time Curve');

%% correct load for drift
uncorrected_load_mN = load_mN;

drift_indices = t_s <= drift_time_val;


% drift_indices = t_s <= 260;

% Calculate raw load/disp again to ensure we don't double-correct if you re-run
load_mN = (relevantData.("Load [N]") - relevantData.("Load [N]")(1)) * 1e3;
disp_um = (relevantData.("Displacement [m]") - relevantData.("Displacement [m]")(1)) * 1e6;

% 3. Fit a linear slope to the Load (mN) vs Time (s) in that region
% polyfit returns [slope, intercept]
p = polyfit(t_s(drift_indices), load_mN(drift_indices), 1);
drift_slope = p(1); % This is in mN/s

fprintf('Drift was extracted from the first %f s of the load data.\n', drift_time_val);
fprintf('Drift rate detected: %.4f mN/s\n', drift_slope);

% 4. Correct the entire load variable
% We subtract the drift contribution (rate * time) from the whole signal
load_mN = load_mN - (drift_slope * t_s);

% 5. Re-zero the load curve (optional, but good practice after correction)
% This ensures the curve starts at 0 after the tilt correction
load_mN = load_mN - load_mN(1);

% 6. Plot the clean load/time curve
figure;
plot(t_s, load_mN, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Load (mN)');
if drift_time_val == 0
    title('Uncorrected Load-Time Curve');
else
    title('Drift-Corrected Load-Time Curve');
end

% 6. Plot the clean result
figure;
plot(disp_um, load_mN, 'LineWidth', 1.5);
grid on;
xlabel('Displacement (\mu m)');
ylabel('Load (mN)');
if drift_time_val == 0
    title('Uncorrected Load-Displacement Curve');
else
    title('Drift-Corrected Load-Displacement Curve');
end

%% 7. PLot displacement over lab time, with DIC, EBSD scans indicated
time_str = regexp(filename, '\d{4}\.\d{2}\.\d{2}-\d{2}\.\d{2}\.\d{2}', 'match');
start_time = datetime(time_str{1}, 'InputFormat', 'yyyy.MM.dd-HH.mm.ss');
actual_timestamps = start_time + seconds(t_s);

disp('First 5 actual timestamps:');
disp(actual_timestamps(1:5));

% 1. Define Directories
time_offset_seconds = 0; % Offset to shift the mechanical data timestamps (in seconds)
DIC_scantime = 65; % Duration of static DIC scans in seconds (set to 0 to auto-estimate)

% Verify directories exist (helps catch disconnected network drives)
if ~isfolder(ebsd_dir)
    warning('Cannot access EBSD directory. Check if network drive N: is connected: %s', ebsd_dir); 
end
if ~isfolder(dic_dir)
    warning('Cannot access DIC directory. Check if network drive N: is connected: %s', dic_dir); 
end
if ~isfolder(dic_vid_dir)
    warning('Cannot access DIC Video directory. Check if network drive N: is connected: %s', dic_vid_dir); 
end

% 2. Extract EBSD Timestamps (.h5oina)
ebsd_files = dir(fullfile(ebsd_dir, '*.h5oina'));
num_ebsd = length(ebsd_files);
ebsd_starts = NaT(num_ebsd, 1);
ebsd_ends   = NaT(num_ebsd, 1);
ebsd_names  = strings(num_ebsd, 1);
ebsd_hvs    = zeros(num_ebsd, 1);

for k = 1:num_ebsd
    fname = fullfile(ebsd_dir, ebsd_files(k).name);
    
    % Save name without extension for the plot label
    ebsd_names(k) = strrep(ebsd_files(k).name, '.h5oina', '');
    
    % Read read-only data from HDF5
    acq_date_raw = h5read(fname, '/1/EBSD/Header/Acquisition Date');
    acq_time_raw = h5read(fname, '/1/EBSD/Header/Acquisition Time');
    acq_hv_raw   = h5read(fname, '/1/EBSD/Header/Beam Voltage');
    
    if iscell(acq_date_raw), acq_date_raw = acq_date_raw{1}; end
    if iscell(acq_time_raw), acq_time_raw = acq_time_raw{1}; end
    if iscell(acq_hv_raw),   acq_hv_raw   = acq_hv_raw{1};   end
    
    start_t = datetime(string(acq_date_raw), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
    duration_s = double(string(acq_time_raw)); 
    
    % Handle HV over-precision (assume already in kV based on Oxford HDF5 standard)
    if isstring(acq_hv_raw) || ischar(acq_hv_raw)
        hv_val = str2double(string(acq_hv_raw));
    else
        hv_val = double(acq_hv_raw);
    end
    ebsd_hvs(k) = round(hv_val);
    
    ebsd_starts(k) = start_t;
    ebsd_ends(k)   = start_t + seconds(duration_s);
end

% 3. Extract DIC Timestamps (.hdr ONLY)
dic_files = dir(fullfile(dic_dir, '*.hdr'));
num_dic = length(dic_files);
dic_starts = NaT(num_dic, 1);
dic_names  = strings(num_dic, 1);
dic_hvs    = zeros(num_dic, 1);

for k = 1:num_dic
    hdr_fname = fullfile(dic_dir, dic_files(k).name);
    
    % Strip the ending to get a clean base name
    base_name = strrep(dic_files(k).name, '-tif.hdr', '');
    % If some files just end in .hdr, this will clean those too
    base_name = strrep(base_name, '.hdr', ''); 
    dic_names(k) = base_name; 
    
    % Read .hdr as plain text
    txt = fileread(hdr_fname);
    
    % Extract Date, Time, and HV
    d_tok  = regexp(txt, 'Date=([^\r\n]+)', 'tokens', 'once');
    t_tok  = regexp(txt, 'Time=([^\r\n]+)', 'tokens', 'once');
    hv_tok = regexp(txt, 'HV=([^\r\n]+)', 'tokens', 'once');
    
    start_t = datetime(strcat(d_tok{1}, " ", t_tok{1}), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    dic_starts(k) = start_t;
    
    % Convert HV (e.g., 20.000e3 -> 20000 V -> 20 kV)
    if ~isempty(hv_tok)
        dic_hvs(k) = round(str2double(hv_tok{1}) / 1000);
    else
        dic_hvs(k) = NaN;
    end
end

% Calculate Estimated DIC Duration
if DIC_scantime > 0
    est_duration = seconds(DIC_scantime);
    fprintf('Using user-defined DIC scan time of %.1f seconds.\n', DIC_scantime);
else
    if num_dic > 1
        % Sort times, find the differences, and pick the smallest gap
        sorted_starts = sort(dic_starts);
        time_gaps = diff(sorted_starts);
        est_duration = min(time_gaps);
        
        warning('DIC scan time was NOT read properly from files. It was estimated using the smallest time gap between consecutive scans (%s).', char(est_duration));
    elseif num_dic == 1
        % Fallback if there is only 1 DIC scan and we cannot calculate a gap
        est_duration = seconds(60); 
        warning('Only 1 DIC scan found. Cannot calculate time gaps. Defaulting to estimated duration of 60 seconds.');
    else
        est_duration = seconds(0);
    end
end

% Apply the estimated duration to find all end times
dic_ends = dic_starts - est_duration;

% 3.5 Extract DIC Video Timestamps (.hdr and .csv)
vid_files = dir(fullfile(dic_vid_dir, '*.hdr'));
num_vid = length(vid_files);
vid_starts = NaT(num_vid, 1);
vid_ends   = NaT(num_vid, 1);
vid_names  = strings(num_vid, 1);
vid_hvs    = zeros(num_vid, 1);

for k = 1:num_vid
    hdr_fname = fullfile(dic_vid_dir, vid_files(k).name);
    
    % Strip the ending to get a clean base name
    base_name = strrep(vid_files(k).name, '.hdr', '');
    vid_names(k) = base_name; 
    
    % Read .hdr as plain text
    txt = fileread(hdr_fname);
    
    % Extract Date, Time, and HV
    d_tok  = regexp(txt, 'Date=([^\r\n]+)', 'tokens', 'once');
    t_tok  = regexp(txt, 'Time=([^\r\n]+)', 'tokens', 'once');
    hv_tok = regexp(txt, 'AccVoltage=([^\r\n]+)', 'tokens', 'once');
    
    start_t = datetime(strcat(d_tok{1}, " ", t_tok{1}), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    vid_starts(k) = start_t;
    
    % Convert HV (e.g., 2.000000e+04 -> 20000 V -> 20 kV)
    if ~isempty(hv_tok)
        vid_hvs(k) = round(str2double(hv_tok{1}) / 1000);
    else
        vid_hvs(k) = NaN;
    end
    
    % Find associated CSV and calculate exact duration
    % Check both naming variations to be safe
    csv_fname1 = fullfile(dic_vid_dir, [vid_files(k).name, '_timestamps.csv']);
    csv_fname2 = fullfile(dic_vid_dir, [base_name, '_timestamps.csv']);
    
    if isfile(csv_fname1)
        target_csv = csv_fname1;
    elseif isfile(csv_fname2)
        target_csv = csv_fname2;
    else
        target_csv = '';
    end
    
    if ~isempty(target_csv)
        try
            vid_data = readtable(target_csv);
            % Duration = last frame rel_time_s - first frame rel_time_s
            duration_s = vid_data.rel_time_s(end) - vid_data.rel_time_s(1);
        catch
            warning('Could not read duration from %s. Defaulting to 10s.', target_csv);
            duration_s = 10;
        end
    else
        warning('Missing CSV for video. Tried %s and %s. Defaulting to 10s.', csv_fname1, csv_fname2);
        duration_s = 10; % Fallback width so the bar is visible
    end
    vid_ends(k) = start_t + seconds(duration_s);
end

% --- Combine and Sort all scans chronologically to detect voltage changes ---
all_starts = [ebsd_starts; dic_starts; vid_starts];
all_ends   = [ebsd_ends; dic_ends; vid_ends];
all_hvs    = [ebsd_hvs; dic_hvs; vid_hvs];

[sorted_starts, sort_idx] = sort(all_starts);
sorted_ends = all_ends(sort_idx);
sorted_hvs  = all_hvs(sort_idx);

% Apply Time Offset to Mechanical Data
if time_offset_seconds ~= 0
    shifted_timestamps = actual_timestamps + seconds(time_offset_seconds);
    warning('A time offset of %.1f seconds has been applied to the mechanical data timestamps.', time_offset_seconds);
else
    shifted_timestamps = actual_timestamps;
end

% 4. Plotting
figure('Position', [100, 100, 1200, 600]);

% --- LEFT Y-AXIS (Background patches and displacement line) ---
yyaxis left
hold on;

% Plot an invisible line first to correctly set the auto-limits for X and Y
plot(shifted_timestamps, disp_um, 'Color', 'none', 'HitTest', 'off');
y_limits_left = ylim;
yl = y_limits_left; 

% 5. Add Acquisition Bars and Text Labels (Drawn in background)
% EBSD Bars (Red, 20% opacity)
for k = 1:num_ebsd
    x_patch = [ebsd_starts(k) ebsd_ends(k) ebsd_ends(k) ebsd_starts(k)];
    y_patch = [yl(1) yl(1) yl(2) yl(2)];
    
    p_ebsd = patch('XData', x_patch, 'YData', y_patch, ...
        'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HitTest', 'off');
    
    % Add text label slightly offset from the start time, near the top
    text(ebsd_starts(k), yl(1) + 0.85 * (yl(2) - yl(1)), ebsd_names(k), ...
        'Rotation', 90, 'Interpreter', 'none', 'FontSize', 8, 'Color', 'r', 'HitTest', 'off');
end

% DIC Bars (Yellow/Orange, 20% opacity)
for k = 1:num_dic
    x_patch = [dic_starts(k) dic_ends(k) dic_ends(k) dic_starts(k)];
    y_patch = [yl(1) yl(1) yl(2) yl(2)];
    
    dic_color = [0.9290, 0.6940, 0.1250]; 
    
    p_dic = patch('XData', x_patch, 'YData', y_patch, ...
        'FaceColor', dic_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HitTest', 'off');
    
    % Add text label slightly offset from the start time, near the bottom
    text(dic_starts(k), yl(1) + 0.15 * (yl(2) - yl(1)), dic_names(k), ...
        'Rotation', 90, 'Interpreter', 'none', 'FontSize', 8, 'Color', dic_color, 'HitTest', 'off');
end

% DIC Video Bars (Green, 20% opacity)
for k = 1:num_vid
    x_patch = [vid_starts(k) vid_ends(k) vid_ends(k) vid_starts(k)];
    y_patch = [yl(1) yl(1) yl(2) yl(2)];
    
    vid_color = [0, 0.8, 0]; % Strong pure green
    
    p_vid = patch('XData', x_patch, 'YData', y_patch, ...
        'FaceColor', vid_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HitTest', 'off');
    
    % Add text label slightly offset from the start time, positioned higher up
    text(vid_starts(k), yl(1) + 0.70 * (yl(2) - yl(1)), vid_names(k), ...
        'Rotation', 90, 'Interpreter', 'none', 'FontSize', 8, 'Color', vid_color, 'HitTest', 'off');
end

% 6. Detect and Plot High Voltage Switches (Drawn in background)
for i = 2:length(sorted_hvs)
    prev_hv = sorted_hvs(i-1);
    curr_hv = sorted_hvs(i);
    
    if ~isnan(prev_hv) && ~isnan(curr_hv) && prev_hv ~= curr_hv
        % Find the midpoint between the previous scan ending and the current scan starting
        switch_time = sorted_ends(i-1) + (sorted_starts(i) - sorted_ends(i-1))/2;
        
        % Plot vertical dashed line
        xline(switch_time, '--k', 'LineWidth', 1.5, 'Alpha', 0.6, 'HitTest', 'off');
        
        % Place text next to the line
        label_str = sprintf(' %d kV to %d kV ', prev_hv, curr_hv);
        text(switch_time, yl(1) + 0.5 * (yl(2) - yl(1)), label_str, ...
            'Rotation', 90, 'Interpreter', 'none', 'FontSize', 10, ...
            'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'HitTest', 'off');
    end
end

% Now draw the actual Displacement curve (visually on top of the background items)
plot(shifted_timestamps, disp_um, '-', 'LineWidth', 1.5, 'Color', '#0072BD');
ylabel('Corrected Displacement (\mu m)');

% --- RIGHT Y-AXIS (Foreground) ---
yyaxis right
% Draw the Load curve (right axis is inherently overlaid on top of the left axis)
plot(shifted_timestamps, load_mN, '-', 'LineWidth', 1.5, 'Color', '#D95319');

if drift_time_val == 0
    ylabel('Uncorrected load (mN)');
else
    ylabel(sprintf('Corrected load (mN) [corrections used the first %.0f seconds of the micromech data]', drift_time_val));
end

hold off;

% 7. Formatting
title('Microcompression Timeline with In-Situ Acquisitions');

if time_offset_seconds ~= 0
    xlabel(sprintf('Lab Time [Shifted by %.1f s]', time_offset_seconds));
else
    xlabel('Lab Time');
end

xtickformat('HH:mm:ss');
grid on;

% Create a robust clean legend
legend_handles = [];
legend_labels = {};

if num_ebsd > 0 && exist('p_ebsd', 'var')
    legend_handles = [legend_handles, p_ebsd];
    legend_labels{end+1} = 'TKD (EBSD) Scans';
end
if num_dic > 0 && exist('p_dic', 'var')
    legend_handles = [legend_handles, p_dic];
    legend_labels{end+1} = 'DIC Scans';
end
if num_vid > 0 && exist('p_vid', 'var')
    legend_handles = [legend_handles, p_vid];
    legend_labels{end+1} = 'DIC Videos';
end

if ~isempty(legend_handles)
    legend(legend_handles, legend_labels, 'Location', 'northeast');
end

% =========================================================================
% 8. EXPORT TIMELINE DATA TO JSON
% =========================================================================
fprintf('\nExporting timeline data to JSON for external inspection...\n');
timelineData = struct();

% Downsample the continuous TDMS curve data by a factor to keep JSON size reasonable
% (Usually TDMS has 10,000+ points, we cap it at ~5000 points for the timeline)
ds_factor = max(1, floor(length(shifted_timestamps) / 5000));
idx_ds = 1:ds_factor:length(shifted_timestamps);

timelineData.tdms_curve = struct();
timelineData.tdms_curve.timestamps = cellstr(string(shifted_timestamps(idx_ds), 'yyyy-MM-dd HH:mm:ss.SSS'));
timelineData.tdms_curve.displacement_um = disp_um(idx_ds);
timelineData.tdms_curve.load_mN = load_mN(idx_ds);

timelineData.events = struct('type', {}, 'name', {}, 'start_time', {}, 'end_time', {}, 'hv_kV', {});
idx = 1;

for i = 1:num_ebsd
    timelineData.events(idx).type = 'EBSD';
    timelineData.events(idx).name = char(ebsd_names(i));
    timelineData.events(idx).start_time = char(string(ebsd_starts(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).end_time = char(string(ebsd_ends(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).hv_kV = ebsd_hvs(i);
    idx = idx + 1;
end

for i = 1:num_dic
    timelineData.events(idx).type = 'DIC_Scan';
    timelineData.events(idx).name = char(dic_names(i));
    timelineData.events(idx).start_time = char(string(dic_starts(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).end_time = char(string(dic_ends(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).hv_kV = dic_hvs(i);
    idx = idx + 1;
end

for i = 1:num_vid
    timelineData.events(idx).type = 'DIC_Video';
    timelineData.events(idx).name = char(vid_names(i));
    timelineData.events(idx).start_time = char(string(vid_starts(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).end_time = char(string(vid_ends(i), 'yyyy-MM-dd HH:mm:ss.SSS'));
    timelineData.events(idx).hv_kV = vid_hvs(i);
    idx = idx + 1;
end

% Save to JSON safely supporting older MATLAB versions
try
    json_text = jsonencode(timelineData, 'PrettyPrint', true);
catch
    json_text = jsonencode(timelineData);
end

[tdms_folder, tdms_base, ~] = fileparts(char(filename));
json_path = fullfile(tdms_folder, [tdms_base, '_timeline.json']);
fid = fopen(json_path, 'w');
if fid ~= -1
    fprintf(fid, '%s', json_text);
    fclose(fid);
    fprintf('  -> Successfully saved timeline JSON: %s\n', json_path);
else
    warning('Failed to open file for JSON export: %s', json_path);
end