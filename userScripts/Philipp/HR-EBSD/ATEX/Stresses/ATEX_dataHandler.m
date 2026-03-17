filename = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0005_Re_7deg_PTP_3\ATEX\Re_7_deg_insitu Specimen 1 4b Map Data 4_AVERAGED\Re_7deg_insitu Specimen 1 4b Map Data 4_AVERAGED_EXPORT_PIXELS.txt';
% requested_columns = {"Stress σ₁₁ in sample frame","Stress σ₁₂ in sample frame",...
%     "Stress σ₁₃ in sample frame","Stress σ₂₂ in sample frame",...
%     "Stress σ₂₃ in sample frame","Stress σ₃₃ in sample frame",... 
%     "Elastic strain ε₁₁ in sample frame","Elastic strain ε₁₂ in sample frame",...
%     "Elastic strain ε₁₃ in sample frame","Elastic strain ε₂₂ in sample frame",...
%     "Elastic strain ε₂₃ in sample frame","Elastic strain ε₃₃ in sample frame",...
%     "Rotation ω₁₃ (/Y) in sample frame ","Indexing-derived rotation ω₁₃ (/Y)",...
%     "Rotation ω₂₁ (/Z) in sample frame ","Indexing-derived rotation ω₂₁ (/Z)",...
%     "Rotation ω₃₂ (/X) in sample frame ","Indexing-derived rotation ω₃₂ (/X)"};
strain_columns = {    "Elastic strain ε₁₁ in sample frame","Elastic strain ε₁₂ in sample frame",...
    "Elastic strain ε₁₃ in sample frame","Elastic strain ε₂₂ in sample frame",...
    "Elastic strain ε₂₃ in sample frame","Elastic strain ε₃₃ in sample frame"};
stress_columns = {"Stress σ₁₁ in sample frame","Stress σ₁₂ in sample frame",...
    "Stress σ₁₃ in sample frame","Stress σ₂₂ in sample frame",...
    "Stress σ₂₃ in sample frame","Stress σ₃₃ in sample frame"};
rotation_columns = {"Rotation ω₁₃ (/Y) in sample frame ","Indexing-derived rotation ω₁₃ (/Y)",...
    "Rotation ω₂₁ (/Z) in sample frame ","Indexing-derived rotation ω₂₁ (/Z)",...
    "Rotation ω₃₂ (/X) in sample frame ","Indexing-derived rotation ω₃₂ (/X)"};

%data_struct = read_atex_data(filename, requested_columns);
ATEX_strains = read_atex_data(filename, strain_columns);
ATEX_stress = read_atex_data(filename, stress_columns);
ATEX_rotation = read_atex_data(filename,rotation_columns);
%fields = fieldnames(data_struct);

% Access data:
%stress = data_struct.Stress_σ₁₁_in_sample_frame;
%strain = data_struct.Elastic_strain_ε₁₁_in_sample_frame;
%%
%what to plot
data_struct = ATEX_stress;

%
fields = fieldnames(data_struct);
% Loop through each field and plot the corresponding matrix
for i = 1:length(fields)
    % Get the field name and data
    field_name = fields{i};
    data = data_struct.(field_name);
    
    % Create a new figure for each plot
    figure;
    
    % Plot using imagesc for a 2D colormap
    imagesc(data);
    
    % Add colorbar to indicate the values of the colormap
    colorbar;
    
    % Set title and labels
    title(field_name);
    xlabel('X');
    ylabel('Y');
    
    % Optionally set the axis to be equal for better visual comparison
    axis equal;
    
    % Optional: adjust the colormap (e.g., 'jet', 'parula', 'hot', etc.)
    colormap('hot');
end
