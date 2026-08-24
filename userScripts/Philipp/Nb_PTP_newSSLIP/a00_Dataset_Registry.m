function [registry, defaults] = a00_Dataset_Registry()
% a00_Dataset_Registry
% Master configuration and dataset mapping registry for Nb_PTP new SSLIP workflows.
%
% Returns:
%   registry - Array of structs defining associated EBSD, SE, and DIC files per step.
%   defaults - Common parameters and root directories for output workflows.

    %% Global Defaults & Directory Paths
    defaults.dir            = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Nb_PTP_newSSLIP';
    defaults.SEImagefolder  = 'DIC_SE_Images';
    defaults.DICfolder      = 'DIC_h5\Processing_subset=25_step=8';
    defaults.sel_points_fname = 'insitu Specimen 1 0 Map Data 1_BW123_20260819_151013_EBSD_alignment_data_EBSD.mat';
    
    % Output folder
    defaults.outDir_undeformed = 'Processed_Output_Undeformed';
    
    % EBSD/TKD processing defaults (from Undeformed_Script01)
    defaults.startx         = 1280;
    defaults.starty         = 1;
    defaults.numx           = 1280; 
    defaults.numy           = 1280; 
    defaults.bandWidth      = 263;
    defaults.tkd_tilt       = -10 * degree;
    defaults.thresh_ci      = 0.01;
    defaults.se_pixel_size  = 5.4687e-3; % um

    fprintf('========================================================================\n');
    fprintf('Loaded Dataset Registry for Nb_PTP_newSSLIP\n');
    fprintf('========================================================================\n');

    %% Dataset Registry
    % Define the 14 increments observed in the original processed output.
    % The undeformed reference EBSD is loaded at Step 01.
    % Note: I guessed the SE filenames based on standard naming. Adjust if needed.
    registry = [ ...
        struct('stepName', 'Step_01', 'ebsd_fname', 'insitu Specimen 1 0 Map Data 1', 'se_fname', '0_10kv_DIC.tif', 'dic_fname', 'B0001.hdf5', 'isReference', true), ...
        struct('stepName', 'Step_02', 'ebsd_fname', '', 'se_fname', '1_10kv_DIC.tif', 'dic_fname', 'B0002.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_03', 'ebsd_fname', '', 'se_fname', '2_10kv_DIC.tif', 'dic_fname', 'B0003.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_04', 'ebsd_fname', '', 'se_fname', '3_10kv_DIC.tif', 'dic_fname', 'B0004.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_05', 'ebsd_fname', '', 'se_fname', '4_10kv_DIC.tif', 'dic_fname', 'B0005.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_06', 'ebsd_fname', '', 'se_fname', '5_10kv_DIC.tif', 'dic_fname', 'B0006.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_07', 'ebsd_fname', '', 'se_fname', '6_10kv_DIC.tif', 'dic_fname', 'B0007.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_08', 'ebsd_fname', '', 'se_fname', '7_10kv_DIC.tif', 'dic_fname', 'B0008.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_09', 'ebsd_fname', '', 'se_fname', '8_10kv_DIC.tif', 'dic_fname', 'B0009.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_10', 'ebsd_fname', '', 'se_fname', '9_10kv_DIC.tif', 'dic_fname', 'B0010.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_11', 'ebsd_fname', '', 'se_fname', '10_10kv_DIC.tif', 'dic_fname', 'B0011.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_12', 'ebsd_fname', '', 'se_fname', '11_10kv_DIC.tif', 'dic_fname', 'B0012.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_13', 'ebsd_fname', '', 'se_fname', '12_10kv_DIC.tif', 'dic_fname', 'B0013.hdf5', 'isReference', false), ...
        struct('stepName', 'Step_14', 'ebsd_fname', '', 'se_fname', '13_10kv_DIC.tif', 'dic_fname', 'B0014.hdf5', 'isReference', false) ...
    ];
end
