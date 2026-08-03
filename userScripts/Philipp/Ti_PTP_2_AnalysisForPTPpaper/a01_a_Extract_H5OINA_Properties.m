%% a01_a_Extract_H5OINA_Properties.m
% Extracts Band Slope, Band Contrast, and synthesizes FSE, BSE, MSE 
% from Kikuchi patterns in Oxford .h5oina files.

clear; close all; clc;

%% 1. Configuration
angDir = 'C:\Users\phkr\OneDrive - empa.ch\MTEX_Data\Philipp\Ti_PTP\Ti_cAxis_PTP_2';
h5Dir = 'N:\EBSD-Analysis-2\Philipp_EBSD_analysis2\0009_Ti_cAxis_PTP_2\H5OINA';
outDir = fullfile(angDir, 'Masks'); % Save alongside masks
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Define Crystal Symmetry for checking EBSD length
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68],'x||a','mineral','Ti');

h5Files = dir(fullfile(h5Dir, '**', '*.h5oina'));
if isempty(h5Files)
    error('No .h5oina files found in %s', h5Dir);
end

chunk_size = 1000; % Number of patterns to read at once (~300 MB per chunk)
verbose = true; % Set to false to hide chunk timing printouts

%% 2. Processing Loop
for i = 1:1%length(h5Files)
    h5Path = fullfile(h5Files(i).folder, h5Files(i).name);
    [~, baseName, ~] = fileparts(h5Files(i).name);
    
    fprintf('\n--- Processing %d/%d: %s ---\n', i, length(h5Files), h5Files(i).name);
    
    % Find corresponding .ang file
    % Assuming the .ang file has the same basename but might have a suffix like _BW158
    angFiles = dir(fullfile(angDir, '**', [baseName, '*.ang']));
    if isempty(angFiles)
        warning('No .ang file found matching %s. Skipping.', baseName);
        continue;
    end
    
    angPath = fullfile(angFiles(1).folder, angFiles(1).name);
    angBaseName = angFiles(1).name(1:end-4);
    outMat = fullfile(outDir, sprintf('%s_h5oina_props.mat', angBaseName));
    
    if exist(outMat, 'file')
        fprintf('Properties already extracted for %s. Skipping...\n', h5Files(i).name);
        continue;
    end
    
    % Load EBSD to check length
    fprintf('Loading EBSD from %s to verify point count...\n', angFiles(1).name);
    ebsd = EBSD.load(angPath, CS, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    num_ebsd_points = length(ebsd);
    
    % Get info from HDF5
    info = h5info(h5Path, '/1/EBSD/Data/Processed Patterns');
    dataSize = info.Dataspace.Size; % MATLAB returns dimensions reversed from C/HDFView
    
    % In HDFView: N x H x W (e.g., 86546 x 512 x 622)
    % In MATLAB: W x H x N (e.g., 622 x 512 x 86546)
    if length(dataSize) ~= 3
        error('Processed Patterns dataset does not have 3 dimensions (size is %s).', mat2str(dataSize));
    end
    
    W = dataSize(1);
    H = dataSize(2);
    N = dataSize(3);
    
    fprintf('HDF5 Pattern Size: %d points, %d (H) x %d (W)\n', N, H, W);
    
    if N ~= num_ebsd_points
        error('Point count mismatch! HDF5 contains %d points, but ANG file contains %d points.', N, num_ebsd_points);
    end
    
    % Read Band Slope and Band Contrast
    fprintf('Reading Band Slope and Band Contrast...\n');
    try
        BS = h5read(h5Path, '/1/EBSD/Data/Band Slope');
        BC = h5read(h5Path, '/1/EBSD/Data/Band Contrast');
    catch ME
        warning(ME.identifier, 'Could not read Band Slope or Band Contrast: %s', ME.message);
        BS = []; BC = [];
    end
    
    % Flatten them to 1D column vectors
    if ~isempty(BS), BS = double(BS(:)); end
    if ~isempty(BC), BC = double(BC(:)); end
    
    % --- NEW PYTHON BRIDGE TO BYPASS MATLAB HDF5 PLUGIN MISSING ---
    fprintf('Reading patterns via Python (bypassing MATLAB HDF5 plugin error)...\n');
    
    pyScriptPath = fullfile(outDir, 'extract_patterns.py');
    fid = fopen(pyScriptPath, 'w');
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'import h5py\n');
    fprintf(fid, 'import hdf5plugin\n');
    fprintf(fid, 'import numpy as np\n');
    fprintf(fid, 'import time\n');
    fprintf(fid, 'from scipy.io import savemat\n\n');
    fprintf(fid, 'def process(h5_path, out_mat_path, N, verbose):\n');
    fprintf(fid, '    verbose = (verbose.lower() == "true")\n');
    fprintf(fid, '    with h5py.File(h5_path, "r") as f:\n');
    fprintf(fid, '        dset = f["/1/EBSD/Data/Processed Patterns"]\n');
    fprintf(fid, '        N = min(N, dset.shape[0])\n');
    fprintf(fid, '        H = dset.shape[1]\n');
    fprintf(fid, '        W = dset.shape[2]\n\n');
    fprintf(fid, '        FSE = np.zeros((N, 1), dtype=np.float64)\n');
    fprintf(fid, '        BSE = np.zeros((N, 1), dtype=np.float64)\n');
    fprintf(fid, '        MSE = np.zeros((N, 1), dtype=np.float64)\n');
    fprintf(fid, '        GradMag = np.zeros((N, 1), dtype=np.float64)\n');
    fprintf(fid, '        LapVar = np.zeros((N, 1), dtype=np.float64)\n\n');
    fprintf(fid, '        y_top_end = H // 3\n');
    fprintf(fid, '        y_mid_end = (2 * H) // 3\n\n');
    fprintf(fid, '        chunk_size = 1000\n');
    fprintf(fid, '        for i in range(0, N, chunk_size):\n');
    fprintf(fid, '            t0 = time.time()\n');
    fprintf(fid, '            end = min(i + chunk_size, N)\n');
    fprintf(fid, '            chunk = dset[i:end, :, :].astype(np.float64)\n');
    fprintf(fid, '            t_io = time.time()\n\n');
    fprintf(fid, '            # Average over W (axis 2) and region of H (axis 1)\n');
    fprintf(fid, '            BSE[i:end, 0] = chunk[:, :y_top_end, :].mean(axis=(1, 2))\n');
    fprintf(fid, '            MSE[i:end, 0] = chunk[:, y_top_end:y_mid_end, :].mean(axis=(1, 2))\n');
    fprintf(fid, '            FSE[i:end, 0] = chunk[:, y_mid_end:, :].mean(axis=(1, 2))\n\n');
    fprintf(fid, '            # Mean Gradient Magnitude\n');
    fprintf(fid, '            grad_x = np.abs(np.diff(chunk, axis=2))\n');
    fprintf(fid, '            grad_y = np.abs(np.diff(chunk, axis=1))\n');
    fprintf(fid, '            GradMag[i:end, 0] = (grad_x.mean(axis=(1, 2)) + grad_y.mean(axis=(1, 2))) / 2.0\n\n');
    fprintf(fid, '            # Variance of Laplacian (Vectorized NumPy Slicing - 100x Faster)\n');
    fprintf(fid, '            lap = (chunk[:, :-2, 1:-1] + chunk[:, 2:, 1:-1] + chunk[:, 1:-1, :-2] + chunk[:, 1:-1, 2:] - 4.0 * chunk[:, 1:-1, 1:-1])\n');
    fprintf(fid, '            LapVar[i:end, 0] = lap.var(axis=(1, 2))\n');
    fprintf(fid, '            t_comp = time.time()\n\n');
    fprintf(fid, '            if verbose:\n');
    fprintf(fid, '                print(f"  Chunk {i//chunk_size + 1}/{int(np.ceil(N/chunk_size))}: I/O {t_io-t0:.2f}s | Compute {t_comp-t_io:.2f}s | Total {t_comp-t0:.2f}s", flush=True)\n\n');
    fprintf(fid, '    savemat(out_mat_path, {"FSE": FSE, "BSE": BSE, "MSE": MSE, "GradMag": GradMag, "LapVar": LapVar})\n\n');
    fprintf(fid, 'if __name__ == "__main__":\n');
    fprintf(fid, '    process(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])\n');
    fclose(fid);
    
    pyOutMat = fullfile(outDir, sprintf('%s_py_patterns.mat', angBaseName));
    
    verb_str = 'false';
    if verbose, verb_str = 'true'; end
    cmd = sprintf('python "%s" "%s" "%s" %d %s', pyScriptPath, h5Path, pyOutMat, N, verb_str);
    
    [status, cmdout] = system(cmd, '-echo');
    if status ~= 0
        error('Python extraction failed:\n%s', cmdout);
    end
    
    pyData = load(pyOutMat);
    FSE = pyData.FSE;
    BSE = pyData.BSE;
    MSE = pyData.MSE;
    GradMag = pyData.GradMag;
    LapVar = pyData.LapVar;
    
    % Cleanup temp python files
    if exist(pyOutMat, 'file'), delete(pyOutMat); end
    if exist(pyScriptPath, 'file'), delete(pyScriptPath); end
    % --------------------------------------------------------------
    
    % Save properties
    fprintf('Saving properties to %s...\n', outMat);
    save(outMat, 'BS', 'BC', 'FSE', 'BSE', 'MSE', 'GradMag', 'LapVar', '-v7.3');
end

fprintf('\n--- Extraction complete! ---\n');
