import os
import glob
import time
import h5py
import hdf5plugin
import psutil
import numpy as np
from scipy.io import savemat
import multiprocessing as mp

def process_chunk(args):
    h5_path, start_idx, end_idx = args
    
    # Each worker opens the file independently in read-only mode to prevent HDF5 locks
    with h5py.File(h5_path, 'r') as f:
        dset = f['/1/EBSD/Data/Processed Patterns']
        H = dset.shape[1]
        
        # Read chunk into RAM
        chunk = dset[start_idx:end_idx, :, :].astype(np.float64)
        
    y_top_end = H // 3
    y_mid_end = (2 * H) // 3
    
    # 1. Intensity Means
    BSE = chunk[:, :y_top_end, :].mean(axis=(1, 2)).reshape(-1, 1)
    MSE = chunk[:, y_top_end:y_mid_end, :].mean(axis=(1, 2)).reshape(-1, 1)
    FSE = chunk[:, y_mid_end:, :].mean(axis=(1, 2)).reshape(-1, 1)
    
    # 2. Gradient Magnitude
    grad_x = np.abs(np.diff(chunk, axis=2))
    grad_y = np.abs(np.diff(chunk, axis=1))
    GradMag = (grad_x.mean(axis=(1, 2)) + grad_y.mean(axis=(1, 2))) / 2.0
    GradMag = GradMag.reshape(-1, 1)
    
    # 3. Laplacian Variance
    lap = (
        chunk[:, :-2, 1:-1] +
        chunk[:, 2:, 1:-1] +
        chunk[:, 1:-1, :-2] +
        chunk[:, 1:-1, 2:] -
        4.0 * chunk[:, 1:-1, 1:-1]
    )
    LapVar = lap.var(axis=(1, 2)).reshape(-1, 1)
    
    return start_idx, end_idx, FSE, BSE, MSE, GradMag, LapVar

def process_file(h5_path, num_cores):
    base_name = os.path.splitext(os.path.basename(h5_path))[0]
    out_mat_path = base_name + '_h5oina_props.mat'
    
    if os.path.exists(out_mat_path):
        print(f"Skipping {os.path.basename(h5_path)}, output {out_mat_path} already exists.")
        return
        
    print(f"\nProcessing {os.path.basename(h5_path)} using {num_cores} CPU cores...")
    t0 = time.time()
    
    with h5py.File(h5_path, 'r') as f:
        dset = f['/1/EBSD/Data/Processed Patterns']
        N = dset.shape[0]
        
        try:
            BS = f['/1/EBSD/Data/Band Slope'][:]
            BC = f['/1/EBSD/Data/Band Contrast'][:]
        except:
            BS = np.array([])
            BC = np.array([])
            
    # --- Dynamic Chunk Size Calculation ---
    # Memory footprint per pattern during Laplacian math (float64 + 5 shifted arrays) ~= 15 MB
    mem_per_pattern_mb = 15.0
    mem_per_pattern_bytes = mem_per_pattern_mb * 1024 * 1024
    
    # We want to use at most 80% of available RAM across all cores
    avail_ram = psutil.virtual_memory().available
    target_ram = avail_ram * 0.8
    
    # Max patterns we can process at once across all cores
    max_patterns_total = int(target_ram / mem_per_pattern_bytes)
    
    # Chunk size per core
    chunk_size = max(100, max_patterns_total // num_cores)
    
    # Cap chunk size to dataset size (removes chunking if we have enough RAM)
    chunk_size = min(chunk_size, N)
    
    print(f"Available RAM: {avail_ram / (1024**3):.1f} GB")
    print(f"Calculated optimal chunk size: {chunk_size} patterns per core (Max Memory Footprint: ~{chunk_size * num_cores * mem_per_pattern_mb / 1024:.1f} GB)")
            
    tasks = []
    for i in range(0, N, chunk_size):
        end = min(i + chunk_size, N)
        tasks.append((h5_path, i, end))
        
    # Preallocate output arrays
    FSE_full = np.zeros((N, 1), dtype=np.float64)
    BSE_full = np.zeros((N, 1), dtype=np.float64)
    MSE_full = np.zeros((N, 1), dtype=np.float64)
    GradMag_full = np.zeros((N, 1), dtype=np.float64)
    LapVar_full = np.zeros((N, 1), dtype=np.float64)
    
    # Spin up 18 parallel workers
    with mp.Pool(processes=num_cores) as pool:
        for i, result in enumerate(pool.imap_unordered(process_chunk, tasks)):
            start_idx, end_idx, FSE, BSE, MSE, GradMag, LapVar = result
            FSE_full[start_idx:end_idx] = FSE
            BSE_full[start_idx:end_idx] = BSE
            MSE_full[start_idx:end_idx] = MSE
            GradMag_full[start_idx:end_idx] = GradMag
            LapVar_full[start_idx:end_idx] = LapVar
            print(f"  [{i+1}/{len(tasks)}] Chunks processed...", flush=True)
            
    # Save to MAT file
    print(f"Saving extracted properties to {out_mat_path}...")
    savemat(out_mat_path, {
        'FSE': FSE_full,
        'BSE': BSE_full,
        'MSE': MSE_full,
        'GradMag': GradMag_full,
        'LapVar': LapVar_full,
        'BS': BS.reshape(-1, 1) if BS.size > 0 else BS,
        'BC': BC.reshape(-1, 1) if BC.size > 0 else BC
    })
    
    print(f"Finished {os.path.basename(h5_path)} in {time.time() - t0:.2f} seconds.")

if __name__ == '__main__':
    # The script processes all .h5oina files in the folder it is run from.
    target_dir = '.' 
    num_cores = 18 
    
    print(f"Looking for .h5oina files in {os.path.abspath(target_dir)}...")
    h5_files = glob.glob(os.path.join(target_dir, '**', '*.h5oina'), recursive=True)
    
    if not h5_files:
        print("No .h5oina files found!")
    else:
        for f in h5_files:
            process_file(f, num_cores)
        print("\nAll files processed successfully!")
