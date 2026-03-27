% MATLAB script to process folders of separate TIFF images
% 1. Scans specified folders for files matching a pattern (e.g., 'det3_frame*.tif').
% 2. Sorts them to ensure correct order.
% 3. Groups them into blocks, averages them, and saves the result.

clear; close all; clc;
%%
% --- CONFIGURATION ---
% Base directory containing the subfolders with images
% Ti
% baseDir = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\insitu\DIC\continual';
% targetFolders = {'Ti_c_Axis_PTP2_FOV6um_DIC2_frames', 'Ti_c_Axis_PTP2_FOV6um_DIC3_frames', 'Ti_c_Axis_PTP2_FOV6um_DIC4_frames', 'Ti_c_Axis_PTP2_FOV6um_DIC5_frames'}; % List of subfolders to process
% outputDirectory = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\insitu\DIC\averaged_continual\averaged5';

% Nb bicrystal
baseDir = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Nb_Bicrystal\Experiment\DIC\videos';
targetFolders = {'Nb insitu 2_frames'}; % List of subfolders to process
outputDirectory = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Nb_Bicrystal\Experiment\DIC\videos\Nb insitu 2_detector3Average_block10';



filePrefix = 'det3_frame'; % The start of the filename (e.g., det3_frame000005.tif)
blockSize = 10;      % Number of images per block
processingMethod = 'average'; % 'average' or 'median'

overwriteAll = false; % Flag to remember if user chose "overwrite all"

% --- PROCESSING ---

if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

for fIdx = 1:length(targetFolders)
    currentFolderName = targetFolders{fIdx};
    currentInputPath = fullfile(baseDir, currentFolderName);
    
    fprintf('Processing folder: %s\n', currentFolderName);
    
    filePattern = fullfile(currentInputPath, [filePrefix, '*.tif']);
    files = dir(filePattern);
    
    if isempty(files)
        warning('No files found in %s with prefix %s', currentInputPath, filePrefix);
        continue;
    end
    
    % Sort files numerically
    [~, reindex] = sort({files.name});
    files = files(reindex);
    
    numFrames = length(files);
    numBlocks = floor(numFrames / blockSize);
    
    fprintf('Found %d frames. Processing %d blocks...\n', numFrames, numBlocks);

    firstInfo = imfinfo(fullfile(currentInputPath, files(1).name));
    imgHeight = firstInfo.Height;
    imgWidth = firstInfo.Width;
    
    for blockIdx = 1:numBlocks
        % 1. Construct the Output Filename FIRST to check existence
        outName = sprintf('%s_%s_%s_%d_%02d.tif', ...
            currentFolderName, filePrefix, processingMethod, blockSize, blockIdx);
        outputFileName = fullfile(outputDirectory, outName);
        
        % --- SAFETY CHECK START ---
        if exist(outputFileName, 'file') && ~overwriteAll
            fprintf('\nWARNING: File already exists:\n%s\n', outputFileName);
            userChoice = input('Do you want to overwrite? (y = yes / n = skip / all = yes to all): ', 's');
            
            if strcmpi(userChoice, 'n')
                fprintf('Skipping file.\n');
                continue; % Skip to next block
            elseif strcmpi(userChoice, 'all')
                overwriteAll = true;
                fprintf('Overwriting ALL future conflicts without asking.\n');
            elseif ~strcmpi(userChoice, 'y')
                fprintf('Invalid input. Skipping file for safety.\n');
                continue;
            end
        end
        % --- SAFETY CHECK END ---

        % 2. Load and Process (Only runs if we passed the check)
        blockImages = zeros(imgHeight, imgWidth, blockSize);
        
        for k = 1:blockSize
            globalFileIdx = (blockIdx - 1) * blockSize + k;
            fullImagePath = fullfile(currentInputPath, files(globalFileIdx).name);
            blockImages(:,:,k) = double(imread(fullImagePath));
        end
        
        if strcmp(processingMethod, 'average')
            processedImage = mean(blockImages, 3);
        elseif strcmp(processingMethod, 'median')
            processedImage = median(blockImages, 3);
        end
        
        processedImage = uint16(processedImage);
        
        % 3. Save
        imwrite(processedImage, outputFileName, 'Compression', 'none');
        fprintf('Saved: %s\n', outName);
    end
end

disp('Processing complete.');