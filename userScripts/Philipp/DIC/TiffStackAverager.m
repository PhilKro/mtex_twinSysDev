% MATLAB script to process a TIFF stack and export separate images
% Import a TIFF stack, divide it into blocks of 20 images, 
% compute the average or median for each block, and save each result as a separate image.

% Specify the input file name and output directory
file = {'01_DIC_Vid_detector3','02_DIC_Vid_detector3','03_DIC_Vid_detector3','04_DIC_Vid_detector3'};
inputFileDir = 'NN:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Titanium\Ti_cAxis_PTP_2\insitu\DIC\continual\Ti_c_Axis_PTP2_FOV6um_DIC1_frames';% Replace with your input file name
outputDirectory = 'N:\Vol1-Th\Analytik\Personal\Philipp Vol1 Analy Personal\Rhenium\20241218_Re_7dc_tkdPTP_3\04_insitu\testarea\DIC_average';  % Directory where separate images will be saved
blockSize = 10;  % Number of images per block


for i = 1:length(file)
    inputFileName = [inputFileDir, file{i}];

    % Choose the processing method ('average' or 'median')
    processingMethod = 'average';  % Change to 'median' if you want the median

    % Read the input TIFF stack information
    info = imfinfo([inputFileName,'.tif']);
    numFrames = length(info);
    blockSize = 20;  % Number of images per block

    % Calculate the number of blocks
    numBlocks = floor(numFrames / blockSize);

    % Create the output directory if it does not exist
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end

    % Initialize a variable to store the image size
    imgHeight = info(1).Height;
    imgWidth = info(1).Width;

    % Loop over blocks
    for blockIdx = 1:numBlocks
        % Initialize a 3D array to store the block of images
        blockImages = zeros(imgHeight, imgWidth, blockSize);

        % Load the images for the current block
        for k = 1:blockSize
            frameIdx = (blockIdx - 1) * blockSize + k;
            blockImages(:,:,k) = double(imread([inputFileName,'.tif'], frameIdx, 'Info', info));
        end

        % Compute the average or median of the block
        if strcmp(processingMethod, 'average')
            processedImage = mean(blockImages, 3);
        elseif strcmp(processingMethod, 'median')
            processedImage = median(blockImages, 3);
        else
            error('Invalid processing method. Choose "average" or "median".');
        end

        % Convert the processed image back to uint8 (or another format if needed)
        processedImage = uint16(processedImage);

        % Create the output file name for this image
        outputFileName = fullfile(outputDirectory, [file{i}, '_', processingMethod, sprintf('_%d_%02d.tif', blockSize, blockIdx)]);

        % Save the processed image as a separate TIFF file
        imwrite(processedImage, outputFileName, 'Compression', 'none');
    end
end

disp('Processing complete. The images have been saved separately.');
