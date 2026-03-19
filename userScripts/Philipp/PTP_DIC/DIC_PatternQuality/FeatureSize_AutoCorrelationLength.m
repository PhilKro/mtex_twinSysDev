clear all
close all
% DIC Pattern Analysis: Feature Size and Correlation Length
% Replace 'dic_pattern.tif' with your actual file name
filename = 'userScripts/Philipp/PTP_DIC/DIC_PatternQuality/Re_PTP_3_speckle.tif';
hdr_fname = 'userScripts/Philipp/PTP_DIC/DIC_PatternQuality/Re_PTP_3_speckle-tif.hdr';

% Read the image and ensure it is grayscale and double precision
img = imread(filename);
if size(img, 3) == 3
    img = rgb2gray(img);
end
img = im2double(img);

%% 0. Image Cropping
% To skip manual cropping in future runs, replace the empty array below
% with the printed crop parameters from your command window.
% Example: cropRect = [150.5, 200.5, 500, 400];
cropRect = [1723.5, 151.5, 368.0, 972.0];

if isempty(cropRect)
    disp('Please draw a rectangle on the image to crop it. Double click inside the rectangle to confirm.');
    figure('Name', 'Interactive Crop');
    
    % imcrop opens an interactive tool and returns the cropped image and the bounding box
    [img_cropped, cropRect] = imcrop(img);
    close; % Close the cropping figure after selection

    % Print the parameters to the command window so you can copy them
    fprintf('\n--- Crop Parameters ---\n');
    fprintf('Copy and paste this line into the script to automate cropping next time:\n');
    fprintf('cropRect = [%.1f, %.1f, %.1f, %.1f];\n\n', cropRect(1), cropRect(2), cropRect(3), cropRect(4));
else
    % If cropRect is not empty, automatically crop using the provided coordinates
    img_cropped = imcrop(img, cropRect);
    fprintf('Automatically cropping using pre-defined parameters: [%.1f, %.1f, %.1f, %.1f]\n', cropRect);
end

% Reassign the cropped image to the main 'img' variable for the analysis steps
img = img_cropped;

%% 0.5 Read Pixel Size from HDR
% Read the hdr file to extract pixel size in meters, then convert to nm
hdr_text = fileread(hdr_fname);
tokens = regexp(hdr_text, 'PixelSizeX=([0-9\.eE\+\-]+)', 'tokens');
if ~isempty(tokens)
    pixelSize_m = str2double(tokens{1}{1});
    pixelSize_nm = pixelSize_m * 1e9; % Convert meters to nanometers
    fprintf('Extracted Pixel Size: %.4f nm/pixel\n', pixelSize_nm);
else
    warning('PixelSizeX not found in the HDR file. Defaulting to 1 nm/pixel.');
    pixelSize_nm = 1;
end

%% 1. Feature Size Estimation
% Use adaptive thresholding to isolate bright dots on a gray background
bw = imbinarize(img, 'adaptive', 'ForegroundPolarity', 'bright');

% Optional: Remove tiny noise artifacts (objects smaller than 4 pixels)
bw = bwareaopen(bw, 4); 

% Calculate the equivalent diameter of each isolated region
props = regionprops(bw, 'EquivDiameter', 'Centroid');
diameters = [props.EquivDiameter];

% Extract the mean feature size and standard deviation in pixels
meanFeatureSize = mean(diameters);
stdFeatureSize = std(diameters);

%% 1.5. Speckle Spacing Estimation
% Extract the (x, y) coordinates for the centroid of each speckle
centroids = cat(1, props.Centroid);

% Calculate the pairwise distances between all speckle centroids
distances = pdist2(centroids, centroids);

% Set the distance from a speckle to itself to infinity
distances(logical(eye(size(distances)))) = inf;

% Find the distance to the nearest neighbor for each speckle
nearestNeighborDists = min(distances, [], 2);

% Calculate the average spacing and standard deviation
meanSpacing = mean(nearestNeighborDists);
stdSpacing = std(nearestNeighborDists);

%% 2. Correlation Length Calculation
% Zero-mean the image to prevent a massive DC offset in the autocorrelation
img_zero_mean = img - mean(img(:));

% Calculate the 2D autocorrelation
autocorr2d = xcorr2(img_zero_mean);

% Normalize the autocorrelation so the peak is exactly 1
autocorr2d = autocorr2d / max(autocorr2d(:));

% Locate the center peak coordinates
[cy, cx] = find(autocorr2d == 1);
cy = cy(1); 
cx = cx(1);

% Extract a 1D slice through the center to fit the Gaussian
% A window of 30 pixels on either side is usually sufficient for fine patterns
window = 30; 
x_data = (-window:window)';
y_data = autocorr2d(cy, cx-window:cx+window)';

% Create a custom Gaussian model with a baseline offset 'd'
% f(x) = a * exp(-((x-b)/c)^2) + d
ft = fittype('a*exp(-((x-b)/c)^2) + d', 'independent', 'x', 'dependent', 'y');

% Set initial guesses to help the solver converge
% StartPoint array: [a (amplitude), b (center), c (width), d (baseline)]
baseline_guess = min(y_data);
amp_guess = max(y_data) - baseline_guess;
opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                  'StartPoint', [amp_guess, 0, 5, baseline_guess]);

% Fit the data using the custom model
fitResult = fit(x_data, y_data, ft, opts);
c_val = fitResult.c;

% Extract 1-sigma (68.27%) confidence intervals to estimate the error on the fit
ci = confint(fitResult, 0.6827);
c_idx = find(strcmp(coeffnames(fitResult), 'c')); % Safely locate parameter 'c'
c_error = (ci(2, c_idx) - ci(1, c_idx)) / 2;

% Calculate correlation length and its error (standard deviation)
correlationLength = c_val / sqrt(2);
stdCorrelationLength = c_error / sqrt(2);

%% 3. Display Results
fprintf('Mean Feature Size: %.2f ± %.2f pixels (%.2f ± %.2f nm)\n', ...
    meanFeatureSize, stdFeatureSize, meanFeatureSize * pixelSize_nm, stdFeatureSize * pixelSize_nm);
fprintf('Mean Speckle Spacing: %.2f ± %.2f pixels (%.2f ± %.2f nm)\n', ...
    meanSpacing, stdSpacing, meanSpacing * pixelSize_nm, stdSpacing * pixelSize_nm);
fprintf('Correlation Length: %.2f ± %.2f pixels (%.2f ± %.2f nm)\n', ...
    correlationLength, stdCorrelationLength, correlationLength * pixelSize_nm, stdCorrelationLength * pixelSize_nm);

% Plotting for visual verification
figure('Name', 'DIC Analysis Results');
subplot(1,3,1);
imshow(img); title('Cropped Image');

subplot(1,3,2);
imshow(bw); title('Binarized Features');

subplot(1,3,3);
plot(x_data, y_data, 'b.', 'MarkerSize', 12); hold on;
plot(fitResult, 'r');
title('Autocorrelation Peak & Fit');
legend('Data', 'Gaussian Fit');
grid on;