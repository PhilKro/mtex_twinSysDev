function average_frames_from_avi(videoPath, framesToAverage)

    if nargin < 2
        disp('Reverted to default averaging of 5 frames.')
        framesToAverage = 5;
    end

    % Extract folder and video name
    [parentFolder, baseName, ~] = fileparts(videoPath);

    % Create output folder next to the video
    outputFolder = fullfile(parentFolder, baseName + "Average");

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Load video
    vid = VideoReader(videoPath);

    buffer = [];
    avgIndex = 0;

    while hasFrame(vid)
        frame = readFrame(vid);
        frame = im2double(frame);

        buffer(:, :, :, end+1) = frame;  %#ok<AGROW>

        % Average when enough frames collected
        if size(buffer, 4) == framesToAverage
            avgFrame = mean(buffer, 4);

            % Save PNG
            outName = fullfile(outputFolder, sprintf('avg_%04d.png', avgIndex));
            imwrite(avgFrame, outName);

            avgIndex = avgIndex + 1;
            buffer = [];  % Reset
        end
    end
end
