% Define test audio and bird images folders
testFolder = "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\TESTITTY";
imageFolder = "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Bird_Images";
placeholderImage = "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Bird_Images\placeholder.jpg";

% List of folders containing reference audio files
referenceFolders = {
    "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Corvus corax";
    "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Dryocopus pileatus";
    "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Tyto alba";
    "C:\Users\moham\Desktop\DSP PROJECT\FINAL EDITED\Zenaida macroura"
};

% Create UI Figure
fig = uifigure('Name', 'Bird Sound Recognition', 'Position', [100 100 1100 700]);

% Title Label
uilabel(fig, 'Text', 'Bird Sound Recognition System', ...
    'Position', [400 650 300 40], 'FontSize', 14, 'FontWeight', 'bold');

% Dropdown for selecting test audio (smaller size)
lblSelect = uilabel(fig, 'Text', 'Select Test Audio:', ...
    'Position', [50 600 120 30], 'FontSize', 12);
audioFiles = list_audio_files(testFolder);
if isempty(audioFiles)
    audioFiles = {'No files found'};
end
dropdownFiles = uidropdown(fig, ...
    'Position', [180 605 150 25], ... % Reduced width to 150, height to 25
    'Items', audioFiles, ...
    'Value', audioFiles{1}, ...
    'ValueChangedFcn', @(dd, event) updateTestFile(fig, dd, testFolder));

% Checkbox for FFT display
chkFFT = uicheckbox(fig, 'Text', 'Show FFT Spectrum', ...
    'Position', [50 560 200 30]);

% Button to run analysis
btnAnalyze = uibutton(fig, 'push', 'Text', 'Analyze Selected File', ...
    'Position', [50 510 200 40], ...
    'ButtonPushedFcn', @(btn, event) processAudio(fig, imageFolder, placeholderImage, referenceFolders));

% Label for Best Match
lblResult = uilabel(fig, 'Text', 'Best Match: Unknown', ...
    'Position', [50 460 300 40], 'FontSize', 12, 'FontWeight', 'bold');

% Play Sound Button
btnPlay = uibutton(fig, 'push', 'Text', 'Play Test Audio', ...
    'Position', [50 410 200 40], ...
    'ButtonPushedFcn', @(btn, event) playAudio(fig));

% Axes for Waveform (left side, shifted down)
axWaveform = uiaxes(fig, 'Position', [300 350 400 250]);
title(axWaveform, 'Waveform');
xlabel(axWaveform, 'Time (s)');
ylabel(axWaveform, 'Amplitude');

% Axes for FFT (right side, side by side with waveform, shifted down)
axFFT = uiaxes(fig, 'Position', [710 350 400 250]);
title(axFFT, 'Frequency Spectrum');
xlabel(axFFT, 'Frequency (Hz)');
ylabel(axFFT, 'Magnitude');

% Axes for Spectrogram (below waveform and FFT, shifted down)
axSpectrogram = uiaxes(fig, 'Position', [300 100 400 250]);
title(axSpectrogram, 'Spectrogram');
xlabel(axSpectrogram, 'Time (s)');
ylabel(axSpectrogram, 'Frequency (Hz)');

% Bird Image Display (moved to right of spectrogram, shifted down)
axImage = uiaxes(fig, 'Position', [710 100 250 250]);
if exist(placeholderImage, 'file')
    imshow(placeholderImage, 'Parent', axImage);
else
    warning('Placeholder image not found!');
end

% Progress Display
txtProgress = uitextarea(fig, 'Position', [50 20 1000 25], 'Editable', 'off');

% Store UI Components in Figure Data
fig.UserData = struct('dropdownFiles', dropdownFiles, ...
                      'chkFFT', chkFFT, ...
                      'lblResult', lblResult, ...
                      'axWaveform', axWaveform, ...
                      'axFFT', axFFT, ...
                      'axSpectrogram', axSpectrogram, ...
                      'axImage', axImage, ...
                      'txtProgress', txtProgress, ...
                      'testFile', fullfile(testFolder, audioFiles{1}));

%% Function Definitions

% Function to list all test audio files
function fileList = list_audio_files(folderPath)
    files = dir(fullfile(folderPath, '*.wav')); % Get all .wav files
    fileList = {files.name}; % Extract file names
    if isempty(fileList)
        fileList = {'No files found'}; % Prevent empty dropdown error
    end
end

% Update test file when user selects from dropdown
function updateTestFile(fig, dd, testFolder)
    if strcmp(dd.Value, 'No files found')
        fig.UserData.testFile = ''; % No valid test file selected
    else
        fig.UserData.testFile = fullfile(testFolder, dd.Value);
    end
end

% Process the selected audio file
function processAudio(fig, imageFolder, placeholderImage, referenceFolders)
    % Ensure a file is selected
    if isempty(fig.UserData.testFile) || strcmp(fig.UserData.testFile, 'No files found')
        uialert(fig, 'No test file selected. Please choose a valid file.', 'Error');
        return;
    end
    testFile = fig.UserData.testFile;

    % Update progress display
    fig.UserData.txtProgress.Value = sprintf('Analyzing: %s...', testFile);

    % Run correlation analysis
    [bestMatch, highestCorr] = run_correlation(testFile, fig, referenceFolders);

    % Update UI Label
    fig.UserData.lblResult.Text = sprintf('Best Match: %s (%.4f)', bestMatch, highestCorr);

    % Load and display bird image
    matchedImage = fullfile(imageFolder, sprintf('%s.jpg', bestMatch));
    if exist(matchedImage, 'file')
        imshow(matchedImage, 'Parent', fig.UserData.axImage);
    else
        imshow(placeholderImage, 'Parent', fig.UserData.axImage);
    end
end

% Play selected test audio file
function playAudio(fig)
    if isfield(fig.UserData, 'testFile') && ~isempty(fig.UserData.testFile)
        [audioData, fs] = audioread(fig.UserData.testFile);
        sound(audioData, fs);
    end
end

% Run correlation for the selected file
function [bestMatch, highestCorr] = run_correlation(testFile, fig, referenceFolders)
    % Load Test Audio
    [test_signal, fs] = audioread(testFile);
    if size(test_signal, 2) > 1
        test_signal = mean(test_signal, 2);
    end
    test_signal = test_signal / max(abs(test_signal));

    % Time vector for plotting
    t = (0:length(test_signal)-1) / fs;

    % FFT Calculation
    N = length(test_signal);
    f_axis = (0:N-1) * (fs/N);
    fft_signal = abs(fft(test_signal));

    % STFT Calculation and Spectrogram
    window = hamming(256); % Window size (adjust as needed)
    noverlap = 128; % Overlap (adjust for smoothness)
    [s, f, t_stft] = stft(test_signal, fs, 'Window', window, 'OverlapLength', noverlap, 'Centered', false);

    % Update Waveform Plot
    plot(fig.UserData.axWaveform, t, test_signal);
    xlabel(fig.UserData.axWaveform, 'Time (s)');
    ylabel(fig.UserData.axWaveform, 'Amplitude');

    % Show FFT if checkbox is checked
    if fig.UserData.chkFFT.Value
        plot(fig.UserData.axFFT, f_axis(1:floor(N/2)), fft_signal(1:floor(N/2)));
        xlabel(fig.UserData.axFFT, 'Frequency (Hz)');
        ylabel(fig.UserData.axFFT, 'Magnitude');
    else
        cla(fig.UserData.axFFT); % Clear FFT plot
    end

    % Show Spectrogram
    imagesc(fig.UserData.axSpectrogram, t_stft, f, 20*log10(abs(s)));
    axis(fig.UserData.axSpectrogram, 'xy');
    colormap(fig.UserData.axSpectrogram, jet);
    colorbar(fig.UserData.axSpectrogram);
    xlabel(fig.UserData.axSpectrogram, 'Time (s)');
    ylabel(fig.UserData.axSpectrogram, 'Frequency (Hz)');
    title(fig.UserData.axSpectrogram, 'Spectrogram');

    % Correlation Analysis
    highestCorr = -Inf;
    bestMatch = '';

    % Loop through each reference folder
    for i = 1:length(referenceFolders)
        folder = referenceFolders{i};
        refFiles = dir(fullfile(folder, '*.wav'));
        speciesName = extractAfter(folder, "FINAL EDITED\");

        % Loop through each reference file
        for j = 1:length(refFiles)
            refPath = fullfile(folder, refFiles(j).name);
            [ref_signal, fs_ref] = audioread(refPath);

            % Ensure same sampling rate
            if fs_ref ~= fs
                ref_signal = resample(ref_signal, fs, fs_ref);
            end

            % Convert to mono and normalize
            if size(ref_signal, 2) > 1
                ref_signal = mean(ref_signal, 2);
            end
            ref_signal = ref_signal / max(abs(ref_signal));

            % Pad signals to the same length
            max_length = max(length(test_signal), length(ref_signal));
            test_signal_padded = pad_signal(test_signal, max_length);
            ref_signal_padded = pad_signal(ref_signal, max_length);

            % Compute cross-correlation
            [xcorr_val, ~] = xcorr(test_signal_padded, ref_signal_padded, 'coeff');
            maxCorr = max(xcorr_val);

            % Update best match
            if maxCorr > highestCorr
                highestCorr = maxCorr;
                bestMatch = speciesName; % Use species name instead of file name
            end
        end
    end

    % Update progress
    fig.UserData.txtProgress.Value = sprintf('Analysis complete. Best Match: %s (%.4f)', bestMatch, highestCorr);
end

% Function to pad a signal to a given length
function padded_signal = pad_signal(signal, target_length)
    padded_signal = zeros(target_length, 1);
    padded_signal(1:length(signal)) = signal;
end