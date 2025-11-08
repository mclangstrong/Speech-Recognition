
% Step 1: Record test audio
fs = 16000; % Sampling rate
nBits = 16; % Bit depth
nChannels = 1; % Mono recording

disp('Recording test audio... Speak now!');
recObj = audiorecorder(fs, nBits, nChannels);
recordblocking(recObj, 2); % Record for 3 seconds
disp('Recording complete.');

% Extract and normalize audio
testAudio = getaudiodata(recObj);
testAudio = testAudio / max(abs(testAudio)); % Normalize amplitude

% Visualize test audio
figure;
plot(testAudio);
title('Test Audio Waveform');
xlabel('Sample Index');
ylabel('Amplitude');

% Step 2: Apply Noise Cancellation to Test Audio
% 2.1 Low-pass filtering (removes high-frequency noise)
fc = 3000; % Cutoff frequency
[b, a] = butter(6, fc / (fs / 2), 'low');
filteredTestAudio = filter(b, a, testAudio);

% 2.2 Spectral Subtraction for noise reduction
% Estimate noise power from the first 0.5 seconds (assuming silence at start)
noiseSample = filteredTestAudio(1:round(0.5 * fs)); % First 0.5 seconds
noisePower = mean(abs(fft(noiseSample)).^2);

% Perform spectral subtraction
N = length(filteredTestAudio);
freqDomain = fft(filteredTestAudio); % Transform to frequency domain
magnitude = abs(freqDomain); % Magnitude of spectrum
phase = angle(freqDomain); % Phase of spectrum
cleanedMagnitude = max(magnitude.^2 - noisePower, 0).^0.5; % Subtract noise
cleanedFreqDomain = cleanedMagnitude .* exp(1j * phase); % Recombine phase
cleanedTestAudio = real(ifft(cleanedFreqDomain)); % Transform back to time domain

% Normalize cleaned test audio
cleanedTestAudio = cleanedTestAudio / max(abs(cleanedTestAudio));

% Visualize the cleaned test audio
figure;
plot(cleanedTestAudio);
title('Cleaned Test Audio Waveform');
xlabel('Sample Index');
ylabel('Amplitude');

% Step 3: Load Reference Audios for Numbers and Words
numberWords = ["one", "two", "three", "four", "five", "six", "seven", ...
    "eight", "nine", "ten", "eleven"];
extraWords = [ "stop", "hi"];
referenceWords = [numberWords, extraWords];

% Folder containing reference audio files
folderPath = 'C:\Users\Marion Laptop\OneDrive\Desktop\Group 1 - Speech Recognition'; % Update this path as needed

% Initialize storage for reference audio and sample rates
referenceAudios = cell(length(referenceWords), 1);
referenceFs = zeros(length(referenceWords), 1);

% Load all reference audio files
disp('Loading reference audio files...');
for i = 1:length(referenceWords)
    filename = fullfile(folderPath, strcat(referenceWords(i), '.wav'));
    try
        disp(['Loading: ', filename]);
        [audio, refFs] = audioread(filename);
        audio = audio / max(abs(audio)); % Normalize amplitude
        referenceAudios{i} = audio;
        referenceFs(i) = refFs;
    catch ME
        disp(['Error loading file: ', filename]);
        disp(ME.message);
        referenceAudios{i} = [];
    end
end

% Resample all reference audios to match the test audio sampling rate
for i = 1:length(referenceWords)
    if ~isempty(referenceAudios{i})
        referenceAudios{i} = resample(referenceAudios{i}, fs, referenceFs(i));
    end
end

% Step 4: Cross-correlation with Reference Audios
disp('Performing cross-correlation with reference audios...');
maxCorrelations = zeros(length(referenceWords), 1); % Store max correlations
timeLags = zeros(length(referenceWords), 1);       % Store time lags

for i = 1:length(referenceWords)
    refAudio = referenceAudios{i};
    if ~isempty(refAudio)
        % Align lengths for correlation
        minLength = min(length(cleanedTestAudio), length(refAudio));
        testSegment = cleanedTestAudio(1:minLength);
        refSegment = refAudio(1:minLength);

        % Perform cross-correlation
        [corr, lags] = xcorr(testSegment, refSegment);
        [maxCorrelations(i), idx] = max(abs(corr));
        timeLags(i) = lags(idx) / fs; % Time lag in seconds
    else
        maxCorrelations(i) = 0; % No correlation if reference is missing
    end
end

% Step 5: Find the Best Match Based on the Highest Correlation
[bestCorr, bestIdx] = max(maxCorrelations);
if bestCorr > 0.7 % Set a suitable threshold
    recognizedWord = referenceWords(bestIdx);
    disp(['Recognized word/number: ', char(recognizedWord)]);
    disp(['Correlation: ', num2str(bestCorr), ', Time Lag: ', num2str(timeLags(bestIdx)), ' seconds']);
else
    disp('No match found.');
end

% Step 6: Display Correlation Results for Debugging
disp('Correlation results for all reference words/numbers:');
for i = 1:length(referenceWords)
    disp([char(referenceWords(i)), ': Correlation = ', num2str(maxCorrelations(i))]);
end

% Step 7: Plot Cross-Correlation of the Best Match
if bestCorr > 0.7
    figure;
    bestRefAudio = referenceAudios{bestIdx};
    [corr, lags] = xcorr(cleanedTestAudio, bestRefAudio);
    plot(lags / fs, corr);
    title(['Cross-Correlation with "', char(referenceWords(bestIdx)), '"']);
    xlabel('Time Lag (seconds)');
    ylabel('Correlation');
else
    disp('No valid match to plot.');
end
