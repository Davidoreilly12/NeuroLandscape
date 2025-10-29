
addpath(genpath('/users/libdor/NeuroLandscape/'));
srate = 256;

% Get SLURM task index
idx = str2double(getenv('SLURM_ARRAY_TASK_ID'));

% Get full list of .csv files (you can refine the path)
basePath = 'D:\Neurolandscape_EEG';  % Update as needed
csvFiles = dir(fullfile(basePath, '**', '*.csv'));

if idx > length(csvFiles)
    fprintf('Index %d exceeds number of files (%d). Skipping.\n', idx, length(csvFiles));
    return;
end

% Load the specific file
csvFilePath = fullfile(csvFiles(idx).folder, csvFiles(idx).name);
data = readtable(csvFilePath, 'ReadVariableNames', true,"NumHeaderLines",0);
data=data(2:end,:);

eeglab('nogui'); % Suppress EEGLAB GUI

% Import EEG data
eeg = [data.RAW_TP9, data.RAW_AF7, data.RAW_AF8, data.RAW_TP10]';
EEG = pop_importdata('dataformat','array','data',eeg,'srate',srate);
EEG.chanlocs = struct( ...
    'labels', {'TP9', 'AF7', 'AF8', 'TP10'}, ...
    'X',      {-72,  -53,   53,   72}, ...
    'Y',      {-17,   43,   43,  -17}, ...
    'Z',      {-28,   36,   36,  -28}, ...
    'type',   {'EEG','EEG','EEG','EEG'} ...
);
EEG = eeg_checkset(EEG);

% Fill missing values
EEG.data = fillmissing(EEG.data, 'linear', 2);

EEG = pop_reref(EEG, [], 'exclude', []);
EEG.ref = 'averef';  % label it as average-referenced
EEG.data=double(EEG.data);
% Bandpass filter: 0.5 - 40 Hz
[b,a] = butter(4, 0.5/(srate/2), 'high');
e = filtfilt(b,a, EEG.data')';
[b,a] = butter(4, 40/(srate/2), 'low');
EEG.data = filtfilt(b,a, e')';

% Run ICA (non-interactive)
EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'off');

% ICLabel and component removal
EEG = pop_iclabel(EEG, 'default');
threshold = 0.9;
labels = EEG.etc.ic_classification.ICLabel.classifications;
toRemove = find(labels(:,2) > threshold | labels(:,3) > threshold | ...
                labels(:,4) > threshold | labels(:,6) > threshold);
EEG = pop_subcomp(EEG, toRemove, 0);
fprintf('Rejected ICA components (artifacts): %s\n', mat2str(toRemove));

chanLabels = {'TP9', 'AF7', 'AF8', 'TP10'};

% Transpose and convert to table with headers
processedData = array2table(EEG.data', 'VariableNames', chanLabels);

% Convert to Excel-compatible path
[folder, fileName, ~] = fileparts(csvFilePath);
excelFilePath = fullfile(folder, [fileName, '.xlsx']);

% Save processed data to a new sheet
writetable(processedData, excelFilePath, 'Sheet', 'ProcessedEEG', 'WriteMode', 'overwrite');
