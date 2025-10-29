BPM = 150;
fs = 256; RMSSD=[];IBIs={};
% Sample MATLAB code to calculate heart rate from PPG signal
for i = 1:size(files,1)
    data = readtable(append(files(i).folder,'\',files(i).name));
    data = fillmissing(table2array(data(:,{'PPG_IR'})),'linear');
    rows_all_zero = find(all(data == 0, 2));
    data(rows_all_zero,:)=[];
    
    ppg=data;

    % Step 1: Preprocess PPG
    ppg = detrend(ppg);
    [b,a] = butter(2, 0.5/(fs/2), 'high');
    ppg_filt = filtfilt(b,a,ppg);
    [b,a] = butter(2, 3/(fs/2), 'low');% Keep cardiac band
    ppg_filt = filtfilt(b,a,ppg_filt);

    % Step 2: First derivative
    dPPG = diff(ppg_filt);

    % Step 3: Zero-crossings (negative→positive slope = onset candidate)
    zc_pos = find(dPPG(1:end-1) < 0 & dPPG(2:end) >= 0);

    % Step 4: Validate onsets (optional, e.g. min distance = 0.4s → 150 bpm max HR)
    minDist = round((60 / BPM) * fs);
    valid_onsets = [];
    for k = 2:length(zc_pos)
        if (zc_pos(k) - zc_pos(k-1)) > minDist
            valid_onsets(end+1) = zc_pos(k); %#ok<SAGROW>
        end
    end

    % Step 5: Inter-beat intervals (IBIs)
    IBI = diff(valid_onsets) / fs;   % in seconds
    HR = 60 ./ IBI;                  % Instantaneous heart rate (bpm)

    % Step 6: HRV (RMSSD example)
    RMSSD = cat(1,RMSSD, sqrt(mean(diff(IBI).^2)));

    IBIs=cat(1,IBIs,IBI);
end


% Preallocate a struct array
RMSSD_by_subject = struct('subject', num2cell(uniqueSubjects), 'RMSSD_values', [], 'RMSSD_mean', []);
for i = 1:length(uniqueSubjects)
subj = uniqueSubjects(i);
% Find indices of files corresponding to this subject
idx = find(subjectIDs == subj);
% Collect all RMSSD values for this subject
subj_RMSSD = [];
for j = 1:length(idx)
if iscell(RMSSD_all)
subj_RMSSD = [subj_RMSSD; RMSSD{idx(j)}(:)];
else
subj_RMSSD = [subj_RMSSD; RMSSD(idx(j))];
end
end
% Store in struct
RMSSD_by_subject(i).RMSSD_values = subj_RMSSD;
RMSSD_by_subject(i).RMSSD_mean = mean(subj_RMSSD);
end
% Display average RMSSD per subject
for i = 1:length(uniqueSubjects)
fprintf('Subject %d: average RMSSD = %.4f\n', RMSSD_by_subject(i).subject, RMSSD_by_subject(i).RMSSD_mean);
end
% Preallocate a struct array
RMSSD_by_subject = struct('subject', num2cell(uniqueSubjects), 'RMSSD_values', [], 'RMSSD_mean', []);
for i = 1:length(uniqueSubjects)
subj = uniqueSubjects(i);
% Find indices of files corresponding to this subject
idx = find(subjectIDs == subj);
% Collect all RMSSD values for this subject
subj_RMSSD = [];
for j = 1:length(idx)
if iscell(RMSSD)
subj_RMSSD = [subj_RMSSD; RMSSD{idx(j)}(:)];
else
subj_RMSSD = [subj_RMSSD; RMSSD(idx(j))];
end
end
% Store in struct
RMSSD_by_subject(i).RMSSD_values = subj_RMSSD;
RMSSD_by_subject(i).RMSSD_mean = mean(subj_RMSSD);
end
% Display average RMSSD per subject
for i = 1:length(uniqueSubjects)
fprintf('Subject %d: average RMSSD = %.4f\n', RMSSD_by_subject(i).subject, RMSSD_by_subject(i).RMSSD_mean);
end



