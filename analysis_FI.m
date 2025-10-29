try
    % Add required paths
    addpath(genpath('/users/libdor/Analysis'));
    addpath(genpath('/users/libdor/packages/'));
    addpath(genpath('/users/libdor/NICF_Final/'));

    % Load input data
    load('/users/libdor/Analysis/IBIs_FI.mat');   % variable: IBIs
    load('/users/libdor/Analysis/Data_FI.mat');   % variable: Data
    Data = permute(Data, [3,1,2]);
    % Initialize output containers
    Rs = [];
    Ss = [];
    UYZs = [];

    % Frequency bands
    delta = [1,4];
    theta = [4,8];
    alpha = [8,13];
    beta  = [13,30];
    gamma = [30,40];
    p_value = 0.95;

    freqz = [delta;theta;alpha;beta;gamma];
    fs = 256;

    combos_f = [nchoosek(1:length(freqz),2); [1:length(freqz);1:length(freqz)]'];
    combos_s = nchoosek(1:4,2);

    iter = 20;

    % Loop over participants
    for i = 1:size(Data,1)
        eeg = Data{i};

        % Loop over frequency band combinations
        for band = 1:size(combos_f,1)
            % Band A
            [b,a] = butter(4, (freqz(combos_f(band,1),2) / (fs/2)), 'low');
            eegA = filtfilt(b,a,eeg);
            [b,a] = butter(4, (freqz(combos_f(band,1),1) / (fs/2)), 'high');
            eegA = filtfilt(b,a,eegA);
            eegA = hilbert(eegA);

            % Band B
            [b,a] = butter(4, (freqz(combos_f(band,2),2) / (fs/2)), 'low');
            eegB = filtfilt(b,a,eeg);
            [b,a] = butter(4, (freqz(combos_f(band,2),1) / (fs/2)), 'high');
            eegB = filtfilt(b,a,eegB);
            eegB = hilbert(eegB);

            % Loop over phase/amplitude (PA) cases
            for PA = 1:4
                switch PA
                    case 1 % amplitude-amplitude
                        EA = interp1(1:length(eegA), abs(eegA), 1:50)';
                        EB = interp1(1:length(eegA), abs(eegB), 1:50)';

                    case 2 % phase-phase (with LEIDA)
                        EA = angle(eegA);
                        EB = angle(eegB);
                        EA = interp1(1:length(EA), unwrap(EA), 1:52, 'cubic');
                        EA = mod(EA + pi, 2*pi) - pi;
                        EB = interp1(1:length(EB), unwrap(EB), 1:52, 'cubic');
                        EB = mod(EB + pi, 2*pi) - pi;

                        [eigenvectorsA,~] = lEIDA(EA);
                        EA = permute(eigenvectorsA(:,1,:), [1,3,2]);
                        [eigenvectorsB,~] = lEIDA(EB);
                        EB = permute(eigenvectorsB(:,1,:), [1,3,2]);

                    case 3 % amplitude-phase
                        EA = interp1(1:length(eegA), abs(eegA), 1:50)';
                        EB = angle(eegB);
                        EB = interp1(1:length(EB), unwrap(EB), 1:52, 'cubic');
                        EB = mod(EB + pi, 2*pi) - pi;
                        [eigenvectorsB,~] = lEIDA(EB);
                        EB = permute(eigenvectorsB(:,1,:), [1,3,2]);

                    case 4 % phase-amplitude
                        EA = angle(eegA);
                        EA = interp1(1:length(EA), unwrap(EA), 1:52, 'cubic');
                        EA = mod(EA + pi, 2*pi) - pi;
                        [eigenvectorsA,~] = lEIDA(EA);
                        EA = permute(eigenvectorsA(:,1,:), [1,3,2]);
                        EB = interp1(1:length(eegA), abs(eegB), 1:50)';

                end

                try
                    % HRV resampling
                    HRV = interp1(1:length(IBIs{i}), IBIs{i}, 1:50, 'linear')';

                    rs=[]; ss=[]; uys=[]; uzs=[];

                    for ii = 1:length(combos_s)
                        [r,s,uy,uz] = Gaussian_PID(EA(combos_s(ii,1),:)', ...
                                                   EB(combos_s(ii,2),:)', ...
                                                   HRV);

                        % Permutation test
                        perms = [];
                        for p = 1:iter
                            [rp,sp,uyp,uzp] = Gaussian_PID( ...
                                EA(combos_s(ii,1), randperm(size(EA,2)))', ...
                                EB(combos_s(ii,2), randperm(size(EB,2)))', ...
                                HRV);
                            perms = [perms; [rp,sp,uyp,uzp]];
                        end

                        mu = mean(perms);
                        sd = std(perms);
                        zscore = norminv(p_value);
                        perm = mu + zscore*sd;

                        if r>perm(1), rs = [rs, r]; else, rs = [rs,0]; end
                        if s>perm(2), ss = [ss, s]; else, ss = [ss,0]; end
                        if uy>perm(3), uys = [uys, uy]; else, uys = [uys,0]; end
                        if uz>perm(4), uzs = [uzs, uz]; else, uzs = [uzs,0]; end
                    end

                    Rs   = [Rs; rs];
                    Ss   = [Ss; ss];
                    UYZs = [UYZs; [uys, uzs]];

                catch innerME
                    warning('Participant %d band %d PA %d failed: %s', i, band, PA, innerME.message);
                    Rs   = [Rs; nan(1,6)];
                    Ss   = [Ss; nan(1,6)];
                    UYZs = [UYZs; nan(1,12)];
                end
            end
        end

        fprintf('Participant %d done \n', i);
    end

    % Save results
    NC = {Rs, Ss, UYZs};
    save('NC_FI.mat','NC');

catch ME
    disp('Fatal error:');
    disp(getReport(ME));
end

exit;
