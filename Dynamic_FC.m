Rs=[];Ss=[];UYZs=[];
delta=[1,4];
theta=[4,8];
alpha=[8,13];
beta=[13,30];
gamma=[30,40];p_value=0.95;
freqz=[delta;theta;alpha;beta;gamma];fs=256;
combos_f=[nchoosek(1:length(freqz),2);[1:length(freqz);1:length(freqz)]'];
iter=20;
combos_s=[nchoosek(1:4,2)];
for i = 1:size(Data,1)
    eeg = Data{i};
    for band=1:length(combos_f)
        [b,a] = butter(4,(freqz(combos_f(band,1),2)/(fs/2)), 'low');
        eegA=filtfilt(b,a,eeg);
        [b,a] = butter(4,(freqz(combos_f(band,1),1)/(fs/2)), 'high');
        eegA=filtfilt(b,a,eegA);
        eegA=hilbert(eegA);

        [b,a] = butter(4,(freqz(combos_f(band,2),2)/(fs/2)), 'low');
        eegA=filtfilt(b,a,eeg);
        [b,a] = butter(4,(freqz(combos_f(band,2),1)/(fs/2)), 'high');
        eegB=filtfilt(b,a,eeg);
        eegB=hilbert(eegB);

        for PA = 1:4
            switch PA
                case 1
                    EA=interp1(1:length(eegA), abs(eegA), 1:50)';
                    EB=interp1(1:length(eegA),abs(eegB), 1:50)';
                case 2
                    EA=angle(eegA);
                    EB=angle(eegB);
                    EA = interp1(1:length(EA),unwrap(EA),1:52,'cubic');
                    EA=mod(EA + pi, 2*pi) - pi;
                    EB = interp1(1:length(EB),unwrap(EB),1:52,'cubic');
                    EB=mod(EB + pi, 2*pi) - pi;

                    [eigenvectorsA,eigenvalues] = lEIDA(EA);
                    EA=permute(eigenvectorsA(:,1,:),[1,3,2]);
                    [eigenvectorsB,eigenvalues] = lEIDA(EB);
                    EB=permute(eigenvectorsB(:,1,:),[1,3,2]);
                case 3
                    EA=interp1(1:length(eegA), abs(eegA), 1:50)';
                    EB=angle(eegB);
                    EB = interp1(1:length(EB),unwrap(EB),1:52,'cubic');
                    EB=mod(EB + pi, 2*pi) - pi;
                    
                    [eigenvectorsB,eigenvalues] = lEIDA(EB);
                    EB=permute(eigenvectorsB(:,1,:),[1,3,2]);
                case 4
                    EA=angle(eegA);
                    EA = interp1(1:length(EA),unwrap(EA),1:52,'cubic');
                    EA=mod(EA + pi, 2*pi) - pi;
                    [eigenvectorsA,eigenvalues] = lEIDA(EA);
                    EA=permute(eigenvectorsA(:,1,:),[1,3,2]);

                    EB=interp1(1:length(eegA),abs(eegB), 1:50)';
            end
            try
                HRV=interp1(1:length(IBIs{i}),IBIs{i},1:50,'linear')';



                rs=[];ss=[];uys=[];uzs=[];
                for ii=1:length(combos_s)
                    [r,s,uy,uz]=Gaussian_PID(EA(combos_s(ii,1),:)', EB(combos_s(ii,2),:)',HRV);
                    perms=[];
                    for p=1:iter
                        [rp,sp,uyp,uzp]=Gaussian_PID(EA(combos_s(ii,1),randperm(size(EA,2)))', EB(combos_s(ii,2),randperm(size(EB,2)))',HRV);
                        perms=[perms;[rp,sp,uyp,uzp]];
                    end
                    mu = mean(perms); sd = std(perms);
                    zscore = norminv(p_value); perm = mu + zscore*sd;
                    if r>perm(1)
                        rs = cat(2,rs,r);
                    else
                        rs = cat(2,rs,0);
                    end
                    if s>perm(2)
                        ss = cat(2,ss,s);
                    else
                        ss = cat(2,ss,0);
                    end
                    if uy>perm(3)
                        uys = cat(2,uys,uy);
                    else
                        uys = cat(2,uys,0);
                    end
                    if uz>perm(4)
                        uzs = cat(2,uzs,uz);
                    else
                        uzs = cat(2,uzs,0);
                    end
                end

                Rs=cat(1,Rs,rs);
                Ss=cat(1,Ss,ss);
                UYZs=cat(1,UYZs,cat(2,uys,uzs));
            catch message
                Rs=cat(1,Rs,[NaN;NaN;NaN;NaN;NaN;NaN]');
                Ss=cat(1,Ss,[NaN;NaN;NaN;NaN;NaN;NaN]');
                UYZs=cat(1,UYZs,cat(2,[NaN;NaN;NaN;NaN;NaN;NaN]',[NaN;NaN;NaN;NaN;NaN;NaN]'));
            end



        end

    end
    NC = [{Rs},{Ss},{UYZs}];
    save('D:\Neurolandscape_EEG\Nova Gorcia_2025\NC.mat','NC');
    fprintf('Participant %d done \n', i)
end