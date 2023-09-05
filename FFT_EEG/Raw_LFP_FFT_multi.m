%% loading a folder of .mat files from biopac generate list of names
testfiledir = uigetdir();
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);
ExportDir = 'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_EEG\EEG_syncope_Jonny\Jonny_EEG_flip_2_atro\MatFiles\TEST\';

%% Huge Loop of all files
for w = 1 : nfiles
    clear data
    clear isi
    clear isi_units
    clear labels
    clear start_sample
    clear units
    clear TrainInfo

    load(fullfile(testfiledir, matfiles(w).name))
    %data1 = data(:,1);
    data2 = data(:,2);
    TTL = data(:,3)>1;

%     if w == 1
%         TTL = data(:,3)>0.01;
%     end

    if strcmp(isi_units, 'ms')
        scale = 0.001;
    else
        scale = 1;
    end
    LFPsr = 1/(isi*scale); %EDITED FOR YANG DAN DATA (1525.9)
    TrainInfo = find_stim_train(TTL,LFPsr,4); % Get train timestamps (custom function)

%      if w == 5
%          TrainInfo(:,1) = TrainInfo(:,2);
%          TrainInfo(:,2) = [];
%      end
    %% Plot FFT power spectrum
    %Get pre-laser Resting EEG window
    nwithinbin = 1;
    WsecS = -240; % window in seconds before stim start EDIT
    WsecE = -30; % window in seconds after stim ends EDIT
    idx_stimS = zeros(1,size(TrainInfo,2));
    idx_stimE = zeros(1,size(TrainInfo,2));
    idx_bsS = zeros(1,size(TrainInfo,2));
    idx_asE = zeros(1,size(TrainInfo,2));
    for n = 1:size(TrainInfo,2)
        idx_stimS(n) = round(TrainInfo(3,n))/nwithinbin; % index of each stim start time
        idx_stimE(n) = round(TrainInfo(4,n))/nwithinbin; % index of each stim end time
        idx_bsS(n) = ((idx_stimS(n)/LFPsr)+(WsecS))*LFPsr/nwithinbin;% index before each stim (window start)
        idx_asE(n) = ((idx_stimE(n)/LFPsr)+(WsecE))*LFPsr/nwithinbin;% index after each stim end (window end)
    end 
    %Get post-laser Resting EEG window
    nwithinbin = 1;
    WsecS = 150; % window in seconds before stim start EDIT
    WsecE = 250; % window in seconds after stim ends EDIT
    idx_stimS = zeros(1,size(TrainInfo,2));
    idx_stimE = zeros(1,size(TrainInfo,2));
    idx_bsSp = zeros(1,size(TrainInfo,2));
    idx_asEp = zeros(1,size(TrainInfo,2));
    for n = 1:size(TrainInfo,2)
        idx_stimS(n) = round(TrainInfo(3,n))/nwithinbin; % index of each stim start time
        idx_stimE(n) = round(TrainInfo(4,n))/nwithinbin; % index of each stim end time
        idx_bsSp(n) = ((idx_stimS(n)/LFPsr)+(WsecS))*LFPsr/nwithinbin;% index before each stim (window start)
        idx_asEp(n) = ((idx_stimE(n)/LFPsr)+(WsecE))*LFPsr/nwithinbin;% index after each stim end (window end)
    end
    
    %FTT figure
    for ch = 2:2 %EDIT FOR DIFFERENT CHAN 
        if ch == 1
            currentData = data1;
        else
            currentData = data2;
        end

        for id = 1:size(TrainInfo,2)
            %calculate PRE-laser FFT
            bsFilt = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'SampleRate',LFPsr); %remove 60Hz line noise
            dataseg = filtfilt(bsFilt,currentData(idx_bsS(id):idx_asE(id))); % get the data segment to calculate FFT using above filter
            %dataseg = currentData(idx_bsS(id):idx_asE(id)); % get the data segment to calculate FFT no filter
            t1 = (length(dataseg)-1)/LFPsr;
            dur = 2; %segment duration in seconds.
            segmentsn = round(t1/dur);
            FFTs = zeros(segmentsn,LFPsr+1);
        
            T = 1/LFPsr; 
            L = LFPsr*dur;        % Length of signal
            t = (0:L-1)*T;        % Time vector
            f = LFPsr*(1:(L/2))/L;
        
            for k = 1:segmentsn
                s = (k*(dur*LFPsr))-(dur*LFPsr)+1;
                e = k*(dur*LFPsr);
                Y = fft(dataseg(s:e,1));
                    P2 = abs(Y/L);
                    P1 = P2(1:L/2+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    FFTs(k,:) = P1;
            end
            power = abs(FFTs);
            FFTmean = mean(power(:,2:2501),1);

            fname = extractBefore(matfiles(w).name, ".");
            drug = contains(fname,'Atro');
            rate = contains(fname,'20Hz');

            DataArray{w,1} = fname; %filename
            DataArray{w,2} = drug; %Drug
            DataArray{w,3} = rate; %Rate
            DataArray{w,4}{ch}{1} = FFTmean;
            DataArray{w,5} = LFPsr;
    
            %semilogx(f,FFTmean, 'k','LineWidth',2)
            %hold;
            %Calculate POST-laser FFT
            bsFilt = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'SampleRate',LFPsr); %remove 60Hz line noise
            dataseg = filtfilt(bsFilt,currentData(idx_bsSp(id):idx_asEp(id))); % get the data segment to calculate FFT using above filter
            %dataseg = currentData(idx_bsSp(id):idx_asEp(id)); % get the data segment to calculate FFT no filter
            t1 = (length(dataseg)-1)/LFPsr;
            dur = 2; %segment duration in seconds.
            segmentsn = round(t1/dur);
            FFTs = zeros(segmentsn,LFPsr+1);
        
            T = 1/LFPsr; 
            L = LFPsr*dur;        % Length of signal
            t = (0:L-1)*T;        % Time vector
            f = LFPsr*(1:(L/2))/L;
        
            for k = 1:segmentsn
                s = (k*(dur*LFPsr))-(dur*LFPsr)+1;
                e = k*(dur*LFPsr);
                Y = fft(dataseg(s:e,1));
                    P2 = abs(Y/L);
                    P1 = P2(1:L/2+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    FFTs(k,:) = P1;
            end
            powerpost = abs(FFTs);
            FFTmeanpost = mean(powerpost(:,2:2501),1);

            fname = extractBefore(matfiles(w).name, ".");
            drug = contains(fname,'Atro');
            rate = contains(fname,'20Hz');

            DataArray{w,1} = fname; %filename
            DataArray{w,2} = drug; %Drug
            DataArray{w,3} = rate; %rate
            DataArray{w,4}{ch}{2} = FFTmeanpost;

        end
    
    RatioPre = FFTmeanpost./FFTmean;

    fname = extractBefore(matfiles(w).name, ".");
    drug = contains(fname,'Atro');
    rate = contains(fname,'20Hz');

    DataArray{w,1} = fname; %filename
    DataArray{w,2} = drug; %Drug
    DataArray{w,3} = rate; %Drug
    DataArray{w,4}{ch}{3} = RatioPre;

    end
end
save(strcat(ExportDir,'Atro_test.mat'),"DataArray")