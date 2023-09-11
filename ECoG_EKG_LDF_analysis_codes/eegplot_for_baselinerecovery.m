bsFilt = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'SampleRate',LFPsr); %remove 60Hz line noise

if g<3
    if j == 1
        id =3;
    else
        id = 1;
    end
    dataseg = filtfilt(bsFilt,EEG_sig(idx_bsS(id):idx_asE(id))); % get the data segment to calculate wavelets using above filter
    fb = cwtfilterbank('SignalLength',numel(dataseg),'SamplingFrequency',LFPsr,'FrequencyLimits',[0.5 120]); % only show frequncies between 0.5-120Hz
    [cfs, frq] = cwt(dataseg,'FilterBank', fb); % cfs = complex number matrix, frq = frequencies
    t = 0:length(dataseg)-1;
    tsec = t/LFPsr;
    power = abs(cfs);
    powerBL = [];
    for f = 1:size(frq,1)
        powerBL(f,:) = power(f,:)/mean(power(f,1:(10*LFPsr)),2); %EDIT TO CHANGE BASELINE PERIOD
    end
    band_EEG = powerBL(4:40,:);
    mean_bandpower = mean(band_EEG,1);
    all_eeg_band(j,g,:)= mean_bandpower(1,1:90*Fs);
else
    for id= 1:2
        dataseg = filtfilt(bsFilt,EEG_sig(idx_bsS(id):idx_asE(id))); % get the data segment to calculate wavelets using above filter
        fb = cwtfilterbank('SignalLength',numel(dataseg),'SamplingFrequency',LFPsr,'FrequencyLimits',[0.5 120]); % only show frequncies between 0.5-120Hz
        [cfs, frq] = cwt(dataseg,'FilterBank', fb); % cfs = complex number matrix, frq = frequencies
        t = 0:length(dataseg)-1;
        tsec = t/LFPsr;
        power = abs(cfs);
        powerBL = [];
        for f = 1:size(frq,1)
            powerBL(f,:) = power(f,:)/mean(power(f,1:(10*LFPsr)),2); %EDIT TO CHANGE BASELINE PERIOD
        end
        band_EEG = powerBL(4:40,:);
        mean_bandpower = mean(band_EEG,1);
        if j <3% check which animal for the swapping of frequency stim
            all_eeg_band(j,id+2,:)= mean_bandpower(1,1:90*Fs);
        elseif j>2 && id == 1
            all_eeg_band(j,4,:)= mean_bandpower(1,1:90*Fs);
        elseif j>2&& id == 2
            all_eeg_band(j,3,:)= mean_bandpower(1,1:90*Fs);    
        end
    end
end

