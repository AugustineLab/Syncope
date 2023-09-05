%% run after opening .mat file from biopac
data1 = data(:,1);
data2 = data(:,2);
TTL = data(:,3)>1;
if strcmp(isi_units, 'ms')
    scale = 0.001;
else
    scale = 1;
end
LFPsr = 1/(isi*scale); %EDITED FOR YANG DAN DATA (1525.9)
TrainInfo = find_stim_train(TTL,LFPsr,4); % Get train timestamps (custom function)
%% !!!only run to get a "no light" idx based on 3min (180s) before first stim train!!!
nwithinbin = 1;
WsecS = -10; % window in seconds before stim start EDIT
WsecE = 10; % window in seconds after stim ends EDIT
idx_stimS = zeros(1,size(TrainInfo,2));
idx_stimE = zeros(1,size(TrainInfo,2));
idx_bsS = zeros(1,size(TrainInfo,2));
idx_asE = zeros(1,size(TrainInfo,2));
for n = 1:size(TrainInfo,2)
    idx_stimS(n) = (round(TrainInfo(3,n))/nwithinbin)-(180*LFPsr); % index of each stim start time
    idx_stimE(n) = (round(TrainInfo(4,n))/nwithinbin)-(180*LFPsr); % index of each stim end time
    idx_bsS(n) = ((idx_stimS(n)/LFPsr)+(WsecS))*LFPsr/nwithinbin;% index before each stim (window start)
    idx_asE(n) = ((idx_stimE(n)/LFPsr)+(WsecE))*LFPsr/nwithinbin;% index after each stim end (window end)
    figure
    plot(data1(idx_bsS(n):idx_asE(n),1));
    hold
    plot(data2(idx_bsS(n):idx_asE(n),1));
end 
%% get the window idx from each stimulation
nwithinbin = 1;
WsecS = -10; % window in seconds before stim start EDIT
WsecE = 10; % window in seconds after stim ends EDIT
idx_stimS = zeros(1,size(TrainInfo,2));
idx_stimE = zeros(1,size(TrainInfo,2));
idx_bsS = zeros(1,size(TrainInfo,2));
idx_asE = zeros(1,size(TrainInfo,2));
for n = 1:size(TrainInfo,2)
    idx_stimS(n) = round(TrainInfo(3,n))/nwithinbin; % index of each stim start time
    idx_stimE(n) = round(TrainInfo(4,n))/nwithinbin; % index of each stim end time
    idx_bsS(n) = ((idx_stimS(n)/LFPsr)+(WsecS))*LFPsr/nwithinbin;% index before each stim (window start)
    idx_asE(n) = ((idx_stimE(n)/LFPsr)+(WsecE))*LFPsr/nwithinbin;% index after each stim end (window end)
    figure
    plot(data1(idx_bsS(n):idx_asE(n),1));
    hold
    plot(data2(idx_bsS(n):idx_asE(n),1));
end 
%% plot LFP power
bsFilt = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'SampleRate',LFPsr); %remove 60Hz line noise
for ch = 1:2 
    if ch == 1
        currentData = data1;
    else
        currentData = data2;
    end
    for id = 1:size(TrainInfo,2) %EDITED FOR YANG DAN DATA (laser,1) DEFAULT (TrainInfo,2)
        figure
        dataseg = filtfilt(bsFilt,currentData(idx_bsS(id):idx_asE(id))); % get the data segment to calculate wavelets using above filter
        fb = cwtfilterbank('SignalLength',numel(dataseg),'SamplingFrequency',LFPsr,'FrequencyLimits',[0.5 120]); % only show frequncies between 0.5-120Hz
        [cfs, frq] = cwt(dataseg,'FilterBank', fb); % cfs = complex number matrix, frq = frequencies
        t = 0:length(dataseg)-1;
        tsec = t/LFPsr;
        power = abs(cfs);
        powerBL = [];
        for f = 1:size(frq,1)
            powerBL(f,:) = power(f,:)/mean(power(f,1:((WsecS*-1)*LFPsr)),2); %BASELINE PERIOD
        end
        pcolor(tsec,frq,powerBL);
        shading interp
        set(gca, 'YScale', 'log') %comment in for log scale visual
        xline((WsecS*-1),'k','LineWidth', 4)
        xline(tsec(end)-WsecE,'k','LineWidth', 4)
        yline(4,'w','LineWidth', 2) % Delta
        yline(8,'w','LineWidth', 2) % Theta
        yline(13,'w','LineWidth', 2) % Alpha
        yline(30,'w','LineWidth', 2) % Beta
        yline(60,'w','LineWidth', 2) % Low Gamma
        colormap(turbo)
        colorbar
        caxis([0 3])
        exportname = strcat('C:\Users\augus\Desktop\Kevin_Test_Analysis\csvs\','128R_20Hz_Veh2',num2str(ch),'.csv');
        writematrix(powerBL,exportname) 
    end
end
clear
clc

