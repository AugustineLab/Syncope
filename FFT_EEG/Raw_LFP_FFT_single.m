%% run after loading single .mat file from biopac
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
ExportDir = 'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_EEG\EEG_syncope_Jonny\PVN_EEG_flip\Matfiles\FFTs\';
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
    figure
    plot(data1(idx_bsS(n):idx_asE(n),1));
    hold
    plot(data2(idx_bsS(n):idx_asE(n),1));
end 
%Get post-laser Resting EEG window
nwithinbin = 1;
WsecS = 40; % window in seconds before stim start EDIT
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
    figure
    plot(data1(idx_bsSp(n):idx_asEp(n),1));
    hold
    plot(data2(idx_bsSp(n):idx_asEp(n),1));
end

%FTT figure
for ch = 1:2
    if ch == 1
        currentData = data1;
    else
        currentData = data2;
    end
    figure
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

        semilogx(f,FFTmean, 'k','LineWidth',2)
        hold;
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
        semilogx(f,FFTmeanpost, 'r','LineWidth',2)
        hold;
    end
    %Graph look / properties
xlim([1 120])
xline(0.5,':k','LineWidth', 2,'LabelHorizontalAlignment','left')
xline(4,':k','LineWidth', 2, 'Label','Delta','LabelHorizontalAlignment','left')
xline(8,':k','LineWidth', 2, 'Label','Theta','LabelHorizontalAlignment','left')
xline(12,':k','LineWidth', 2,'Label','Alpha','LabelHorizontalAlignment','left')
xline(30,':k','LineWidth', 2,'Label','Beta','LabelHorizontalAlignment','left')
xline(60,':k','LineWidth', 2,'Label','LowGamma','LabelHorizontalAlignment','left')
xline(120,':k','LineWidth', 2,'Label','HighGamma','LabelHorizontalAlignment','left')
xlabel('Frequency(Hz)') 
ylabel('Power (mV^2)')
    if ch == 1
        title('Temporal')
    else
        title('Frontal')
    end
xticks([0.5 4 8 12 30 60 120])
hold off
legend('Pre','Post','Location','northwest')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
grid on

RatioPre = FFTmeanpost./FFTmean;
figure
semilogx(f,RatioPre,'b','LineWidth', 2)
hold on

xlim([1 120])
ylim([0 2])
yline(1,'--k','LineWidth',2)
xline(0.5,':k','LineWidth', 2,'LabelHorizontalAlignment','left')
xline(4,':k','LineWidth', 2, 'Label','Delta','LabelHorizontalAlignment','left')
xline(8,':k','LineWidth', 2, 'Label','Theta','LabelHorizontalAlignment','left')
xline(12,':k','LineWidth', 2,'Label','Alpha','LabelHorizontalAlignment','left')
xline(30,':k','LineWidth', 2,'Label','Beta','LabelHorizontalAlignment','left')
xline(60,':k','LineWidth', 2,'Label','LowGamma','LabelHorizontalAlignment','left')
xline(120,':k','LineWidth', 2,'Label','HighGamma','LabelHorizontalAlignment','left')
xlabel('Frequency(Hz)') 
ylabel('Post/Pre') 
    if ch == 1
        title('Temporal')
    else
        title('Frontal')
    end
xticks([0.5 4 8 12 30 60 120])
hold off
legend('Ratio','Location','northwest')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
grid on
end

distFig('Rows',2,'Columns',3)
