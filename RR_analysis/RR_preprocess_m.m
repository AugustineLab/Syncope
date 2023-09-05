
function [Rrms_final,stim_name,tss,Amps] = RR_preprocess_m (test,path,sti_number,plotting,nstart,nend,smoothing)
Rrms_20HZ= [];
Amps = [];
for n = nstart:nend
[cage, mouse, ~, group, fr]  = readcagemouse(n); % Read the mouse information (edited in readmouse.m)
filename = [ mouse '_' cage '.mat']; % Get the name of mat file downloaded from Acknowledge
Data = load([path '\' test '\' group '\' filename]); % Load dataData.isi

% Import TTL and RR
for m = 4:size(Data.labels,1) % In Acknowledge the smoothed (3000 sample mean smoothing) Breathing raw data should be named as RR
    label = Data.labels (m,:);
    label_logic = (label (1:2)== 'RR' ); % find the channel number of RR
    if  label_logic(2) ==1 
        RR = m;
    end 
end 
TTL = Data.data (:,2); % Load TTL named as TTL
Respiration_smooth = Data.data (:,RR); % Load RR named as Respiration_smooth
plot(Respiration_smooth);
% Make sure all fr is 0.5
if fr == 0.4 % some data is recorded in different SR
originalFs = 2500;
desiredFs = 2000;
[p,q] = rat(desiredFs / originalFs);
TTL = resample(TTL,p,q);
Respiration_smooth = resample(Respiration_smooth,p,q);
end
% Make sure all fr is 0.5
if fr == 0.25 % some data is recorded in different SR
TTL = downsample(TTL,2);
Respiration_smooth = downsample(Respiration_smooth,2);
end
% Make sure all fr is 0.5
if fr == 0.1 % some data is recorded in different SR
TTL = downsample(TTL,5);
Respiration_smooth = downsample(Respiration_smooth,5);
end 
% do the smoothing to avoid additional peaks
if smoothing ==1 
Respiration_smooth = movmean(Respiration_smooth,2000);
disp(Respiration_smooth)
end 
plot(Respiration_smooth);

Trains = find_stim_train(TTL, 2000, 4); % use the function 'find_stim_train2.mat' to find the onset of TTL
% Find the onset of stimulus (TTL on)
sti = Trains(3,sti_number)/2000; % unit: Secs
stim_name = [num2str(Trains(2,sti_number)) 'HZ']; % find the stimulus name 
be = sti-60 ; % 60 Secs before the Stim 
ed = sti + 90; % 30s stim + 60s PS
TTL_20HZ = TTL((be*1000*2):(ed*1000*2)); % Extract TTL from specific time window (2000 is sampling rate)
Respiration_smooth_20HZ = Respiration_smooth((be*1000*2):(ed*1000*2)); % Extract the RR from specific time window
tss = [be:0.0005:ed]; % Produce the time-window (Secs)
tss = tss - sti; % Move the original point to the onset of Stim
tss_test = tss; % Rename it
if length(tss)>length(Respiration_smooth_20HZ)
    tss = tss(:,1:end-1);
end 

% Peak detection
[pks_all,locs_all]  = findpeaks(Respiration_smooth_20HZ); % Find the peaks of RR raw data (pks_all: peak amplitude; locs_all: the idx/location of detected peaks)
locs_all_positive = locs_all(find (pks_all>0)); % filter all the positive peaks
pks_all_positive_order = sort(Respiration_smooth_20HZ(locs_all_positive)); % Sort the peaks from larger amplitude to smaller amplitude
if n == 13 % n =13 is a highly sparse RR, it should eliminate less noises
    histo_thre = 0.99;
else 
    histo_thre = 0.97;
end
pks_thre = pks_all_positive_order(round(histo_thre*length(pks_all_positive_order))); % get the threshold based on histo_thre
pks_positive = Respiration_smooth_20HZ(locs_all_positive); % get the idx of positive peaks in raw data 

% Plotting 
if plotting ==1 
figure;
subplot(2,1,1)
histogram(Respiration_smooth_20HZ(locs_all_positive),100)
ylim([0 200])
hold on 
line([pks_thre pks_thre],[0 200],'Color','r')
ylabel('Count');xlabel('Peak Amplitude')
x= locs_all(find(pks_all>pks_thre));
y = pks_all(find(pks_all>pks_thre));
subplot(2,1,2)
plot(tss,Respiration_smooth_20HZ)
hold on 
plot(tss(x), y,'*');
ylabel('Volts');xlabel('Time(S)')
legend('Respiration','Peak')
ylim([-0.2,0.2])
sgtitle([mouse cage]) 
end 

% Moving thresholds
movemax = movmax(Respiration_smooth_20HZ,2000);% Find the moving max every 2000 samples
peak_group = []; 
peak_group_amp = [];
for  m = 1:150
    [pks,locs] = findpeaks(Respiration_smooth_20HZ((2000*(m-1)+1):2000*m)); % get all peaks 
    amp = Respiration_smooth_20HZ(locs+2000*(m-1)); % get the amp of all peaks 
    peak = find (amp> 0.5*max(movemax(2000*(m-1)+1:2000*m))); % filter the peaks that larger than half of the moving max 
    peak = locs (peak)+2000*(m-1); % each run add 2000*(m-1)
    peak_group = [peak_group peak']; % group of peak locus 
    peak_amp = Respiration_smooth_20HZ(peak);
    peak_group_amp = [peak_group_amp peak_amp']; % group of peak amp
end 
idx = find(peak_group_amp<pks_thre/5); % detele some local max-based peaks with too small peak amplitude (that may be the ECG noise) by using pks_thre
peak_group_amp(:,idx) =[];
peak_group(:,idx) = [];

% Fix a sliding window to calculate the mean RR across tw
tw = 3; % 3s *2*1000
Rrs = [];
Rrms = [];
for s = 1:(150/tw)
fpeaks = find(peak_group>=(1+2000*tw*(s-1)) & peak_group<=(2000*tw*s));
if fpeaks == 0 
    Rr = 0;
    Rrm = 0;
else 
    Rr = length(fpeaks)/tw; % breaths per s
    Rrm = (length(fpeaks)/tw)*60; % breath per min 
end 
Rrs = [Rrs repelem(Rr,2000*tw)];% breaths per secs
Rrms = [Rrms repelem(Rrm,2000*tw)]; % breath per mins
end 

thres = [];
for m = 1:150
    thre =0.5*max(movemax(2000*(m-1)+1:2000*m));
    thres = [thres repelem(thre,2000)];
end 
% Plotting of final peak detection and breath per minutes 
if plotting ==1
    figure;
    subplot(2,1,1)
    plot(tss,Respiration_smooth_20HZ)
    hold on 
    plot(tss(peak_group), peak_group_amp,'*');
    hold on 
    if length(tss_test)>length(Respiration_smooth_20HZ)
    plot(tss,thres)
    else
    plot(tss(1:end-1),thres)
    end
    hold on 
    line([tss(1) tss(end)],[pks_thre/5 pks_thre/5],'Color','g')
    title([mouse cage]) 
    ylabel('Volts');xlabel('Time(S)')
    legend('Respiration','Peak','MovThres','Thres')
    xlim([-60 90]);ylim([-0.3 0.5])
    subplot(2,1,2)
    if length(tss) == length(Rrms)
        plot(tss,Rrms)
    else 
        plot(tss(1:end-1),Rrms)
    end 
    ylabel('BPM');xlabel('Time(S)')
    ylim([-10 200])
    xlim([-60 90])
set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[30 25])
% print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\' num2str(tw) '\RR' mouse cage]);
end 

% Self-adjusted sliding time-window to calculate the mean RR across a
% specific time-window
b_Rr = 1/mean(Rrs(1:2000*6)); % b_Rr:the time(Sec) needed for one breath [calculate the baseline BPS(breath per sec) by using the first 6 Secs as a baseline]
% this Rrs is from primiary tw=3
tw = 1.5*b_Rr; % 1.5 time of b_Rr. 1.5 is decided based on multiple tests and its performance

peaks_dig = repelem(0,length(Rrs));
Rrs = []; % new Breath per Secs
Rrms = []; % new breath per mins 
peaks_dig(peak_group) = 1;
Rrs = movmean(peaks_dig,2000*tw)*2000;
Rrms = movmean(peaks_dig,2000*tw)*2000*60;
Rrms_mean = movmean(Rrms,3000);

Rrs = []; % new Breath per Secs
Rrms = []; % new breath per mins 
tw = round(1.5*b_Rr,2); %EDITED TO FIX
if 150/tw == fix(150/tw) % if 150 is divisible by tw
    for s = 1:150/tw
        fpeaks = find(peak_group>=(1+2000*tw*(s-1)) & peak_group<=(2000*tw*s)); % select the peak within a bin 
        if fpeaks == 0 
            Rr = 0;
            Rrm = 0;
        else 
            Rr = length(fpeaks)/tw; % breaths per s
            Rrm = (length(fpeaks)/tw)*60; % breath per min 
        end 
        Rrs = [Rrs repelem(Rr,2000*tw)]; % repeat the value based on the bins
        Rrms = [Rrms repelem(Rrm,2000*tw)];
    end 
else 
    for s = 1:fix(150/tw)  % if 150 is not divisible by tw
        fpeaks = find(peak_group>=(1+2000*tw*(s-1)) & peak_group<=(2000*tw*s));
        if fpeaks == 0 
            Rr = 0;
            Rrm = 0;
        else 
            Rr = length(fpeaks)/tw; % breaths per s
            Rrm = (length(fpeaks)/tw)*60; % breath per min 
        end 
        Rrs = [Rrs repelem(Rr,2000*tw)];
        Rrms = [Rrms repelem(Rrm,2000*tw)];
    end 
    fpeaks = find(peak_group>=(1+2000*tw*fix(150/tw)));
    if fpeaks == 0 
        Rr = 0;
        Rrm = 0;
    else 
        Rr = length(fpeaks)/tw; % breaths per s
        Rrm = (length(fpeaks)/tw)*60; % breath per min 
    end 
        Rrs = [Rrs repelem(Rr,round(2000*tw*(150/tw - fix(150/tw))))];
        Rrms = [Rrms repelem(Rrm,round(2000*tw*(150/tw - fix(150/tw))))];
end 

% Normalize the amplitude
tw = 1.5*b_Rr; % create a sliding window
amp_dig = repelem(0,length(Rrs)); % create a empty vector
amp_dig2 =amp_dig; 
amp_dig(peak_group) = peak_group_amp; % fill the peak locus by using the amplitude value
amp_dig2(peak_group) = 1;
movamp = movmean(amp_dig,2000*tw)*2000*tw; % the moving mean of peak amplitude
movamp2 = round(movmean(amp_dig2,2000*tw)*2000*tw); % the moving count of peaks
movamp_final = movamp./movamp2;
% T = tw*2000;                 % Window Size
% K = length(Rrs) - T;   % Number of repetition
% movamp_tws=peakcounts(K,amp_dig,T,movamp);
% peak_counts_end= length(find(amp_dig(end-2000*tw+1:end)>0));
% movamp_2 = movamp(end-round(2000*tw)+1:end)./repelem(peak_counts_end,round(2000*tw));
% movamp_final = [movamp_tws movamp_2];

%% Plotting 
if plotting ==2
    figure;
    % plot the peak detection 
    subplot(4,1,1);
    plot(tss,Respiration_smooth_20HZ)
    hold on 
    plot(tss(peak_group), peak_group_amp,'*');
    hold on 
        if length(tss)>length(thres)
        plot(tss(1:end-1),thres)
        else
        plot(tss,thres)
        end
    hold on 
    line([tss(1) tss(end)],[pks_thre/5 pks_thre/5],'Color','g')
    title([mouse cage]) 
    ylabel('Volts');xlabel('Time(S)')
    legend('Respiration','Peak','MovThres','Thres')
    xlim([-60 90]);ylim([-0.3 0.5])
    % plot the RR 
    subplot(4,1,2);
    if length(tss) == length(Rrms)
        plot(tss,Rrms)
    else 
        plot(tss(1:end-1),Rrms)
    end 
    ylabel('BPM');xlabel('Time(S)')
    ylim([-10 200])
    xlim([-60 90])
    title('RR with fixed tw')
    % plot the moving RR
    subplot(4,1,3)
    if length(tss) == length(Rrms_mean)
        plot(tss,Rrms_mean)
    else 
        plot(tss(1:end-1),Rrms_mean)
    end 
    ylabel('BPM');xlabel('Time(S)')
    ylim([-10 250])
    xlim([-60 90])
    title('RR with sliding tw')
    % plot the amplitude
    subplot(4,1,4)
    movamp_final2 = movmean(movamp_final,2000);
    if length(tss) == length(movamp_final2)
        plot(tss,movamp_final2)
    else 
        plot(tss(1:end-1),movamp_final2)
    end 
    ylim([0 0.5]);xlim([-60 90])
    ylabel('Volts');xlabel('Time(S)')
     title('Breathing amplitude with sliding tw')
    % save the whole figure
    sgtitle([mouse cage ' self-adjusted tw3' ' ' stim_name])
    set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[25 20])
    %print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\' stim_name '\' mouse cage group]);
    %close all
end 
Rrms_20HZ= [Rrms_20HZ; Rrms_mean]; % final results
% movamp_final2(isnan(movamp_final2))=0;
% Amps =[Amps;movamp_final2] ;% save the amplitude
end 
Rrms_final = Rrms_20HZ;
end 
% tempdir  = [savepath stim_name '\'];
% save([tempdir 'Rrms_' stim_name '.mat'],['Rrms_' stim_name])