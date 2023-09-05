
function [HR_group_ps,stim_name,tss] =HR_optimized (test,path,sti_number,nstart,nend)
HR_group_ps = [];
for n =nstart:nend
[cage, mouse, day,group,fr]  = readcagemouse(n); % Read the mouse information (edited in readmouse.m)
filename = [ mouse '_' cage '.mat']; % Get the name of mat file downloaded from Acknowledge
Data = load([path '\' test '\' group '\' filename]); % Load data
% Import TTL and HR
for m= 4:size(Data.labels,1) % In Acknowledge the smoothed (3000 sample mean smoothing) Breathing raw data should be named as RR
    label = Data.labels (m,:);
    label_logic = (label (1:2)== 'HR' ); % find the channel number of HR
    if  label_logic(2) ==1 
        HR = m;
    end 
end 
TTL = Data.data (:,2); % Load TTL named as TTL
HR_smooth = Data.data (:,HR); % Load processed HR named as HR_smooth
if fr == 0.25 % some data is recorded in different SR
TTL = downsample(TTL,2);
HR_smooth = downsample(HR_smooth,2);
end 
Trains = find_stim_train2(TTL, 2000, 5); % use the function 'find_stim_train2.mat' to find the onset of TTL
% Find the onset of stimulus (TTL on)
sti = Trains(3,sti_number)/2000; % unit: Secs
stim_name = [num2str(Trains(2,sti_number)) 'HZ']; % find the stimulus name 
be = sti-60 ; % 60 Secs before the Stim 
ed = sti + 90; % 30s stim + 60s PS
TTL_20HZ = TTL((be*1000*2):(ed*1000*2)); % Extract TTL from specific time window
HR_smooth_20HZ = HR_smooth((be*1000*2):(ed*1000*2)); % Extract the HR from specific time window
tss = [be:0.0005:ed]; % Produce the time-window (Secs)
tss = tss - sti; % Move the original point to the onset of Stim
tss_test = tss; % Rename it
if length(tss)>length(HR_smooth_20HZ)
    tss = tss(:,1:end-1);
end 

bs = HR_smooth_20HZ(1:2000*60);
ss = HR_smooth_20HZ(2000*60+1:2000*90);
ps = HR_smooth_20HZ(2000*90+1:end);
mean_bs = mean(bs);
bs_HR = bs./mean_bs;
ss_HR = ss./mean_bs; 
ps_HR = ps./mean_bs; 
HR_group_p = [bs_HR' ss_HR' ps_HR']; %% the /baseline data in whole group
HR_group_p = HR_group_p*100;
HR_group_ps =[HR_group_ps;HR_group_p];
end 
end 