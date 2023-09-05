%% 0. Create a cell to save experiment
HBECG= {};
%% 1.1 inmport data
path='D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST'; % Open the folder
%% 1.2 save mouse information
for n = 1:15
 [cage, mouse, day,group,fr]  = readcagemouse(n);
HBECG{n,1} = mouse;
HBECG{n,2} = cage;
HBECG{n,3} = group;
end 
save('HBECG.mat','HBECG');
%% 2.1 Import and analyze BP data from csv files 
test = 'BPwithHR';
nummiceBP = 12; % ther number of mice tested
for n = 1:nummiceBP 
%% 2.2 Load BP Data
[cage, mouse, day,group,fr]  = readcagemouse(n); % read the mouse information (edited in readmouse.m)
Data = readtable([path '\' test '\BP\' group '\' mouse '_' cage '.csv']);
Data = Data(:,[2 7 10]);
Data = table2cell(Data);
Run = Data(:,1);
Run = [Run{:}];
Accepted = Data(:,2);
Accepted = [Accepted{:}];
Mean = Data(:,3); % import the BP value
Mean = [Mean{:}];
if n == 6
ss =10; % the onset of stimulus among run
else 
ss = 11;
end 
bs_range = [(ss-5):(ss-1)];% baseline start 
ss_range = [ss:(ss+9)]; % stimulus
ps_range = [(ss+10):(ss+19)]; %post stimulus 
bs = mean(Mean(bs_range)); % mean baseline 
% normalization by baseline
bs_c = 100*Mean(bs_range)/bs; % propotional change during baseline 
ss_c = 100*Mean(ss_range)/bs; % propotional change during stimulus 
ps_c = 100*Mean(ps_range)/bs; % propotional change post stimulus 
HBECG{n,4} = bs_c;
HBECG{n,5} = ss_c;
HBECG{n,6} = ps_c;
end 
bs_matrixs  = [];
ss_matrixs = [];
ps_matrixs = [];
for n = 1:nummiceBP 
bs_matrixs = [bs_matrixs; HBECG{n,4}];
ss_matrixs = [ss_matrixs;HBECG{n,5}]; 
ps_matrixs = [ps_matrixs;HBECG{n,6}];
end 
bs_matrixs= bs_matrixs';
ss_matrixs =ss_matrixs';
ps_matrixs = ps_matrixs';
matrix=[bs_matrixs;ss_matrixs;ps_matrixs];  % final matrix of all mice and normalized data,used for Prism

% % Additonal group for prism 
% Ai32_mean = mean (matrix(:,1:8),2);
% Ai9_mean = mean(matrix(:,9:12),2);
% Ai32_std = std(matrix(:,1:8)')';
% Ai9_std = std(matrix(:,9:12)')';
% fir = mean(ss_matrixs(1:3,:)); % the first 3 after ss 
% fir_sd = std(ss_matrixs(1:3,:));
%% 2.3 Calculate the first 3 min during Stim in BP
ssc_matrixs = [];
for n = 1:12
ssc_matrixs =HBECG{n,5}(1:3); % calculation for change (100/base)
HBECG{n,7} =ssc_matrixs ;
end 
save('HBECG.mat','HBECG');
save ('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\rank\HBECG.mat',HBECG)
%% 3. LOAD SKIN TEMPERATURE 
test = 'ThermalwithHR';
tem_bin_group = [];
%%
for n = 1:7
[cage, mouse, day,group]  = readcagemouse(n);
after_ss = 60;
nwithinbin = 0;
[tem_bins,tem_ts]= read_skin(path,cage ,test , mouse, group, 0, 260, 0,after_ss,nwithinbin);
tem_bin_group =[tem_bin_group tem_bins];
end 
for n = 9:12
[cage, mouse, day,group]  = readcagemouse(n);
after_ss = 60;
nwithinbin = 0;
[tem_bins,tem_ts]= read_skin(path,cage ,test , mouse, group, 0, 260, 0,after_ss,nwithinbin);
tem_bin_group =[tem_bin_group tem_bins];
end 
mean_bins = [];
for m = 1:5400/60
mean_bin = mean(tem_bin_group (1:60*m,:)); % 2s
mean_bins = [mean_bins; mean_bin];
end 
mean_bins_t = [2:2:5400/30];

%% 4.1 Load ECG data
path='D:\Dropbox\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST'; % Open the folder
test = 'Breathing_Lightgradient'; % open the test folder
savepath = ''; % CHOOSE A FOLDER TO SAVE THE RESULTS
TTL_sum = 3;
plotting =0 ; % if plotting = 1, plot the figures; if not, no plotting
nstart = 1;
nend = 15;
for s = 1:TTL_sum % call for 3 stimulus: 5HZ, 10HZ, and 20HZ
sti_number = s; % the number of stimulus
[HR_group_ps,stim_name,tss] =HR_optimized (test,path,sti_number,nstart,nend);
tempdir  = [savepath stim_name '\'];
save([tempdir 'ECG_' stim_name '.mat'],['HR_group_ps'])
end 
% Files are named as ECG_20HZ.mat ECG_10HZ.mat and ECG_5HZ.mat .... They
% are normalized by the baseline (unit: change%)
%% 4.2 HR Plotting
ECG_20HZ = load('ECG_20HZ.mat');% name HR_group_p
ECG_10HZ = load('ECG_10HZ.mat');
ECG_5HZ = load('ECG_5HZ.mat');
load tss
HR_group_p = ECG_20HZ.HR_group_p; % Need to change 20, 10 and 5
figure; 
for n = 1:15
[cage, mouse, day,group,fr]  = readcagemouse(n);
subplot(3,5,n)
plot(tss(1:end-1),HR_group_p(n,:))
title([mouse ' ' cage ' ' group]);
ylim([-0.5 2])
end
%% 4.3 Save HR for plotting in Prism 
HR_group_p = ECG_20HZ.HR_group_p ; % need to change 20, 10 and 5 
Ai32 = [];
Ai9 = [];
for n =1:15
[cage, mouse, day,group,fr] = readcagemouse(n);
if str2double(group(end)) == 2
    Ai32 = [Ai32 n];
else 
    Ai9 = [Ai9 n];
end 
end 
HR_Ai9_p = HR_group_p(Ai9,:);
HR_Ai32_p=HR_group_p(Ai32,:);
HR_Ai32_meanp = mean(HR_Ai32_p);
HR_Ai9_meanp = mean(HR_Ai9_p);
HR_Ai9_p_d =downsample(HR_Ai9_p,1000);
HR_Ai32_p_d = downsample(HR_Ai32_p,1000);
tss_d = downsample(tss,1000);
if length(tss_d) > length(HR_Ai32_p_d)
    HR_group_d_p = [tss_d(1:end-1) HR_Ai32_p_d HR_Ai9_p_d];
else 
    HR_group_d_p = [tss_d HR_Ai32_p_d HR_Ai9_p_d]; % group data with time, downsampled Ai32 and Ai9
end 
%% 4.4 Save 20HZ ECG to a big cell HBECG
for n = 1:15
HBECG{n,10} = ECG_20HZ.HR_group_p (n,(60*2000+1):90*2000); % during stimulation respiration 20 HZ
end 
save('HBECG.mat','HBECG');