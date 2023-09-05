%% 1.1 Import path
% path='D:\Dropbox\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST'; % % Open the folder in Hanbing's computer 
path='D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_ECG\npy2r_branch_opto\RR_analysis\trunk'; % % Open the folder in GPU station
test = 'Matfiles'; % Enter the name of experiment folder
savepath = 'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_ECG\npy2r_branch_opto\RR_analysis\trunk\Matfiles\Results\'; % enter the folder to save results
%% 1.2 Enter the parameters
prompt = {'Enter # of the number of TTL trials:','Enter the start mouse number:','Enter the last mouse number:','Enter the plotting type','Enter whether smoothing(1) or not(0)'};
dlg_title = 'Settings';
num_lines = 1;
defaultans = {'3','1','12','2','0'}; % change these values to alter default input!!!
settingsUI = inputdlg(prompt,dlg_title,num_lines,defaultans);
TTL_sum = str2double(settingsUI{1});
nstart = str2double(settingsUI{2});
nend = str2double(settingsUI{3});
plotting = str2double(settingsUI{4}); % if plotting = 1, plot the figures; if not, no plotting; if plot = 2, only the final plot
smoothing = str2double(settingsUI{5}); 
%% 2. Calculate and save Rrms for each frequency Stimus (if done can skip this step)
for s = 1:TTL_sum % call for 3 stimulus: 5HZ, 10HZ, and 20HZ
sti_number = s; % the number of stimulus
[Rrms_final,stim_name,tss,Amps] = RR_preprocess_m (test,path,sti_number,plotting,nstart,nend,smoothing); % big function that calculate the Breath per minutes for each mouse
% Rrms format: 15 * 30000 (15 mice * 30000 samples) with SR = 2000. Rrms is
% not normalized, and is the real value of breathing rate per minute. 
tempdir  = [savepath stim_name '\'];
save([tempdir 'Rrms_' stim_name '.mat'],['Rrms_final']) % save Rrms
save ([tempdir 'Amps' stim_name '.mat'],['Amps']) % save the amplitude 
end 
%% 3. Load Rrms and Amps
ouputfolder = 'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_ECG\npy2r_branch_opto\RR_analysis\trunk\Matfiles\Results\';
load(strcat(ouputfolder, '20HZ\Rrms_20HZ.mat'));
Rrms_20HZ = Rrms_final;
% load(strcat(ouputfolder, '10HZ\Rrms_10HZ.mat'));
% Rrms_10HZ = Rrms_final;
% load(strcat(ouputfolder, '5HZ\Rrms_5HZ.mat'));
% Rrms_5HZ = Rrms_final;
%Amps_20HZ = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\20HZ\Amps20HZ.mat');
%Amps_10HZ = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\10HZ\Amps10HZ.mat');
%Amps_5HZ = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\5HZ\Amps5HZ.mat');
load tss 
tss = tss(1:end-1);
%% 4. Post process (baseline-based normalization & sorted for Prism plotting)
pp_group_prism = {};
for s = 1:TTL_sum 
if s == 1 
    test = Rrms_20HZ;
    name = '20HZ';
end 
if s ==2 
    test = Rrms_10HZ;
    name = '10HZ';
end 
if s ==3
   test = Rrms_5HZ;
   name = '5HZ';
end 

bs = test(:,1:2000*60);
ss = test(:,2000*60+1:2000*90);
ps = test(:,2000*90+1:end);

mean_bs = mean(bs,2);
bs_p = bs./mean_bs;
ss_p = ss./mean_bs; 
ps_p = ps./mean_bs; 
pp = [bs_p ss_p ps_p];

Ai32 = [];
Ai9 = [];
for n = nstart:nend % need to be optimized
[cage, mouse, ~, group, fr]  = readcagemouse(n);

%if str2double(group(end)) == 2
%    Ai32 = [Ai32 n];
%else 
%    Ai9 = [Ai9 n];
%end 
    if strcmp(group,'trunk')
        Ai32 = [Ai32 n];
    else
        Ai9 = [Ai9 n];
    end
end 

pp_Ai32 = pp(Ai32,:);
pp_Ai9 = pp(Ai9,:);

figure;
for n = 1:7
    subplot(5,2,n) 
plot(tss(1:end),pp_Ai32(n,:))
ylabel('BPM(/Baseline)');xlabel('Time(S)')
ylim([0 4])
xlim([-60 90])
[cage, mouse, day,group,fr] = readcagemouse(Ai32(n));
title([mouse ' ' cage ' ' 'Ai32' ' ' name])
end
set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[30 25])
%print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\' name '\Ai32 group']);
figure;
%         for n = 1:6
%                 subplot(5,2,n) 
%             plot(tss(1:end),pp_Ai9(n,:))
%             ylabel('BPM(/Baseline)');xlabel('Time(S)')
%             ylim([0 4])
%             xlim([-60 90])
%             [cage, mouse, day,group,fr] = readcagemouse(Ai9(n));
%             mouse
%             title([mouse ' ' cage ' ' 'Ai9' ' '  name])
%         end
set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[30 25])
%print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\' name '\Ai9 group']);

pp_Ai32_d = pp_Ai32(:,1:2000:end);
        %pp_Ai9_d = pp_Ai9(:,1:2000:end);
tss_d = tss(1:2000:end-1);
pp_group_prism{s} = [tss_d' pp_Ai32_d'];%[tss_d' pp_Ai32_d' pp_Ai9_d'];
end 
%close all;
%% 5.1 Rearange dataset for Google sheet 
% Make the mouse order same as the order in Google_sheet
order_google= [1 3 13 6 2 7 5 4 8 14 9 10 11 12 15] ;
Rrms_20HZ= Rrms_20HZr.Rrms_final(order_google,:);
Rrms_10HZ = Rrms_10HZr.Rrms_final(order_google,:);
Rrms_5HZ = Rrms_5HZr.Rrms_final(order_google,:);
Amps_20HZ= Amps_20HZr.Amps(order_google,:);
Amps_10HZ = Amps_10HZr.Amps(order_google,:);
Amps_5HZ = Amps_5HZr.Amps(order_google,:);
%% 5.2 Produce a normalized RR matrix [ For Prism ]
RR_group={};
TTL_sum  =3;
for s = 1:TTL_sum 
if s == 1 
    test = Rrms_20HZ;
    name = '20HZ';
end 
if s ==2 
    test = Rrms_10HZ;
    name = '10HZ';
end 
if s ==3
   test = Rrms_5HZ;
   name = '5HZ';
end 

bs = test(:,1:2000*60); % Stim()
ss = test(:,2000*60+1:2000*90); % Stim 
ps = test(:,2000*90+1:end); % Post 

mean_bs = mean(bs,2);
bs_p = bs./mean_bs;
ss_p = ss./mean_bs; 
ps_p = ps./mean_bs; 

RR_all = [bs_p(:,1:2000:end)'; ss_p(:,1:2000:end)'; ps_p(:,1:2000:end)'];
RR_group{s,7} = [tss(1:2000:end-1)' RR_all];
end 
save([ouputfolder 'final_RR.mat'],'RR_group')
%% 5.3 Produce a normalized RR matrix [For google sheet]
for s = 1:TTL_sum 
    if s == 1 
        test = Rrms_20HZ;
        name = '20HZ';
    end 
    if s ==2 
        test = Rrms_10HZ;
        name = '10HZ';
    end 
    if s ==3
       test = Rrms_5HZ;
       name = '5HZ';
    end 

bs = test(:,1:2000*60); % Stim()
ss = test(:,2000*66+1:2000*75); % Stim (6s - 15s)
ps = test(:,2000*96+1:2000*105); % Post (36s - 45s)

mean_bs = mean(bs,2);
mean_ss = mean(ss,2);
mean_ps = mean(ps,2);
RR_group{s,1}=mean_bs;
RR_group{s,2}=mean_ss;
RR_group{s,3}=mean_ps;
RR_group{s,4}=((mean_ss./mean_bs) -1)*100;
RR_group{s,5}=((mean_ps./mean_bs) -1)*100;
RR_group{s,6} = [mean_bs mean_ss (((mean_ss./mean_bs) -1)*100) mean_ps (((mean_ps./mean_bs) -1)*100)]; % this is pasted in Google sheet
RR_group{s,8} = test./mean_bs*100-100;
end 
%% 5.4 Amplitude analysis [for Prism]
Amps_group = {};
TTL_sum =3;
for s = 1:TTL_sum 
if s == 1 
    test = Amps_20HZ;
    name = '20HZ';
end 
if s ==2 
    test = Amps_10HZ;
    name = '10HZ';
end 
if s ==3
   test = Amps_5HZ;
   name = '5HZ';
end 

bs = test(:,1:2000*60); % Stim()
ss = test(:,2000*60+1:2000*90); % Stim (6s - 15s)
ps = test(:,2000*90+1:end); % Post (36s - 45s)

mean_bs = mean(bs,2); 
if  s ==3
    mean_bs (13,:) = 0.000001;% non-zero
end 
bs_p = bs./mean_bs;
ss_p = ss./mean_bs; 
ps_p = ps./mean_bs; 
 
Amps_all = [bs_p(:,1:2000:end)'; ss_p(:,1:2000:end)'; ps_p(:,1:2000:end)'];
Amps_group{s,7} = [tss(1:2000:end)' Amps_all*100-100];
end 
%% 6.1 Apnea calculation 
apnea_Ai32_group = [];
apnea_Ai9_group = [];
for s = 1:3
if s == 1 
    test = Rrms_20HZ;
    name = '20HZ';
end 
if s ==2 
    test = Rrms_10HZ;
    name = '10HZ';
end 
if s ==3
   test = Rrms_5HZ;
   name = '5HZ';
end 

bs = test(:,1:2000*60);
ss = test(:,2000*60+1:2000*90);
ps = test(:,2000*90+1:end);

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

ss_Ai32 = ss(Ai32,:);
ss_Ai9 = ss(Ai9,:);
apnea_Ai32 =[];
apnea_Ai9 = [];
for m = 1:10
apnea_Ai32 = [apnea_Ai32 length(find (ss_Ai32(m,:) == 0))];
end 

for m = 1:5
apnea_Ai9 = [apnea_Ai9 length(find(ss_Ai9(m,:) == 0))];
end 

apnea_Ai32_group = [apnea_Ai32_group apnea_Ai32'];
apnea_Ai9_group = [apnea_Ai9_group apnea_Ai9'];
end 
apnea_Ai32_group =apnea_Ai32_group';
apnea_Ai9_group = apnea_Ai9_group';

apnea_group = [apnea_Ai32_group apnea_Ai9_group];
apnea_group =apnea_group/2000;
%% 6.2 Apnea calculation NEW
apneas =[];
for s = 1:3
if s == 1 
    test = Rrms_20HZ;
    name = '20HZ';
end 
if s ==2 
    test = Rrms_10HZ;
    name = '10HZ';
end 
if s ==3
   test = Rrms_5HZ;
   name = '5HZ';
end 

ss = test(:,2000*60+1:2000*90);
apnea =[];
for m = 1:15
apnea = [apnea length(find(ss(m,:) == 0))];
end 
apneas = [apneas;apnea];
end 
apneas =apneas/2000;