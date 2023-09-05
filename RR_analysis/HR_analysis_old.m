%% Path 
path='D:\Dropbox\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST'; % Open the folder
test = 'Breathing_Lightgradient'; % open the test folder
%% HR_analysis
% ECG should first be processed by Acknowledge and transformed to HR, and
% same frequency stimulus should be put in the same file and export as Mat
% file. 
inte = 20; % stimulation fre
filename = [num2str(inte) 'HZ_group_HR.mat'];
Data_ECG = load([path '\' test '\' 'HR_analysis\' filename]);
labels = Data_ECG.labels;
mark_20HZ = [];
for n = 1:15
[cage, mouse, day,group,fr]  = readcagemouse(n);
cname = [num2str(inte) 'Hz_' mouse '_' cage];

for m= 1:size(labels,1)
    label = labels (m,:);
    if inte == 10 || inte == 20
        label_logic = (label (1:13)== cname(1:13));
        if  length(find(label_logic ==1))>12
        mark = m;
        else
        end 
    else 
        label_logic = (label (1:12)== cname(1:12));
        if  length(find(label_logic ==1))>11
        mark = m;
        else
        end 
    end
end
if n ==7 
    if inte == 20 
    mark = 19;
    end 
    if inte == 5 || inte == 10
        mark =14;
    end 
end
mark_20HZ = [mark_20HZ mark]; % the channel number of each mouse
end 

HR = Data_ECG.data;
HR = downsample(HR,2);
HR = HR(2:end,:);
HR =HR';

HR_group =[];
for n = 1:15
[cage, mouse, day,group,fr]  = readcagemouse(n);
mark = mark_20HZ(n);
HR_mouse = HR(mark,:);
HR_group =[HR_group;HR_mouse];
end 

bs = HR_group(:,1:2000*60);
ss = HR_group(:,2000*60+1:2000*90);
ps = HR_group(:,2000*90+1:end);
mean_bs = mean(bs,2);
bs_HR = bs./mean_bs;
ss_HR = ss./mean_bs; 
ps_HR = ps./mean_bs; 
HR_group_p = [bs_HR ss_HR ps_HR]; %% the /baseline data in whole group
tempdir = ['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\HR_analysis\']
save([tempdir 'ECG_' num2str(inte) 'HZ.mat'],'HR_group_p')
%% analysis 
load breathing_20HZ.mat % name pp
ECG_20HZ = load('ECG_20HZ.mat');% name HR_group_p
ECG_10HZ = load('ECG_10HZ.mat');
ECG_5HZ = load('ECG_5HZ.mat');
load tss % load time-series
tss = tss(1:end-1);
%% cross-correlation of RR change and HR change
figure
for n = 1:15
subplot(3,5,n)
[cage, mouse, day,group,fr]  = readcagemouse(n);
[r, lags] = xcorr(pp(n,:),HR_group_p(n,:));
stem(lags,r,'MarkerSize',1)
xlabel('lags'); ylabel('cor idx')
title([mouse ' ' cage ' ' group]);
end
%% Plot HR
figure; 
for n = 1:15
[cage, mouse, day,group,fr]  = readcagemouse(n);
subplot(3,5,n)
plot(tss(1:end),HR_group_p(n,:))
title([mouse ' ' cage ' ' group]);
ylim([-0.5 2])
end
suptitle('20HZ HR (/baseline)')
%% plot HR and RR together
figure; 
for n = 1:15
[cage, mouse, day,group,fr]  = readcagemouse(n);
subplot(5,3,n)
yyaxis left
plot(tss(1:end),pp(n,:))
ylabel('RR(/bs)')
ylim([-0.5 10])
yyaxis right
plot(tss(1:end),HR_group_p(n,:))
ylabel('HR(/bs)')
title([mouse ' ' cage ' ' group]);
ylim([-0.5 2])
xlim([-60 90])
end
suptitle('20HZ HR & RR')
set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[30 25])
print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\correlation test\20HZ group']);

%% Save HR data for Prism
HR_group_p = ECG_20HZ.HR_group_p ;
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
HR_Ai9_p_d =downsample(HR_Ai9_p',1000);
HR_Ai32_p_d = downsample(HR_Ai32_p',1000);
tss_d = downsample(tss',1000);
HR_group_d_p = [tss_d HR_Ai32_p_d HR_Ai9_p_d]; % group data with time, downsampled Ai32 and Ai9

%%
pp_group_prism = pp_group_prism(:,2:end);
HR_group_d_p =HR_group_d_p(:,2:end);
%% Linear correlation
figure
for n = 1:15
    
[rho,pval] = corr(pp_group_prism(60:90,n),downsample(HR_group_d_p(120:180,n),2));
if pval<0.05
   n 
   [cage, mouse] = readcagemouse(n)
end 
hold on 
scatter(pp_group_prism(60:90,n), downsample(HR_group_d_p(120:180,n),2))
end

figure
for n = 1:15
[rho,pval] = corr(pp_group_prism(60:90,n),downsample(HR_group_d_p(120:180,n),2));
subplot(5,3,n)
if pval<0.05
   n 
   [cage, mouse] = readcagemouse(n)
   scatter(pp_group_prism(60:90,n), downsample(HR_group_d_p(120:180,n),2),'filled','MarkerFaceColor','r')
   text(1,1,['R =' num2str(rho) ' P value = ' num2str(pval)]) 
else 
    scatter(pp_group_prism(60:90,n), downsample(HR_group_d_p(120:180,n),2),'filled')
    text(1,1,['R =' num2str(rho) ' P value = ' num2str(pval)]) 
end
xlabel('RR');ylabel('HR');
xlim([0 3]); ylim([0 3])
title([mouse cage])
end

for n = 1:15
[rho,pval] = corr(pp_group_prism(1:60,n),downsample(HR_group_d_p(1:120,n),2));
scatter(pp_group_prism(1:60,n), downsample(HR_group_d_p(1:120,n),2))
if pval<0.05
   n 
   [cage, mouse] = readcagemouse(n)
   figure
end 
end
%% Save ECG to HBECG
for n = 1:15
HBECG{n,10} = ECG_20HZ.HR_group_p (n,(60*2000+1):90*2000); % during stimulation respiration 20 HZ
end 
save('HBECG.mat','HBECG');