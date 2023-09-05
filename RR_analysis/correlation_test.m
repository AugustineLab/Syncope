load HBECG 
%% Give the group idx (Ai32 or Ai9)
 groups = [];
 for n = 1:15
 [cage, mouse, day,group,fr] = readcagemouse(n);
groups = [groups group];
end
%% fill the empty data as NaN
empties = cellfun('isempty',HBECG);
HBECG(empties) = {NaN}; % fill the empty region as NAN
%% Import all parameters 
dBPs =[];
dRRs =[];
dapneas =[];
dHRs = [];
dRTPPs = [];
dflips = [];
for n = 1:15
dBP = mean(HBECG{n,5}); % the first 3min in BP sti || 20HZ
dRR = HBECG{n,8}(1:1/3*(length(HBECG{n,8}))); % the first 10s of RR 
dapnea = HBECG{n,9}; % apnea duration(s) 20HZ
dHR = HBECG{n,10}(1:2/3*(length(HBECG{n,8}))); % the first 20s of HR
dRTPP = HBECG{n,11};%RTPP mean(day3-5) - mean(day 1-2) || 20HZ 
dflip = HBECG{n,12}; %Flip duration (s) || 20HZ control 
% mean or other 
dRR = mean(dRR);
dHR = mean(dHR);
dapnea = -dapnea;
dflip = -dflip;
%save 
dBPs =[dBPs; dBP];
dRRs =[dRRs; dRR];
dapneas =[dapneas; dapnea];
dHRs = [dHRs; dHR];
dRTPPs = [dRTPPs; dRTPP];
dflips = [dflips; dflip];
end 
NPY2Rmatrix = [dBPs dRRs dapneas dHRs dRTPPs dflips];
%%  Z score normalization of each parameter set
NPY2Rmatrix_z=[];
for m = 1: 6
p = find(~isnan(NPY2Rmatrix(:,m)));
[Z,mu,sigma] = zscore(NPY2Rmatrix(p,m));
for n = 1:length(p)
NPY2Rmatrix_z (p(n),m) = Z(n);
end 
end 
%% correlation matrix 
figure
clims = ([-2 2]);
imagesc(NPY2Rmatrix_z,clims)
% xticks = 1:5:26;  %adjust as appropriate, positive integers only
% xlabels = {'dBPs', 'dRRs', 'dapneas', 'dHRs', 'dRTPPs', 'dflips'};  %time labels
% set(gca, 'XTick',  1:5:26, 'XTickLabel', xlabels);
xt = get(gca, 'XTick');                                             % Original 'XTick' Values
yt = get(gca, 'YTick');                                             % Original 'XTick' Values
yt = 1:15;
xtlbl = linspace(-150.36, 265.8773, numel(xt));                     % New 'XTickLabel' Vector
xlabels = {'dBPs', 'dRRs', 'dapneas', 'dHRs', 'dRTPPs', 'dflips'};
set(gca, 'XTick',xt, 'XTickLabel',xlabels , 'XTickLabelRotation',30)   % Label Ticks
ytlbl = linspace(-150.36, 265.8773, numel(yt));                     % New 'XTickLabel' Vector
%ylabels = {'11N','3R', '6L', '32N' ,'18N', '1N', '16N', '34L', '11R', '12L' ,'15L' ,'16N' ,'5R' ,'21R', '20R'};
ylabels = {'NPY2RxAi32','NPY2RxAi32', 'NPY2RxAi32', 'NPY2RxAi32','NPY2RxAi32','NPY2RxAi32', 'NPY2RxAi32','NPY2RxAi32','NPY2RxAi9', 'NPY2RxAi9','NPY2RxAi9','NPY2RxAi9','NPY2RxAi32','NPY2RxAi32' ,'NPY2RxAi9'};
set(gca, 'YTick',yt, 'YTickLabel',ylabels , 'YTickLabelRotation',30)   % Label Ticks
colorbar ('Position',[0.92 0.4 0.02 0.33],'Location','east',...
        'Ticks',[-2 -1 1 2],...
    'TickLabels',{'High','','','Low'});
%% Aligned correlation matrix
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

NPY2Rmatrix_z_align=[NPY2Rmatrix_z(Ai32,:); NPY2Rmatrix_z(Ai9,:)] ;
[sA,index] = sort(mean(NPY2Rmatrix_z_align,2));
NPY2Rmatrix_z_align_sort = NPY2Rmatrix_z_align(index,:);
figure
clims = ([-2 2]);
imagesc(NPY2Rmatrix_z_align_sort ,clims)
% xticks = 1:5:26;  %adjust as appropriate, positive integers only
% xlabels = {'dBPs', 'dRRs', 'dapneas', 'dHRs', 'dRTPPs', 'dflips'};  %time labels
% set(gca, 'XTick',  1:5:26, 'XTickLabel', xlabels);
xt = get(gca, 'XTick');                                             % Original 'XTick' Values
yt = get(gca, 'YTick');                                             % Original 'XTick' Values
yt = 1:15;
xtlbl = linspace(-150.36, 265.8773, numel(xt));                     % New 'XTickLabel' Vector
xlabels = {'dBPs', 'dRRs', 'dapneas', 'dHRs', 'dRTPPs', 'dflips'};
set(gca, 'XTick',xt, 'XTickLabel',xlabels , 'XTickLabelRotation',30)   % Label Ticks
ytlbl = linspace(-150.36, 265.8773, numel(yt));                     % New 'XTickLabel' Vector
%ylabels = {'11N','3R', '6L', '32N' ,'18N', '1N', '16N', '34L', '11R', '12L' ,'15L' ,'16N' ,'5R' ,'21R', '20R'};
ylabels = {'NPY2RxAi32','NPY2RxAi32', 'NPY2RxAi32', 'NPY2RxAi32','NPY2RxAi32','NPY2RxAi32', 'NPY2RxAi32','NPY2RxAi32','NPY2RxAi32', 'NPY2RxAi32','NPY2RxAi9','NPY2RxAi9','NPY2RxAi9','NPY2RxAi9' ,'NPY2RxAi9'};
ylabels = ylabels(index);
set(gca, 'YTick',yt, 'YTickLabel',ylabels , 'YTickLabelRotation',30)   % Label Ticks
colorbar ('Position',[0.92 0.4 0.02 0.33],'Location','east',...
        'Ticks',[-2 -1 1 2],...
    'TickLabels',{'1','','','-1'});
colormapeditor
set (gcf,'PaperPosition',[-1,10,25,20],'PaperSize',[30 25])
print(gcf,'-dtiff','-r300',['D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\rank\parameter']);
%% cross-correlation of RR change and HR change
ECG_20HZ = load('ECG_20HZ.mat');% name HR_group_p
HR_group_p = ECG_20HZ.HR_group_p;
load Rrms_20HZ
figure
for n = 1
subplot(3,5,n)
[cage, mouse, day,group,fr]  = readcagemouse(n);
[r, lags] = xcorr(Rrms_20HZ(n,:),HR_group_p(n,:));
stem(lags,r,'MarkerSize',1)
xlabel('lags'); ylabel('cor idx')
title([mouse ' ' cage ' ' group]);
end
%% GLM (apnea)
load apneas % 3*15 matrix (3 rows:20HZ 10HZ 5HZ)
load Rrms_final % cell
load Amps_group % cell 
% feature matrix 
% predictors 1. AI32 or Ai9  2. BS RR 3.BS HR 4. 20HZ iso HR reduction
% binary: apnea 
Genotype = ['Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai32';'Ai9 ';'Ai9 ';'Ai9 ';'Ai9 ';'Ai9 '];
RRbs = RR_group{3,1};
HRbs = HR(:,1);
dHR= HR(:,2);
apnea = apneas(1,:)';
apnea2 = apnea;
apnea2(find(apnea>0))=1;
feature = table(Genotype,RRbs,HRbs,dHR,apnea,apnea2);

%  initcoeff = glmfit(feature(:, 1:(w-1)), [feature(:, w) N], ...
%    'binomial', 'link', 'logit');
% 
% mdl = fitlm(hospital,'interactions','ResponseVar','Weight',...
%     'PredictorVars',{'Sex','Age','Smoker'},...
%     'CategoricalVar',{'Sex','Smoker'})

model = fitlm(feature,'ResponseVar','apnea',...
    'PredictorVars',{'Genotype','RRbs','HRbs','dHR'},...
    'CategoricalVar',{'Genotype'});
model.Coefficients

modelspec = 'apnea2 ~ Genotype+dHR+RRbs+HRbs';
mdl = fitglm(feature,modelspec,'Distribution','binomial')
%% Load data
ECG_20HZ = load('ECG_20HZ.mat');% name HR_group_p
ECG_10HZ = load('ECG_10HZ.mat');
ECG_5HZ = load('ECG_5HZ.mat');
Rrms_20HZr = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\20HZ\Rrms_20HZ.mat');
Rrms_10HZr = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\10HZ\Rrms_10HZ.mat');
Rrms_5HZr = load('D:\Dropbox (Scripps Research)\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST\Breathing_Lightgradient\Breathing_results\sw\5HZ\Rrms_5HZ.mat');
order_google= [1 3 13 6 2 7 5 4 8 14 9 10 11 12 15] ;
ECG_20HZ = ECG_20HZ.HR_group_p(order_google,:) ;% name HR_group_p
ECG_10HZ = ECG_10HZ.HR_group_p(order_google,:);
ECG_5HZ = ECG_5HZ.HR_group_p(order_google,:);
Rrms_20HZ= Rrms_20HZr.Rrms_final(order_google,:);
Rrms_10HZ = Rrms_10HZr.Rrms_final(order_google,:);
Rrms_5HZ = Rrms_5HZr.Rrms_final(order_google,:);
%% 
cc = {};
cc{1,1} = ECG_20HZ;
cc{2,1} = ECG_10HZ;
cc{3,1} = ECG_5HZ;
cc{1,2} = Rrms_20HZ;
cc{2,2} = Rrms_10HZ;
cc{3,2} = Rrms_5HZ ;
%% Cross-correlation
%Cross-correlation measures the similarity between a vector x and shifted 
%(lagged) copies of a vector y as a function of the lag. If lags >0, it
%means x is prior to the y.
% y = x+ lag
for s = 1:TTL_sum 
    if s == 1
    name = '20HZ';
    end 
    if s ==2 
        name = '10HZ';
    end 
    if s ==3
       name = '5HZ';
    end 
    figure
    for n = 1:15
        subplot(3,5,n)
        [c,lags] = xcorr(cc{s,1}(n,:),cc{s,2}(n,:));
        stem(lags,c)
        lag = lags(find(c ==max(c)));
        hold on 
        plot(lag, max(c),'Color','r')
        title(['Lag' ' = ' num2str(lag)])
    end 
            sgtitle(name)
end 
%% Plot the HR(%) and RR together, Label the apnea onset and its corresponding HR%
for s = 1:3
    if s == 1 
        name = '20HZ';
    end 
    if s ==2 
        name = '10HZ';
    end 
    if s ==3
       name = '5HZ';
    end 
    figure
    for n =1:15
        subplot(3,5,n)
        Rrm=cc{s,2}(n,:);
        HR =(cc{s,1}(n,:)-1).*100;
        idx = find(Rrm== 0);
        yyaxis left
        plot(tss(1:end),Rrm)
        ylabel('RR(BPM)')
        ylim([-100 300])
        xlim([-60 90])
        yyaxis right
        plot(tss(1:end),HR)
        ylabel('HR(%)')
        ylim([-100 100])
        xlim([-60 90])
        if isempty(idx)
            idx1 = 0;
        else 
            idx(find(idx<2000*60)) = [];
            idx1 = idx(1);
            HR1 = HR(idx1);
            hold on
            line([idx1/120000 idx1/120000],[-100 HR1],'color','r')
            text (idx1/120000,-50,['HR% =' num2str(HR1)])
        end 
    end 
end 
