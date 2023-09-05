%% Data settings Dialog Boxes
prompt = {'Enter # of Channels:','Enter Sampling Rate (Hz):','Enter # of Wavelets:','Enter Start time (ms)','Enter End time (ms)','Enter # of Monte Carlo Permutations (at least 1000 recomended)','# of rows to skip'};
dlg_title = 'Settings';
num_lines = 1;
defaultans = {'1','2500','80','0','50000','5000','0'}; % change these values to alter default input!!!
settingsUI = inputdlg(prompt,dlg_title,num_lines,defaultans);
channels = str2double(settingsUI{1});
samplingrate = str2double(settingsUI{2});
HzBin = str2double(settingsUI{3});
StartTimeUI = str2double(settingsUI{4});
EndTimeUI = str2double(settingsUI{5});
perms = str2double(settingsUI{6});
rowskip = str2double(settingsUI{7});
%% importdata
testfiledir = uigetdir();
matfiles = dir(fullfile(testfiledir, '*.csv'));
nfiles = length(matfiles);
data  = cell(1,nfiles);
for w = 1 : nfiles
   fid = fullfile(testfiledir, matfiles(w).name);
   data{w} = readmatrix(fid);
end
load Freq.mat
%% Data settings Dialog Boxes
loopnum = 0;
inttrue = false;
while inttrue==false
testint = (((EndTimeUI-loopnum)-StartTimeUI)/1000)*samplingrate;
    if floor(testint)==testint
        inttrue = true;
        EndTimeUI = EndTimeUI-loopnum;
        disp(strcat('NewEndTime:',num2str(EndTimeUI)))
    elseif EndTimeUI-loopnum == 0
        disp('not possible') 
        break
    else
        loopnum = loopnum+1;
        inttrue = false;
    end
end
timebinoptions1 = divisors(testint);
binopIND = 1;
timebinoptions2 = {length(timebinoptions1)};
for binop = 1 : length(timebinoptions1)
    if timebinoptions1(binop) >= 50 && timebinoptions1(binop) <=600
        binopIND = binopIND+1;
        timebinoptions2{binopIND} = strcat(num2str(timebinoptions1(binop)),',');
    end
end
timebinoptions2{1} = 'How many time bins? (recomended entered below):';
prompt = strjoin(timebinoptions2);
dlg_title = 'Settings';
num_lines = 1;
[~,recbins] = min(abs(timebinoptions1-300));
defaultans = {num2str(timebinoptions1(recbins))};
settingsUI = inputdlg(prompt,dlg_title,num_lines,defaultans);
timebin = str2double(settingsUI{1});
samplesbin = (((EndTimeUI-StartTimeUI)/1000)*samplingrate)/timebin;
startsample = floor((StartTimeUI/1000)*samplingrate);
list = {matfiles(:).name};
[group1] = listdlg('PromptString','Ctrl+click to Select control group files','ListString',list);

%% create ROI from data, and reduce matrix size
binData = zeros(HzBin,timebin,channels,nfiles);
for w = 1 : nfiles
    current = cell2mat(data(1,w));
    for chan = 1:channels
        for c = 1:HzBin
            for j = 1 : timebin
                s = (startsample+(j*samplesbin))-(samplesbin-1);
                e = startsample+(j*samplesbin);
                binData(c,j,chan,w) = mean(current(((HzBin*chan)-HzBin)+c,s:e));
            end
        end
    end
end
binData = flipud(binData);
%% Create Real Group Data
group1n = 0;
group2n = 0;
group1Chan = zeros(HzBin,timebin,channels,size(group1,2));
group2Chan = zeros(HzBin,timebin,channels,nfiles-size(group1,2));
for m = 1:size(binData,4)
    if ismember(m,group1)
        group1n = group1n+1;
        group1Chan(:,:,:,group1n) = binData(:,:,:,m);
    else
        group2n = group2n+1;
        group2Chan(:,:,:,group2n) = binData(:,:,:,m);
    end
end
group1means = mean(group1Chan,4);
group2means = mean(group2Chan,4);
tChan = zeros(HzBin,timebin,channels);
for chan = 1:channels
    for cc = 1:timebin
        for rr = 1:HzBin
            [~,~,~,tempstatAC] = ttest2(group1Chan(rr,cc,chan,:),group2Chan(rr,cc,chan,:));
            tChan(rr,cc,chan) = [tempstatAC.tstat];
        end
    end
end
%% Monte Carlo permutation: Tperm,POS and NEG ClusterSizes
group1size = size(group1,2);
group2size = nfiles-size(group1,2);
bothgroupsize = group1size+group2size;
parfor permu = 1:perms
    for chan = 1:channels
        [TpermChan{permu,chan},ClustSizePOSpermChan{permu,chan},ClustSizeNEGpermChan{permu,chan}] = mymontecarlo(binData(:,:,chan,:),group1n,group2n);
    end
end
TpermChanmat = zeros(HzBin,timebin,channels,perms);
for chan = 1:channels
    for mm = 1:perms
        TpermChanmat(:,:,chan,mm) = cell2mat(TpermChan(mm,chan));
    end
end
%% Monte Carlo Cluster Distribution
allClustSizePOSchan = cell(2,1);
allClustSizeNEGchan = cell(2,1);
for chan = 1:channels
        allClustSizePOSchan{chan} = cat(2,ClustSizePOSpermChan{:,chan});
        allClustSizeNEGchan{chan} = cat(2,ClustSizeNEGpermChan{:,chan});
end
%% Find signifncant clusters in real group
POSpvalChan = zeros(HzBin,timebin,channels);
NEGpvalChan = zeros(HzBin,timebin,channels);
for chan = 1:channels
    for cc = 1:timebin
        for rr = 1:HzBin
            POSpvalChan(rr,cc,chan) = sum(TpermChanmat(rr,cc,chan,:) > tChan(rr,cc,chan))/perms;
            NEGpvalChan(rr,cc,chan) = sum(TpermChanmat(rr,cc,chan,:) < tChan(rr,cc,chan))/perms;
        end
    end
end
sigPOSsizeChan = zeros(2,1);
sigNEGsizeChan = zeros(2,1);
for chan = 1:channels
    [BsigPOSchan{chan},LsigPOSchan{chan},NsigPOSchan{chan}] = bwboundaries(POSpvalChan(:,:,chan)<0.025);
    [BsigNEGchan{chan},LsigNEGchan{chan},NsigNEGchan{chan}] = bwboundaries(NEGpvalChan(:,:,chan)<0.025);
    sigPOSsizeChan(chan) = prctile(allClustSizePOSchan{chan},97.5);
    sigNEGsizeChan(chan) = prctile(allClustSizeNEGchan{chan},97.5);
end
for chan = 1:channels
    for ii = 1:max(max(LsigPOSchan{chan}))
        if sum(sum(LsigPOSchan{chan} == ii))>=sigPOSsizeChan(chan)
            BsigPOSchan{chan}{ii,2} = true;
        else
            BsigPOSchan{chan}{ii,2} = false;
        end
    end
end
for chan = 1:channels
    for ii = 1:max(max(LsigNEGchan{chan}))
        if sum(sum(LsigNEGchan{chan} == ii))>=sigNEGsizeChan(chan)
            BsigNEGchan{chan}{ii,2} = true;
        else
            BsigNEGchan{chan}{ii,2} = false;
        end
    end
end
%% Draw Figures
prompt = {'Control group name:','Experimental group name:','color bounds(+/-):','Time unit (0 = ms, 1 = sec)','X-axis start(ms):','X-axis step(ms):','Font Size:','Draw single groups? (1=yes):','color bounds single(max):','Seperate Windows? (1=yes)','Font Style:'};
dlg_title = 'Figure Settings';
num_lines = 1;
defaultans = {'Baseline','Light','1','1','0','5000','14','1','2','0','Arial'}; %change these for defaults, order in "prompt" 3 lines above
settingsUI = inputdlg(prompt,dlg_title,num_lines,defaultans);
ctrlname = settingsUI{1};
expname = settingsUI{2};
clrbound = str2double(settingsUI{3});
timeunit = str2double(settingsUI{4});
xstart = str2double(settingsUI{5});
xstep = str2double(settingsUI{6});
fntsize = str2double(settingsUI{7});
drawsingle = str2double(settingsUI{8});
clrbound2 = str2double(settingsUI{9});
sepplot = str2double(settingsUI{10});
fntstyle = settingsUI{11};

if timeunit == 0
    xlbl = 'Time (ms)';
else
    xlbl = 'Time (s)';
end
Chanfigs = group2means-group1means;
if drawsingle ~= 1
    for chan = 1:channels
        figure('Name',strcat('Channel',num2str(chan)));
        clims = [(-1*clrbound) clrbound];
        xms = (xstep/1000)/((1/samplingrate)*samplesbin);
        if xstart < 0
            x = [0,(-xstart/xstep)*xms:xms:timebin];
            if timeunit == 0
                xtimelbl = [xstart,0:xstep:((EndTimeUI-StartTimeUI)-xstart)];
            else
                xtimelbl = [(xstart/1000),0:(xstep/1000):(((EndTimeUI/1000)-(StartTimeUI/1000))-(xstart/1000))];
            end
        else
            x = 0:xms:timebin;
            if timeunit == 0
                xtimelbl = xstart:xstep:((EndTimeUI-StartTimeUI)+xstart);
            else
                xtimelbl = (xstart/1000):(xstep/1000):(((EndTimeUI/1000)-(StartTimeUI/1000))+(xstart/1000));
            end
        end
        y = flip(frq); %0:(HzBin/10):HzBin; % loaded y axis from another script
        imagesc([0 timebin],y,Chanfigs(:,:,chan), clims)
        set(gcf,'Position',[(20*chan) (450-(20*chan)) 600 400])
        %set(gcf,'InnerPosition',[~ ~ ~ 390])
        set(gca, 'YScale', 'log')
        set(gca,'FontName',fntstyle,'fontsize',fntsize,'YDir','Normal')
        ylabel('Frequency (Hz)')
        xlabel(xlbl)
        title(strcat(expname, ' -',{' '}, ctrlname))
        xticks(x)
        xticklabels(string(xtimelbl))
        colormap jet
        colorbar
            xline((10)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            xline((40)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((30)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((150)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            yline(4,'w','LineWidth', 5) % Delta
            yline(8,'w','LineWidth', 5) % Theta
            yline(13,'w','LineWidth', 5) % Alpha
            yline(30,'w','LineWidth', 5) % Beta
            yline(60,'w','LineWidth', 5) % Low Gamma
            yline(4,'k','LineWidth', 1) % Delta
            yline(8,'k','LineWidth', 1) % Theta
            yline(13,'k','LineWidth', 1) % Alpha
            yline(30,'k','LineWidth', 1) % Beta
            yline(60,'k','LineWidth', 1) % Low Gamma
            yticks([0.5 4 8 13 30 60 120])
        hold on
        %flipval = HzBin+1;
        for k = 1:size(BsigPOSchan{chan},1)
            if BsigPOSchan{chan}{k,2} == true
                boundary = BsigPOSchan{chan}{k,1};
                yy = [];
                for h = 1:size(boundary,1)
                    yy(h,1) = y(boundary(h,1),1);
                end
                plot(boundary(:,2), yy(:,1), 'w', 'LineWidth', 5, 'LineStyle','-')
                plot(boundary(:,2), yy(:,1), 'k', 'LineWidth', 2, 'LineStyle','-')
            end
        end
        for k = 1:size(BsigNEGchan{chan},1)
            if BsigNEGchan{chan}{k,2} == true
                boundary = BsigNEGchan{chan}{k,1};
                yy = [];
                for h = 1:size(boundary,1)
                    yy(h,1) = y(boundary(h,1),1);
                end
                plot(boundary(:,2), yy(:,1), 'w', 'LineWidth', 4, 'LineStyle','-')
                plot(boundary(:,2), yy(:,1), 'k', 'LineWidth', 3, 'LineStyle',':')
            end
        end
    end
else
   for chan = 1:channels
        if sepplot ~= 1
            figure('Name',strcat('Channel',num2str(chan)));
            set(gcf,'Position',[(20*chan) (450-(20*chan)) 1800 400])
        else
            figure('Name',strcat('Channel',num2str(chan),'_',expname));
        end
        %set(gcf,'InnerPosition',[~ ~ ~ 390])
        clims = [(-1*clrbound) clrbound];
        clims2 = [0 clrbound2];
        xms = (xstep/1000)/((1/samplingrate)*samplesbin);
        if xstart < 0
            x = [0,(-xstart/xstep)*xms:xms:timebin];
            if timeunit == 0
                xtimelbl = [xstart,0:xstep:((EndTimeUI-StartTimeUI)-xstart)];
            else
                xtimelbl = [(xstart/1000),0:(xstep/1000):(((EndTimeUI/1000)-(StartTimeUI/1000))-(xstart/1000))];
            end
        else
            x = 0:xms:timebin;
            if timeunit == 0
                xtimelbl = xstart:xstep:((EndTimeUI-StartTimeUI)+xstart);
            else
                xtimelbl = (xstart/1000):(xstep/1000):(((EndTimeUI/1000)-(StartTimeUI/1000))+(xstart/1000));
            end
        end
        y = flip(frq); %0:(HzBin/10):HzBin; % loaded y axis from another script
        if sepplot ~= 1
            subplot(1,3,2);
        end
        imagesc([0 timebin],y,group2means(:,:,chan), clims2)
        set(gca,'FontName',fntstyle,'fontsize',fntsize,'YDir','Normal')
        set(gca, 'YScale', 'log')
        ylabel('Frequency (Hz)')
        xlabel(xlbl)
        title(expname)
        xticks(x)
        xticklabels(string(xtimelbl))
        colormap jet
        colorbar
            xline((10)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            xline((40)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((30)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((150)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            yline(4,'w','LineWidth', 5) % Delta
            yline(8,'w','LineWidth', 5) % Theta
            yline(13,'w','LineWidth', 5) % Alpha
            yline(30,'w','LineWidth', 5) % Beta
            yline(60,'w','LineWidth', 5) % Low Gamma
            yline(4,'k','LineWidth', 1) % Delta
            yline(8,'k','LineWidth', 1) % Theta
            yline(13,'k','LineWidth', 1) % Alpha
            yline(30,'k','LineWidth', 1) % Beta
            yline(60,'k','LineWidth', 1) % Low Gamma
                %yline(100,'k','LineWidth', 3) % 100Hz
            yticks([0.5 4 8 13 30 60 120])
        
        if sepplot ~= 1
            subplot(1,3,1);
        else
            figure('Name',strcat('Channel',num2str(chan),'_',ctrlname));
        end
        imagesc([0 timebin],y,group1means(:,:,chan), clims2)
        set(gca,'FontName',fntstyle,'fontsize',fntsize,'YDir','Normal')
        set(gca, 'YScale', 'log')
        ylabel('Frequency (Hz)')
        xlabel(xlbl)
        title(ctrlname)
        xticks(x)
        xticklabels(string(xtimelbl))
        colormap jet
        colorbar
            xline((10)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            xline((40)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((30)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%             xline((150)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            yline(4,'w','LineWidth', 5) % Delta
            yline(8,'w','LineWidth', 5) % Theta
            yline(13,'w','LineWidth', 5) % Alpha
            yline(30,'w','LineWidth', 5) % Beta
            yline(60,'w','LineWidth', 5) % Low Gamma
            yline(4,'k','LineWidth', 1) % Delta
            yline(8,'k','LineWidth', 1) % Theta
            yline(13,'k','LineWidth', 1) % Alpha
            yline(30,'k','LineWidth', 1) % Beta
            yline(60,'k','LineWidth', 1) % Low Gamma
            yticks([0.5 4 8 13 30 60 120])
        if sepplot ~= 1
            subplot(1,3,3);
        else
            figure('Name',strcat('Channel',num2str(chan),'_',strcat(expname, '_-','_', ctrlname)));
        end
        imagesc([0 timebin],y,Chanfigs(:,:,chan), clims)
        set(gca,'FontName',fntstyle,'fontsize',fntsize,'YDir','Normal')
        set(gca, 'YScale', 'log')
        ylabel('Frequency (Hz)')
        xlabel(xlbl)
        title(strcat(expname, ' -',{' '}, ctrlname))
        xticks(x)
        xticklabels(string(xtimelbl))
        colormap jet
        colorbar
            xline((10)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            xline((40)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%              xline((30)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
%              xline((150)/((1/samplingrate)*samplesbin),'k','LineWidth', 4)
            yline(4,'w','LineWidth', 5) % Delta
            yline(8,'w','LineWidth', 5) % Theta
            yline(13,'w','LineWidth', 5) % Alpha
            yline(30,'w','LineWidth', 5) % Beta
            yline(60,'w','LineWidth', 5) % Low Gamma
            yline(4,'k','LineWidth', 1) % Delta
            yline(8,'k','LineWidth', 1) % Theta
            yline(13,'k','LineWidth', 1) % Alpha
            yline(30,'k','LineWidth', 1) % Beta
            yline(60,'k','LineWidth', 1) % Low Gamma
            yticks([0.5 4 8 13 30 60 120])
        hold on
        %flipval = HzBin+1;
        for k = 1:size(BsigPOSchan{chan},1)
            if BsigPOSchan{chan}{k,2} == true
                boundary = BsigPOSchan{chan}{k,1};
                yy = [];
                for h = 1:size(boundary,1)
                    yy(h,1) = y(boundary(h,1),1);
                end
                plot(boundary(:,2), yy(:,1), 'w', 'LineWidth', 5, 'LineStyle','-')
                plot(boundary(:,2), yy(:,1), 'k', 'LineWidth', 2, 'LineStyle','-')
            end
        end
        for k = 1:size(BsigNEGchan{chan},1)
            if BsigNEGchan{chan}{k,2} == true
                boundary = BsigNEGchan{chan}{k,1};
                yy = [];
                for h = 1:size(boundary,1)
                    yy(h,1) = y(boundary(h,1),1);
                end
                plot(boundary(:,2), yy(:,1), 'w', 'LineWidth', 5, 'LineStyle','-')
                plot(boundary(:,2), yy(:,1), 'k', 'LineWidth', 2, 'LineStyle',':')
            end
        end
    end 
end
%% Monte Carlo Function for use with parfor
function [Tperm, ClustSizesPOSmc, ClustSizesNEGmc] = mymontecarlo(bindata,groupsize1,groupsize2)
    randorder = randperm(size(bindata,4));
    bothgroupsize = groupsize1+groupsize2; 
    rowsN = size(bindata,1);
    columnN = size(bindata,2);
    randGroup1 = zeros(rowsN,columnN,groupsize1);
    randGroup2 = zeros(rowsN,columnN,groupsize2);   
    for m = 1:bothgroupsize
        if m <= groupsize1
            randGroup1(:,:,m) = bindata(:,:,randorder(1,m));
        elseif m > groupsize1
            randGroup2(:,:,m-groupsize1) = bindata(:,:,randorder(1,m));
        else
            disp('file has no group id!!!')
            stop
        end
    end
    Tperm = buncha_ttest2s(randGroup1,randGroup2);

    dF = bothgroupsize-2;
    TsigPOSmc = tcdf(Tperm,dF)>=0.975;
    TsigNEGmc = tcdf(Tperm,dF)<=0.025;
    [BsigPOSmc,LsigPOSmc] = bwboundaries(TsigPOSmc);
    [BsigNEGmc,LsigNEGmc] = bwboundaries(TsigNEGmc);
    ClustSizesPOSmc = zeros(1,size(BsigPOSmc,1));
    for ii = 1:size(BsigPOSmc,1)
        ClustSizesPOSmc(ii) = sum(sum(LsigPOSmc == ii));
    end
    ClustSizesNEGmc = zeros(1,size(BsigNEGmc,1));
    for ii = 1:size(BsigNEGmc,1)
        ClustSizesNEGmc(ii) = sum(sum(LsigNEGmc == ii));
    end
end
%% t-test func
function [tstats] = buncha_ttest2s(g1,g2)

tstats = (mean(g1,3)-mean(g2,3))./sqrt( ( (std(g1,0,3).^2)./size(g1,3)) + ( (std(g2,0,3).^2)./size(g2,3)) );
end
%% custom divisors 
 function d = divisors( n )
%DIVISORS(N) Returns array of divisors of n
if ~isscalar(n)
    error('n must be a scalar');
end
if n < 1
    error('n must be positive integer');
end
if n == 1
    d = 1;
    return;
end
f = factor(n);
pf = unique(f);
for i = 1:length(pf)
    m(i) = sum(f == pf(i));
end
mi = zeros(size(m));
d = zeros(1,prod(m+1));
i = 1;
carry = 0;
while ~carry
    d(i) = prod(pf.^mi);
    i = i + 1;
    if mi(1) < m(1)
        mi(1) = mi(1) + 1;
    else
        mi(1) = 0;
        carry = 1;
    end
    for j = 2:length(m)
        if carry
            if mi(j) < m(j)
                mi(j) = mi(j) + 1;
                carry = 0;
            else
                mi(j) = 0;
                carry = 1;
            end
        end
    end
end
d = sort(d);
end