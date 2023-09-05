%% Data settings Dialog Boxes 1
prompt = {'Enter # of Channels:','Enter Sampling Rate (Hz):','Enter # of Wavelets:','Enter Start time (ms)','Enter End time (ms)','Enter # of Monte Carlo Permutations (at least 1000 recomended)','# of rows to skip'};
dlg_title = 'Settings';
num_lines = 1;
defaultans = {'1','2500','80','0','50000','1000','0'}; % change these values to alter default input!!!
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

%% calculate mean lowHz and mean hiHz, find syncope thresh (gamma 50% drop, 80% rebound)
for i = 1:nfiles
    t = 0:length(data{1,i})-1;
    tsec = t/samplingrate;
    lowHz{i} = mean(data{1,i}(47:70,:),1); % 47:70 = 5Hz-1Hz
    hiHz{i} = mean(data{1,i}(4:40,:),1); % 4:40 = 100Hz-8Hz
    hiHz2{i} = mean(data{1,i}(4:40,:),1); % 11:21 = 60Hz-30Hz
    lowHzS{i} = movmean(lowHz{i},6250,2); %3125
    hiHzS{i} = movmean(hiHz{i},3125,2); %3125
    hiHz2S{i} = movmean(hiHz2{i},3125,2); %3125
    figure;
    plot(tsec,hiHzS{i});
    hold;
    %plot(tsec,lowHzS{i});
    %plot(tsec,hiHz2S{i});
    xline(10,'k','LineWidth', 2) %edit Dr Dan
    xline(40,'k','LineWidth', 2)
    yline(1,'k','LineWidth', 2,'LineStyle',':')
    yline(0.5,'r','LineWidth', 2)
    yline(0.80,'b','LineWidth', 2)
    ylim([0 2])
    xlim([0 60])

        %xline(16.368+10,'y','LineWidth', 2) %eyeroll

    b = true;
    loop = 0;

    while b
        loop = loop+1;
        if loop == 1
            tempthresh = find(hiHzS{i}(1,1:end)<0.5);
                if size(tempthresh,2) == 0
                    b = false;
                    thresh(loop,i) = 0;
                    thresh2(loop,i) = 0;
                    break
                end
            thresh(loop,i) = (tempthresh(1)/samplingrate)-10;
            threshsamp(loop,i) = tempthresh(1);
            xline(thresh(loop,i)+10,'r','LineWidth', 2)
        else
            tempthresh = find(hiHzS{i}(1,threshsamp2(loop-1,i):end)<0.5);
                if size(tempthresh,2) == 0
                    b = false;
                    thresh(loop,i) = 0;
                    thresh2(loop,i) = 0;
                    break
                end
            thresh(loop,i) = ((threshsamp2(loop-1,i)+tempthresh(1))/samplingrate)-10;
            threshsamp(loop,i) = threshsamp2(loop-1,i)+tempthresh(1);
            xline(thresh(loop,i)+10,'r','LineWidth', 2)
        end    
        tempthresh2 = find(hiHzS{i}(1,threshsamp(loop,i):end)>0.80);
        if size(tempthresh2,2) == 0
            b = false;
            break
        end
        thresh2(loop,i) = ((threshsamp(loop,i)+tempthresh2(1))/samplingrate)-10;
        threshsamp2(loop,i) = threshsamp(loop,i)+tempthresh2(1);
        xline(thresh2(loop,i)+10,'b','LineWidth', 2)
    end
end

%% plot power with syncope thresh
for i = 1:nfiles
    figure
    t = 0:length(data{1,i})-1;
    %tsec = t/samplingrate;
    %t = downsample(t,250);
    tsec = t/samplingrate;
    %datap = downsample(data{1,i}',250);
    %datap = datap';
    %pcolor(tsec,frq,datap);    
    imagesc(tsec,frq,data{1,i});
    %imagesc(tsec,frq,datap);
    %pcolor(tsec,frq,data{1,i});
    shading interp
    set(gca, 'YScale', 'log','Ydir','normal') %comment in for log scale visual
    xline(10,'k','LineWidth', 4) %edit Dr Dan
    xline(40,'k','LineWidth', 4)
    if exist('thresh','var') == 1
        for ii = 1:nnz(thresh(1:size(thresh,1),i))
            xline(thresh(ii,i)+10,'r','LineWidth', 5)
            xline(thresh2(ii,i)+10,'w','LineWidth', 5)
            xline(thresh(ii,i)+10,'k','LineWidth', 2, 'LineStyle',':')
            xline(thresh2(ii,i)+10,'k','LineWidth', 2,'LineStyle',':')
        end
    end
    yline(4,'w','LineWidth', 2) % Delta
    yline(8,'w','LineWidth', 2) % Theta
    yline(13,'w','LineWidth', 2) % Alpha
    yline(30,'w','LineWidth', 2) % Beta
    yline(60,'w','LineWidth', 2) % Low Gamma
    yticks([0.5 4 8 13 30 60 120])
    colormap(turbo)
    colorbar
    caxis([0 3])
    name = matfiles(i).name;
    name = name(1:end-4);
    print(gcf, '-dtiff', strcat(name,'.tiff')); %creates .tiff image
end
avgLat = mean(thresh(1,1:end));
% semLat = avgLat/sqrt(nfiles-8);
% minlat = min(thresh(1,15:end));
% maxlat = max(thresh(1,15:end));

%% plot power based on syncope latency
for i = 1:nfiles
    tstart =  (threshsamp(1,i))-(15*samplingrate);
    tend =  (threshsamp(1,i))+(30*samplingrate);
    datathresh = data{1,i}(:,tstart:tend);
    datathreshM(i,:,:) = datathresh(:,:);
    figure
    t = 0:length(datathresh(1,:))-1;
    tsec = t/samplingrate;
    pcolor(tsec,frq,datathresh);
    shading interp
    set(gca, 'YScale', 'log') %comment in for log scale visual
    xline(15-thresh(1,i),'k','LineWidth', 4)
    xline(15,'k','LineWidth', 3,'LineStyle',':')
    xline((15-thresh(1,i))+30,'k','LineWidth', 4)
    yline(4,'w','LineWidth', 2) % Delta
    yline(8,'w','LineWidth', 2) % Theta
    yline(13,'w','LineWidth', 2) % Alpha
    yline(30,'w','LineWidth', 2) % Beta
    yline(60,'w','LineWidth', 2) % Low Gamma
    yticks([0.5 4 8 13 30 60 120])
    colormap(turbo)
    colorbar
    caxis([0 3])
   name = matfiles(i).name;
   name = name(1:end-4);
   print(gcf, '-dtiff', strcat(name,'.tiff'));
end

%% calculate mean difference
threshplot = squeeze(mean(datathreshM(1:end,:,:),1));
pcolor(tsec,frq,threshplot);
    shading interp
    set(gca, 'YScale', 'log')
    xline(15,'k','LineWidth', 3,'LineStyle',':')
    yline(4,'w','LineWidth', 2) % Delta
    yline(8,'w','LineWidth', 2) % Theta
    yline(13,'w','LineWidth', 2) % Alpha
    yline(30,'w','LineWidth', 2) % Beta
    yline(60,'w','LineWidth', 2) % Low Gamma
    yticks([0.5 4 8 13 30 60 120])
    colormap(turbo)
    colorbar
    caxis([0 3])

%% save file
 savename = '20Hz_syncope_bouts_Temporal_CNO_Veh.mat';
 save(savename,'thresh','thresh2','matfiles','-v7.3') % save only required info