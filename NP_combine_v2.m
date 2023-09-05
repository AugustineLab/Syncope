%% make file list from exported .mat files in single directory
testfiledir = uigetdir();
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);
spikes_save_all = [];

%% loops for spikes_save_all
for w = 1 : nfiles
    fid = fullfile(testfiledir, matfiles(w).name);
    load(fid)
    APsr = sp.sample_rate; % samping rate (HZ) for AP channels
    NIsr = str2double(NImeta.niSampRate); %sampling rate for NI channels
    ProbeEntireDur = sp.st(end); % length of file in Secs. (Time of last recorded spike)
    tt = 0:1/APsr:ProbeEntireDur; % set the time scale (x-axis) for spikes

    trainIndex = 2; % 4th train = 20Hz 30sec, ect

    if w == 21 || w ==22
        trainIndex = 1; % exceptions to fix lack of "no light" recordings
    end
    
    % get the window idx from each stimulation
    nwithinbin = 3000; % set samples per bin 100ms resolution currently
    WsecS = -30; % window in seconds before stim start
    WsecE = 60; % window in seconds after stim ends
    idx_stimS = zeros(1,size(TrainInfo,2));
    idx_stimE = zeros(1,size(TrainInfo,2));
    idx_bsS = zeros(1,size(TrainInfo,2));
    idx_asE = zeros(1,size(TrainInfo,2));
    for n = 1:size(TrainInfo,2)
        idx_stimS(n) = ((round(TrainInfo(3,n)/NIsr))*APsr)/nwithinbin; % index of each stim start time
        idx_stimE(n) = ((round(TrainInfo(4,n)/NIsr))*APsr)/nwithinbin; % index of each stim end time
        idx_bsS(n) = idx_stimS(n)+(WsecS*APsr)/nwithinbin;% index before each stim (window start)
        idx_asE(n) = idx_stimE(n)+(WsecE*APsr)/nwithinbin;% index after each stim end (window end)
    end

    % Spikes post-processing, ifr calculation (only 1 train)
    for n = 0:(size(sp.cids,2)-1)
        m = find(sp.clu == n);
        timeofspike = unique(sp.st(m)); % (Secs)
        cluster_group = sp.cgs(n+1);
        try
            [ifr] = instantfr(timeofspike,tt);% Instantaneous firing rate of a spike train
        catch
            warning(['Sample points must be unique and sorted in ascending order. with Cluter id' num2str(n)]);
            ifr=0;
        end 
        nbin = fix(length(ifr)/nwithinbin);
        ifr_bins = [];
        for k = 1:nbin
            ifr_bin = mean(ifr(nwithinbin*(k-1)+1):ifr(nwithinbin*k));
            ifr_bin(isnan(ifr_bin)) = 0; % TEST
            ifr_bins = [ifr_bins ifr_bin];
        end
        spikes_save{n+1,1} = 0;% placeholder for global rowID    
        spikes_save{n+1,2} = n;% local cluster id
        spikes_save{n+1,3} = cluster_group; %  1: mua, 2:good
        spikes_save{n+1,4} = ifr_bins; %temp calculation, replaced by windowed Z-score
        spikes_save{n+1,5} = cell2mat(cluster_region(2,n+1)); %allen atlas ID
        spikes_save{n+1,6} = string(cluster_region(3,n+1)); %allen atlas Abrev
        spikes_save{n+1,7} = matfiles(w).name; %file name per recordings
    end
    
    % Z-score(BL and SD)
    id = trainIndex;
    for n = 1:size(sp.cids,2) % calculate Z-score
        current_spikes = spikes_save{n,4};
        blmean = mean(current_spikes(1,idx_bsS(id):idx_stimS(id)));
        blsdev = std(current_spikes(1,idx_bsS(id):idx_stimS(id)));
        spikes_save{n,4} = ((current_spikes(1,idx_bsS(id):idx_asE(id)-1))-blmean)/blsdev; %overwrite ifr_bins with BL Z-score
        spikes_save{n,8} = blmean; % mean spiking rate before laser stimulation for thresholding purposes
    end
    spikes_save_all = [spikes_save_all; spikes_save];
    clearvars -except w spikes_save_all nfiles matfiles testfiledir
end
for k = 1:size(spikes_save_all,1)
        spikes_save_all{k,1} = k; % make global row IDs
end

%% make single plotable matrix out of "spikes_save_all"
BigMatrix = zeros(size(spikes_save_all,1),120);
for n = 1:size(spikes_save_all,1)
    leng = size(spikes_save_all{n,4},2);
    BigMatrix(n,1:leng) = spikes_save_all{n,4}; %make raw data array for use with filters
end
load('RWB_colormap.mat') %load custom colormap settings saved to disk

%% Sorting Figures, Z-score(BL) mean change, select unit, no regions
% sort by mean change during stim compared to pre-stim window, threshold active units.

unitype = 3; % EDIT 1 = mua, 2 = good, 3 = all
smoothfactor = 1; % EDIT to change smooth samples
smoothdim = 2;% EDIT to change smooth dim (1=cluster unit,2=time)
clims = ([-3 3]); % EDIT to change color scale
thresh = 0.1; % EDIT to change min firing rate during Basline

if unitype == 3
    unifilt = ones(size(spikes_save_all,1),1);
else
    unifilt = cell2mat(spikes_save_all(:,3))==unitype;
end
figure;
threfilt = cell2mat(spikes_save_all(:,8))>thresh; %threshold baseline spikerate
%create baseline Z-scores
finalfilt = and(threfilt,unifilt);
filt_plot = BigMatrix(finalfilt,:);
mean_change = [];
    for n = 1:size(filt_plot,1) % calculate mean-change
            meanstimZ = mean(filt_plot(n,301:600));
            mean_change = [mean_change; meanstimZ];
    end
[B,I] = sort(mean_change);
smooth_plot = movmean(filt_plot,smoothfactor,smoothdim);
imagesc(smooth_plot(I,:),clims);
xline(301,'k','LineWidth', 2) % draw line for stim start
xline(600,'k','LineWidth', 2) % draw line for stim end
colorbar
colormap(Red_White_Blue)
ylabel('Unit')

%% grand average across all units, top-half vs bottom half
sort_plot = smooth_plot(I,:);
half = round(size(sort_plot,1)/2);
BigMean1 = mean(sort_plot(1:half,:),1);
BigMean2 = mean(sort_plot(half+1:end,:),1);
figure
plot(BigMean1)
ylim([-0.75 5]); 
hold;
plot(BigMean2)
xline(301,'k','LineWidth', 2) % draw line for stim start
xline(600,'k','LineWidth', 2) % draw line for stim end

%% Make Region classes/categories based on allen atlas
mainregions = {'ACA','AI','BLA','BMA','BST','DORpm','DORsm','HPC','ILA','LA','LS','LZ','MBmot','MBsta','MEZ','MOp','MOs','OLF','ORB','p-mot','PA','PALd','PALm','PALv','PIR','PL','PVR','PVZ','RSP','SC','SSp','SSs','STRd','STRv','sAMY','VISC','PTLp','VS','Fiber','error'};
region={};
region{1} = {'ACAd1','ACAd2/3', 'ACAd5', 'ADAd6a','ACAv1','ACAv2/3','ACAv5','ACAv6a','ACAv6b', 'ACAd6a'}; %ACA
region{2} = {'AId6a','AIp6a', 'AIp6b', 'AIv6a','CLA','GU6a'}; %AI
region{3} = {'BLAa','BLAp', 'CTXsp'}; %BLA
region{4} = {'BMAa','BMAp'}; %BMA
region{5} = {'BST'}; %BST
region{6} = {'AD','AMd','AMv','AV','CL','CM','Eth','IAM','IMD','LD','LH','LP','MD','MH','PCN','PF','PO','PT','PVT','RE','RH','RT','SMT','SubG','TH','Xi','IAD','SGN'}; %DORpm (add IAD,SGN)
region{7} = {'LGd-co','LGd-ip','LGd-sh','LGv','VAL','VM','VPL','VPLpc','VPM','VPMpc','PoT','SPFp'}; %DORsm (add PoT,SPFp)
region{8} = {'CA1','CA2','CA3','DG-mo','DG-po','DG-sg','HPF','ProS','SUB','IG'}; %HPC (add IG)
region{9} = {'ILA1','ILA2/3','ILA5','ILA6a','ILA6b'}; %ILA
region{10} = {'LA'}; %LA
region{11} = {'LSc','LSr','LSv','SF'}; %LS
region{12} = {'FF','HY','LHA','LPO','PSTN','ZI'}; %LZ
region{13} = {'APN','MB','MRN','MT','NOT','PAG','PPT','PR','RN','RR','SNr','VTA'}; %MBmot
region{14} = {'PPN','SNc'}; %MBsta
region{15} = {'MPN','PH','PMv','PMd','TMv','TU','VMH','AHN'}; %MEZ (add AHN)
region{16} = {'MOp2/3','MOp5','MOp6a','MOp6b'}; %MOp
region{17} = {'MOs1','MOs2/3','MOs5','MOp6a','MOp6b'}; %MOs
region{18} = {'AON','COApl','COApm','DP','EPd','EPv','NLOT3','OLF','PAA','TTd','TR'}; %OLF (add TR)
region{19} = {'FRP5','FRP6a','ORBl1','ORB12/3','ORBl5','ORBl6a','ORBm1','ORBm2/3','ORBm5','ORBvl6a'}; %ORB
region{20} = {'PG','PRNc','PRNr','TRN'}; %p-mot
region{21} = {'PA'}; %PA
region{22} = {'GPe','GPi','PAL'}; %PALd
region{23} = {'MS','NDB','TRS'}; %PALm
region{24} = {'Sl'}; %PALv
region{25} = {'PIR'}; %PIR
region{26} = {'PL1','PL2/3','PL5','PL6a','PL6b'}; %PL
region{27} = {'DMH','SBPV'}; %PVR
region{28} = {'PVH','PVHd','PVi'}; %PVZ
region{29} = {'RSPv2/3','RSPv5','RSPv6a'}; %RSP
region{30} = {'SCdg','SCdw','SCig','SCiw','SCop'}; %SC
region{31} = {'SSp-n4','SSp-n5','SSp-n6a','SSp-ul6a','SSP-ul6b','SSp-bfd6a','SSp-bfd6b','SSp-ll2/3','SSp-ll4','SSp-ll5','SSp-ll6a','SSp-ll6b','SSp-n1','SSp-n2/3','SSp-tr1','SSp-ul1','SSp-ul2/3','SSp-ul4','SSp-ul5','SSp-ul6b','SSp-un2/3','SSp-un4','SSp-un5','SSp-un6a','SSp-un6b'}; %SSp (add SSp-bfd6a->)
region{32} = {'SSs5','SSs6a','SSs6b'}; %SSs
region{33} = {'CP','FS','STR'}; %STRd
region{34} = {'ACB','OT'}; %STRv
region{35} = {'AAA','CEac','CEAl','CEAm','IA','MEA'}; %sAMY
region{36} = {'VISC6a','VISC6b'}; %VISC
region{37} = {'VISa1','VISa2/3','VISa4','VISa5','VISa6a','VISa6b'}; %PTLp

region{38} = {'SEZ','V3','VL'}; %VS- ventricular system
region{39} = {'aco','act','alv','amc','bsc','ccb','ccg','ccs','cing','cpd','ec','ee','fa','fi','fiber tracts','em','ml','or','vhc','scwm'}; %Fiber tracts
region{40} = {'outside_brain','root'}; % channel not in brain


% assign sub-regions to main regions
for n = 1:size(region,2)
    Rindex = find(ismember(string(spikes_save_all(:,6)),region{n}));
    spikes_save_all(Rindex,9) = mainregions(n);
end

% find unassigned regions (used for troubleshooting missing subregions)
needasignID = find(cellfun(@isempty,spikes_save_all(:,9)));
needasign = unique(string(spikes_save_all(needasignID,6))); % subregions that do not have a main region assignment

%% Sorting Figures, Z-score(BL) mean change, region, select unit
% sort by mean change during stim compared to pre-stim window, threshold
% active units, main regionwise

unitype = 3; % 1 = mua, 2 = good, 3 = all
smoothfactor = 5; % EDIT to change smooth samples
smoothdim = 2;% EDIT to change smooth dim (1=unit,2=time)
clims = ([-3 3]); % EDIT to change color scale
thresh = 0.5; % EDIT to change min firing rate during Basline
minmax = 0; %0 = mean change, 1 = min lat, 2 = max lat (sorting method within each main region)

if unitype == 3
    unifilt = ones(size(spikes_save_all,1),1);
else
    unifilt = cell2mat(spikes_save_all(:,3))==unitype;
end

threfilt = cell2mat(spikes_save_all(:,8))>thresh; %thresh baseline spikerate
%create baseline Z-scores
mean_change = [];
for n = 1:size(BigMatrix,1) % calculate mean-change
        meanstimZ = mean(BigMatrix(n,301:600)); % edit to change "mean change" sorting window
        mean_change = [mean_change; meanstimZ];
end
finalfilt = and(threfilt,unifilt);
[activeregions,~,ic] = unique(string(spikes_save_all(finalfilt,9)),'stable'); %row identifiers of probeinfo regions
unit_counts = accumarray(ic,1);
value_counts = [activeregions, unit_counts];
deleteidx = [];
for kk = 1:size(activeregions,1) % exclude regions with fewer than X units, remove specigic main regions
    if str2double(value_counts(kk,2)) < 20 || strcmp(value_counts(kk,1),'Fiber')||strcmp(value_counts(kk,1),'error')||strcmp(value_counts(kk,1),'VS')
        deleteidx = [deleteidx; kk];
    end
end
activeregions([deleteidx],:) = [];
regionmean = zeros(1,size(activeregions,1));
for jj = 1:size(activeregions,1)
    regionIDX = string(spikes_save_all(:,9)) == activeregions(jj);
    regionIDX2 = and(regionIDX, finalfilt);
    regionmean(jj) = mean(BigMatrix(regionIDX2,360:430),'all'); % area to calculate mainregion-wise sorting
end
[~,RI] = sort(regionmean,'ascend'); % direction of sort: 'descend' or 'ascend'
for jj = 1:size(activeregions,1)
        regionIDX = string(spikes_save_all(:,9)) == activeregions(RI(jj));
        regionIDX2 = and(regionIDX, finalfilt);
        if minmax == 0
            [~,I] = sort(mean_change(regionIDX2));
        elseif minmax == 1
            [~,lat] = min(BigMatrix(regionIDX2,301:600),[],2);
            [~,I] = sort(lat);
        else
            [~,lat] = max(BigMatrix(regionIDX2,301:600),[],2);
            [~,I] = sort(lat);
        end
        regionmatrix = movmean(BigMatrix(regionIDX2,:),smoothfactor,smoothdim);
        if jj == 1
            y(jj) = 1;
            finalplot = regionmatrix(I,:);
            ytick(jj) = round(size(regionmatrix(I,:),1)/2);
        else
            y(jj) = size(finalplot, 1);
            ytick(jj) = round((size(regionmatrix(I,:),1)/2)+size(finalplot,1));
            finalplot = cat(1,finalplot, regionmatrix(I,:));
        end
end
figure;
imagesc(finalplot,clims)
for ii = 1:size(y,2)
    yline(y(ii),'k', 'LineWidth', 1)
end
yticks(ytick);
yticklabels(activeregions(RI));
xline(301,'k','LineWidth', 2)
xline(600,'k','LineWidth', 2)
colormap(Red_White_Blue)
colorbar
clear y
clear ytick
