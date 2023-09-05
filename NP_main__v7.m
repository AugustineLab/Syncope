%% Import Data
% The script will ask for a Folder where the AP and LFP bin files are
% located, and then do the rest
load('RWB_colormap.mat') %load custom colormap settings saved to disk

% get Imec Probe Data after Kilosort
disp('find the Folder where the .ap.bin folder is for the probe')
ProbePath = uigetdir();
sp = loadKSdir(ProbePath);% Steinmetz function to load alot of basic data info into "sp" datastructure
APsr = sp.sample_rate; % samping rate (HZ) for AP channels
ProbeEntireDur = sp.st(end); % length of file in Secs. (Time of last recorded spike)
tt = 0:1/APsr:ProbeEntireDur; % set the time scale (x-axis) for spikes
[~, ~, sp.templateDepthsT] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps); %get depth from templates
sp.templateDepthsT = sp.templateDepthsT.';

% Import LFP channels info
lfpD = dir(fullfile(ProbePath, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(ProbePath, lfpD.name);
LFPmeta = ReadMeta(lfpD.name, ProbePath);
LFPsr = str2double(LFPmeta.imSampRate);
LFPsamples = str2double(LFPmeta.fileTimeSecs)*LFPsr;

% Find NI binary file based on previous path
idcs   = strfind(ProbePath,'\'); % find all folder markers
NIpath = ProbePath(1:idcs(end)-1); % go back 1 folder from Probe path
NIbinInfo = dir(fullfile(NIpath, '*.bin')); % get file name of .bin in new folder (NI .bin file)
clear idcs

% Parse the corresponding NI metafile
NImeta = ReadMeta(NIbinInfo.name, NIpath);
NIsr = str2double(NImeta.niSampRate); %sampling rate for NI channels

% Import NI data
NInSamp = floor(str2double(NImeta.fileTimeSecs)*NIsr); %Total NI samples
NIdataArray = ReadBin(0, NInSamp, NImeta, NIbinInfo.name, NIpath); % dimension: [nChan,nSamp] (1-8:analog channel; 9:digital channel)
clear NIbinName
clear NIpath

% Gain Correct analog channels (1-based for MATLAB).
for ch = 1:str2double(NImeta.nSavedChans)-1
    if strcmp(NImeta.typeThis, 'imec')
        NIdataArray = GainCorrectIM(NIdataArray, [ch], NImeta);
    else
        NIdataArray = GainCorrectNI(NIdataArray, [ch], NImeta);
    end
end
clear ch

% Import digital channels from digital word in channel 9
dw=1;
dLineList = [0,1,2]; %Read these lines in dw (0-based).
digArray = ExtractDigital(NIdataArray, NImeta, dw, dLineList);
clear dw
clear dLineList
clear ch

%% Get SHARPtrack depth boundries from file and assign region index to clusters (not required to generate plots)
disp('find depth_table.mat in SHARPtrack "processed" folder')
[DepthFile, DepthPath] = uigetfile(); %find depth_table.mat in SHARPtrack "processed" folder after running custom JL script
load(strcat(DepthPath,DepthFile));
prompt = {'Enter Probe #'};
dlgtitle = 'Depth file Loaded';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
probenum = str2num(answer{1});
tip_regionconvert = (cell2mat(ProbeInfo{probenum}(:,1:2))-max(cell2mat(ProbeInfo{probenum}(:,2))))*-1; % depth calculation from tip in SHARPtrack
cluster_region = cell(4,size(sp.cids,2));
for j = 1:size(tip_regionconvert,1)
    for k = 1:size(sp.cids,2)
        if (sp.templateDepthsT(1,k)>=tip_regionconvert(j,2) && (sp.templateDepthsT(1,k)<tip_regionconvert(j,1)))
            cluster_region{1,k} = j;
            cluster_region{2,k} = cell2mat(ProbeInfo{probenum}(j,5));
            cluster_region{3,k} = ProbeInfo{probenum}(j,3);
            cluster_region{4,k} = ProbeInfo{probenum}(j,4);
        else
            if sp.templateDepthsT(1,k)>tip_regionconvert(1,1)
            cluster_region{1,k} = size(tip_regionconvert,1)+1;
            cluster_region{2,k} = size(tip_regionconvert,1)+1;
            cluster_region{3,k} = 'outside_brain';
            cluster_region{4,k} = 'outside_brain';
            end
        end
    end
end
%% Find TTL stim trains and get time stamps (required)
lowestHzStim = 4; % lowest stim rate in Hz
TrainChan = 3; %Stim channel # in digitalArray
TrainInfo = find_stim_train(digArray(TrainChan,:),NIsr,lowestHzStim); % Get train timestamps (custom function)
clear lowestHzStim

%% Spikes post-processing, ifr calculation (required for spikes) (OPTIMIZE, Use par for, need to optimize still)
for n = 0:(size(sp.cids,2)-1)
    m = find(sp.clu == n);
    numspikes = length(m); % number of spikes in each cluster(a single neuron)
    timeofspike = unique(sp.st(m)); % (Secs)
    cluster_group = sp.cgs(n+1);
    try
        [ifr] = instantfr(timeofspike,tt);% Instantaneous firing rate of a spike train
    catch
        warning(['Sample points must be unique and sorted in ascending order. with Cluter id' num2str(n)]);
        ifr=0;
    end 
    % save
    nwithinbin = 30000; % set samples per bin
    nbin = fix(length(ifr)/nwithinbin);
    ifr_bins = [];
    for k = 1:nbin
        ifr_bin = mean(ifr(nwithinbin*(k-1)+1):ifr(nwithinbin*k));
        ifr_bin(isnan(ifr_bin)) = 0; % TEST
        ifr_bins = [ifr_bins ifr_bin];
    end
    spikes_save{n+1,1} = n;% cluster id
    spikes_save{n+1,2} = numspikes;% number of spikes in each cluster(a single neuron)
    spikes_save{n+1,3} = sp.templateDepthsT(n+1);
    spikes_save{n+1,4} = cluster_group; %  1: mua, 2:good
    spikes_save{n+1,5} = ifr_bins;
    %spikes_save{n+1,6} = cell2mat(cluster_region(1,n+1)); % !!!ccomment out if no SHARPtrack !!!
end
clear cluster_group
clear ifr
clear ifr_bin
clear ifr_bins
clear k
clear m
clear nbin
clear numspikes

%% get the window idx from each stimulation (spikes) (required)
nwithinbin = 30000; % set samples per bin
WsecS = -300; % window in seconds before stim start EDIT
WsecE = 900; % window in seconds after stim ends EDIT
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
clear WsecS
clear WsecE

%% Sorting Figures, spike rate mean change, select unit (option)
% sort by mean change during stim compared to pre-stim window, threshold active units

unitype = 2; % 1 = mua, 2 = good, 3 = all
smoothfactor = 3; % EDIT to change smooth samples
clims = ([0 3]); % EDIT to change color scale
thresh = 0.05; % EDIT to change min firing rate during Basline

if unitype == 3
    current_spikes = spikes_save;
else
    current_spikes = spikes_save(cell2mat(spikes_save(:,4))==unitype,:); % get spike_save from cluster_groupid
end
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
    current_spikes2 = [];
        for n = 1:size(current_spikes,1) % only include active units during baseline
            current_spikes2 = current_spikes{n,5};
            if mean(current_spikes2(1,idx_bsS(id):idx_stimS(id)))>thresh % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    mean_change = [];
    current_spikes_plot = [];
        for n = 1:size(current_spikes,1) % calculate mean-change
                mean1 = mean(current_spikes{n,5}(idx_stimS(id):idx_stimE(id)));% data{1,1}(15:25,1)
                mean2 = mean(current_spikes{n,5}(idx_bsS(id):idx_stimS(id)));
                mean_change = [mean_change; mean1/mean2];
                current_spikes_plot = [current_spikes_plot; current_spikes{n,5}];
        end
    [~,I] = sort(mean_change(thre2));
    current_spikes_plot = current_spikes_plot(thre2,:);
    smooth_plot = movmean(current_spikes_plot,smoothfactor,2);
    imagesc(smooth_plot(I,idx_bsS(id):idx_asE(id)),clims)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim start
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim end
    colorbar
    colormap(Red_White_Blue)
    ylabel('Unit')
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change
%% Sorting Figures, Z-score mean change, select unit (option)
% sort by mean change during stim compared to pre-stim window, threshold active units.

unitype = 3; % 1 = mua, 2 = good, 3 = all
smoothfactor = 5; % EDIT to change smooth samples
clims = ([-2 2]); % EDIT to change color scale
thresh = 0.05; % EDIT to change min firing rate during Basline

if unitype == 3
    current_spikes = spikes_save;
else
    current_spikes = spikes_save(cell2mat(spikes_save(:,4))==unitype,:); % get spike_save from cluster_groupid
end
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
    current_spikes2 = [];
        for n = 1:size(current_spikes,1) % only include active units during baseline
            current_spikes2 = current_spikes{n,5};
            current_spikes_plot(n,:) = current_spikes{n,5};% make full plot matrix
            if mean(current_spikes2(1,idx_bsS(id):idx_stimS(id)))>thresh % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    current_spikes_plotZ = zscore(current_spikes_plot(thre2,idx_bsS(id):idx_asE(id)-1),1,2);
    smooth_plot = movmean(current_spikes_plotZ,smoothfactor,2);
    winshift = idx_bsS(id)-1;
    mean_change = [];
        for n = 1:size(smooth_plot,1) % calculate mean-change
                mean1 = mean(smooth_plot(n,idx_stimS(id)-winshift:idx_stimE(id)-winshift-1));
                mean2 = mean(smooth_plot(n,idx_bsS(id)-winshift:idx_stimS(id)-winshift-1));
                mean_change = [mean_change; mean1-mean2];
        end
    [~,I] = sort(mean_change);
    imagesc(smooth_plot(I,:),clims);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim start
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim end
    colorbar
    colormap(Red_White_Blue)
    ylabel('Unit')
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change

%% Sorting Figures, Z-score(BL) mean change, select unit (option)
% sort by mean change during stim compared to pre-stim window, threshold active units.

unitype = 3; % 1 = mua, 2 = good, 3 = all
smoothfactor = 1; % EDIT to change smooth samples
clims = ([-3 3]); % EDIT to change color scale
thresh = 0.05; % EDIT to change min firing rate during Basline

if unitype == 3
    current_spikes = spikes_save;
else
    current_spikes = spikes_save(cell2mat(spikes_save(:,4))==unitype,:); % get spike_save from cluster_groupid
end
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
    current_spikes2 = [];
        for n = 1:size(current_spikes,1) % only include active units during baseline
            current_spikes2 = current_spikes{n,5};
            current_spikes_plot(n,:) = current_spikes{n,5};% make full plot matrix
            blmean(n) = mean(current_spikes2(1,idx_bsS(id):idx_stimS(id)));
            blsdev(n) = std(current_spikes2(1,idx_bsS(id):idx_stimS(id)));
            if blmean(n) > thresh % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    % create baseline Z-scores
    current_spikes_plotZ = [];
        for n = 1:size(current_spikes_plot,1)
            current_spikes_plotZ(n,:) = ((current_spikes_plot(n,idx_bsS(id):idx_asE(id)-1))-blmean(n))/blsdev(n);
        end
    smooth_plot = movmean(current_spikes_plotZ(thre2,:),smoothfactor,2);
    winshift = idx_bsS(id)-1;
    mean_change = [];
        for n = 1:size(smooth_plot,1) % calculate mean-change
                mean1 = mean(smooth_plot(n,idx_stimS(id)-winshift:idx_stimE(id)-winshift-1));
                mean2 = mean(smooth_plot(n,idx_bsS(id)-winshift:idx_stimS(id)-winshift-1));
                mean_change = [mean_change; mean1-mean2];
        end
    [B,I] = sort(mean_change);
    imagesc(smooth_plot(I,:),clims);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim start
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim end
    %ytemp = find(B>1);
    %yline(ytemp(1),'k','Linewidth',2)
    %ytemp = find(B<-1);
    %yline(ytemp(end),'k','Linewidth',2)
    colorbar
    colormap(Red_White_Blue)
    ylabel('Unit')
end
%clear id
%clear thre
%clear n
%clear B
%clear I
%clear mean_change
%% Sorting Figures, Z-score latency, select unit (option)
% sort by mean change during stim compared to pre-stim window, threshold active units.
unitype = 3; % 1 = mua, 2 = good, 3 = all
smoothfactor = 5; % EDIT to change smooth samples
clims = ([-1 1]); % EDIT to change color scale
thresh = 0.05; % EDIT to change min firing rate during Basline
minmax = 1; % EDIT to change min=0 max=1

if unitype == 3
    current_spikes = spikes_save;
else
    current_spikes = spikes_save(cell2mat(spikes_save(:,4))==unitype,:); % get spike_save from cluster_groupid
end
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
    current_spikes2 = [];
        for n = 1:size(current_spikes,1) % only include active units during baseline
            current_spikes2 = current_spikes{n,5};
            current_spikes_plot(n,:) = current_spikes{n,5};% make full plot matrix
            if mean(current_spikes2(1,idx_bsS(id):idx_stimS(id)))>thresh % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    current_spikes_plotZ = zscore(current_spikes_plot(thre2,idx_bsS(id):idx_asE(id)-1),1,2);
    smooth_plot = movmean(current_spikes_plotZ,smoothfactor,2);
    if minmax == 0
        [peak, peak_latency] = min(smooth_plot(:,idx_stimS(id)-idx_bsS(id):idx_asE(id)-idx_bsS(id)),[],2);
    else
        [peak, peak_latency] = max(smooth_plot(:,idx_stimS(id)-idx_bsS(id):idx_asE(id)-idx_bsS(id)),[],2);
    end
    [~,I] = sort(peak_latency);
    imagesc(smooth_plot(I,:),clims)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim start
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim end
    colorbar
    colormap(Red_White_Blue)
    ylabel('Unit')
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change

%% Sorting Figures, Z-score(BL) latency, select unit (option)
% sort by mean change during stim compared to pre-stim window, threshold active units.

unitype = 3; % 1 = mua, 2 = good, 3 = all
smoothfactor = 5; % EDIT to change smooth samples
clims = ([-2 2]); % EDIT to change color scale
thresh = 0.05; % EDIT to change min firing rate during Basline
minmax = 1; % EDIT to change min=0 max=1

if unitype == 3
    current_spikes = spikes_save;
else
    current_spikes = spikes_save(cell2mat(spikes_save(:,4))==unitype,:); % get spike_save from cluster_groupid
end
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
    current_spikes2 = [];
        for n = 1:size(current_spikes,1) % only include active units during baseline
            current_spikes2 = current_spikes{n,5};
            current_spikes_plot(n,:) = current_spikes{n,5};% make full plot matrix
            blmean(n) = mean(current_spikes2(1,idx_bsS(id):idx_stimS(id)));
            blsdev(n) = std(current_spikes2(1,idx_bsS(id):idx_stimS(id)));
            if blmean(n) > thresh % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    % create baseline Z-scores
    current_spikes_plotZ = [];
        for n = 1:size(current_spikes_plot,1)
            current_spikes_plotZ(n,:) = ((current_spikes_plot(n,idx_bsS(id):idx_asE(id)-1))-blmean(n))/blsdev(n);
        end
    smooth_plot = movmean(current_spikes_plotZ(thre2,:),smoothfactor,2);
    if minmax == 0
        [peak, peak_latency] = min(smooth_plot(:,idx_stimS(id)-idx_bsS(id):idx_asE(id)-idx_bsS(id)),[],2);
    else
        [peak, peak_latency] = max(smooth_plot(:,idx_stimS(id)-idx_bsS(id):idx_asE(id)-idx_bsS(id)),[],2);
    end
    [~,I] = sort(peak_latency);
    imagesc(smooth_plot(I,:),clims);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim start
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2) % draw line for stim end
    colorbar
    colormap(Red_White_Blue)
    ylabel('Unit')
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change
%% Sorting Figures, Z-score, REGION, select unit (labels fixed), DOES NOT WORK - UPDATE TO NEW DATA ARRAYS
% sort by mean change during stim compared to pre-stim window, all
% units, threshold active units.
unitype = spikes_matrix; %spikes_matrix, spikes_matrix_good, spikes_matrix_mua
finalplot = [];
y = [];
ytick = [];
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
        for n = 1:size(unitype,1) % only include active units during baseline
            if mean(unitype(n,idx_bsS(id):idx_stimS(id)),2)>0.1 % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    activeregions = unique(cell2mat(spikes_save(thre2,6))); %row identifiers of probeinfo regions that have active units
    activelabels = ProbeInfo{probenum}(activeregions,3); %pull out string label for active regions
    for jj = 1:size(activeregions,1)
            regionIDX = cell2mat(spikes_save(:,6)) == activeregions(jj);
            regionIDX2 = regionIDX & thre2;
            spikes_matrixZ = zscore(unitype(regionIDX2,idx_bsS(id):idx_asE(id)),1,2); % z-score thresholded clusters within time window
            mean_changeZ = mean(spikes_matrixZ(:,idx_stimS(id)-idx_bsS(id):idx_stimE(id)-idx_bsS(id)),2);% Z-scored mean_change during stim compared to baseline
            [~,I] = sort(mean_changeZ);
            if jj == 1
                y(jj) = 1;
                finalplot = spikes_matrixZ(I,:);
                ytick(jj) = round(size(spikes_matrixZ(I,:),1)/2);
            else
                y(jj) = size(finalplot, 1);
                ytick(jj) = round((size(spikes_matrixZ(I,:),1)/2)+size(finalplot,1));
                finalplot = cat(1,finalplot, spikes_matrixZ(I,:));
            end
    end
    clims = ([-1 +1]); %z score limits
    imagesc(finalplot,clims)
    clear clims
    for ii = 1:size(y,2)
        yline(y(ii), 'LineWidth', 1)
    end
    yticks(ytick(1,:));
    yticklabels(activelabels);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2)
    colormap(Red_White_Blue)
    colorbar
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change

%% Sorting Figures, Z-score smooth latency, REGION, select unit (labels fixed), DOES NOT WORK - UPDATE TO NEW DATA ARRAYS 
% sort by mean change during stim compared to pre-stim window, all
% units, threshold active units.
unitype = spikes_matrix; %spikes_matrix, spikes_matrix_good, spikes_matrix_mua
finalplot = [];
y = [];
ytick = [];
figure;
for id = 1:(size(TrainInfo,2)) % id = window #
    subplot(1,size(TrainInfo,2),id)
    thre2 = [];
        for n = 1:size(unitype,1) % only include active units during baseline
            if mean(unitype(n,idx_bsS(id):idx_stimS(id)),2)>0.1 % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    activeregions = unique(cell2mat(spikes_save(thre2,6))); %row identifiers of probeinfo regions that have active units
    activelabels = ProbeInfo{probenum}(activeregions,3); %pull out string label for active regions
    for jj = 1:size(activeregions,1)
            regionIDX = cell2mat(spikes_save(:,6)) == activeregions(jj);
            regionIDX2 = regionIDX & thre2;
            spikes_matrixZ = zscore(unitype(regionIDX2,idx_bsS(id):idx_asE(id)),1,2); % z-score thresholded clusters within time window
            smoothZ = movmean(spikes_matrixZ,5,2);
            [min_peak, min_peak_latency] = min(smoothZ(:,idx_stimS(id)-idx_bsS(id):idx_asE(id)-idx_bsS(id)),[],2);
            [B,I] = sort(min_peak_latency);
            if jj == 1
                y(jj) = 1;
                finalplot = smoothZ(I,:);
                ytick(jj) = round(size(smoothZ(I,:),1)/2);
            else
                y(jj) = size(finalplot, 1);
                ytick(jj) = round((size(smoothZ(I,:),1)/2)+size(finalplot,1));
                finalplot = cat(1,finalplot, smoothZ(I,:));
            end
    end
    clims = ([-1 +1]); %z score limits
    imagesc(finalplot,clims)
    clear clims
    for ii = 1:size(y,2)
        yline(y(ii), 'LineWidth', 2)
    end
    yticks(ytick(1,:));
    yticklabels(activelabels);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2)
    colormap(Red_White_Blue)
    colorbar
end
clear id
clear thre
clear n
clear B
clear I
clear mean_change

%% Load .CSV files from DLC (For Pupil Diameter)
[pupilcsv, pupilfolder] = uigetfile(); %18N_d2_pupil-10082021140812-0000DLC_resnet50_Jonny_Pupil_Test_9_2_2021Sep2shuffle2_550000.csv
pupildata = readmatrix(strcat(pupilfolder, pupilcsv));

% Calculate Pupil Diameter
plds = [];
pbs = [];
for n = 1:size(pupildata,1)
pl = [pupildata(n,8) pupildata(n,9)];
pr = [pupildata(n,11) pupildata(n,12)];
pt = [pupildata(n,2) pupildata(n,3)];
pb = [pupildata(n,5) pupildata(n,6)];
pld = norm(pl-pr);
plds = [plds pld];
pbs = [pbs; pb];
end 
plds = movmedian(plds,120);
pbsxy = pbs(:,2);
pbsxy2 = movmedian(pbsxy,10);
pbsxysmooth = movmean(pbsxy2,10);
pbsxyZ = zscore(pbsxysmooth); % try to z-score within window?
plot(pbsxyZ)
figure; 
plot(plds)
figure;
plot(pbs)
%% Plot the pupil stimulus window (I THINK THERE IS SYNC PROBLEM HERE)
dig_2_fo = find(diff(digArray(2,:))==1); %find all frame TTL
lag = dig_2_fo(1); %find first from TTL
lag_s = round(lag/NIsr); % seconds PART OF SYNC ISSUE?
plds_align = repelem(plds,1000); % 30fps to 30,000 sampling rate
plds_d = downsample(plds,30);
pbs_align = repelem(pbsxyZ,1000);
pbs_d = downsample(pbsxyZ,30);
frameloss = size(dig_2_fo,2) - size(plds,2);
timeloss = frameloss/30; % offsets should not be much different from this value
disp(strcat('number of frames lost = ', num2str(frameloss)));
disp(strcat('time lost = ', num2str(timeloss)));
offset = [-30 -30 -30 -30, -30]; % window offsets (s) to allign pupildata, one offset per plot, why?

figure
for id = 1:size(TrainInfo,2)
    subplot(1,size(TrainInfo,2),id)
    plot(plds_d(idx_bsS(id)+lag_s+offset(id):idx_asE(id)+lag_s+offset(id)));
    %plot(pbs_d(idx_bsS(id)+lag_s+offset(id):idx_asE(id)+lag_s+offset(id)));
    xticks(0:15:idx_asE(id)-idx_bsS(id));
    xlim([0 idx_asE(id)-idx_bsS(id)]);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
    pyscale(id,1) = min(plds_d(idx_bsS(id)+lag_s+offset(id):idx_asE(id)+lag_s+offset(id)))-10; %find global min y-scale 10% buffer
    pyscale(id,2) = max(plds_d(idx_bsS(id)+lag_s+offset(id):idx_asE(id)+lag_s+offset(id)))+10; %find global max y-scale 10% buffer
    ylim([pyscale(id,1) pyscale(id,2)]); % large window to capture everything on first pass
end
%% Sorted figures, spike rate, and pupil, select units DOES NOT WORK - UPDATE TO NEW DATA ARRAYS
unitype = spikes_matrix_good; %spikes_matrix, spikes_matrix_good, spikes_matrix_mua
plotnum = size(TrainInfo,2);
figure
for id = 1:plotnum
    subplot(6,plotnum,id)
        plot(plds_d(idx_bsS(id)+lag_s+offset(id) :idx_asE(id)+lag_s+offset(id)));
        xticks(0:15:idx_asE(id)-idx_bsS(id));
        ylim([min(pyscale(:,1)) max(pyscale(:,2))]); 
        xlim([0 idx_asE(id)-idx_bsS(id)])
        xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
        xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
end
for id = 1:plotnum
    subplot(6,plotnum,[id+(plotnum) id+(plotnum*2) id+(plotnum*3) id+(plotnum*4) id+(plotnum*5)])
    thre = [];
        for n = 1:size(unitype,1)
        if mean(unitype(n,idx_bsS(id):idx_stimS(id)),2)>0.01
           thre = [thre n];
        else 
        end
    end 
    mean_change = mean(unitype(thre,idx_stimS(id):idx_stimE(id)),2)./mean(unitype(thre,idx_bsS(id):idx_stimS(id)),2);
[B,I] = sort(mean_change);
clims = ([0 1]);
imagesc(unitype(I,idx_bsS(id):idx_asE(id)),clims)
ylabel('Unit')
xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
%xticks(0:15:120);
end 
colormap(Red_White_Blue)
colorbar('Position',[1 0.9 0.02 0.33],'Location','east')

%% Sorted figures, z-scored spike rate, and pupil, select units DOES NOT WORK - UPDATE TO NEW DATA ARRAYS
unitype = spikes_matrix_good; %spikes_matrix, spikes_matrix_good, spikes_matrix_mua
plotnum = size(TrainInfo,2);
figure
for id = 1:plotnum
    subplot(6,plotnum,id)
        plot(plds_d(idx_bsS(id)+lag_s+offset(id) :idx_asE(id)+lag_s+offset(id)));
        xticks(0:15:idx_asE(id)-idx_bsS(id));
        ylim([min(pyscale(:,1)) max(pyscale(:,2))]); 
        xlim([0 idx_asE(id)-idx_bsS(id)])
        xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
        xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
end
for id = 1:plotnum % id = window #
    subplot(6,plotnum,[id+(plotnum) id+(plotnum*2) id+(plotnum*3) id+(plotnum*4) id+(plotnum*5)])
    thre = [];
        for n = 1:size(unitype,1) % only include units with spikes during baseline
            if mean(unitype(n,idx_bsS(id):idx_stimS(id)),2)>0.1
                thre = [thre n];
            else 
            end
        end 
    spikes_matrixZ = zscore(unitype(thre,idx_bsS(id):idx_asE(id)),1,2);
    mean_change = mean(spikes_matrixZ(:,31:61),2);
    [B,I] = sort(mean_change);
    clims = ([-1.5 +1.5]); %z score limits
    imagesc(spikes_matrixZ(I,:),clims)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2)
    colormap(Red_White_Blue)
    %colorbar
    ylabel('Unit')
end
%% Sorted figures, z-scored spike rate, region, and pupil, select units (BEST) (only works on all units right now) DOES NOT WORK - UPDATE TO NEW DATA ARRAYS
unitype = spikes_matrix; %spikes_matrix, spikes_matrix_good, spikes_matrix_mua
plotnum = size(TrainInfo,2);
figure
for id = 1:plotnum
    subplot(6,plotnum,id)
        plot(plds_d(idx_bsS(id)+lag_s+offset(id) :idx_asE(id)+lag_s+offset(id)));
        xticks(0:15:idx_asE(id)-idx_bsS(id));
        ylim([min(pyscale(:,1)) max(pyscale(:,2))]); 
        xlim([0 idx_asE(id)-idx_bsS(id)])
        xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
        xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
end
finalplot = [];
y = [];
for id = 1:plotnum % id = window #
    subplot(6,plotnum,[id+(plotnum) id+(plotnum*2) id+(plotnum*3) id+(plotnum*4) id+(plotnum*5)])
    thre2 = [];
        for n = 1:size(unitype,1) % only include active units during baseline
            if mean(unitype(n,idx_bsS(id):idx_stimS(id)),2)>0.02 % change this number for thresh: >x
                thre2 = [thre2 1];
            else
                thre2 = [thre2 0];
            end
        end
    thre2 = logical(thre2)';
    activeregions = unique(cell2mat(spikes_save(thre2,6))); %row identifiers of probeinfo regions that have active units
    activelabels = ProbeInfo{probenum}(activeregions,3); %pull out string label for active regions
    for jj = 1:size(activeregions,1)
            regionIDX = cell2mat(spikes_save(:,6)) == activeregions(jj);
            regionIDX2 = regionIDX & thre2;
            spikes_matrixZ = zscore(unitype(regionIDX2,idx_bsS(id):idx_asE(id)),1,2); % z-score thresholded clusters within time window
            mean_changeZ = mean(spikes_matrixZ(:,idx_stimS(id)-idx_bsS(id):idx_stimE(id)-idx_bsS(id)),2);% Z-scored mean_change during stim compared to baseline
            [~,I] = sort(mean_changeZ);
            if jj == 1
                y(jj) = 1;
                finalplot = spikes_matrixZ(I,:);
                ytick(jj) = round(size(spikes_matrixZ(I,:),1)/2);
            else
                y(jj) = size(finalplot, 1);
                ytick(jj) = round((size(spikes_matrixZ(I,:),1)/2)+size(finalplot,1));
                finalplot = cat(1,finalplot, spikes_matrixZ(I,:)) ;
            end
    end
    %clims = ([-(max(max(spikes_matrix_good))/150) max(max(spikes_matrix_good))/150]);
    clims = ([-1 +1]); %z score limits
    imagesc(finalplot,clims)
    %imagesc(1:1:size(spikes_matrixZ,2),B,spikes_matrixZ(I,:),clims)
    clear clims
    for ii = 1:size(y,2)
        yline(y(ii), 'LineWidth', 2)
    end
    yticks(ytick(1,:));
    yticklabels(activelabels);
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 2)
    colormap(Red_White_Blue)
    %colorbar
    ylabel('Unit')
end




%% Basic LFP data import

nChansInFile = str2double(LFPmeta.nSavedChans);
nSamps = lfpD.bytes/2/nChansInFile;
mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
LFPca1 = mean(mmf.Data.x,1); %Specific Channels in x

%% get the window idx from each stimulation (LFP)
WsecS = -60; % window in seconds before stim start
WsecE = 60; % window in seconds after stim ends
nwithinbin = 10;
idx_stimS = zeros(1,size(TrainInfo,2));
idx_stimE = zeros(1,size(TrainInfo,2));
idx_bsS = zeros(1,size(TrainInfo,2));
idx_asE = zeros(1,size(TrainInfo,2));
for n = 1:size(TrainInfo,2)
    idx_stimS(n) = round(TrainInfo(3,n)/NIsr*LFPsr)/nwithinbin; % index of each stim start time
    idx_stimE(n) = round(TrainInfo(4,n)/NIsr*LFPsr)/nwithinbin; % index of each stim end time
    idx_bsS(n) = idx_stimS(n)+round(WsecS*LFPsr)/nwithinbin;% index before each stim (window start)
    idx_asE(n) = idx_stimE(n)+round(WsecE*LFPsr)/nwithinbin;% index after each stim end (window end)
end 

%% Simple LFP plot

%Plot all stim Hz, mean LFP across entire probe
figure
LFPmeanDS = downsample(LFPmean,nwithinbin);
plot(LFPmeanDS);
for id = 1:size(TrainInfo,2)
    subplot(1,size(TrainInfo,2),id)
    dataseg = LFPmeanDS(idx_bsS(id):idx_asE(id));
    plot(dataseg);
    %periodogram(dataseg,rectwin(length(dataseg)),length(dataseg),LFPsr)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
    ylim([-500 500])
end

%Plot 20Hz 30sec stim mean LFP across entire probe
for id = 1:size(TrainInfo, 2)
    figure
    dataseg = LFPmeanDS(idx_bsS(id):idx_asE(id));
    plot(dataseg);
    %periodogram(dataseg,rectwin(length(dataseg)),length(dataseg),LFPsr)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 3)
%     xline(((3+6.4788)*LFPsr)/nwithinbin,'r','LineWidth', 2,'LineStyle',':')
%     xline((13+6.4788)*LFPsr,'k','LineWidth', 2,'LineStyle',':')
%     xline((6.4788)*LFPsr,'r','LineWidth', 2)
%     xline((16.4788)*LFPsr,'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 3)
    ylim([-500 500])
end

%% Morse Wavelet analysis in Time X Freq domain
bsFilt = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'SampleRate',LFPsr); %remove 60Hz line noise
for id = 1:size(TrainInfo, 2)
    figure
    dataseg = filtfilt(bsFilt,LFPmean(idx_bsS(id):idx_asE(id))); % get the data segment to calculate wavelets using above filter
    fb = cwtfilterbank('SignalLength',numel(dataseg),'SamplingFrequency',LFPsr,'FrequencyLimits',[0.5 120]); % only show frequncies between 0.5-120Hz
    [cfs, frq] = cwt(dataseg,'FilterBank', fb); % cfs = complex number matrix, frq = frequencies
    t = 0:length(dataseg)-1;
    tsec = t/LFPsr;
    power = abs(cfs);
    powerBL = [];
    for f = 1:size(frq,1)
        powerBL(f,:) = power(f,:)/mean(power(f,1:abs(WsecS*LFPsr)),2);
    end
    pcolor(tsec,frq,powerBL);
    shading interp
    set(gca, 'YScale', 'log') %comment in for log scale visual
    xline(abs(WsecS),'k','LineWidth', 4)
    xline(tsec(end)-WsecE,'k','LineWidth', 4)
    yline(4,'w','LineWidth', 2) % Delta
    yline(8,'w','LineWidth', 2) % Theta
    yline(13,'w','LineWidth', 2) % Alpha
    yline(30,'w','LineWidth', 2) % Beta
    yline(60,'w','LineWidth', 2) % Low Gamma
    colormap(turbo)
    colorbar
    caxis([0 3])
    exportname = strcat('D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\Wavelet_Analysis\',num2str(id),'_',lfpD.name,'.csv');
    writematrix(powerBL,exportname) 
end

%% ECG
ECGraw = NIdataArray(1,:);
bsFilt = designfilt('bandstopiir','FilterOrder',10,'HalfPowerFrequency1',50,'HalfPowerFrequency2',100,'SampleRate',NIsr);
ECGraw = filtfilt(bsFilt,ECGraw);
hpFilt = designfilt('highpassiir','FilterOrder',8,'PassbandFrequency',2,'SampleRate',NIsr);
ECGraw = filtfilt(hpFilt,ECGraw);
% get the window idx from each stimulation (ECG)
WsecS = -10; % window in seconds before stim start
WsecE = 3; % window in seconds after stim ends
nwithinbin = 10;
idx_stimS = zeros(1,size(TrainInfo,2));
idx_stimE = zeros(1,size(TrainInfo,2));
idx_bsS = zeros(1,size(TrainInfo,2));
idx_asE = zeros(1,size(TrainInfo,2));
for n = 1:size(TrainInfo,2)
    idx_stimS(n) = round(TrainInfo(3,n))/nwithinbin; % index of each stim start time
    idx_stimE(n) = round(TrainInfo(4,n))/nwithinbin; % index of each stim end time
    idx_bsS(n) = idx_stimS(n)+round(WsecS)*NIsr/nwithinbin;% index before each stim (window start)
    idx_asE(n) = idx_stimE(n)+round(WsecE)*NIsr/nwithinbin;% index after each stim end (window end)
end 
figure
ECGrawDS = downsample(ECGraw,nwithinbin);
for id = 1:size(TrainInfo,2)
    subplot(1,size(TrainInfo,2),id)
    dataseg = ECGrawDS(idx_bsS(id):idx_asE(id));
    plot(dataseg);
    %periodogram(dataseg,rectwin(length(dataseg)),length(dataseg),LFPsr)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 1)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 1)
    ylim([-3 5])
end

%Plot stim lines and syncope on figure
for id = 1:size(TrainInfo, 2)
    figure
    dataseg = ECGrawDS(idx_bsS(id):idx_asE(id));
    plot(dataseg);
    %periodogram(dataseg,rectwin(length(dataseg)),length(dataseg),LFPsr)
    xline(idx_stimS(id)-idx_bsS(id),'k','LineWidth', 3)
    xline(((10+6.4788)*NIsr)/nwithinbin,'r','LineWidth', 3,'LineStyle',':')
    %xline((13+6.4788)*NIsr,'k','LineWidth', 2,'LineStyle',':')
    %xline((6.4788)*NIsr,'r','LineWidth', 2)
    %xline((16.4788)*NIsr,'k','LineWidth', 2)
    xline(idx_stimE(id)-idx_bsS(id),'k','LineWidth', 3)
    ylim([-2 3])
    xlim([0 size(dataseg,2)])
end
%% Plot entire set
bin_label = 5; % what measure to plot from "spikes_save(:,bin_label)" [5=ifr_bins]
spikes_matrix = [];
for n = 1:length(spikes_save) % 1:number of total clusters
spikes_plot = spikes_save{n,bin_label};
spikes_matrix = [spikes_matrix; spikes_plot];
end 
clims = ([0 2]);
figure;
imagesc(spikes_matrix,clims)
colorbar
colormap(Red_White_Blue)
clear clims
% Good clusters
spikes_matrix_good = [];
for n = 1:length(spikes_save)
    spikes_plot = spikes_save{n,bin_label};
    if spikes_save{n,4} == 2 & spikes_save{n,2}>0
        spikes_matrix_good = [spikes_matrix_good; spikes_plot];
    else 
    end
end 
clims = ([0 2]);
figure;
imagesc(spikes_matrix_good,clims)
colorbar
colormap(Red_White_Blue)
clear clims
% Mua clusters
spikes_matrix_mua = [];
for n=1:length(spikes_save)
spikes_plot=spikes_save{n,bin_label};
    if spikes_save{n,4}==1 & spikes_save{n,2}>0
    spikes_matrix_mua = [spikes_matrix_mua; spikes_plot];
    else 
    end
end 
clims = ([0 2]);
figure;
imagesc(spikes_matrix_mua,clims)
colormap(Red_White_Blue)
colorbar

clear bin_label
clear clims