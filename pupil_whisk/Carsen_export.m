%% Import Data
%Data organization script for exporting data with selected paramter

% The script will get all folders listed below and create individual files
% for each probe recording
ProbePath = {
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\16N_D1_Insular_g0\16N_D1_Insular_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\22R_D1_SFO_MidBrain_g0\22R_D1_SFO_MidBrain_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\22R_D1_SFO_MidBrain_g0\22R_D1_SFO_MidBrain_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\22R_D2_SFO_MidBrain_g0\22R_D2_SFO_MidBrain_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\22R_D2_SFO_MidBrain_g0\22R_D2_SFO_MidBrain_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\47L_D1_Amyg_ACC_g1\47L_D1_Amyg_ACC_g1_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\47L_D1_Amyg_ACC_g1\47L_D1_Amyg_ACC_g1_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\47L_D2_Amyg_ACC_g0\47L_D2_Amyg_ACC_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\47L_D2_Amyg_ACC_g0\47L_D2_Amyg_ACC_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\57L_mPFC_Insular_D1_g0\57L_mPFC_Insular_D1_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\57L_mPFC_Insular_D1_g0\57L_mPFC_Insular_D1_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\57L_mPFC_Insular_D2_g0\57L_mPFC_Insular_D2_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\57L_mPFC_Insular_D2_g0\57L_mPFC_Insular_D2_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\59L_D1_Hippocampus_g0\59L_D1_Hippocampus_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\59L_D2_Hippocampus_g0\59L_D2_Hippocampus_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\69N_D1_mPFC_Insular_g0\69N_D1_mPFC_Insular_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\69N_D2_mPFC_Insular_g0\69N_D2_mPFC_Insular_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\69N_D2_mPFC_Insular_g0\69N_D2_mPFC_Insular_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\74L_D1_PVN_SN_g0\74L_D1_PVN_SN_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\74L_D1_PVN_SN_g0\74L_D1_PVN_SN_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\74L_D2_PVN_SN_g0\74L_D2_PVN_SN_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\74L_D2_PVN_SN_g0\74L_D2_PVN_SN_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_18N_Day1_g0\NPY2R_Ai32_18N_Day1_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_18N_Day1_g0\NPY2R_Ai32_18N_Day1_g0_imec1',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_32N_Day1_g0\NPY2R_Ai32_32N_Day1_g0_imec0',
    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_32N_Day2_g0\NPY2R_Ai32_32N_Day2_g0_imec0'
    }; % all valid recording above, bellow are "sort error" or no SHARPtrack

%    'D:\Dropbox (Scripps Research)\Augustine-2_g0\NPY2R_Ai32_18N_DLab\NPY2R_NP\69N_D1_mPFC_Insular_g0\69N_D1_mPFC_Insular_g0_imec0', 
%    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_18N_Day2_g0\NPY2R_Ai32_18N_Day2_g0_imec0',
%    'D:\Dropbox (Scripps Research)\Augustine-Lab\NPY2R_NP\NPY2R_Ai32_18N_Dayay2_g0_imec1',

probenum = [1,1,3,2,4,1,3,2,4,1,3,2,4,1,2,3,2,4,1,3,2,4,1,2,1,2]; % probe num from sharptrack depth_table in order for each file listed above
for i = 1:size(ProbePath,1)
% get Imec Probe Data after Kilosort
ProbePathCur = ProbePath{i};
sp = loadKSdir(ProbePathCur);% Steinmetz function to load alot of basic data info into "sp" datastructure
APsr = sp.sample_rate; % samping rate (HZ) for AP channels
ProbeEntireDur = sp.st(end); % length of file in Secs. (Time of last recorded spike)
tt = 0:1/APsr:ProbeEntireDur; % set the time scale (x-axis) for spikes
[~, ~, sp.templateDepthsT] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps); %get depth from templates
sp.templateDepthsT = sp.templateDepthsT.';

% Import LFP channels info
%lfpD = dir(fullfile(ProbePath, '*.lf.bin')); % LFP file from spikeGLX specifically
%lfpFilename = fullfile(ProbePath, lfpD.name);
%LFPmeta = ReadMeta(lfpD.name, ProbePath);
%LFPsr = str2double(LFPmeta.imSampRate);
%LFPsamples = str2double(LFPmeta.fileTimeSecs)*LFPsr;

% Find NI binary file based on previous path
idcs   = strfind(ProbePathCur,'\'); % find all folder markers
NIpath = ProbePathCur(1:idcs(end)-1); % go back 1 folder from Probe path
NIbinInfo = dir(fullfile(NIpath, '*.bin')); % get file name of .bin in new folder (NI .bin file)

% Parse the corresponding NI metafile
NImeta = ReadMeta(NIbinInfo.name, NIpath);
NIsr = str2double(NImeta.niSampRate); %sampling rate for NI channels

% Import NI data
NInSamp = floor(str2double(NImeta.fileTimeSecs)*NIsr); %Total NI samples
NIdataArray = ReadBin(0, NInSamp, NImeta, NIbinInfo.name, NIpath); % dimension: [nChan,nSamp] (1-8:analog channel; 9:digital channel)

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

% Get SHARPtrack depth boundries from file and assign region index to clusters (not required to generate plots)
%auto script to get depth_table.mat
depthfile = fullfile(NIpath, 'depth_table.mat'); % saved file with all probe information combined
load(depthfile);
probenumcur = probenum(i); %get SHARPtrack for current data
tip_regionconvert = (cell2mat(ProbeInfo{probenumcur}(:,1:2))-max(cell2mat(ProbeInfo{probenumcur}(:,2))))*-1; % depth calculation from tip in SHARPtrack
cluster_region = cell(4,size(sp.cids,2));
for j = 1:size(tip_regionconvert,1)
    for k = 1:size(sp.cids,2)
        if (sp.templateDepthsT(1,k)>=tip_regionconvert(j,2) && (sp.templateDepthsT(1,k)<tip_regionconvert(j,1)))
            cluster_region{1,k} = j;
            cluster_region{2,k} = cell2mat(ProbeInfo{probenumcur}(j,5));
            cluster_region{3,k} = ProbeInfo{probenumcur}(j,3);
            cluster_region{4,k} = ProbeInfo{probenumcur}(j,4);
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
% Find TTL stim trains and get time stamps
lowestHzStim = 4; % lowest stim rate in Hz
TrainChan = 3; %Stim channel # in digitalArray
TrainInfo = find_stim_train(digArray(TrainChan,:),NIsr,lowestHzStim); % Get train timestamps (custom function)
clear lowestHzStim

% Load .CSV files from DLC (For Pupil Diameter) in predermined location
if or(i==4,i==5) % exclude probe data with no pupil video
    pupildata = [];
    plds = [];
else
    dlcpath = strcat(NIpath,'\Video\'); % get dlc path based on NIpath
    dlcinfo = dir(fullfile(dlcpath, '*.csv')); % get file name of .csv in Video folder
    pupildata = readmatrix(strcat(dlcpath, dlcinfo.name)); % create the pupildata matrix from file

% Calculate Pupil Diameter
    plds = [];
    for n = 1:size(pupildata,1)
        pl = [pupildata(n,8) pupildata(n,9)];
        pr = [pupildata(n,11) pupildata(n,12)];
        pt = [pupildata(n,2) pupildata(n,3)];
        pb = [pupildata(n,5) pupildata(n,6)];
        pld = norm(pl-pr);
        plds = [plds pld];
    end 
end
% Name and Export
savename = strcat('D:\Dropbox (Scripps Research)\Neuropixels_syncope\Neuropixels\Data2\',ProbePathCur(idcs(5)+1:end),'_data.mat');
save(savename,'NIdataArray','NImeta','TrainInfo','digArray','sp','cluster_region','pupildata','plds','-v7.3') % save only required info
end