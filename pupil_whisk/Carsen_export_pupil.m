%% Import Data
%Data organization script for carsen with all parameters

% The script will get all files needed and export data
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
for i = 1:size(ProbePath,1)
ProbePathCur = ProbePath{i};

% Find NI binary file based on previous path
idcs   = strfind(ProbePathCur,'\'); % find all folder markers
NIpath = ProbePathCur(1:idcs(end)-1); % go back 1 folder from Probe path


% Load .CSV files from DLC (For Pupil Diameter)
if or(i==4,i==5) % exclude probe data with no pupil video
    pupildata = [];
    plds = [];
    disp('flagged')
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
savename = strcat('D:\Dropbox (Scripps Research)\Neuropixels_syncope\Neuropixels\Data\Pupil\',ProbePathCur(idcs(5)+1:end),'_pupil.mat');
save(savename,'pupildata','plds','-v7.3')
end