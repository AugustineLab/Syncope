%% combine pupil data with other data
testfiledir = 'D:\Dropbox (Scripps Research)\Neuropixels_syncope\Neuropixels\Data';
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);

% get pupil data directory
pupilpath = 'D:\Dropbox (Scripps Research)\Neuropixels_syncope\Neuropixels\Data\Pupil';
matfiles2 = dir(fullfile(pupilpath, '*.mat'));