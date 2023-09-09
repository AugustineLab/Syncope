%% add TTL to PBG breathing data timestamp
% get all .mat files
testfiledir = uigetdir();
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);
data  = cell(1,nfiles);
for w = 1 : nfiles
   fid = fullfile(testfiledir, matfiles(w).name);
end

% manual timestamp times in minutes
ts(1) = 14.98;
ts(2) = 6.036;
ts(3) = 5.465;
ts(4) = 5.675;
ts(5) = 25.616;
ts(6) = 15.366;
ts(7) = 10.349;
ts(8) = 11.727;
ts(9) = 8.098;
ts(10) = 20.571;
ts(11) = 12.223;
ts(12) = 10.257;
ts(13) = 10.228;
ts(14) = 11.183;
ts(15) = 8.240;
ts(16) = 15.320;

%% loop .mat files and add timestamps
for i = 1:nfiles
    load(fullfile(testfiledir, matfiles(i).name))
    if strcmp(isi_units, 'ms')
        scale = 0.001;
    else
        scale = 1;
    end
    sr = 1/(isi*scale);
    sample = ts(i)*(60*sr); %convert min timestamp to sample
    
    % generate 20Hz train marker for findstimtrain function
    T = 0.1; % 100 milliseconds timeperiod, 1 square wave in 100 milliseconds
    t = linspace(0, T*10, 1000);
    y = (square(t/T*2*pi)+1)*2.5; %scale from 0-5V
    train = repelem(y,10)';% expland to 10kHz sr
    
    % write TTL train to timestamp location
    data(sample:sample+9999,2) = train;

    % make new file with TTL channel with new marker
    filename = strcat(testfiledir,'\TTL_Add\',matfiles(i).name);
    save(filename,'data','units','labels','isi','isi_units','start_sample');
end

