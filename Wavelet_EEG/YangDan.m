%% run after opening .mat files from Yang Dan
data1 = EEG(1,:)';
LFPsr = 1525.9;

%% !!!only run to get a "no light" marker based on 3min before stim train!!!
WsecS = -30; % window in seconds before stim start EDIT
WsecE = 30; % window in seconds after stim ends EDIT
idx_stimS = zeros(1,size(laser,1));
idx_stimE = zeros(1,size(laser,1));
idx_bsS = zeros(1,size(laser,1));
idx_asE = zeros(1,size(laser,1));
for n = 1:size(laser,1)
    idx_stimS(n) = round(laser(n,1))-(200); % index of each stim start time
    idx_stimE(n) = round(laser(n,2))-(200); % index of each stim end time
    idx_bsS(n) = (idx_stimS(n)+WsecS)*LFPsr;% index before each stim (window start)
    idx_asE(n) = (idx_stimE(n)+WsecE)*LFPsr;% index after each stim end (window end)
    figure
    plot(data1(idx_bsS(n):idx_asE(n),1));
end 
%% Yang Dan Data, get the window idx from each stimulation
WsecS = -30; % window in seconds before stim start EDIT
WsecE = 30; % window in seconds after stim ends EDIT
idx_stimS = zeros(1,size(laser,1));
idx_stimE = zeros(1,size(laser,1));
idx_bsS = zeros(1,size(laser,1));
idx_asE = zeros(1,size(laser,1));
for n = 1:size(laser,1)
    idx_stimS(n) = round(laser(n,1)); % index of each stim start time
    idx_stimE(n) = round(laser(n,2)); % index of each stim end time
    idx_bsS(n) = (idx_stimS(n)+WsecS)*LFPsr;% index before each stim (window start)
    idx_asE(n) = (idx_stimE(n)+WsecE)*LFPsr;% index after each stim end (window end)
    figure
    plot(data1(idx_bsS(n):idx_asE(n),1));
end 
