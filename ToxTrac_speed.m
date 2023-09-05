fid = 'Instant_Speed.txt';
data = readmatrix(fid);
[B,ia,ic] = unique(data(:,2));
a_counts = accumarray(ic,1);

%% Atro group stim time video
time(1) = 309;
time(2) = 297;
time(3) = 310;
time(4) = 298;
time(5) = 309;
time(6) = 118;
time(7) = 301;
time(8) = 303;
time(9) = 299;
time(10) = 385;
time(11) = 301;
time(12) = 295;
time(13) = 300;
time(14) = 298;
time(15) = 300;
time(16) = 327;
time(17) = 301;
time(18) = 297;
time(19) = 301;
time(20) = 300;
time(21) = 298;
time(22) = 304;
time(23) = 302;
time(24) = 302;

%% Hm4di group stim time video
time(1) = 298;
time(2) = 297;
time(3) = 296;
time(4) = 311;
time(5) = 300;
time(6) = 296;
time(7) = 299;
time(8) = 295;
% time(9) = 293; 198RL
% time(10) = 297; 198RL
time(9) = 298;
time(10) = 302;
time(11) = 296;
time(12) = 297;
time(13) = 297;
time(14) = 297;

%% Hm3dq group stim time video
time(1) = 299;
time(2) = 290;
time(3) = 301;
time(4) = 264;
time(5) = 290;
time(6) = 289;
time(7) = 294;
time(8) = 268;
time(9) = 297;
time(10) = 342;
time(11) = 296;
time(12) = 295;
time(13) = 302;
time(14) = 285;
time(15) = 299;
time(16) = 319;
time(17) = 299;
time(18) = 306;
time(19) = 299;
time(20) = 336;

%% stim alligned traces
for i = 1:size(B)
    idx = (data(:,2))==B(i);
    data2{i} = data(idx,5);
%     if i == 6
%     data2{i} = NaN(9601,1);
%         continue
%     end
    s = (time(i)*20)-4800; % 4800 = 4min before stim (4min*60s*20fps)
    e = (time(i)*20)+4800; % 4800 = 4min after stim (4min*60s*20fps)
    data2{i} = movmean(data2{i}(s:e),1);
    data2mean(1,i) = mean(data2{i}(1:3600)); % 20FPS, 3600/20 = 180s = 3min
    data2mean(2,i) = mean(data2{i}(end-3600:end));
end

%% avergae trace2
data3 =[];
for i = 1:size(B)
    idx = (data(:,2))==B(i);
    datatemp = movmean(data2{i},800,'omitnan');
    data3 = [data3 datatemp];
end
meanCNO = mean(data3(:,1:2:end), 2,'omitnan');
meanVEH = mean(data3(:,2:2:end), 2,'omitnan');
figure
plot(meanVEH,'k','LineWidth', 2);
hold
plot(meanCNO,'r','LineWidth', 2);
xlim([0 9600])
ylim([0 25])

%% Prism fix order
tracesDrug = data3(:,1:2:end);
tracesVeh = data3(:,2:2:end);

meanDrug = data2mean(:,1:2:end);
meanVeh = data2mean(:,2:2:end);