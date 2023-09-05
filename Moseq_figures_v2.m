%% Import .csv file from Moseq output
Data = uigetfile();
T = readtable(Data);

% Organize table and sylables
SylArray = table2array(T(:,12));
%SylArray = table2array(T(:,29));
syl={};
syl{1} = [2 6 8 15 21 39 35 37 64 65 66 67 71 72 73]; %freeze
syl{2} = [9 20 25 50 52]; % move(slow)
syl{3} = [0 3 7 10 12 22 23 31 39 40 42];% move(norm)
syl{4} = [19 24 27 28 32 48 74];%move(fast)
syl{5} = [6 11 13 43 44 45 46 53 55 58 59 60 62]; % rear
syl{6} = [1 26 29 41 49]; % sniff
syl{7}= [4 14 30 38 47 51 54 56 61 63 70]; %groom
syl{8}= [16 17 33 57 75 76]; % Artifact
syl{9} = [18 36 68 69]; % Hunched
maxsyl = max(table2array(T(:,12)));
syl{10} = [77:maxsyl]; % unassaugned

% assign sylable to category
for n = 1:size(syl,2)
    index = find(ismember(SylArray(:,1),syl{n})==1);
    SylArray(index,2)= n;
end 
T.SylCat = SylArray(:,2);% add sylable category to original table
%% fix session labels
T.SessionName(strcmp(T.SessionName,'No_light')) = {'No_Light'};
T.SubjectName(strcmp(T.SubjectName,'NPY2R_AI9_C1451_20R')) = {'NPY2R_AI9_C1451A_20R'};

%% Stim Videos
SylMatrix = [];
SylCol = size(T,2);

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1341A_11R');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1341_12L');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1389B_16N');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1451A_20R');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C2038A_25L');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
FFF = [table2array(T(FF((1:77377),1),SylCol))' zeros(1,80000-77377)];
FFF(1,77377:80000)= FFF(1,77377:80000)+8;
SylMatrix = [SylMatrix; FFF]; %short vid
clear FFF

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1369_21R');% actually Ai32 mouse
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI32_C1422_34L');
F(:,2) = strcmp(T.SessionName,'Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

FF = find(strcmp(T.SubjectName,'6L-stim-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];
FF = find(strcmp(T.SubjectName,'3R-stim-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];
FF = find(strcmp(T.SubjectName,'1N-stim-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];
FF = find(strcmp(T.SubjectName,'11N-stim-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:80000),1),SylCol))'];

figure
load('MoseqColorMap_10v2.mat')
colormap(Moseq10v2)
%colormap(turbo)
SylMatrixAdj = SylMatrix + 0.5;
clims = [1 size(syl,2)+1];
subplot(11,10,[1:8 11:18 21:28 31:38 41:48 51:58 61:68 71:78 81:88 91:98 101:108])
imagesc(SylMatrixAdj, clims);
%imagesc(SylMatrix, clims);
colorbar

BLframes = 30*15*60; % 30fps, 15 min, 60 sec
xline(BLframes,'--w','LineWidth', 5)

BLtotals = zeros(size(SylMatrix,1),size(syl,2));
STIMtotals = zeros(size(SylMatrix,1),size(syl,2));
RATIOtotals = zeros(size(SylMatrix,1),size(syl,2));
for n = 1:size(SylMatrix,1)
    for k = 1:size(syl,2)
        BLtotals(n,k) = sum(SylMatrix(n,1:BLframes)==k);
        STIMtotals(n,k) = sum(SylMatrix(n,BLframes+1:end)==k);
        RATIOtotals(n,k) = (STIMtotals(n,k)/size(SylMatrix(n,BLframes+1:end),2))/(BLtotals(n,k)/size(SylMatrix(n,1:BLframes),2));
    end
    subplot(11,10,(n*10)-1)
    labels = cell(1,size(BLtotals(n,:),2));
    labels(:) = {''};
    pie(BLtotals(n,:),labels)
    subplot(11,10,n*10)
    labels = cell(1,size(BLtotals(n,:),2));
    labels(:) = {''};
    pie(STIMtotals(n,:), labels)
end

%% Baseline Videos
SylMatrix = [];
SylCol = size(T,2);

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1341A_11R');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1341_12L');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1389B_16N');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1451A_20R');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C2038A_25L');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];




F(:,1) = strcmp(T.SubjectName,'NPY2R_AI9_C1369_21R');% error no "Light"
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

F(:,1) = strcmp(T.SubjectName,'NPY2R_AI32_C1422_34L');
F(:,2) = strcmp(T.SessionName,'No_Light');
FF = find(F(:,1)&F(:,2));
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),SylCol))'];

FF = find(strcmp(T.SubjectName,'6L-baseline-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),16))'];
FF = find(strcmp(T.SubjectName,'3R-baseline-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),16))'];
FF = find(strcmp(T.SubjectName,'1N-baseline-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),16))'];
FF = find(strcmp(T.SubjectName,'11N-baseline-day1   45min' ));%11N-stim-day1   45min
SylMatrix = [SylMatrix; table2array(T(FF((1:79000),1),16))'];

figure
load('MoseqColorMap_10v2.mat')
colormap(Moseq10v2)
%colormap(turbo)
SylMatrixAdj = SylMatrix + 0.5;
clims = [1 size(syl,2)+1];
subplot(11,10,[1:8 11:18 21:28 31:38 41:48 51:58 61:68 71:78 81:88 91:98 101:108])
imagesc(SylMatrixAdj, clims);
%imagesc(SylMatrix, clims);
colorbar

BLframes = 30*15*60; % 30fps, 15 min, 60 sec
xline(BLframes,'--w','LineWidth', 5)

BLtotals = zeros(size(SylMatrix,1),size(syl,2));
STIMtotals = zeros(size(SylMatrix,1),size(syl,2));
RATIOtotals = zeros(size(SylMatrix,1),size(syl,2));
for n = 1:size(SylMatrix,1)
    for k = 1:size(syl,2)
        BLtotals(n,k) = sum(SylMatrix(n,1:end)==k);
        STIMtotals(n,k) = sum(SylMatrix(n,1:end)==k);
        RATIOtotals(n,k) = (STIMtotals(n,k)/size(SylMatrix(n,BLframes+1:end),2))/(BLtotals(n,k)/size(SylMatrix(n,1:BLframes),2));
    end
    subplot(11,10,(n*10)-1)
    labels = cell(1,size(BLtotals(n,:),2));
    labels(:) = {''};
    pie(BLtotals(n,:),labels)
    subplot(11,10,n*10)
    labels = cell(1,size(BLtotals(n,:),2));
    labels(:) = {''};
    pie(STIMtotals(n,:), labels)
end

%%