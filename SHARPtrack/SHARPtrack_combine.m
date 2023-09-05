%% manually input "processed" folder directories

dirlist{1} = 'D:\SHARPTrack\16N_12-28\processed';
dirlist{2} = 'D:\SHARPTrack\18N_10-12\processed';
dirlist{3} = 'D:\SHARPTrack\32N_10-19\processed';
dirlist{4} = 'D:\SHARPTrack\47L_12-30\processed';
dirlist{5} = 'D:\SHARPTrack\57L_01-28\processed';
dirlist{6} = 'D:\SHARPTrack\59L_01-07\processed';
dirlist{7} = 'D:\SHARPTrack\69N_01-20\processed';
dirlist{8} = 'D:\SHARPTrack\74L_02-11\processed';
dirlist{9} = 'D:\SHARPTrack\22R_03-10\processed';

%% compile all probe tracks into "master"
pointsListMaster.pointList = [];
pointsListMaster.pointHands = [];
pointlist2 = [];
pointhands2 = [];
for k = 1:size(dirlist,2)
    curmouse = strcat(dirlist{k},'\probe_pointselectrodetrack1.mat');
    load(curmouse);
    pointlist2 = [pointlist2; pointList.pointList];
    pointhands2 = [pointhands2; pointList.pointHands];
end

%% Save all probe trajectories into file that "Display_Probe_Track_ALL.m can read
pointsListMaster.pointList = pointlist2;
pointsListMaster.pointHands = pointhands2;
pointList = pointsListMaster;
save('D:\SHARPTrack\probe_pointselectrodetrack1.mat','pointList');
