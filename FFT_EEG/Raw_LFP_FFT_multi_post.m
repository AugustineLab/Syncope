%% run after loading .mat file from "Raw_LFP_multi" saved file
%row1 = name, row2 = Drug (CNO=1, Veh=0), row3 = Hz (20Hz=1, 10Hz=0),
%row4 = data{ch}{pre,post,ratio}

% Set Filter Settings
IDXfilt1 = cell2mat(DataArray(:,2)) == 1; % find Veh=0 CNO=1
IDXfilt2 = cell2mat(DataArray(:,3)) == 1; % find 10Hz=0 20Hz=1
IDX = IDXfilt1&IDXfilt2;
DataFilt1 = DataArray(IDX,4);
labels = DataArray(IDX,1);

    L = 2500*2; % Length of signal
    f = 2500*(1:(L/2))/L; %Frequency list

% Make custom array
clear Power
for i = 1:size(DataFilt1,1)
    tempdata(i,:) = cell2mat(DataFilt1{i}{2}(3)); % {mouse}{ch(1,2)}(1(pre),2(post),3(ratio))
    Power(1,i) = mean(tempdata(i,find(f==0.5):find(f==4))); 
    Power(2,i) = mean(tempdata(i,find(f==4.5):find(f==8))); 
    Power(3,i) = mean(tempdata(i,find(f==8.5):find(f==12))); 
    Power(4,i) = mean(tempdata(i,find(f==12.5):find(f==30))); 
    Power(5,i) = mean(tempdata(i,find(f==31.5):find(f==59)));
    Power(6,i) = mean(tempdata(i,find(f==61):find(f==120)));
end
tempdata(:,find(f==59):find(f==61)) = NaN;
tempmean1 = mean(tempdata,1);
tempmean1(:,find(f==59):find(f==61)) = NaN;

%% draw figure with multiple lines
clear Power
for i = 1:size(DataFilt1,1)
    tempdata(i,:) = cell2mat(DataFilt1{i}{2}(1)); % {mouse}{ch(1,2)}(1(pre),2(post),3(ratio))
    Power(1,i) = mean(tempdata(i,find(f==0.5):find(f==4))); 
    Power(2,i) = mean(tempdata(i,find(f==4.5):find(f==8))); 
    Power(3,i) = mean(tempdata(i,find(f==8.5):find(f==12))); 
    Power(4,i) = mean(tempdata(i,find(f==12.5):find(f==30))); 
    Power(5,i) = mean(tempdata(i,find(f==31.5):find(f==59)));
    Power(6,i) = mean(tempdata(i,find(f==61):find(f==120)));
end
tempmean2 = mean(tempdata,1);
tempmean2(:,find(f==59):find(f==61)) = NaN;

%% draw figure
% semilogx(f,tempmean1,'k','LineWidth', 3)
% hold on
% semilogx(f,tempmean2,'r','LineWidth', 3)

loglog(f,tempmean1,'b','LineWidth', 2)
hold on
loglog(f,tempmean2,'g','LineWidth', 2)

xlim([1 120])
%ylim([0.001 0.03])
ylim([0.001 0.03])
%yline(1,'--k','LineWidth',2)
xline(0.5,':k','LineWidth', 2,'LabelHorizontalAlignment','left')
xline(4,':k','LineWidth', 2, 'Label','Delta','LabelHorizontalAlignment','left')
xline(8,':k','LineWidth', 2, 'Label','Theta','LabelHorizontalAlignment','left')
xline(12,':k','LineWidth', 2,'Label','Alpha','LabelHorizontalAlignment','left')
xline(30,':k','LineWidth', 2,'Label','Beta','LabelHorizontalAlignment','left')
xline(60,':k','LineWidth', 2,'Label','LowGamma','LabelHorizontalAlignment','left')
xline(120,':k','LineWidth', 2,'Label','HighGamma','LabelHorizontalAlignment','left')
xlabel('Frequency(Hz)') 
ylabel('Post/Pre') 
%
xticks([0.5 4 8 12 30 60 120])
hold off
legend('Ratio','Location','northwest')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
grid on
