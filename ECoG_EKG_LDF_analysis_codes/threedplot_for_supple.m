

eeg20veh = squeeze(all_eeg_band(:,1,:));
ldf20veh = squeeze(all_Ldf(:,:,1));
hr20veh = squeeze(nHr_all(:,:,1));
p_gamma20veh = mean(eeg20veh,1);
mean_ldf20veh = mean(ldf20veh,2);
Hr20veh = mean(hr20veh,2);
x1 = movmean(p_gamma20veh,10000);%1s filter
x2 = medfilt1(mean_ldf20veh,10000);% 1s filter
x3 = movmean(Hr20veh,1000);% 100 msec filter

eeg20atro = squeeze(all_eeg_band(:,2,:));
ldf20atro = squeeze(all_Ldf(:,:,2));
hr20atro = squeeze(nHr_all(:,:,2));
p_gamma20atro = mean(eeg20atro,1);
mean_ldf20atro = mean(ldf20atro,2);
Hr20atro = mean(hr20atro,2);
y1 = movmean(p_gamma20atro,10000);
y2 = medfilt1(mean_ldf20atro,10000);
y3 = movmean(Hr20atro,1000);

eeg5veh = squeeze(all_eeg_band(:,3,:));
ldf5veh = squeeze(all_Ldf(:,:,3));
hr5veh = squeeze(nHr_5all(:,:,3));
p_gamma5veh = mean(eeg5veh,1);
mean_ldf5veh = mean(ldf5veh,2);
Hr5veh = mean(hr5veh,2);
z1 = movmean(p_gamma5veh,10000);
z2 = medfilt1(mean_ldf5veh,10000);
z3 = movmean(Hr5veh,1000);

eeg10veh = squeeze(all_eeg_band(:,4,:));
ldf10veh = squeeze(all_Ldf(:,:,4));
hr10veh = squeeze(nHr_10all(:,:,4));
p_gamma10veh = mean(eeg10veh,1);
mean_ldf10veh = mean(ldf10veh,2);
Hr10veh = mean(hr10veh,2);
w1 = movmean(p_gamma10veh,10000);
w2 = medfilt1(mean_ldf10veh,10000);
w3 = movmean(Hr10veh,1000);

st = 10;
vidObj = VideoWriter('Trajectories_color_post.AVI');
vidObj.FrameRate = 20;
open(vidObj);
Color_arr = [0.5 0.5 0.5;...
        1 0 0;...
        1 0.8 0.7;...
        0.9 0.4 0.4];
for i = st*Fs:0.5*Fs:90*Fs
    rot_fac = i/(10*Fs);
    figure(8)
    hold on;
    if i<30*Fs
        plot3(x1(st*Fs:i),x2(st*Fs:i),x3(st*Fs:i),'LineStyle',':','Color','k','LineWidth',3); hold on;
        plot3(y1(st*Fs:i),y2(st*Fs:i),y3(st*Fs:i),'LineStyle',':','Color','k','LineWidth',3); hold on;
        plot3(z1(st*Fs:i),z2(st*Fs:i),z3(st*Fs:i),'LineStyle',':','Color','k','LineWidth',3); hold on;
        plot3(w1(st*Fs:i),w2(st*Fs:i),w3(st*Fs:i),'LineStyle',':','Color','k','LineWidth',3); hold on;
        
    else
        title(strcat('Stim:',num2str((i/Fs)-30,'%0.1f'),'s'))
        plot3(x1(st*Fs:i),x2(st*Fs:i),x3(st*Fs:i),'Color',Color_arr(1,:),'LineWidth',3); hold on;
        plot3(y1(st*Fs:i),y2(st*Fs:i),y3(st*Fs:i),'Color',Color_arr(2,:),'LineWidth',3); hold on;
        plot3(z1(st*Fs:i),z2(st*Fs:i),z3(st*Fs:i),'Color',Color_arr(3,:),'LineWidth',3); hold on;
        plot3(w1(st*Fs:i),w2(st*Fs:i),w3(st*Fs:i),'Color',Color_arr(4,:),'LineWidth',3); hold on;
    end
    xlabel('\DeltaPower');
    ylabel('\DeltaLDF')
    zlabel('\DeltaHR')
    % hold off;
    set(gca,'FontSize', 32)
    grid on;
    xlim([0.4 1.6])
    ylim([-0.7 0.4])
    zlim([-1 0.2])
    set(gcf,'position',[0,0,1800,1300])
    view([-30/rot_fac -rot_fac 0.8*(0.5*rot_fac)]);% for specific view
    drawnow
   frame = getframe(gcf) ;
%     legend('','','','','20Hz (Vehicle)','20Hz')   
    writeVideo(vidObj, frame);
end
% close the object
close(vidObj);
winopen('Trajectories_color_post.avi')