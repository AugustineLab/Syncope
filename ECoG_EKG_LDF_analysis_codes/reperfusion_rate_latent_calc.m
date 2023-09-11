
% Sampling rate = 10kHz
% ECG = Ch1
% EEG = Ch2
% Laser = Ch3
%
% Piezo = Ch4
%
% LDF = Ch5
% HR BPM = Ch6
% Cam TTL = Ch7
clearvars
close all
clc

%------------------------Section1: Loading section-----------------------------------
tic
namesveh20hz = {'216N_awake_headfix_5Hz_10Hz_20Hz_Vehicle',...
    '203L_20Hz_take2_Veh',...
    '231R_20Hz_Veh',...
    '236N_20Hz_Veh',...
    '237R_20Hz_Veh',...
    '239RL_20Hz_Veh',...
    '240RR_20Hz_Veh';};

namesatro20hz = {'216N_awake_headfix_5Hz_10Hz_20Hz_Atropine',...
    '203L_20Hz_take2_Atro',...
    '231R_20Hz_Atro',...
    '236N_20Hz_Atro',...
    '237R_20Hz_Atro',...
    '239RL_20Hz_Atro',...
    '240RR_20Hz_Atro';};

namesveh5_10hz = {'216N_awake_headfix_5Hz_10Hz_20Hz_Vehicle',...
    '203L_10Hz_5Hz_take2_Veh',...
    '231R_10Hz_5Hz_Veh',...
    '236N_10Hz_5Hz_Veh',...
    '237R_10Hz_5Hz_Veh',...
    '239RL_10Hz_5Hz_Veh',...
    '240RR_10Hz_5Hz_Veh';};

namesveh_HR = {'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Vehicle_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_20Hz_Veh_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_20Hz_Veh_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_20Hz_Veh_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_20Hz_Veh_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_20Hz_Veh_HRonly'};

namesatro_HR = {'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Atro_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_20Hz_Atro_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_20Hz_Atro_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_20Hz_Atro_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_20Hz_Atro_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_20Hz_Atro_HRonly_FIXED'};

namesveh_5_10hz_HR = {'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Vehicle_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_5Hz_10Hz_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_10Hz_5Hz_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_10Hz_5Hz_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_10Hz_5Hz_HRonly',...
    'D:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_10Hz_5Hz_HRonly'};

load latencies

%% --------------------------------------Section2: Process and individual plots--------------------
Color_arr = jet(7);
fac = 1;% moving avg fac
Fs = 10000;
win_tim = 270;% seconds
bef = 10*Fs;
aft = 20*Fs;
groups = {'namesveh20hz','namesatro20hz','namesveh5_10hz'};
for g = 1:length(groups)
    g
    close all
    grp_name = eval(groups{g});
    % switch case for variable for plotting...
    if g == 1 % for vehicle data
        latent = laten_veh;
        names_HR = namesveh_HR;
    elseif g == 2 % for atropine data
        latent = laten_atro;
        names_HR = namesatro_HR;
    else
        names_HR = namesveh_5_10hz_HR;
    end
    counter =1;
    % get data for all and store the data
    
    for j = 1:length(grp_name)
        load(grp_name{j});
        ldf_sig = data(:,5);
        EEG_sig = data(:,2);
        Stim = data(:,3);
        ECG = data(:,2);
        HR = data(:,7);% data 7 for atleast 216N from observation not from the info given
        %convert stim to square pulse and get start and end times
        %         [stim_start, stim_end]= getStims(Stim,Fs); my stim calc
        TTL = Stim>1;
        LFPsr = 10000;
        TrainInfo = find_stim_train(TTL,LFPsr,4); % Get train timestamps (custom function)
        % get the window idx from each stimulation
        nwithinbin = 1;
        WsecS = -30; % window in seconds before stim start EDIT
        WsecE = 30; % window in seconds after stim ends EDIT
        idx_stimS = zeros(1,size(TrainInfo,2));
        idx_stimE = zeros(1,size(TrainInfo,2));
        idx_bsS = zeros(1,size(TrainInfo,2));
        idx_asE = zeros(1,size(TrainInfo,2));
        for n = 1:size(TrainInfo,2)
            idx_stimS(n) = round(TrainInfo(3,n))/nwithinbin; % index of each stim start time
            idx_stimE(n) = round(TrainInfo(4,n))/nwithinbin; % index of each stim end time
            idx_bsS(n) = ((idx_stimS(n)/LFPsr)+(WsecS))*LFPsr/nwithinbin;% index before each stim (window start)
            idx_asE(n) = ((idx_stimE(n)/LFPsr)+(WsecE))*LFPsr/nwithinbin;% index after each stim end (window end)
        end
        % processing as per groups
        
        if g<3
            % preprocess ldf
            ldf_sig1 = detrend(ldf_sig);
            ldf_sig2 = ldf_sig1-mean(ldf_sig1);
            lp_ldf = lowpass(ldf_sig2,1,Fs);
            %         ldf_plots;  % for different freq....this gets the normalized ldfs with specific time windows.. look into the code for window spec
            ldfplots_alln;
            all_Ldf(:,j,g)= nLDF_20hz(1:270*Fs);% taking 90 seconds by force to keep variable size incheck\
            loc_min(j,g)=min(nLDF_20hz(30:60*Fs));
        else
            ldfplots_alln_5n10hz;
            all_Ldf(:,j,g)= nLDF_5hz(1:270*Fs);% taking 90 seconds by force to keep variable size incheck\
            all_Ldf(:,j,g+1)= nLDF_10hz(1:270*Fs);% taking 90 seconds by force to keep variable size incheck\
            loc_min(j,g)=min(nLDF_5hz(30:60*Fs));
            loc_min(j,g+1)=min(nLDF_10hz(30:60*Fs));
        end
       
        toc
    end
end
%% --- reperfusion rate caluculation post stim---
for g = 1:3
    g
    grp_name = eval(groups{g});
    for j = 1:length(grp_name)
        if g<3
            nLDF_20hz=all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
            poststim_ep = movmean(nLDF_20hz(60*Fs:270*Fs),100000);% 100k is random smoothing factor
            slp_post = diff(movmean(poststim_ep,50000));
            inflx_ids = find(slp_post<0); %points slope stops to increase as in reperfusion stops
            pt_inflx = inflx_ids(1);
            std_poststim_ep = 1*std(nLDF_20hz(1:30*Fs),[],1);% taking 2*std to see when the reperfusion reaches 2*std of baseline
            idx2 = find(poststim_ep(pt_inflx:end)<=std_poststim_ep);% time at which the slope stops to increase
            late_ldf = nLDF_20hz(end-100*Fs:end); % taking last 100s
            mean_late_ldf(j,g) = mean(late_ldf(:));
            if (idx2)
                rest_point(j,g) = (idx2(1)+pt_inflx);
            else
                rest_point(j,g)= 210*Fs;
            end
            
            figure, plot(poststim_ep); hold on; xline(rest_point(j,g));
            hold off;
            pause
        else
            nLDF_5hz= all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
            nLDF_10hz=all_Ldf(:,j,g+1);% taking 90 seconds by force to keep variable size incheck\
            poststim_ep5 = movmean(nLDF_5hz(60*Fs:270*Fs),10000);
            poststim_ep10 = movmean(nLDF_10hz(60*Fs:270*Fs),10000) ;
            
            slp_post5 = diff(movmean(poststim_ep5,50000));
            inflx_ids5 = find(slp_post5<0); %points slope stops to increase as in reperfusion stops
            pt_inflx5 = inflx_ids5(1);
            std_poststim_ep5 = 1*std(nLDF_5hz(1:30*Fs),[],1);% taking 2*std to see when the reperfusion reaches 2*std of baseline
            idx25 = find(poststim_ep5(pt_inflx5:end)<=std_poststim_ep5);
            late_ldf5 = nLDF_5hz(end-100*Fs:end); % taking last 100s
            mean_late_ldf(j,g) = mean(late_ldf5(:));
            if (idx25)
                rest_point(j,g) = (idx25(1)+pt_inflx5);
            else
              rest_point(j,g)= 210*Fs;  
            end
            figure, plot(poststim_ep5); hold on; xline(rest_point(j,g));
            hold off;
            pause
            
            slp_post10 = diff(movmean(poststim_ep10,50000));
            inflx_ids10 = find(slp_post10<0); %points slope stops to increase as in reperfusion stops
            pt_inflx10 = inflx_ids10(1);
            std_poststim_ep10 = 1*std(nLDF_10hz(1:30*Fs),[],1);% taking 2*std to see when the reperfusion reaches 2*std of baseline
            idx210 = find(poststim_ep10(pt_inflx10:end)<=std_poststim_ep10);
            late_ldf10 = nLDF_10hz(end-100*Fs:end); % taking last 100s
            mean_late_ldf(j,g+1) = mean(late_ldf10(:));
            if (idx210)
                rest_point(j,g+1) = (idx210(1)+pt_inflx10);
            else
                rest_point(j,g+1)= 210*Fs;
            end
            
            figure, plot(poststim_ep10); hold on; xline(rest_point(j,g+1));
            hold off;
            pause
%             rep_time_all(j,g)=rep_time5./Fs;
%             rep_time_all(j,g+1)=rep_time10./Fs;
%             rep_rate_all(j,g)=mean(slp_post5(1:rep_time5));
%             rep_rate_all(j,g+1)=mean(slp_post10(1:rep_time10));
        end
    end
end
%% --------------------------------------Section7: Avg line plots: plot smoothened mean for different groups and different stimulation frequencies
    Color_arr = [0.33 0.33 0.34;...
        1 0 0;...
        0.95 0.8 0.8;...
        0.9 0.4 0.4];
    time2 = [1:size(all_Ldf,1)]./Fs;
   time = time2-60;
    time1 = downsample(time,900);
    xticks = [-60 -30 0 30 60 90 120 150 180 210];
    for g = 1:4;%size(all_Ldf,3)
        
        samp_all = all_Ldf(:,:,g);
        mean_ldf = mean(samp_all,2);
        sem_ldf = std(samp_all,[],2)./sqrt(7);
        up_ldf = mean_ldf + sem_ldf;
        down_ldf = mean_ldf -sem_ldf;
        sm_mean_ldf = downsample(medfilt1(mean_ldf,5000),900);
        sm_up = downsample(medfilt1(up_ldf,5000),900);
        sm_down = downsample(medfilt1(down_ldf,5000),900);
        figure(1);
        fill([time1 fliplr(time1)], [sm_up' fliplr(sm_down')], Color_arr(g,:), 'linestyle', 'none');
        alpha(0.8)
        toc
        hold on
        plot(time1,sm_mean_ldf','Color',[0 0 0])
        %     xlim([min(time1) max(time1)])
        % ylim([-50 70])
        set(gca,'FontSize', 24);
        toc
        xlabel('Time (s)')
        ylabel('\DeltaLDF')
        box off
        title('')
        xlim([-60 210])
   ax = gca;
%     ax.YMinorTick = 'on';
    ax.XAxis.TickValues = xticks;
    end
