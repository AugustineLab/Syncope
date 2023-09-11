%% Main code for analysing doppler-updated to the nof 7 till jan 11th 2023
%--Run as sections for different plots and goals (includes a section with chronux implementation for EEG analysis)
% Section 1: Loading section (load from the address or just load workspace....)
% Section 2: Process (normalization) and individual plots
% Section 2: NOTE: loops exist to compensate for animal specific variations
% inthe data
% Section 3: Syncope triggered avg plot
% Section 4: Stats for LDF parameters
% Section 5: Video for syncope triggered avg plot : 3d rotation plot
% Section 6: 3d plot and Video for stimulus triggered avg plot: for all 4 parameters(3
% diff stims and atropine):3d rotation plot
% Section 7: Avg LDF trace plotting section
% Section 8: ECOG processing using chronux (as per DK)



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

namesveh_HR = {'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Vehicle_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_20Hz_Veh_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_20Hz_Veh_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_20Hz_Veh_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_20Hz_Veh_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_20Hz_Veh_HRonly'};

namesatro_HR = {'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Atro_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_20Hz_Atro_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_20Hz_Atro_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_20Hz_Atro_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_20Hz_Atro_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_20Hz_Atro_HRonly_FIXED'};

namesveh_5_10hz_HR = {'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\216N_Vehicle_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\231R_5Hz_10Hz_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\236N_10Hz_5Hz_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\237R_10Hz_5Hz_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\239RL_10Hz_5Hz_HRonly',...
    'E:\Vinny_Jonny_doppler_collaboration\Matfiles_HR\240RR_10Hz_5Hz_HRonly'};

load latencies
%% -------------------------------------------------Section1: OR----------------
tic
load workspace
%% --------------------------------------Section2: Process and individual plots--------------------
Color_arr = jet(7);
fac = 1;% moving avg fac
Fs = 10000;
win_tim = 200;% seconds
bef = 10*Fs;
aft = 20*Fs;
groups = {'namesveh20hz','namesatro20hz','namesveh5_10hz'};
for g = 3:length(groups)
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
            all_Ldf(:,j,g)= nLDF_20hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\
            loc_min(j,g)=min(nLDF_20hz(30:60*Fs));
        else
            ldfplots_alln_5n10hz;
            all_Ldf(:,j,g)= nLDF_5hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\
            all_Ldf(:,j,g+1)= nLDF_10hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\
            loc_min(j,g)=min(nLDF_5hz(30:60*Fs));
            loc_min(j,g+1)=min(nLDF_10hz(30:60*Fs));
        end
        eeg_allplot;
        % for 20hz only use the latencies to shift
        if g<3
            temp_LDF = all_Ldf(:,j,g);
            temp_EEG = all_eeg_band(j,g,:);
            
            if (latent(j).latency)
                lat_zer = latent(j).latency;
            end
            shift_stim = (30+lat_zer)*Fs;% adding to the baseline 30 sec taken from the window....
            shift_LDF(:,j,g) = temp_LDF(shift_stim-bef:shift_stim+aft);
            shift_EEG(j,g,:) = temp_EEG(shift_stim-bef:shift_stim+aft);
        end
        if g<3
            if j~=2
                load(names_HR{counter});
                HR = data(:,1);
                if j==1
                    HR_20hz= HR(idx_stimS(3)-base_len:idx_stimE(3)+postim_len);
                else
                    HR_20hz= HR(idx_stimS(1)-base_len:idx_stimE(1)+postim_len);
                end
                base_hr = HR_20hz(1:base_len);
                nhr = (HR_20hz - mean(base_hr(:)))./(mean(base_hr(:)));
                shift_HR(:,counter,g) = nhr(shift_stim-bef:shift_stim+aft);
                
                Hr_all(:,counter)=HR_20hz(1:90*Fs);
                nHr_all(:,counter,g)=nhr(1:90*Fs);% normalsed HR data
                
                counter = counter+1;
            end
        else
            load(names_HR{counter});
            HR = data(:,1);
            if j<3
                HR_5hz= HR(idx_stimS(1)-base_len:idx_stimE(1)+postim_len);
                HR_10hz= HR(idx_stimS(2)-base_len:idx_stimE(2)+postim_len);
            else
                HR_5hz= HR(idx_stimS(2)-base_len:idx_stimE(2)+postim_len);
                HR_10hz= HR(idx_stimS(1)-base_len:idx_stimE(1)+postim_len);
            end
            base_5hr = HR_5hz(1:base_len);
            base_10hr = HR_10hz(1:base_len);
            
            nhr5 = (HR_5hz - mean(base_5hr(:)))./(mean(base_5hr(:)));
            nhr10 = (HR_10hz - mean(base_10hr(:)))./(mean(base_10hr(:)));
            
            Hr_5all(:,counter)=HR_5hz(1:90*Fs);
            nHr_5all(:,counter,g)=nhr5(1:90*Fs);% normalsed HR data
            Hr_10all(:,counter)=HR_10hz(1:90*Fs);
            nHr_10all(:,counter,g+1)=nhr10(1:90*Fs);% normalsed HR data
            
            counter = counter+1;
        end
        toc
    end
end
%% -------------------------------Section 3:for plots of syncope trig avg for comparing atropine and vehicle
tic
syncope_trig_avg;
%% -------------------------------- Section 3.1: Smoothen individual LDF animal trs

for i = 1:size(all_Ldf,3)
    for j = 1:size(all_Ldf,2)
        temp_ldf = all_Ldf(:,j,i);
        sm_temp = medfilt1(temp_ldf,2*Fs);
        figure(1);
        plot(sm_temp);
        hold on;
        %     pause
        sm_Ldf_2s_win(:,j,i)=sm_temp;
    end
end
%%
for i = 3:size(shift_LDF,3)
    for j = 1:size(shift_LDF,2)
        temp_ldf = shift_LDF(:,j,i);
        sm_temp = medfilt1(temp_ldf,1*Fs);
        figure(1);
        plot(sm_temp);
        hold on;
        %     pause
        syncope_Ldf_1s_win(:,j,i)=sm_temp;
        
        temp_eeg = squeeze(shift_EEG(j,i,:));
        sm_temp1 = medfilt1(temp_eeg,1*Fs);
        figure(2);
        plot(sm_temp1);
        hold on;
        %     pause
        syncope_EEG_1s_win(j,i,:)=sm_temp1;
        
        temp_hr = shift_HR(:,j,i);
        sm_temp2 = medfilt1(temp_hr,1*Fs);
        figure(3);
        plot(sm_temp2);
        hold on;
        %     pause
        syncope_HR_1s_win(:,j,i)=sm_temp2;
    end
end
%% --------------------------Section 4: Stats------------------------------numbers for stats for all group comparisons
for g = 1:3
    g
    %     close all
    grp_name = eval(groups{g});
    for j = 1:length(grp_name)
        if g<3
            nLDF_20hz=all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
            stim_ep = movmean(nLDF_20hz(30*Fs:60*Fs),5000);
            int_tim = 30*Fs:60*Fs;
            
            all_min = find(stim_ep == min(stim_ep));
            loc_min(j,g)= all_min(1);
            val_min(j,g)= min(stim_ep);
            
            base = mean(nLDF_20hz(1:30*Fs));
            % Find the half max value.
            thresh = 0.5*val_min(j,g);%-2.*std(nLDF_20hz(1:30*Fs));
            z_stim = thresh.*ones(size(stim_ep,1),size(stim_ep,2));
            z_stim(stim_ep<thresh)=stim_ep(stim_ep<thresh);
            in_zstim = -z_stim;% so that the area is positive
            auc(j,g) = trapz(int_tim,in_zstim);% area under the curve
            idx = find(in_zstim>-thresh);
            if (idx)
                tt = numel(idx)./Fs;
                mean_transit_time(j,g)=tt;
            end            
            [pt_inflex, drop_rate, rep_rate]=slpparcal(z_stim);
            drop_rate_all(j,g)=drop_rate;
            rep_rate_all(j,g) = rep_rate;
            rate_coeff(j,g)=drop_rate/rep_rate;            
            figure()
            plot(stim_ep); hold on;
            %             xline(loc_min(j,g)); yline(val_min(j,g));
            plot(z_stim);
            yline(thresh)
            if (idx)
                xline(idx(1)); xline(idx(end))
            end
            hold off;
            %             pause
        else
            nLDF_5hz= all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
            nLDF_10hz=all_Ldf(:,j,g+1);% taking 90 seconds by force to keep variable size incheck\
            stim_ep5 = movmean(nLDF_5hz(30*Fs:60*Fs),5000);
            stim_ep10 = movmean(nLDF_10hz(30*Fs:60*Fs),5000) ;
            
            all_min1 = find(stim_ep5 == min(stim_ep5));
            all_min2 = find(stim_ep10 == min(stim_ep10));
            
            loc_min(j,g)= all_min1(1);
            loc_min(j,g+1)= all_min2(1);
            
            val_min(j,g)= min(stim_ep5);
            val_min(j,g+1)= min(stim_ep10);
            
            int_tim = 30*Fs:60*Fs;
            base5 = mean(nLDF_5hz(1:30*Fs));
            thresh5 = 0.5*val_min(j,g);%-2*std(nLDF_5hz(1:30*Fs));%(min(in_zstim5) + max(in_zstim5)) / 2
            z_stim5 = thresh5.*ones(size(stim_ep5,1),size(stim_ep5,2));
            z_stim5(stim_ep5<thresh5)=stim_ep5(stim_ep5<thresh5);
            in_zstim5 = -z_stim5;
            auc(j,g) = trapz(int_tim,in_zstim5);% area under the curve
            idx5 = find(in_zstim5>-thresh5);
            
            if (idx5)
                tt5 = numel(idx5)./Fs;
                mean_transit_time(j,g)=tt5;
            end
            base10 = mean(nLDF_10hz(1:30*Fs));
            thresh10 = 0.5*val_min(j,g+1);%-2*std(nLDF_10hz(1:30*Fs));%(min(in_zstim10) + max(in_zstim10)) / 2
            z_stim10 = thresh10.*ones(size(stim_ep10,1),size(stim_ep10,2));
            z_stim10(stim_ep10<thresh10)=stim_ep10(stim_ep10<thresh10);
            in_zstim10 = -z_stim10;
            auc(j,g+1) = trapz(int_tim,in_zstim10);% area under the curve
            idx10 = find(in_zstim10>-thresh10);
            
            if (idx10)
                tt10 = numel(idx10)./Fs;
                mean_transit_time(j,g+1)=tt10;
            end
                [pt_inflex5, drop_rate5, rep_rate5]=slpparcal(z_stim5);
                drop_rate_all(j,g)=drop_rate5;
                rep_rate_all(j,g) = rep_rate5;
                rate_coeff(j,g)=drop_rate5/rep_rate5;
                [pt_inflex10, drop_rate10, rep_rate10]=slpparcal(z_stim10);
                drop_rate_all(j,g+1)=drop_rate10;
                rep_rate_all(j,g+1) = rep_rate10;
                rate_coeff(j,g+1)=drop_rate10/rep_rate10;
                              
                figure()                              
                subplot(211)
                plot(stim_ep5); hold on;
                %             plot(loc_min(j,g), val_min(j,g),'*');
                plot(z_stim5)
                yline(thresh5)
                if (idx5)
                    xline(idx5(1)); xline(idx5(end))
                end
                hold off;
                subplot(212)
                plot(stim_ep10); hold on;
                %             plot(loc_min(j,g+1), val_min(j,g+1),'*');
                plot(z_stim10)
                yline(thresh10)
                if (idx10)
                    xline(idx10(1)); xline(idx10(end))
                end
                hold off;
                %             pause
            end
            
        end
    end
    %%
    % close all
    figure();
    
%     plot([0:5],4.*ones(1,6),'--','Color',[0 0 0])
%     hold on;
    s1= scatter(val_min(:,1),reperfusion_time(:,1),'filled');
    s1.SizeData=100;
    s1.MarkerFaceColor = [0.5 0.5 0.5];
    hold on;
    
    s2 = scatter(val_min(:,2),reperfusion_time(:,2),'filled');
    s2.SizeData=100;
    s2.MarkerFaceColor = [1 0 0];
    
    s3 = scatter(val_min(:,3),reperfusion_time(:,3),'filled');
    s3.SizeData=100;
    s3.MarkerFaceColor = [1 0.8 0.7];%[1 0.8 0.65];
    
    s4 = scatter(val_min(:,4),reperfusion_time(:,4),'filled');
    s4.SizeData=100;
    s4.MarkerFaceColor = [0.9 0.4 0.4];%[0.9 0.4 0.4];
    
    set(gca,'FontSize',30)
    xlabel('Max decrease in LDF (a.u.)')
    ylabel('LDF Recovery (s)')
    box on;
    %%
    %for stats...
    % plot graph with error bars
    % fractional drop in LDF
    x1 = [1:4];
    y2 = mean(val_min,1);
    y1a = val_min(:,1);
    y2a =  val_min(:,2);
    y3a = val_min(:,3) ;
    y4a = val_min(:,4) ;
    
    figure(1);
    scatter(x1,y2,200,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
    hold on;
    scatter(x1(1),y1a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0])
    scatter(x1(2),y2a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0]);
    scatter(x1(3),y3a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1]);
    scatter(x1(4),y3a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.8 0 0.8]);
    set(gca,'xlim',[0.75 4.25],'FontSize',28);
    box off;
    xticks([1 2 3 4])
    xticklabels({'20Hz(Veh.)', '20Hz(Atro.)','5Hz(Veh.)','10Hz(Veh.)'})
    ylabel('Mean \DeltaLDF')
    
    % latency to minimal drop during stim in LDF
    x1 = [1:4];
    y2 = mean(loc_min./Fs,1);
    y1a = loc_min(:,1)./Fs;
    y2a =  loc_min(:,2)./Fs;
    y3a = loc_min(:,3)./Fs ;
    y4a = loc_min(:,4)./Fs ;
    
    figure(2);
    scatter(x1,y2,200,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
    hold on;
    scatter(x1(1),y1a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0])
    scatter(x1(2),y2a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0]);
    scatter(x1(3),y3a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1]);
    scatter(x1(4),y3a,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.8 0 0.8]);
    set(gca,'xlim',[0.75 4.25],'FontSize',26);
    box off;
    xticks([1 2 3 4])
    xticklabels({'20Hz(Veh.)', '20Hz(Atro.)','5Hz(Veh.)','10Hz(Veh.)'})
    ylabel('Mean Duration (s)')
    title('Time to reach LDF minima')
    %% ------------------Section 5: VIDEO:for SYNCOPE TRIGGERED avg 3d line plot for ECOg, LDF and HR-------------------
    %
    videoObj = VideoWriter('Syncope_trig_veh_vs_atro.AVI');
    videoObj.FrameRate = 2;
    open(videoObj);
    for ang = 1:20
        increment = 5*Fs; %sec for the color gradient change on the 3dplot
        for grp = 1:2% two grps veh and atropine
            eeg = squeeze(shift_EEG(:,grp,:));
            ldf = squeeze(shift_LDF(:,:,grp));
            hr = squeeze(shift_HR(:,:,grp));
            colorcode = grp;
            [~, ~, ~] = my3dplot(eeg,ldf,hr,colorcode,increment,Fs,ang);
        end
        drawnow
        frame = getframe(gcf) ;
        %     legend('','','','','20Hz (Vehicle)','20Hz')
        writeVideo(videoObj, frame);
    end
    close(videoObj);
    
    %% ----------------------------Section 6: 3dplot and VIDEO: for Stimulus TRIGGERED avg 3d line plot for ECOg, LDF and HR
    %
    % increment = 10*Fs; %sec for the color gradient change on the 3dplot
    dwlen =900;
    increment = 10*Fs; %sec for the color gradient change on the 3dplot
    count = 1;
    for grp = 1:4% 4 grps veh and atropine
        eeg = squeeze(all_eeg_band(:,grp,:));
        ldf = squeeze(all_Ldf(:,:,grp));
        if grp<3
            hr = squeeze(nHr_all(:,:,grp));
        elseif grp ==3
            hr = squeeze(nHr_5all(:,:,grp));
        else
            hr = squeeze(nHr_10all(:,:,grp));
        end
        
        
        
        colorcode = grp;
        [~, ~, ~] = my3dplot(eeg,ldf,hr,count,increment,Fs);
        %     pause
        count = count+1;
        
    end
    
    %%
    threedplot_for_supple;
    %% --------------------------------------Section7: Avg line plots: plot smoothened mean for different groups and different stimulation frequencies
    Color_arr = [0.33 0.33 0.34;...
        1 0 0;...
        0.95 0.8 0.8;...
        0.9 0.4 0.4];
    time1 = downsample(time,900);
    for g = 1:size(all_Ldf,3)
        
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
    end
    %% ---------------------------------------Section 8:EEG processing using Chronux for David...
    
    % remove DC and normalize and plot and save
    %     spectrumplots;
    %     spectrogram_plots;
    Sall_j{j}=Sallstim;
    fall_j{j}=fallstim;
    EEG_20hz= EEG_sig(idx_stimS(3)-base_len:stim_end(3)+postim_len);
    
    params_EEG;
    for grp = 2:length(EEG_20hz)
        h(grp)= EEG_20hz(grp)-EEG_20hz(grp-1);
    end
    
    [Seeg,t5,feeg]=mtspecgramc(h,window,paramsSpecgramEEG);% spectrogram for whole signal
    figure(18);
    subplot(211)
    imagesc(t5,feeg,log10(Seeg)'); axis xy;
    colorbar;
    %     c.Limits = [-3.5 -0.5];
    colormap parula
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(strcat('LFP Power spectrum'))
    set(gca,'FontSize', 24);
    temp = LDF_20hz;
    base = temp(1:10*Fs);
    mbase = mean(temp(:));
    ntemp = (temp-mbase)./mbase;
    
    subplot(212)
    plot(time,ntemp);
    xlim([0 50])
    xlabel('Time (s)')
    ylabel('\DeltaLDF')
    set(gca,'FontSize', 24);
    % extract ECOG band for the 3d plot
    band= Seeg(:,feeg>39&feeg<81);
    p_gamma = sum(band,2);
    
    HR_20hz= HR(stim_start(3)-base_len:stim_end(3)+postim_len);
    
    % end
    toc
    % %% 3d line plot for ECOg, LDF and HR
    % % need to downsample and shorten the duration for the 3d plot and color
    % % code artificially
    %
    % increment = 5*Fs; %sec
    % nonstim = 10;% sec for pre and post stim taken for the 3d plot
    % stim_dur = 30;% sec
    % dwsamp = 1000;
    % % p_gamma = mean(shift_EEG,1);
    % % mean_ldf20 = mean(shift_LDF,2);
    % % nLDF_20hz = medfilt1(mean_ldf20,5000);% med filt smooth
    %
    % p_gamma = mean(all_eeg_band,1);
    % mean_ldf20 = mean(all_Ldf,2);
    % nLDF_20hz = medfilt1(mean_ldf20,5000);% med filt smooth
    %
    % Hr = mean(Hr_all,2);
    %
    % % n=7;
    % % p_gamma = all_eeg_band(n,:);
    % %  nLDF_20hz = all_Ldf20(:,n);
    % %  Hr = Hr_all(:,n);
    % %  x1 = downsample(p_gamma, floor(length(p_gamma)/dwsamp));
    % %  x2 = downsample(nLDF_20hz,floor(length(nLDF_20hz)/dwsamp));
    % %  x3 = downsample(Hr,floor(length(Hr)/dwsamp));
    %
    %
    % x1 = movmean(p_gamma,10000);
    % x2 = movmean(nLDF_20hz,1000);
    % x3 = movmean(Hr,1000);
    % %
    % %  Color_arr = [0.7 0 0;...
    % %      0.9 0.3 0.1;...
    % %      0.9 0.9 0;...
    % %      0 0.8 0;...
    % %      0 0 1;...
    % %      0.5 0.2 0.5];
    %
    %
    % Color_arr = [0 0.2 0;...
    %     0 0.4 0;...
    %     0 0.6 0;...
    %     0 0.8 0];
    %
    % grp = 10*Fs;
    % figure(8)
    % hold on;
    % plot3(x1(1*Fs:10*Fs),x2(1*Fs:10*Fs),x3(1*Fs:10*Fs),'LineStyle',':','Color','k','LineWidth',3); hold on;
    % count = 1;
    % while grp <30*Fs
    %     plot3(x1(grp:grp+increment),x2(grp:grp+increment),x3(grp:grp+increment),'Color',Color_arr(count,:),'LineWidth',3); hold on;
    %     grp = grp+increment;
    %     count = count+1;
    % end
    % % plot3(x1(60*Fs:70*Fs),x2(60*Fs:70*Fs),x3(60*Fs:70*Fs),'LineStyle',':','Color','k','LineWidth',3); hold off;
    %  plot3(x1(10*Fs),x2(10*Fs),x3(10*Fs),'bo','LineWidth',10)
    % xlabel('Avg EEG power');
    % ylabel('\DeltaLDF')
    % zlabel('Heart rate (bpm)')
    % % hold off;
    % set(gca,'FontSize', 32)
    %
    % %%
