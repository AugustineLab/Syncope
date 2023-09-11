
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

%%
Color_arr = jet(7);
fac = 1;% moving avg fac
Fs = 10000;
win_tim = 200;% seconds
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
        all_Ldf(:,j,g)= nLDF_20hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\
        loc_min(j,g)=min(nLDF_20hz(30:60*Fs));
    else
        ldfplots_alln_5n10hz;
        all_Ldf(:,j,g)= nLDF_5hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\        
        all_Ldf(:,j,g+1)= nLDF_10hz(1:200*Fs);% taking 90 seconds by force to keep variable size incheck\
        loc_min(j,g)=min(nLDF_5hz(30:60*Fs));
        loc_min(j,g+1)=min(nLDF_10hz(30:60*Fs));
    end
    % for 20hz only use the latencies to shift
    eegplot_for_baselinerecovery;

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
                % shift_HR(:,counter,g) = nhr(shift_stim-bef:shift_stim+aft);
                
                Hr_all(:,counter)=HR_20hz(1:200*Fs);
                nHr_all(:,counter,g)=nhr(1:200*Fs);% normalsed HR data
                
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
            
            Hr_5all(:,counter)=HR_5hz(1:200*Fs);
            nHr_5all(:,counter,g)=nhr5(1:200*Fs);% normalsed HR data
            Hr_10all(:,counter)=HR_10hz(1:200*Fs);
            nHr_10all(:,counter,g+1)=nhr10(1:200*Fs);% normalsed HR data
            
            counter = counter+1;
        end
    toc
end
end
%%
for g = 1:3
     g
    close all 
    grp_name = eval(groups{g});  
    for j = 1:length(grp_name)
        if g<3
            nLDF_20hz=all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
%             poststim_ep = movmean(nLDF_20hz(30*Fs:60*Fs),5000);
            poststim_ep = movmean(nLDF_20hz(60*Fs:150*Fs),5000);
           
            % Find the half max value.
            thresh = -1.*std(nLDF_20hz(1:30*Fs));
            z_stim = thresh.*ones(size(poststim_ep,1),size(poststim_ep,2));
            z_stim(poststim_ep<thresh)=poststim_ep(poststim_ep<thresh);
            in_zstim = -z_stim;% so that the area is positive
            auc(j,g) = trapz(int_tim,in_zstim);% area under the curve
            idx = find(in_zstim>-thresh);
            if (idx)   
                tt = numel(idx)./Fs;
                mean_transit_time(j,g)=tt;
            end
            all_min = find(poststim_ep == min(poststim_ep));
            loc_min(j,g)= all_min(1);
            val_min(j,g)= min(poststim_ep);
          
            
            
            figure(1)
            plot(poststim_ep); hold on; 
%             xline(loc_min(j,g)); yline(val_min(j,g)); 
            plot(z_stim); 
            yline(thresh)
            if (idx)
                xline(idx(1)); xline(idx(end))
            end
            hold off;
            pause
        else           
            nLDF_5hz= all_Ldf(:,j,g);% taking 90 seconds by force to keep variable size incheck\
            nLDF_10hz=all_Ldf(:,j,g+1);% taking 90 seconds by force to keep variable size incheck\
            stim_ep5 = movmean(nLDF_5hz(30*Fs:60*Fs),5000);
            stim_ep10 = movmean(nLDF_10hz(30*Fs:60*Fs),5000) ;
            int_tim = 30*Fs:60*Fs;
            base5 = mean(nLDF_5hz(1:30*Fs));
            thresh5 = -1*std(nLDF_5hz(1:30*Fs));%(min(in_zstim5) + max(in_zstim5)) / 2
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
            thresh10 = -1*std(nLDF_10hz(1:30*Fs));%(min(in_zstim10) + max(in_zstim10)) / 2
            z_stim10 = thresh10.*ones(size(stim_ep10,1),size(stim_ep10,2));
            z_stim10(stim_ep10<thresh10)=stim_ep10(stim_ep10<thresh10);
            in_zstim10 = -z_stim10;
            auc(j,g+1) = trapz(int_tim,in_zstim10);% area under the curve 
            idx10 = find(in_zstim10>-thresh10);

            if (idx10)
                tt10 = numel(idx10)./Fs;
                mean_transit_time(j,g+1)=tt10;
            end
            
            all_min1 = find(stim_ep5 == min(stim_ep5));
            all_min2 = find(stim_ep10 == min(stim_ep10));

            loc_min(j,g)= all_min1(1);
            loc_min(j,g+1)= all_min2(1);

            val_min(j,g)= min(stim_ep5);
            val_min(j,g+1)= min(stim_ep10);
           
            figure(2)
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
            pause
        end
        
    end
end