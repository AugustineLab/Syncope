%% make file list from exported .mat files in single directory
testfiledir = uigetdir();
matfiles = dir(fullfile(testfiledir, '*.mat'));
nfiles = length(matfiles);

%% loops for pupil and whisking
% raw data
for w = 1 : nfiles
    fid = fullfile(testfiledir, matfiles(w).name);
    load(fid)
    B = 1910-size(pupil_area,2);
    C = zeros(1,B);
    pupil_area = [pupil_area C];
    whisk = [whisk C];
    pupil(w,:) = pupil_area(1,:);
    whiskall(w,:) = whisk(1,:);
end
% figure
% pupilmean = mean(pupil,1,'omitnan');
% plot(pupilmean)
% figure
% whiskmean = mean(whiskall,1,'omitnan');
% plot(whiskmean)
% figure
% pupilsmooth = movmean(pupilmean,10,2,'omitnan');
% plot(pupilsmooth)
% figure
% whisksmooth = movmean(whiskmean,15,2,'omitnan');
% plot(whisksmooth)
%% Baseline normalization
pupilBL = mean(pupil(:,1:900),[1 2],'omitnan'); %normalize to group mean
whiskBL = mean(whiskall(:,1:900),[1 2],'omitnan');
for w = 1 : size(pupil,1)
    if isnan(any(pupil(w)))==0        
        %pupilBL = mean(pupil(w,1:900),2,'omitnan'); %normalize to single
        pupiln(w,:) = pupil(w,:)./pupilBL;
    end
    %whiskBL = mean(whiskall(w,1:900),2,'omitnan');
    whisks = movmean(whiskall(w,:),15,2,'omitnan');
    whiskn(w,:) = whisks./whiskBL;
    %whiskn(w,:) = whiskall(w,:)./whiskBL;
end
figure
pupilmeann = mean(pupiln,1,'omitnan');
plot(pupilmeann)
xlim([870 960])
figure
whiskmeann = mean(whiskn,1,'omitnan');
plot(whiskmeann)
xlim([870 960])
figure
pupilsmoothn = movmean(pupilmeann,10,2,'omitnan');
plot(pupilsmoothn)
xlim([870 960])
figure
whisksmoothn = movmean(whiskmeann,15,2,'omitnan');
plot(whisksmoothn)
xlim([870 960])