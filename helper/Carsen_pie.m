%% pie charts from Carsen data
    labels = ['decreased'; 'all other'];
    regionlabels = {'decreased','all other'}; %EDIT?
    %decreasedN = [21,129,15,0,197,8,5,82,302,76,51,157,124,260,490,322,97,31,210,320,24,29,10,69,399,450]; %EDITNEW
    decreasedN = [30,172,20,1,225,19,6,143,346,81,50,201,105,261,494,306,96,33,305,428,24,52,10,226,448,579]; %EDITNEW
    %otherN = [17,131,9,12,77,18,13,86,123,43,11,146,43,25,124,65,31,19,192,251,9,50,21,211,184,404];
    otherN = [8,88,4,11,49,7,12,25,79,38,12,102,62,24,120,81,32,17,97,143,8,27,21,54,135,275];
    totalN = [decreasedN; otherN];
    figure
    subplot(size(totalN,2),1,1)
    labels = {'',''};
    for i = 1:size(totalN,2)
        subplot(size(totalN,2),1,i)
        pie(totalN(:,i),labels);
    end
