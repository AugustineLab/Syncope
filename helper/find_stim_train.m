function Trains = find_stim_train(digArray, srate, LowHz)
if isa(digArray, 'logical')
else
    digArray = digArray>1;
end
digArrayON = find(diff(digArray)==1)+1; % find the oneset of each '1' after '0'. like 0'1'110'1'11 (rising edge)
digArrayONs = digArrayON./srate; % convert to seconds
intervalSec = round(diff(digArrayONs),3); % round to nearest 1ms
intervalSamp = diff(digArrayON); % sample between each stim
HzStim = round(1./intervalSec); % calculate HZ of each stim
Trains = [0 find(HzStim<LowHz)']; % row1 = find the first stim of each train
Trains = [Trains; HzStim(Trains+1)']; % row2 = Hz of each train
Trains = [Trains; digArrayON(Trains(1,:)+1)']; % row3 = time (samp) for each train start
Trains = [Trains; (digArrayON(Trains(1,2:size(Trains,2)))+intervalSamp(Trains(1,2:size(Trains,2))-1))' digArrayON(end)+intervalSamp(end)']; % row4 = time (samp) for each train end
Trains = [Trains; ((Trains(4,:)-Trains(3,:))./srate)]; % row5 = duration of trains (sec)
end
