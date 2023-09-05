Hello Everyone! Here is a brief usage of each function or matfile. :) 
Feel free to contact me if you have any question! 

----------------------------------------------------
List: 

readcagemouse.m 
RR_main.m
RR_preprocess.m 
find_stim_train2.m
load_data2.m
HR_analysis_old.m
HR_optimized.m 
correlation_test.m
---------------------------------------------------
Description: 

readcagemouse.m : this is a function that can call and trace the mouse information 
[cage, mouse, day,group,fr] = readcagemouse(n)

RR_main.m: the main script to call for other functions and postprocess breathing data

RR_preprocess.mï¼? key function to preprocess the breathing data, including peak detection and RR calculations. 

find_stim_train2.m: adapted from Jonny, it can automatically identify the frequency, strart point and duration of each stimulus

load_data2.m: load_data2.mat is a main script. 1. create a cell to save all features. 2. Load blood pressure data. 3. Load skin temperature data. 4. Load ECG data.

HR_analysis_old.m: the old version to extract HR. Works in previous data. Before using this function, you should first preprocessed the ECG signals in the Acknowledge, including FIR filtering, find rate, and smoothing. Second, you should name the ouput channel as 'HR'.  Third, you need to find the first cycle and mark the onset of TTL. Finally, you need to paste each cycle in a seperate window that contains all mice data.

HR_optimized.m: Optimized function to extract HR from the Acknowledge. You can call this function in the load_data2. mat. Before using this function, you should first preprocessed the ECG signals in the Acknowledge, including FIR filtering, find rate, and smoothing. Second, you should name the ouput channel as 'HR'.  There is no need to find the first cycle and mark the onset of TTL. 

correlation_test.m: First, this script is used to incoprate all parameters between NPY2R Ai32 and Ai9. In this script, it will form a matrix, we can see a series of measures of each mouse, such as BP change, HR change, apnea ..... Second, this script can be used to see the cross-correlation between breathing and heart rate change. 

--------------------------------------------------------
IMPORTANT: 

***How to name my Acq file? 

	Because in MATLAB, we have  the following codes to find our files as an example 
	path='D:\Dropbox\Neuropixels-Jonny_Hanbing\WITH_HANBING_TEST'; % Open the parent folder
	test = 'Breathing_Lightgradient'; % Open the test folder within the parent folder
	filename = [ mouse '_' cage '.mat']; % Get the name of mat file downloaded from Acknowledge, like 11N_C1294
	Data = load([path '\' test '\' group '\' filename]); % Load data, group is the genotype 
----------------------------------------------------------
RR: Add smoothing part here! 
add more calculation. 
WaveletStats_ASSR WINDOW!
Physizoo to automate ECG
qrs=mhrv.ecg.wjqrs(ECG',niSampRate,0.5,0.03,10); % Peak detection. unit:n samples.
z-score 