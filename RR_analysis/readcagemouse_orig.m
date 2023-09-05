
function [cage, mouse, day,group,fr] = readcagemouse(n)
% read experiments to be analyzed
if n == 1 
mouse = '11N';
cage = 'C1294';
day = '071221';
group = 'NPY2R_Ai32';
fr = 0.25; % ISI
end 

if n  == 2 
mouse = '3R';
cage = 'C1278';
day = '071421';
group = 'NPY2R_Ai32';
fr = 0.25; 
end 

if n == 3
mouse = '6L';
cage = 'C1279A';
day = '071421';
group = 'NPY2R_Ai32';
fr = 0.25;
end 

if n == 4 
mouse = '32N';
cage = 'C1422A';
day = '090921';
group = 'NPY2R_Ai32';
fr = 0.5;
end 

if n == 5 
mouse = '18N';
cage = 'C1339';
day = '090921';
group = 'NPY2R_Ai32';
fr = 0.5;
end 

if n == 6
mouse = '1N';
cage = 'C1278A';
day = '091321';
group = 'NPY2R_Ai32';
fr = 0.25; 
end 

if n ==7
mouse = '16N';
cage = 'C1337';
day = '091321';
group = 'NPY2R_Ai32';
fr = 0.5; 
end 

if n == 8
mouse = '34L';
cage = 'C1422';
day = '091321';
group = 'NPY2R_Ai32';
fr = 0.5; 
end 

if n == 9 % dead
mouse = '11R';
cage = 'C1341A';
group = 'NPY2R_Ai9';
day = '';
fr = 0.5; 
end 

if n == 10 % dead 
mouse = '12L';
cage = 'C1341';
group = 'NPY2R_Ai9';
day = '';
fr = 0.5; 
end 

if n == 11 
mouse = '15L';
cage = 'C1388';
group = 'NPY2R_Ai9';
day = '';
fr = 0.5; 
end 

if n == 12
mouse = '16N';
cage = 'C1389B';
group = 'NPY2R_Ai9';
day = '';
fr = 0.5; 
end 

if n == 13
mouse = '5R';
cage = 'C1279B';
group = 'NPY2R_Ai32';
day = '';
fr = 0.25; 
end 

if n == 14
mouse = '21R';
cage = 'C1369';
group = 'NPY2R_Ai32';
day = '';
fr = 0.5; 
end 

if n == 15
mouse = '20R';
cage = 'C1451A';
group = 'NPY2R_Ai9';
day = '';
fr = 0.5; 
end 

if n == 16
mouse = '7N';
cage = 'C1293';
group = 'NPY2R_NTS';
day = '';
fr = 0.5; 
end 

if n == 17
mouse = '100N';
cage = 'C8888';
group = 'NPY2R_Ai32';
day = '';
fr = 0.5; 
end 

end 
 
