%combine all PCR detection data for M gene control experiments
%Geoff Zath
%4-5-21 v1.1 remove Cy% channel

%4-16-21
%M gene std curve 100 cpd in 50 um drops (10^7 copies/uL)

%rename import file names, 'Cycles' variable, .mat file name

%0.70V FAM gain

clear; clc

%Cycles = [1 10 13 16 19 22 25 28 40]; %UPDATE
Cycles = [1 10 13 16 19 22 28 40]; %outlier removed


filename = 'detection_data_std_041621_outlier.mat'; %UPDATE

%% Load Data

A{1} = importdata('std_1_peaks_100kHz_17.txt');
A{2} = importdata('std_2_peaks_100kHz_15.txt');
A{3} = importdata('std_3_peaks_100kHz_12.txt');
A{4} = importdata('std_4_peaks_100kHz_10.txt');
A{5} = importdata('std_5_peaks_100kHz_9.txt');
A{6} = importdata('std_6_peaks_100kHz_7.txt');
%A{7} = importdata('std_7_peaks_100kHz_6.txt');
A{7} = importdata('std_8_peaks_100kHz_5.txt');
A{8} = importdata('std_40_peaks_100kHz_4.txt');
%A{10} = importdata('19_peaks_100kHz_4.txt');
% A{11} = importdata('19_peaks_100kHz_5.txt');
% A{12} = importdata('22_peaks_100kHz_1.txt');
% A{13} = importdata('22_peaks_100kHz_2.txt');
% A{14} = importdata('22_peaks_100kHz_3.txt');
% A{15} = importdata('22_peaks_100kHz_4.txt');
% A{16} = importdata('22_peaks_100kHz_5.txt');
% A{17} = importdata('40_peaks_100kHz.txt');

%% Combine Data

L = length(A);

for i = 1 : L
    
        temp_cell = A{i};
        temp_data = temp_cell.data;

        Time_Inside{i} = temp_data(:,13);
        %CB{i} = temp_data(:,1);
        FAM{i} = temp_data(:,2);
        ROX{i} = temp_data(:,3);
        %CY5{i} = temp_data(:,4);
        
end

README = ['Cell Structure = (cycle)'];

%save(filename,'Time_Inside','FAM','ROX','CY5','Cycles','README')
save(filename,'Time_Inside','FAM','ROX','Cycles','README')