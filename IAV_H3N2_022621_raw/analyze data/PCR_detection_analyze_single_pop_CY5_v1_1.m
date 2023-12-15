%Convert single cycle PCR detection data with std curve from PCR efficiency model
%7-2-20 v1.1 %added median instead of mean
%Geoff Zath

%2-26-21
%CY5 (B actin)

%use ROX normalized data only (what std curve was made with)

%add ability to measure multiple cycles

clear; clc

%% Inputs

cycle = 25; %cycle data to convert
nbins = 20;
split_ratio = 1/8; %split ratio of EVO chip

edges_lin = linspace(0,100,nbins);
edges_log = linspace(0,3,nbins);

%% Load Data
%curve_data = load('stdcurve_delRnRn0_detection_data_std_070320.mat');
%PCR_data = load('processed_delRnRn0_detection_data_unk_1T_070320.mat');
%curve_data = load('stdcurve_delRn_detection_data_std_070320.mat');
curve_data = load('eff_CY5_stdcurve_delRn_detection_data_std_022621.mat');
PCR_data = load('processed_delRn_detection_data_IAV_022621_v7.mat');
%PCR_data = load('processed_delRn_detection_data_std_121420.mat');


PCR_cycles = PCR_data.cycle;
PCR_cycle_data = PCR_data.delRn_CY5_FAM_FINAL;

cycle_loc = find(PCR_cycles == cycle);

PCR_convert = PCR_cycle_data{cycle_loc};

Xmodel_data = curve_data.model_scaled;
Xmodel = Xmodel_data(cycle,:);
Ymodel = curve_data.conc_stdc;


%% Process Data

%find location of PCR_convert data in Xmodel and use matching Ymodel data
L = length(PCR_convert);

for i = 1 : L
    
    k(i) = dsearchn(Xmodel',PCR_convert(i));
    conc_convert(i) = Ymodel(k(i));
    
end



%stats
conc_avg = mean(conc_convert/split_ratio);
conc_log10_avg = mean(log10(conc_convert/split_ratio));
conc_std = std(conc_convert/split_ratio);
conc_median = median(conc_convert/split_ratio);
conc_CV = conc_std/conc_avg

%% Figures

blue = linspecer('blue');
green = linspecer('green');

%converted data (log scale)
figure(1); clf(1)

min_log_diff = conc_avg - 3*conc_std;
if min_log_diff < 0
    min_log = 0;
else
    min_log = log10(min_log_diff);
end

max_log = log10(conc_avg + 3*conc_std);

%edges_log = linspace(min_log,max_log,nbins);

N = length(conc_convert);

histogram(log10(conc_convert/split_ratio),edges_log,...
    'normalization','pdf','facecolor',blue(96,:),'facealpha',0.5)
%axis([1e1 1e4 0 inf])
box on
%axis([1 2.5 0 inf])
% title(['log scale, mean = ',num2str(conc_avg,'%1.2E'),...
%     ' +/- ',num2str(conc_std,'%1.2E')])
% title(['mean = ',num2str(conc_avg,'%1.2E'),...
%     ' +/- ',num2str(conc_std,'%1.2E')])
title(['median = ',num2str(conc_median,'%1.2E'),' cpd'])
xlabel('log10(cpd)')
%ylabel('Number of Drops')
ylabel('PDF')
legend(['N = ',num2str(N),' Drops'])
%title(['Converted Data Before Split (Cycle ',num2str(cycle,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on')






%converted data (linear scale)
figure(2); clf(2)

min_lin = conc_avg - 3*conc_std;
max_lin = conc_avg + 3*conc_std;

%edges_lin = linspace(min_lin,max_lin,nbins);
%edges_lin = linspace(0,5,nbins);

N = length(conc_convert);

histogram(conc_convert/split_ratio,edges_lin,...
    'normalization','pdf','facecolor',green(96,:),'facealpha',0.5)
%axis([1e1 1e4 0 inf])
box on
%axis([20 120 0 inf])
% title(['linear scale, mean = ',num2str(conc_avg,'%1.2E'),...
%     ' +/- ',num2str(conc_std,'%1.2E')])
title(['median = ',num2str(conc_median,'%1.2E'),' cpd'])
xlabel('cpd')
%ylabel('Number of Drops')
ylabel('PDF')
legend(['N = ',num2str(N),' Drops'])
%title(['Converted Data Before Split (Cycle ',num2str(cycle,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1)










