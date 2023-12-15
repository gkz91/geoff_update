%Convert single cycle PCR detection data with std curve from PCR efficiency model
%3-4-21 v1.2 %add gini coeff and lorentz curve, #cells per fraction of RNA
%graph, and Shaprio-Wilk normality test
%Geoff Zath

%2-26-21
%FAM (M gene)

%use ROX normalized data only (what std curve was made with)

%add ability to measure multiple cycles

clear; clc

%% Inputs

cycle = 25; %cycle data to convert
nbins = 30;
split_ratio = 1/8; %split ratio of EVO chip

edges_lin = linspace(0,200,nbins);
edges_log = linspace(0,6,nbins);

nbins_gini = 100;
edges_gini = logspace(0,6,nbins_gini);

%% Load Data
%curve_data = load('stdcurve_delRnRn0_detection_data_std_070320.mat');
%PCR_data = load('processed_delRnRn0_detection_data_unk_1T_070320.mat');
%curve_data = load('stdcurve_delRn_detection_data_std_070320.mat');
curve_data = load('eff_FAM_stdcurve_delRn_detection_data_std_022621.mat');
PCR_data = load('processed_delRn_detection_data_IAV_022621_v7.mat');
%PCR_data = load('processed_delRn_detection_data_std_121420.mat');


PCR_cycles = PCR_data.cycle;
PCR_cycle_data = PCR_data.delRn_FAM_FINAL;

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


%Gini coefficient
%BS_data = log10(conc_convert/split_ratio);
BS_data = conc_convert/split_ratio;
[N_BS edges_BS] = histcounts(BS_data,edges_gini);
[G,Lor] = gini(N_BS,edges_BS(1:end-1));
%gini(N_BS,edges_BS(1:end-1),true);



%Shaprio-Wilk normality test
BS_data_log = log10(BS_data);
[H, p_SW, SWstatistic]  = swtest(BS_data_log); %if H = 1, not normal




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
% title(['mean = ',num2str(conc_avg,'%1.2E'),...
%     ' +/- ',num2str(conc_std,'%1.2E')])
title(['median = ',num2str(conc_median,'%1.2E'),' cpd'])
xlabel('cpd')
%ylabel('Number of Drops')
ylabel('PDF')
legend(['N = ',num2str(N),' Drops'])
%title(['Converted Data Before Split (Cycle ',num2str(cycle,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1)





%Lorentz curve
figure(3); clf(3)

gini(N_BS,edges_BS(1:end-1),true);

xlabel('fraction of drops (low to high)')
ylabel('fraction of RNA')
set(gca,'fontsize',12,'linewidth',1)




%bar graph of cpd vs cell in order
figure(4); clf(4)

BS_data_sort = sort(BS_data);

bar(BS_data_sort,'facecolor',green(100,:))

xlabel('drop #')
ylabel('RNA per drop')
axis([0 inf 1e1 1e5])
set(gca,'fontsize',14,'linewidth',1,'yscale','log')



%histogram of #cells vs RNA
figure(5); clf(5)

% edges_f5 = linspace(0,0.6,20);
% 
% BS_data_norm = BS_data/sum(BS_data);
% histogram(BS_data_norm,edges_f5,'facecolor',blue(96,:),'facealpha',0.5)
% 
% xlabel('fraction of RNA')
% ylabel('number of cells')
% axis([0 inf 0 inf])
% set(gca,'fontsize',14,'linewidth',1,'yscale','lin')

edges_f5 = linspace(0,4e4,20);

N = length(conc_convert);

histogram(BS_data,edges_f5,'facecolor',blue(96,:),'facealpha',0.5)
%axis([1e1 1e4 0 inf])
box on
%axis([20 120 0 inf])
% title(['mean = ',num2str(conc_avg,'%1.2E'),...
%     ' +/- ',num2str(conc_std,'%1.2E')])
%title(['median = ',num2str(conc_median,'%1.2E'),' cpd'])
xlabel('RNA per drop')
%ylabel('Number of Drops')
ylabel('# of drops')
legend(['N = ',num2str(N),' cells'])
%title(['Converted Data Before Split (Cycle ',num2str(cycle,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1)



%cumulative fraction plot
figure(6); clf(6)


Lor_inv = 1 - Lor;
% L_BS_data = length(BS_data);


plot(Lor_inv(:,1),Lor_inv(:,2),'-k','linewidth',2)
xlabel('fraction of drops (high to low)')
ylabel('cumulative fraction of RNA')
title(['Gini coefficient = ',num2str(G)]);
set(gca,'fontsize',14,'linewidth',1,'yscale','lin')






%fit normal dist to log10 data
figure(7); clf(7)

h = histfit(BS_data_log,10,'normal');

h(1).FaceColor = blue(96,:);
h(1).FaceAlpha = 0.5;

%title(['median = ',num2str(conc_median,'%1.2E'),' cpd'])

if p_SW < 0.05
    
    title('Shapiro-Wilk Test p < 0.05 ')
    
else
    
    title(['Shapiro-Wilk Test p = ',num2str(p_SW,'%1.2E')])

end

xlabel('log10(RNA per drop)')
%ylabel('Number of Drops')
ylabel('# of drops')
legend(['N = ',num2str(N),' Drops'])
%title(['Converted Data Before Split (Cycle ',num2str(cycle,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on')









