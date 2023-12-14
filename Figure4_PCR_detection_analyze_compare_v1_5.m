%plot pooled burst size data from all reps/cycles and compare H1N1 to H3N2
%1-19-22 v1.5 update graphs for paper
%Figure 4 for Manuscript
%Geoff Zath

%based on PCR_detection_analyze_single_pop_FAM_v1_4_pool

%FAM (M gene)

%H1N1 + H3N2 pooled


clear; clc

%% Inputs

N_sample = 100; %number of data points for resampling


nbins = 20;
split_ratio = 1/8; %split ratio of EVO chip

edges_lin = linspace(0,200,nbins);
edges_log = linspace(0,5,nbins);

nbins_gini = 100;
edges_gini = logspace(0,6,nbins_gini);

%% Load Data

%H1N1 Data
H1N1{1} = load('Mgene_filter_021221.mat');
H1N1{2} = load('Mgene_filter_021921.mat');
H1N1{3} = load('Mgene_filter_022421.mat');

%H3N2 Data
H3N2{1} = load('Mgene_filter_022621.mat');
H3N2{2} = load('Mgene_filter_031721.mat');
H3N2{3} = load('Mgene_filter_032421.mat');



%% Process Data

%H1N1
L_H1N1 = length(H1N1);

for i = 1 : L_H1N1
    
    temp = H1N1{i}.Mgene_data_filter;
    
    H1N1_avg(i) = mean(temp);
    H1N1_min(i) = min(temp);
    H1N1_max(i) = max(temp);
    H1N1_IQR(i) = iqr(temp);
    H1N1_IQR_range(i,1) = prctile(temp,25);
    H1N1_IQR_range(i,2) = prctile(temp,75);
    
    H1N1_conc_convert_cell{i} = temp;
    
end

BS_data_H1N1_pooled = [H1N1_conc_convert_cell{:}];
H1N1_pooled_avg = mean(BS_data_H1N1_pooled);
H1N1_pooled_IQR_range(1) = prctile(BS_data_H1N1_pooled,25);
H1N1_pooled_IQR_range(2) = prctile(BS_data_H1N1_pooled,75);



%H3N2
L_H3N2 = length(H3N2);

for i = 1 : L_H3N2
    
    temp = H3N2{i}.Mgene_data_filter;
    
    H3N2_avg(i) = mean(temp);
    H3N2_min(i) = min(temp);
    H3N2_max(i) = max(temp);
    H3N2_IQR(i) = iqr(temp);
    H3N2_IQR_range(i,1) = prctile(temp,25);
    H3N2_IQR_range(i,2) = prctile(temp,75);
    
    H3N2_conc_convert_cell{i} = temp;
    
end

BS_data_H3N2_pooled = [H3N2_conc_convert_cell{:}];
H3N2_pooled_avg = mean(BS_data_H3N2_pooled);
H3N2_pooled_IQR_range(1) = prctile(BS_data_H3N2_pooled,25);
H3N2_pooled_IQR_range(2) = prctile(BS_data_H3N2_pooled,75);



%Gini coefficient

%H1N1
%BS_data_H1N1 = H1N1_conc_convert;
[N_BS_H1N1 edges_BS_H1N1] = histcounts(BS_data_H1N1_pooled,edges_gini);
[G_H1N1,Lor_H1N1] = gini(N_BS_H1N1,edges_BS_H1N1(1:end-1));

%H3N2
%BS_data_H3N2 = H3N2_conc_convert;
[N_BS_H3N2 edges_BS_H3N2] = histcounts(BS_data_H3N2_pooled,edges_gini);
[G_H3N2,Lor_H3N2] = gini(N_BS_H3N2,edges_BS_H3N2(1:end-1));

%update for both datasets

%Shaprio-Wilk normality test
BS_data_log_H1N1 = log10(BS_data_H1N1_pooled);
[H_SW_H1N1, p_SW_H1N1, SWstatistic_H1N1]  = swtest(BS_data_log_H1N1 ); %if H = 1, not normal

BS_data_log_H3N2 = log10(BS_data_H3N2_pooled);
[H_SW_H3N2, p_SW_H3N2, SWstatistic_H3N2]  = swtest(BS_data_log_H3N2 ); %if H = 1, not normal

%Kolmogorov-Smirnov test
avg_KS_H1N1 = mean(BS_data_log_H1N1);
std_KS_H1N1 = std(BS_data_log_H1N1);
test_KS_H1N1 = (BS_data_log_H1N1 - avg_KS_H1N1)/std_KS_H1N1;
[H_KS_H1N1 p_KS_H1N1 kstat_H1N1] = kstest(test_KS_H1N1);

avg_KS_H3N2 = mean(BS_data_log_H3N2);
std_KS_H3N2 = std(BS_data_log_H3N2);
test_KS_H3N2 = (BS_data_log_H3N2 - avg_KS_H3N2)/std_KS_H3N2;
[H_KS_H3N2 p_KS_H3N2 kstat_H3N2] = kstest(test_KS_H3N2);

%Lilliefors test
[H_L_exp p_L_exp] = lillietest(BS_data_H1N1_pooled,'Distribution','exponential');
[H_L_ev p_L_ev] = lillietest(log(BS_data_H1N1_pooled),'Distribution','extreme value');

%Jarque-Bera test
[H_JB P_JB] = jbtest(BS_data_log_H1N1);


%two sample Kolmogorov-Smirnov test
[H_KS2 p_KS2 KS2stat] = kstest2(test_KS_H1N1,test_KS_H3N2,'tail','larger');




%randomly sample data for bar graphs
s = RandStream('mlfg6331_64'); 

H1N1_resample = datasample(s,BS_data_H1N1_pooled,N_sample,'Replace',true);
H3N2_resample = datasample(s,BS_data_H3N2_pooled,N_sample,'Replace',true);



%% Figures

blue = linspecer('blue');
green = linspecer('green');
red = linspecer('red');
gray = linspecer('gray');



%Lorentz curve
figure(1); clf(1)

hold on

gini(N_BS_H1N1,edges_BS_H1N1(1:end-1),true);
gini(N_BS_H3N2,edges_BS_H3N2(1:end-1),true);

hold off

xlabel('fraction of drops (low to high)')
ylabel('fraction of RNA')
set(gca,'fontsize',12,'linewidth',1)




%bar graph of cpd vs cell in order (H1N1)
figure(2); clf(2)



BS_data_sort_H1N1 = sort(BS_data_H1N1_pooled);

bar(BS_data_sort_H1N1,'facecolor',red(96,:))


xlabel('drop #')
ylabel('RNA per drop')
title('H1N1')
axis([0 inf 1e0 1e7])
set(gca,'fontsize',14,'linewidth',1,'yscale','log')





%bar graph of cpd vs cell in order (H3N2)
figure(3); clf(3)



BS_data_sort_H3N2 = sort(BS_data_H3N2_pooled);


bar(BS_data_sort_H3N2,'facecolor',blue(96,:))


xlabel('drop #')
ylabel('RNA per drop')
title('H3N2')
axis([0 inf 1e0 1e7])
set(gca,'fontsize',14,'linewidth',1,'yscale','log')





%histogram of #cells vs RNA
figure(4); clf(4)


edges_f4 = linspace(0,1e4,20);

N_H1N1 = length(BS_data_H1N1_pooled);
N_H3N2 = length(BS_data_H3N2_pooled);

hold on

histogram(BS_data_H3N2_pooled,edges_f4,'normalization','probability',...
    'facecolor',blue(96,:),'facealpha',0.33)
histogram(BS_data_H1N1_pooled,edges_f4,'normalization','probability',...
    'facecolor',red(96,:),'facealpha',0.33)


hold off

box on

xlabel('RNA per drop')
ylabel('Fraction')
legend(['H3N2, N = ',num2str(N_H1N1),' drops'],['H1N1, N = ',num2str(N_H3N2),' drops'],...
    'fontsize',10)

set(gca,'fontsize',14,'linewidth',1)



%cumulative fraction plot (lorenz plot)
figure(5); clf(5)


Lor_inv_H1N1 = 1 - Lor_H1N1;
Lor_inv_H3N2 = 1 - Lor_H3N2;

hold on

h1 = plot(Lor_inv_H3N2(:,1),Lor_inv_H3N2(:,2),':','linewidth',1.5,...
    'color',blue(96,:));

h2 = plot(Lor_inv_H1N1(:,1),Lor_inv_H1N1(:,2),'--','linewidth',1.5,...
    'color',red(96,:));


hold off

box on
xlabel('Fraction of drops')
ylabel('Fraction of vRNA')
legend([h1 h2],{['H3N2, {\itG}_{pool} = ',num2str(G_H3N2,'%0.3f')],['H1N1, {\itG}_{pool} = ',num2str(G_H1N1,'%0.3f')]},...
    'fontsize',14,'location','se')
set(gca,'fontsize',16,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf,'Position',  [100, 100, 650, 425]);

%print -painters -depsc Fig4K_pooled_Lorenz.eps




%fit normal dist to log10 data
figure(6); clf(6)

hold on

h_H3N2 = histfit(BS_data_log_H3N2,nbins,'gamma');

h_H1N1 = histfit(BS_data_log_H1N1,nbins,'gamma');

hold off

h_H1N1(1).FaceColor = red(96,:);
h_H1N1(1).FaceAlpha = 0.5;
h_H1N1(2).Color = [.2 .2 .2];

h_H3N2(1).FaceColor = blue(96,:);
h_H3N2(1).FaceAlpha = 0.5;
h_H3N2(2).Color = [.2 .2 .2];


box on

xlabel('log10(RNA per drop)')
ylabel('# of drops')

axis([0 7 0 inf])

legend([h_H3N2(1) h_H1N1(1)],{['H3N2, N = ',num2str(N_H3N2),' Drops'],...
    ['H1N1, N = ',num2str(N_H1N1),' Drops']},'fontsize',10)
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on')






%log10 histogram
figure(7); clf(7)

hold on

histogram(BS_data_log_H3N2,edges_log,'normalization','probability','facecolor',blue(96,:),...
    'facealpha',0.5)

histogram(BS_data_log_H1N1,edges_log,'normalization','probability','facecolor',red(96,:),...
    'facealpha',0.5)

hold off



box on

xlabel('log10(cpd)')
ylabel('Fraction')

axis([min(edges_log) max(edges_log) 0 0.25])

legend(['H3N2, N = ',num2str(N_H3N2),' drops'],...
    ['H1N1, N = ',num2str(N_H1N1),' drops'],'location','nw','fontsize',14)
set(gca,'fontsize',16,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf,'Position',  [100, 100, 650, 425]);

%print -painters -depsc Fig4J_pooled_histogram.eps


%bar graph of cpd vs drop in order (10,100,1000 cpd combined)
figure(8); clf(8)

H1N1_resample_sort = sort(H1N1_resample);
H3N2_resample_sort = sort(H3N2_resample);

hold on

bar(H3N2_resample_sort,'facecolor',blue(96,:),'facealpha',0.5,'edgecolor','none')
bar(H1N1_resample_sort,'facecolor',red(96,:),'facealpha',0.5,'edgecolor','none')

hold off

box on


legend('H3N2','H1N1','fontsize',14,'location','nw')
axis([0 inf 1e0 1e5])
set(gca,'fontsize',14,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')

xlabel('Drop index','fontsize',16)
ylabel('cpd','fontsize',16)

set(gcf,'Position',  [100, 100, 650, 425]);

%print -painters -depsc Fig4I_pooled_bar.eps


