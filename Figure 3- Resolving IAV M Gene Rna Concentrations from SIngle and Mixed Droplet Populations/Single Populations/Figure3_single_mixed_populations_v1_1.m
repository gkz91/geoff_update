%plot single and mized mgene populations with the same graphs as burst size
%data
%3-30-22 v1.1 add fold change calc
%Figure 3 for manuscript
%Geoff Zath

clear; clc

%% Inputs

cpd_min = 1; %set threshold above theoretical detection limir

filename{1} = 'conv_10cpd_c22-25_2020exp.mat';
filename{2} = 'conv_100cpd_c19-22_2020exp.mat';
filename{3} = 'conv_1000cpd_c16-22_2020exp.mat';

N_sample = 100; %number of data points for resampling
nbins = 50;

nbins_gini = 100;
edges_gini_10cpd = logspace(1,2,nbins_gini);
edges_gini_100cpd = logspace(1.5,3,nbins_gini);
edges_gini_1000cpd = logspace(2,4.5,nbins_gini);

nbins_combo = 30;
edges_combo = linspace(0.5,4.5,nbins_combo); %bin edges for combo graph

conc_expected = [1.71e1 1.71e2 1.71e3]; %cpd



%% Process data

%load data


file_10cpd = load(filename{1});
data_10cpd = file_10cpd.conc_convert;
data_10cpd = data_10cpd(data_10cpd<500);

file_100cpd = load(filename{2});
data_100cpd = file_100cpd.conc_convert;

file_1000cpd = load(filename{3});
data_1000cpd = file_1000cpd.conc_convert;
data_1000cpd = data_1000cpd(data_1000cpd>100);



%randomly sample data for bar graphs
s = RandStream('mlfg6331_64'); 

data_10cpd_resample = datasample(s,data_10cpd,N_sample,'Replace',true);
data_100cpd_resample = datasample(s,data_100cpd,N_sample,'Replace',true);
data_1000cpd_resample = datasample(s,data_1000cpd,N_sample,'Replace',true);




%Gini coefficient


%10cpd
[N_10cpd edges_10cpd] = histcounts(data_10cpd,edges_gini_10cpd);
[G_10cpd,Lor_10cpd] = gini(N_10cpd,edges_10cpd(1:end-1));
CV_10cpd = std(data_10cpd)/mean(data_10cpd);
avg_10cpd = mean(data_10cpd);

%100cpd
[N_100cpd edges_100cpd] = histcounts(data_100cpd,edges_gini_100cpd);
[G_100cpd,Lor_100cpd] = gini(N_100cpd,edges_100cpd(1:end-1));
CV_100cpd = std(data_100cpd)/mean(data_100cpd);
avg_100cpd = mean(data_100cpd);

%1000cpd
[N_1000cpd edges_1000cpd] = histcounts(data_1000cpd,edges_gini_1000cpd);
[G_1000cpd,Lor_1000cpd] = gini(N_1000cpd,edges_1000cpd(1:end-1));
CV_1000cpd = std(data_1000cpd)/mean(data_1000cpd);
avg_1000cpd = mean(data_1000cpd);

%measure %diff from expected

%using the average value
P_diff_10cpd = (mean(data_10cpd) - conc_expected(1))/conc_expected(1);
P_diff_100cpd = (mean(data_100cpd) - conc_expected(2))/conc_expected(2);
P_diff_1000cpd = (mean(data_1000cpd) - conc_expected(3))/conc_expected(3);

P_diff_mean = (P_diff_10cpd + P_diff_100cpd + P_diff_1000cpd)/3
P_diff_std = std([P_diff_10cpd, P_diff_100cpd, P_diff_1000cpd])

%10cpd
if P_diff_10cpd >= 0
    
    Fold_10cpd = P_diff_10cpd + 1;
    
else
    
    Fold_10cpd = 1/(P_diff_10cpd + 1);
    
end

%100cpd
if P_diff_100cpd >= 0
    
    Fold_100cpd = P_diff_100cpd + 1;
    
else
    
    Fold_100cpd = 1/(P_diff_100cpd + 1);
    
end

%1000cpd
if P_diff_1000cpd >= 0
    
    
    Fold_1000cpd = P_diff_1000cpd + 1;
    
else
    
    Fold_1000cpd = 1/(P_diff_1000cpd + 1);
    
end

avg_fold = mean([Fold_10cpd Fold_100cpd Fold_1000cpd]);
std_fold = std([Fold_10cpd Fold_100cpd Fold_1000cpd]);
CV_fold = std_fold/avg_fold;



%using all values
P_diff_10cpd_all = (data_10cpd - conc_expected(1))/conc_expected(1);
P_diff_100cpd_all = (data_100cpd - conc_expected(2))/conc_expected(2);
P_diff_1000cpd_all = (data_1000cpd - conc_expected(3))/conc_expected(3);


%10cpd
L_10cpd = length(P_diff_10cpd_all);

for i = 1 : L_10cpd
    
    if P_diff_10cpd_all(i) >= 0

        Fold_10cpd_all(i) = P_diff_10cpd_all(i) + 1;

    else

        Fold_10cpd_all(i) = 1/(P_diff_10cpd + 1);

    end
end

Fold_10cpd_all = Fold_10cpd_all(Fold_10cpd_all<10);

%100cpd
L_100cpd = length(P_diff_100cpd_all);

for i = 1 : L_100cpd
    
    if P_diff_100cpd_all(i) >= 0

        Fold_100cpd_all(i) = P_diff_100cpd_all(i) + 1;

    else

        Fold_100cpd_all(i) = 1/(P_diff_100cpd + 1);

    end
end

Fold_100cpd_all = Fold_100cpd_all(Fold_100cpd_all<10);

%1000cpd
L_1000cpd = length(P_diff_1000cpd_all);

for i = 1 : L_1000cpd
    
    if P_diff_1000cpd_all(i) >= 0

        Fold_1000cpd_all(i) = P_diff_1000cpd_all(i) + 1;

    else

        Fold_1000cpd_all(i) = 1/(P_diff_1000cpd + 1);

    end
end

Fold_1000cpd_all = Fold_1000cpd_all(Fold_1000cpd_all<10);

avg_fold_10cpd_all = mean(Fold_10cpd_all);
std_fold_10cpd_all = std(Fold_10cpd_all);

avg_fold_100cpd_all = mean(Fold_100cpd_all);
std_fold_100cpd_all = std(Fold_100cpd_all);

avg_fold_1000cpd_all = mean(Fold_1000cpd_all);
std_fold_1000cpd_all = std(Fold_1000cpd_all);

avg_fold_all = mean([avg_fold_10cpd_all avg_fold_100cpd_all avg_fold_1000cpd_all]);
std_fold_all = std([avg_fold_10cpd_all avg_fold_100cpd_all avg_fold_1000cpd_all]);



%% Figures

%colors
red = linspecer('red');
green = linspecer('green');
blue = linspecer('blue');
gray = linspecer('gray');





%histogram of 10cpd m gene drop data
figure(2); clf(2)

histogram(log10(data_10cpd),floor(nbins/3),'normalization','probability','facecolor',red(64,:),'facealpha',0.5,...
    'linewidth',0.5)


title('10 cpd data')
xlabel('log_{10}(cpd)')
ylabel('Fraction')
axis([0 2 0 inf])

box on
set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf, 'Position',  [100, 100, 420,250])





%histogram of 100cpd m gene drop data
figure(3); clf(3)



histogram(log10(data_100cpd),floor(nbins/3),'normalization','probability','facecolor',blue(64,:),'facealpha',0.5,...
    'linewidth',0.5)


title('100 cpd data')
xlabel('log_{10}(cpd)')
ylabel('Fraction')
axis([1.5 3 0 inf])
box on
set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf, 'Position',  [100, 100, 420,250])






%histogram of 1000cpd m gene drop data
figure(4); clf(4)


histogram(log10(data_1000cpd),floor(nbins/3),'normalization','probability','facecolor',green(96,:),'facealpha',0.5,...
    'linewidth',0.5)

%title('1000 cpd data')
xlabel('log_{10}(cpd)')
ylabel('Fraction')
axis([2 4 0 inf])
box on
set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')


set(gcf, 'Position',  [100, 100, 405,230])



%print -painters -depsc Fig2E_1000cpd_cycle22_2020exp.eps




%bar graph of cpd vs drop in order (10,100,1000 cpd combined)
figure(6); clf(6)

data_10cpd_resample_sort = sort(data_10cpd_resample);
data_100cpd_resample_sort = sort(data_100cpd_resample);
data_1000cpd_resample_sort = sort(data_1000cpd_resample);

hold on

bar(data_1000cpd_resample_sort,'facecolor',green(96,:),'facealpha',0.5,'edgecolor','none')
bar(data_100cpd_resample_sort,'facecolor',blue(96,:),'facealpha',0.5,'edgecolor','none')
bar(data_10cpd_resample_sort,'facecolor',red(96,:),'facealpha',0.5,'edgecolor','none')

hold off

box on
xlabel('Drop index')
ylabel('cpd')
axis([0 inf 1e0 1e4])
set(gca,'fontsize',10,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')

%set(gcf, 'Position',  [100, 100, 420,250])
%set(gcf, 'Position',  [100, 100, 300,250])
set(gcf, 'Position',  [100, 100, 405,230])


%print -painters -depsc Fig2I_cpd_bar_graph.eps




%bar graph of cpd vs drop in order (10 cpd combined)
figure(7); clf(7)


bar(data_10cpd_resample_sort,'facecolor',red(96,:),'facealpha',0.5,'edgecolor','none')


box on
xlabel('drop #')
ylabel('RNA per drop')
axis([0 inf 1e0 1e4])
set(gca,'fontsize',8,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf, 'Position',  [100, 100, 420,250])





%bar graph of cpd vs drop in order (100 cpd combined)
figure(8); clf(8)


bar(data_100cpd_resample_sort,'facecolor',blue(96,:),'facealpha',0.5,'edgecolor','none')


box on
xlabel('drop #')
ylabel('RNA per drop')
axis([0 inf 1e0 1e4])
set(gca,'fontsize',8,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf, 'Position',  [100, 100, 420,250])






%bar graph of cpd vs drop in order (1090 cpd combined)
figure(9); clf(9)


bar(data_1000cpd_resample_sort,'facecolor',green(96,:),'facealpha',0.5,'edgecolor','none')


box on
xlabel('drop #')
ylabel('RNA per drop')
axis([0 inf 1e0 1e4])
set(gca,'fontsize',8,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.03 1],'layer','top')

set(gcf, 'Position',  [100, 100, 420,250])




%inverse Lorenz plot
figure(10); clf(10)


Lor_10cpd_inv = 1 - Lor_10cpd;
Lor_100cpd_inv = 1 - Lor_100cpd;
Lor_1000cpd_inv = 1 - Lor_1000cpd;

hold on

plot([0 1],[0 1],'--k','linewidth',1)
plot(Lor_10cpd_inv(:,1),Lor_10cpd_inv(:,2),'--k','linewidth',1,'color',red(64,:))
plot(Lor_100cpd_inv(:,1),Lor_100cpd_inv(:,2),':k','linewidth',1.5,'color',blue(64,:))
plot(Lor_1000cpd_inv(:,1),Lor_1000cpd_inv(:,2),'-k','linewidth',1,'color',green(64,:))



hold off


box on
xlabel('Fraction of drops')
ylabel('Fraction of RNA')


legend(['uniform, {\itG} = 0'],...
    ['10 cpd, {\itG} = ',num2str(G_10cpd,'%0.3f')],...
    ['100 cpd, {\itG} = ',num2str(G_100cpd,'%0.3f')],...
    ['1000 cpd, {\itG} = ',num2str(G_1000cpd,'%0.3f')],...
    'fontsize',8,'location','se')
    

set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')

set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc Fig2F.eps


%histogram of all individual measurements (10,100,1000 cpd)
figure(11); clf(11)

hold on

histogram(log10(data_10cpd),edges_combo,'normalization','probability','facecolor',red(96,:),'facealpha',0.5,...
    'linewidth',0.5)
histogram(log10(data_100cpd),edges_combo,'normalization','probability','facecolor',blue(96,:),'facealpha',0.5,...
    'linewidth',0.5)
histogram(log10(data_1000cpd),edges_combo,'normalization','probability','facecolor',green(96,:),'facealpha',0.5,...
    'linewidth',0.5)

hold off

xlabel('log_{10}(cpd)')
ylabel('Fraction')
axis([0.5 4.5 0 0.4])
box on
set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on','layer','top','ticklength',[0.015 1])


set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc Fig2_pooled.eps

