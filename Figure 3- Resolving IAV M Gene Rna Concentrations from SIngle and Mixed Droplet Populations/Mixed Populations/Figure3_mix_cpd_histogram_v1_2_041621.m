%plot histograms of mixed cpd populations
%4-6-21 v1.2_041621 add pooled histogram
%Figure 3 for Manuscript
%Geoff Zath

clear; clc

%% Inputs

nbins = 30;
edges_delRn = linspace(-0.2,3,nbins);
edges_cpd = linspace(0.5,4.5,nbins);

cpd_min = 1; %set threshold above theoretical detection limit

nbins_gini = 100;
edges_gini = logspace(1,5,nbins_gini);

N_sample = 100; %number of data points for resampling

%% Load data

%mixed M gene (fluorescence)
A = load('conv_mixture_pooled_filter_041621.mat');
data_pooled = A.conc_convert;



%% Process data


%pool cpd data
cpd_mix_pool = data_pooled;



%Gini coefficient
[N_BS edges_BS] = histcounts(cpd_mix_pool,edges_gini);
[G_mix,Lor_mix] = gini(N_BS,edges_BS(1:end-1));  


%randomly sample data for bar graphs
s = RandStream('mlfg6331_64'); 

data_mix_resample = datasample(s,cpd_mix_pool,N_sample,'Replace',true);



%% Figures

%colors
blue = linspecer('blue');
red = linspecer('red');
green = linspecer('green');
gray = linspecer('gray');



%histogram of pooled mixture
figure(3); clf(3)

hold on

histogram(log10(cpd_mix_pool),edges_cpd,'normalization','probability','facecolor',gray(96,:),'facealpha',0.5,...
    'linewidth',0.5)

hold off

N = length(cpd_mix_pool);

xlabel('log_{10}(cpd)')
ylabel('Fraction')
axis([min(edges_cpd) max(edges_cpd) 0 0.15])
box on
set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on','layer','top','ticklength',[0.015 1])


set(gcf, 'Position',  [100, 100, 325,250])




%inverse Lorenz plot
figure(4); clf(4)


Lor_mix_inv = 1 - Lor_mix;


hold on

plot([0 1],[0 1],'--k','linewidth',1)
plot(Lor_mix_inv(:,1),Lor_mix_inv(:,2),'-k','linewidth',1)


hold off


box on
xlabel('Fraction of drops')
ylabel('Fraction of RNA')

legend(['uniform, {\itG} = 0'],...
    ['mix, {\itG} = ',num2str(G_mix,'%0.3f')],...
    'fontsize',8,'location','se')


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')


set(gcf, 'Position',  [100, 100, 325,250])

%print -painters -depsc Fig2J_Lorenz_plot.eps





%bar graph of cpd vs drop in order (mixed)
figure(5); clf(5)

data_mix_resample_sort = sort(data_mix_resample);

bar(data_mix_resample_sort,'facecolor',gray(96,:),'facealpha',0.5,'edgecolor','none')

box on
xlabel('Drop index')
ylabel('cpd')
axis([0 inf 1e0 1e5])
set(gca,'fontsize',10,'linewidth',0.5,'yscale','log','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')

set(gcf, 'Position',  [100, 100, 325,250])





