%graph multiplex PCR data for FAM and CY5 probes in drops
%2-15-21 v1.2
%Geoff Zath

%2-12-21 burst size experiment

clear; clc

%% Inputs

split_ratio = 1/8; %EVO chip split ratio

%heatmap
Nx = 50; %x bin
Ny = 50; %y bin

xmin = 0;
xmax = 2000;
ymin = 0;
ymax = 500;

xmin_log = 2;
xmax_log = 4;
ymin_log = 1;
ymax_log = 3;



%% Load data


%converted CPD data (linear)
data = load('C19_conc_Bactin.mat');
Bactin_C19 = data.conc_convert/split_ratio;

data = load('C22_conc_Bactin.mat');
Bactin_C22 = data.conc_convert/split_ratio;

data = load('C25_conc_Bactin.mat');
Bactin_C25 = data.conc_convert/split_ratio;

data = load('C28_conc_Bactin.mat');
Bactin_C28 = data.conc_convert/split_ratio;


data = load('C19_conc_Mgene.mat');
Mgene_C19 = data.conc_convert/split_ratio;

data = load('C22_conc_Mgene.mat');
Mgene_C22 = data.conc_convert/split_ratio;

data = load('C25_conc_Mgene.mat');
Mgene_C25 = data.conc_convert/split_ratio;

data = load('C28_conc_Mgene.mat');
Mgene_C28 = data.conc_convert/split_ratio;



% %converted CPD data (log)
% Bactin_C22_log = log10(Bactin_C22);
% Bactin_C25_log = log10(Bactin_C25);
% 
% Mgene_C22_log = log10(Mgene_C22);
% Mgene_C25_log = log10(Mgene_C25);

%% Process data





%% Figures

%colors
map = colormap(linspecer);
%map(1,:) = 1; %background white = 1 or black = 0
blue = linspecer('blue');
green = linspecer('green');





%heatmap of Bactin CPD vs Mgene CPD at cycle 19
figure(1); clf(1)

%Make 2D histogram
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);
[N_2D Xedge Yedge] = histcounts2(Mgene_C19,Bactin_C19,x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 19')
xlabel('M gene (CPD)')
ylabel('\beta actin (CPD)')
axis([xmin xmax ymin ymax])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])






%heatmap of Bactin CPD vs Mgene CPD at cycle 22
figure(2); clf(2)

%Make 2D histogram
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);
[N_2D Xedge Yedge] = histcounts2(Mgene_C22,Bactin_C22,x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 22')
xlabel('M gene (CPD)')
ylabel('\beta actin (CPD)')
axis([xmin xmax ymin ymax])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])







%heatmap of Bactin CPD vs Mgene CPD at cycle 25
figure(3); clf(3)

%Make 2D histogram
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);
[N_2D Xedge Yedge] = histcounts2(Mgene_C25,Bactin_C25,x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 25')
xlabel('M gene (CPD)')
ylabel('\beta actin (CPD)')
axis([xmin xmax ymin ymax])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])






%heatmap of Bactin CPD vs Mgene CPD at cycle 28
figure(4); clf(4)

%Make 2D histogram
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);
[N_2D Xedge Yedge] = histcounts2(Mgene_C28,Bactin_C28,x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 28')
xlabel('M gene (CPD)')
ylabel('\beta actin (CPD)')
axis([xmin xmax ymin ymax])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])







%heatmap of Bactin CPD vs Mgene CPD at cycle 22 (log)
figure(5); clf(5)

%Make 2D histogram
x = linspace(xmin_log, xmax_log, Nx);
y = linspace(ymin_log, ymax_log, Ny);
[N_2D Xedge Yedge] = histcounts2(log10(Mgene_C22),log10(Bactin_C22),x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 22 log)')
xlabel('M gene (log10(CPD))')
ylabel('\beta actin (log10(CPD))')
axis([xmin_log xmax_log ymin_log, ymax_log])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])





%heatmap of Bactin CPD vs Mgene CPD at cycle 25 (log)
figure(6); clf(6)

%Make 2D histogram
x = linspace(xmin_log, xmax_log, Nx);
y = linspace(ymin_log, ymax_log, Ny);
[N_2D Xedge Yedge] = histcounts2(log10(Mgene_C25),log10(Bactin_C25),x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title('Cycle 25 log)')
xlabel('M gene (log10(CPD))')
ylabel('\beta actin (log10(CPD))')
axis([xmin_log xmax_log ymin_log, ymax_log])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])






