%plot reference curves with respective unk data
%11-14-22 v1.6_041621 add efficiency to graph
%Figure 2 - dqPCR Manuscript
%Geoff Zath

clear; clc

%% Inputs

cycle40 = 1:40;
N_span_10 = 10; %number of reference curve
N_span_1000 = 1000;
T_delRn_C1_max_1000cpd = 0.1697; %min threshold on fluorescence from histogram script
cycle_std = 18; %cycle for standard curve 
N_cycle = 5; %cycle index, 5=cycle18 6=cycle19


%% Load data

%unknown data

%10 cpd unk
filename = 'processed_delRn_detection_data_10^6_unk_041621_nofilter.mat';
PCR_curve = load(filename);
cycle_10cpd = PCR_curve.cycle;
C_10cpd = length(cycle_10cpd);
delRn_10cpd_avg = PCR_curve.delRn_FAM_avg_FINAL;
delRn_10cpd_std = PCR_curve.delRn_FAM_std_FINAL;

%100 cpd unk
filename = 'processed_delRn_detection_data_10^7_unk_041621_nofilter.mat';
PCR_curve = load(filename);
cycle_100cpd = PCR_curve.cycle;
C_100cpd = length(cycle_100cpd);
delRn_100cpd_avg = PCR_curve.delRn_FAM_avg_FINAL;
delRn_100cpd_std = PCR_curve.delRn_FAM_std_FINAL;

%10 cpd unk
filename = 'processed_delRn_detection_data_10^8_unk_041621_nofilter.mat';
PCR_curve = load(filename);
cycle_1000cpd = PCR_curve.cycle;
C_1000cpd = length(cycle_100cpd);
delRn_1000cpd_avg = PCR_curve.delRn_FAM_avg_FINAL;
delRn_1000cpd_std = PCR_curve.delRn_FAM_std_FINAL;


%reference curves

%10 cpd ref
A = load('model_cycle_outlier_041621.mat');
curve_cycle_10cpd = A.cycle_curve;
A = load('model_curve_outlier_041621.mat');
curve_data_10cpd = A.model_best;
A = load('model_avg_outlier_041621.mat');
curve_avg_10cpd = A.D_avg_curve;
A = load('model_std_outlier_041621.mat');
curve_std_10cpd = A.D_std_curve;

%100 cpd ref
A = load('model_cycle_outlier_041621.mat');
curve_cycle_100cpd = A.cycle_curve;
A = load('model_curve_outlier_041621.mat');
curve_data_100cpd = A.model_best;
A = load('model_avg_outlier_041621.mat');
curve_avg_100cpd = A.D_avg_curve;
A = load('model_std_outlier_041621.mat');
curve_std_100cpd = A.D_std_curve;

%1000 cpd ref
A = load('model_cycle_outlier_041621.mat');
curve_cycle_1000cpd = A.cycle_curve;
A = load('model_curve_outlier_041621.mat');
curve_data_1000cpd = A.model_best;
A = load('model_avg_outlier_041621.mat');
curve_avg_1000cpd = A.D_avg_curve;
A = load('model_std_outlier_041621.mat');
curve_std_1000cpd = A.D_std_curve;



%efficiency data from fit
B = load('model_eff_041621.mat');
E_best_rev = B.E_best_rev;

%R^2 model fit
C = load('model_R2_outlier_041621.mat');
SCFPE_R2 = C.R_sq_eff_best;



%amplification curve library

%1000 cpd (N = 10)
A = load('eff_FAM_stdcurve_delRn_detection_data_std_041621_outlier_N10.mat');
conc_1000cpd_N10 = A.conc_stdc;
amp_library_1000cpd_N10 = A.model_scaled;

%1000 cpd (N = 1000)
A = load('eff_FAM_stdcurve_delRn_detection_data_std_041621_outlier_N1000.mat');
conc_1000cpd_N1000 = A.conc_stdc;
amp_library_1000cpd_N1000 = A.model_scaled;



%% Process data

%find threshold where X% drops amplify at cycle 1

%1000 cpd
x_min = T_delRn_C1_max_1000cpd ;
x_max = curve_avg_1000cpd(end) - 1*curve_std_1000cpd(end); %copies/uL



%% Figures

%colors
blue = linspecer('blue');
red = linspecer('red');
green = linspecer('green');
gray = linspecer('gray');
color_span_10 = linspecer(N_span_10);
color_span_1000 = linspecer(N_span_1000);
N_gray_10 = 1;
gray_span_10 = linspecer(10 + N_gray_10,'gray');
N_gray_1000 = 100;
gray_span_1000 = linspecer(1000 + N_gray_1000,'gray');

%reference curves plotted with repsective unk cpd data
figure(1); clf(1)

fig = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

hold on

%fluorescence data
yyaxis right

%ref curve
%h1 = plot(cycle40,curve_data_10cpd,'-','linewidth',1,'color',gray(128,:)); %model fit


errorbar(curve_cycle_1000cpd,curve_avg_1000cpd,curve_std_1000cpd,'.','color',gray(128,:),'linewidth',0.5,...
    'markersize',10) %reference data
plot(cycle40,curve_data_1000cpd,'-k','linewidth',1)%,'color',green(96,:)) %model fit
%plot(cycle40-3.33,curve_data_1000cpd,':','linewidth',1.5,'color',green(96,:)) %shifted model fit


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])


ylabel('\DeltaR_N (a.u.)')
ylim([-0.2 3.5])
%ylim([-0.2 4])


%efficiency data
yyaxis left
h5 = plot(cycle40,E_best_rev,'--','linewidth',1,'color',gray(128,:));
ylabel('PCR Efficiency')
ylim([0 1.05])
%ylim([0 1.55])

hold off

box on


xlabel('Cycle')

xlim([0 41])

legend('PCR Eff.','1.71e+02 cpd','SCF-E','fontsize',8,'location','se')
set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 375,230])

%print -painters -depsc Fig2B.eps



%1000 cpd unk data and shifted red curve
figure(2); clf(2)



hold on




%1000 cpd data
plot(cycle40,curve_data_1000cpd,'-k','linewidth',1.5)%,'color',green(96,:)) %model fit

errorbar(curve_cycle_1000cpd,curve_avg_1000cpd,curve_std_1000cpd,'.','color',gray(128,:),'linewidth',0.5,...
    'markersize',7) %reference data

errorbar(cycle_1000cpd(3:end-3),delRn_1000cpd_avg(3:end-3),delRn_1000cpd_std(3:end-3),'.','markersize',7,...
    'color',green(96,:)) %unknown data


hold off

box on

xlabel('Cycle')
ylabel('\DeltaR_N (a.u.)')
legend('Reference curve','1.71e+02 cpd reference','1.71e+03 cpd unknown','fontsize',8,'location','nw')

ylim([-0.2 5])
xlim([0 41])


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 405,230])


%print -painters -depsc Fig2B_1000cpd_unk.eps



%10 cpd unk data and shifted ref curve
figure(3); clf(3)


hold on

%10 cpd data
plot(cycle40,curve_data_10cpd,'-k','linewidth',1.5)%,'color',red(96,:)) %model fit
errorbar(curve_cycle_10cpd,curve_avg_10cpd,curve_std_10cpd,'.','color',gray(128,:),'linewidth',0.5,...
    'markersize',7) %reference data

errorbar(cycle_10cpd(3:end-3),delRn_10cpd_avg(3:end-3),delRn_10cpd_std(3:end-3),'.','markersize',7,...
    'color',red(96,:)) %unknown data




hold off

box on

xlabel('Cycle')
ylabel('\DeltaR_N (a.u.)')
legend('Reference curve','1.71e+02 cpd reference','1.71e+01 cpd unknown','fontsize',8,'location','nw')

ylim([-0.2 5])
xlim([0 41])


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
%set(gcf, 'Position',  [100, 100, 370,230])
set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc Fig2D_10cpd_unk.eps



%100 cpd unk data and shifted red curve
figure(4); clf(4)


hold on


%100 cpd data
plot(cycle40,curve_data_100cpd,'-k','linewidth',1.5)%,'color',blue(96,:)) %model fit
errorbar(curve_cycle_100cpd,curve_avg_100cpd,curve_std_100cpd,'.','color',gray(128,:),'linewidth',0.5,...
    'markersize',7) %reference data

errorbar(cycle_100cpd(3:end-3),delRn_100cpd_avg(3:end-3),delRn_100cpd_std(3:end-3),'.','markersize',7,...
    'color',blue(96,:)) %unknown data


hold off

box on

xlabel('Cycle')
ylabel('\DeltaR_N (a.u.)')
legend('Reference curve','1.71e+02 cpd reference','1.71e+02 cpd unknown','fontsize',8,'location','nw')

ylim([-0.2 5])
%ylim([1e-2 10])
xlim([0 41])


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc Fig2C_100cpd_unk.eps





%1000 cpd with amp library N = 10
figure(5); clf(5)


hold on



%amp library
for i = 1 : N_span_10
    
    
    h(i) = plot(cycle40,amp_library_1000cpd_N10(1:40,i),'--','linewidth',1,'color',gray_span_10(i + N_gray_10,:));
    
end

colormap gray
colorbar('ticks',[0 0.5 1],'ticklabels',{'low cpd' 'int. cpd' 'high cpd'})

%ref curve
h1 = plot(cycle40,curve_data_10cpd,'-','linewidth',1,'color',gray(128,:)); %model fit

%10 cpd data
h2 = plot(cycle_10cpd(N_cycle),delRn_10cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',red(96,:)); %unknown data cycle 19=6 cycle 18=5

%100 cpd data
h3 = plot(cycle_100cpd(N_cycle),delRn_100cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',blue(96,:)); %unknown data cycle 19=6 cycle 18=5


%1000 cpd data
h4 = plot(cycle_1000cpd(N_cycle),delRn_1000cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',green(96,:)); %unknown data cycle 19=6 cycle 18=5






hold off


box on

xlabel('Cycle')
ylabel('\DeltaR_N (a.u.)')

ylim([-0.2 3.5])
xlim([0 41])


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 405,230])


%print -painters -depsc Fig2C.eps



figure(6); clf(6)


hold on

%amp library
for i = 1 :2: N_span_1000
    
    
    plot(cycle40,amp_library_1000cpd_N1000(1:40,i),'-','linewidth',1,'color',...
        gray_span_1000(i+N_gray_1000,:))
    
end

xline(cycle_std,':k','linewidth',1.5);

%ref curve
h1 = plot(cycle40,curve_data_10cpd,'-','linewidth',1,'color',gray(128,:)); %model fit

%10 cpd data
h2 = plot(cycle_10cpd(N_cycle),delRn_10cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',red(96,:)); %unknown data cycle 19=6 cycle 18=5

%100 cpd data
h3 = plot(cycle_100cpd(N_cycle),delRn_100cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',blue(96,:)); %unknown data cycle 19=6 cycle 18=5


%1000 cpd data
h4 = plot(cycle_1000cpd(N_cycle),delRn_1000cpd_avg(N_cycle),...
    '.','markersize',10,...
    'color',green(96,:)); %unknown data cycle 19=6 cycle 18=5





hold off


box on


xlabel('Cycle')
ylabel('\DeltaR_N (a.u.)')

ylim([-0.2 3.5])
xlim([0 41])


set(gca,'fontsize',10,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc Fig2F_amp_curve_library.eps



figure(7); clf(7)

hold on

semilogy(amp_library_1000cpd_N1000(cycle_std,:),conc_1000cpd_N1000,':k','linewidth',1.5)




%plot min and max of valid region at cycle 20

x_shade = [x_min x_max x_max x_min ];
y_shade = [1e-3 1e-3 1e7 1e7];
patch(x_shade,y_shade,gray(64,:),'facealpha',0.2,'edgecolor','none');%,'facecolor',gray(64,:));

xline(x_min,'-k','linewidth',0.5);
xline(x_max,'-k','linewidth',0.5);


hold off

box on



axis([1e-1 3.5 1e0 1e6])

set(gca,'fontsize',10,'linewidth',0.5,'yscale','log','xscale','log','xminortick','on','yminortick','on','layer','top','ticklength',[0.015 1])
ylabel('cpd','fontsize',10)
xlabel('\DeltaR_N (a.u.)','fontsize',10)
set(gcf, 'Position',  [100, 100, 405,230])


%print -painters -depsc Fig2D.eps



