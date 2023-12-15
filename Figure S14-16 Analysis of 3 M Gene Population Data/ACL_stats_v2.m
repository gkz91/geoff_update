%Stats on Summmer 2020 ACL data
%11-28-22 v2 Add LME model
%Geoff Zath

%Summer 2020 data

clear; clc
close all

%% Inputs

N_sample = 500; %random sample of N drops for ANOVA

cpd_expected = [17.1 171 1710]; %epxected concentration (cpd)
cycles = [22 25; 19 22; 16 22]; %cycles measured


%% Load Data

%10 cpd data
file_10cpd{1} = load('conv_10cpd_c22');
file_10cpd{2} = load('conv_10cpd_c23');
file_10cpd{3} = load('conv_10cpd_c24');
file_10cpd{4} = load('conv_10cpd_c25');

%100 cpd data
file_100cpd{1} = load('conv_100cpd_c19');
file_100cpd{2} = load('conv_100cpd_c20');
file_100cpd{3} = load('conv_100cpd_c21');
file_100cpd{4} = load('conv_100cpd_c22');

%1000 cpd data
file_1000cpd{1} = load('conv_1000cpd_c16');
file_1000cpd{2} = load('conv_1000cpd_c17');
file_1000cpd{3} = load('conv_1000cpd_c18');
file_1000cpd{4} = load('conv_1000cpd_c19');
file_1000cpd{5} = load('conv_1000cpd_c20');
file_1000cpd{6} = load('conv_1000cpd_c21');
file_1000cpd{7} = load('conv_1000cpd_c22');








%% Process Data

%10 cpd stats

L = length(file_10cpd);
data_10cpd_pool = [];
count = 0;

for i = 1 : L

    data_10cpd{i} = file_10cpd{i}.conc_convert;
    avg_10cpd(i) = mean(data_10cpd{i});
    std_10cpd(i) = std(data_10cpd{i});
    CV_10cpd(i) = std_10cpd(i)/avg_10cpd(i)*100;
    N_10cpd(i) = length(data_10cpd{i});
    
    %pool
    data_10cpd_pool = [data_10cpd_pool data_10cpd{i}];
    
end

%10cpd pool
avg_10cpd_pool = mean(data_10cpd_pool);
std_10cpd_pool = std(data_10cpd_pool);
CV_10cpd_pool = std_10cpd_pool/avg_10cpd_pool*100;
N_10cpd_pool = length(data_10cpd_pool);



%100 cpd data

L = length(file_100cpd);
data_100cpd_pool = [];
count = 0;

for i = 1 : L

    data_100cpd{i} = file_100cpd{i}.conc_convert;
    avg_100cpd(i) = mean(data_100cpd{i});
    std_100cpd(i) = std(data_100cpd{i});
    CV_100cpd(i) = std_100cpd(i)/avg_100cpd(i)*100;
    N_100cpd(i) = length(data_100cpd{i});
    
    %pool
    data_100cpd_pool = [data_100cpd_pool data_100cpd{i}];
    
end

%100cpd pool
avg_100cpd_pool = mean(data_100cpd_pool);
std_100cpd_pool = std(data_100cpd_pool);
CV_100cpd_pool = std_100cpd_pool/avg_100cpd_pool*100;
N_100cpd_pool = length(data_100cpd_pool);



%1000 cpd data

L = length(file_1000cpd);
data_1000cpd_pool = [];
count = 0;

for i = 1 : L

    data_1000cpd{i} = file_1000cpd{i}.conc_convert;
    avg_1000cpd(i) = mean(data_1000cpd{i});
    std_1000cpd(i) = std(data_1000cpd{i});
    CV_1000cpd(i) = std_1000cpd(i)/avg_1000cpd(i)*100;
    N_1000cpd(i) = length(data_1000cpd{i});
    
    %pool
    data_1000cpd_pool = [data_1000cpd_pool data_1000cpd{i}];
    
end

%100cpd pool
avg_1000cpd_pool = mean(data_1000cpd_pool);
std_1000cpd_pool = std(data_1000cpd_pool);
CV_1000cpd_pool = std_1000cpd_pool/avg_1000cpd_pool*100;
N_1000cpd_pool = length(data_1000cpd_pool);

%% Stats

%calculate % difference

%prep arrays for LME
cycle_record = []; %cycle for each data point
cpd_record = []; %converted cpd of each data point

%10 cpd
num = 1;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;






for i = 1 : length(cycles_all)
    
    temp = data_10cpd{i};
    
    for j = 1 : N_10cpd(i)
        
        perc_diff_10cpd_temp(j) = (temp(j) - cpd_expected(num))/cpd_expected(num) * 100;
        cycle_temp(j) = cycles_all(i);
        
    end
    
    perc_diff_10cpd{i} = perc_diff_10cpd_temp;
    perc_diff_10cpd_avg(i) = mean(perc_diff_10cpd_temp);
    perc_diff_10cpd_std(i) = std(perc_diff_10cpd_temp);
    
    %cycle_temp = datasample(cycle_temp, N_sample);
    cycle_record = [cycle_record cycle_temp];
    
    %perc_diff_10cpd_temp = datasample(perc_diff_10cpd_temp, N_sample);
    cpd_record = [cpd_record perc_diff_10cpd_temp];
    
    clear cycle_temp
    clear perc_diff_10cpd_temp;
    
end

L1 = length(cycle_record);
exp_record(1:L1) = 10; %expected cpd for each data point


%100 cpd
num = 2;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;

for i = 1 : length(cycles_all)
    
    temp = data_100cpd{i};
    
    for j = 1 : N_100cpd(i)
        
        perc_diff_100cpd_temp(j) = (temp(j) - cpd_expected(num))/cpd_expected(num) * 100;
        cycle_temp(j) = cycles_all(i);
        
    end
    
    perc_diff_100cpd{i} = perc_diff_100cpd_temp;
    perc_diff_100cpd_avg(i) = mean(perc_diff_100cpd_temp);
    perc_diff_100cpd_std(i) = std(perc_diff_100cpd_temp);
    
    %cycle_temp = datasample(cycle_temp, N_sample);
    cycle_record = [cycle_record cycle_temp];
    
    %perc_diff_100cpd_temp = datasample(perc_diff_100cpd_temp, N_sample);
    cpd_record = [cpd_record perc_diff_100cpd_temp];
    
    clear cycle_temp
    clear perc_diff_100cpd_temp;
    
end

L2 = length(cycle_record);
exp_record(L1+1:L2) = 100; %expected cpd for each data point


%1000 cpd
num = 3;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;

for i = 1 : length(cycles_all)
    
    temp = data_1000cpd{i};
    
    for j = 1 : N_1000cpd(i)
        
        perc_diff_1000cpd_temp(j) = (temp(j) - cpd_expected(num))/cpd_expected(num) * 100;
        cycle_temp(j) = cycles_all(i);
        
    end
    
    perc_diff_1000cpd{i} = perc_diff_1000cpd_temp;
    perc_diff_1000cpd_avg(i) = mean(perc_diff_1000cpd_temp);
    perc_diff_1000cpd_std(i) = std(perc_diff_1000cpd_temp);
    
    %cycle_temp = datasample(cycle_temp, N_sample);
    cycle_record = [cycle_record cycle_temp];
    
    %perc_diff_1000cpd_temp = datasample(perc_diff_1000cpd_temp, N_sample);
    cpd_record = [cpd_record perc_diff_1000cpd_temp];
    
    clear cycle_temp
    clear perc_diff_1000cpd_temp;
    
end

L3 = length(cycle_record);
exp_record(L2+1:L3) = 1000; %expected cpd for each data point


%% ANOVA

%10 cpd
L = length(file_10cpd);

for i = 1 : L
    
    ANOVA_10cpd(:,i) = datasample(data_10cpd{i}, N_sample);
    perc_diff_sample_10cpd(:,i) = datasample(perc_diff_10cpd{i}, N_sample);
    
end

[p_10cpd,tbl_10cpd,stats_10cpd] = anova1(ANOVA_10cpd);



%100 cpd
L = length(file_100cpd);

for i = 1 : L%-1
    
    %ANOVA_100cpd(:,i) = datasample(data_100cpd{i+1}, N_sample);
    ANOVA_100cpd(:,i) = datasample(data_100cpd{i}, N_sample);
    perc_diff_sample_100cpd(:,i) = datasample(perc_diff_100cpd{i}, N_sample);
    
end

[p_100cpd,tbl_100cpd,stats_100cpd] = anova1(ANOVA_100cpd);



%1000 cpd
L = length(file_1000cpd);

for i = 1 : L
    
    ANOVA_1000cpd(:,i) = datasample(data_1000cpd{i}, N_sample);
    perc_diff_sample_1000cpd(:,i) = datasample(perc_diff_1000cpd{i}, N_sample);
    
end

[p_1000cpd,tbl_1000cpd,stats_1000cpd] = anova1(ANOVA_1000cpd);
    
%combine % diff data
perc_diff_sampe_combine = [perc_diff_sample_10cpd perc_diff_sample_100cpd perc_diff_sample_1000cpd];



%% LME

%save data in table
tbl_LME = table(cpd_record',cycle_record',exp_record','VariableNames',{'Measured_CPD','cycle','Expected_CPD'});
%tbl_LME = table(cpd_record(L1+1:end)',cycle_record(L1+1:end)',exp_record(L1+1:end)','VariableNames',{'Measured_CPD','cycle','Expected_CPD'});

%run LME
%lme = fitlme(tbl_LME,'Measured_CPD~cycle+(1|Expected_CPD)')
lme = fitlme(tbl_LME,'Measured_CPD~Expected_CPD+(1|cycle)')

%grab coefficients
lme_coeff = lme.Coefficients;
intercept_mixed = double(lme_coeff(1,2));
slope_mixed = double(lme_coeff(2,2));

intercept_mixed_upper = double(lme_coeff(1,8));
slope_mixed_upper = double(lme_coeff(2,8));

intercept_mixed_lower = double(lme_coeff(1,7));
slope_mixed_lower = double(lme_coeff(2,7));







%% Figures



%percent diffs 10 cpd
figure(7); clf(7)

num = 1;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;

cycles_combine = [];
cycles_combine = [cycles_combine cycles_all];

boxplot(perc_diff_sample_10cpd,cycles_all)

box on
xlabel('Cycle')
ylabel('Percent Difference from Expected')
title('10 CPD (Summer 2020)')

axis([0 inf -100 350])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')



%percent diffs 100 cpd
figure(8); clf(8)

num = 2;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;
cycles_combine = [cycles_combine cycles_all];

boxplot(perc_diff_sample_100cpd,cycles_all)

box on
xlabel('Cycle')
ylabel('Percent Difference from Expected')
title('100 CPD (Summer 2020)')

axis([0 inf -100 110])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')


%percent diffs 1000 cpd
figure(9); clf(9)

num = 3;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;
cycles_combine = [cycles_combine cycles_all];

boxplot(perc_diff_sample_1000cpd,cycles_all)

box on
xlabel('Cycle')
ylabel('Percent Difference from Expected')
title('1000 CPD (Summer 2020)')

axis([0 inf -100 200])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')


%percent diffs 1000 cpd (first 4 cycles)
figure(10); clf(10)

num = 3;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;

boxplot(perc_diff_sample_1000cpd(:,1:4),cycles_all(1:4))

box on
xlabel('Cycle')
ylabel('Percent Difference from Expected')
title('1000 CPD (Summer 2020, First 4 Cycles)')

axis([0 inf -100 110])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')




%percent diffs combined CPD
figure(11); clf(11)

cycles_combine_string = string(cycles_combine);
conc_string = {'10^1 cpd'; '10^1 cpd'; '10^1 cpd'; '10^1 cpd';
    '10^2 cpd'; '10^2 cpd'; '10^2 cpd'; '10^2 cpd';
    '10^3 cpd'; '10^3 cpd'; '10^3 cpd'; '10^3 cpd';
    '10^3 cpd'; '10^3 cpd'; '10^3 cpd'}';

for i = 1 : length(conc_string)
    
    cycle_label{i} = append(cycles_combine_string{i}," (",conc_string{i},")");
    
end

boxplot(perc_diff_sampe_combine,cycle_label)

box on
xlabel('Cycle + Concentration')
ylabel('Percent Difference from Expected')
title('Combined Data (Summer 2020)')

xtickangle(45)

axis([0 inf -100 350])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')



%percent diff combined cpd (errorbars)
figure(12); clf(12)

hold on

%10 cpd
num = 1;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;
errorbar(cycles_all,perc_diff_10cpd_avg,perc_diff_10cpd_std,'o','markersize',9,'linewidth',0.5,'color','k')

%100 cpd
num = 2;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;
errorbar(cycles_all,perc_diff_100cpd_avg,perc_diff_100cpd_std,'x','markersize',9,'linewidth',0.5,'color','b')

%1000 cpd
num = 3;
start = cycles(num,1);
stop = cycles(num,2);
cycles_all = start:stop;
errorbar(cycles_all,perc_diff_1000cpd_avg,perc_diff_1000cpd_std,'sq','markersize',9,'linewidth',0.5,'color','r')

yline(100,'--k')


%LME model
x = 15:26;
y = slope_mixed*x + intercept_mixed;
plot(x,y,':k','linewidth',0.75)

y_lower = slope_mixed_lower*x + intercept_mixed_lower;
plot(x,y_lower,':k','linewidth',0.75)

y_upper = slope_mixed_upper*x + intercept_mixed_upper;
plot(x,y_upper,':k','linewidth',0.75)

yline(-50,'--k')

hold off


box on
xlabel('Cycle')
ylabel('Percent Difference from Expected')
title('Combined Data (Summer 2020)')

legend({'10 cpd'; '100 cpd'; '1000 cpd'; '+/- 2-fold'; 'LME w/ 95% CI'},'location','nw')


axis([15 26 -90 110])

set(gca,'fontsize',12,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on',...
    'ticklength',[0.015 1],'layer','top')



























