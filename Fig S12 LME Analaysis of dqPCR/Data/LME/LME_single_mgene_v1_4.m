%compile converted drop data from single M gene population experiments and
%run LME model
%1-13-21 v1.4 clean up for paper data
%Geoff Zath

clear; clc

%% Inputs

filename_tech = 'compiled drop data v1.txt';
filename_tech_days = 'compiled drop data v3.txt'; %day1 cycle22 removed
cycles = [19 20 21 22]; %cycles detected
cycle_start = 18;
reps = 3; %technical replicates
days = 3; %experiments run
cpd_exp = 171; %expected CPD



%% Process data

%load data
file_tech = importdata(filename_tech);
[x1 y1] = size(file_tech);

file_tech_days = importdata(filename_tech_days);
[x2 y2] = size(file_tech_days);

array_LME_tech = [];
array_string_tech = " ";
count = 1;

%tech rep data only

%save data for LME (tech rep only)
for i = 1 : y1
    
    cycle = file_tech(1,i);
    tube = file_tech(2,i);
    temp1 = file_tech(3:end,i);
    temp2 = temp1(~isnan(temp1));
    L = length(temp2);
    array_LME_tech(count:count+L-1,2) = cycle;
    array_LME_tech(count:count+L-1,3) = tube;
    
    %change tube into string
    if tube == 1
        
        array_string_tech(count:count+L-1) = "A";
        
    elseif tube == 2 
        
        array_string_tech(count:count+L-1) = "B";
        
    else
        
        array_string_tech(count:count+L-1) = "C";
        
    end
        
    array_LME_tech(count:count+L-1,1) = temp2;
    
    count = count + L;
    
    data_tech{i} = temp1(~isnan(temp1)); %data of technical replicates
    
   
    
end

array_string_tech = array_string_tech'; %for LME

%save data for bar graph
C = length(cycles);
R = reps;
for i = 1 : C %cycles
    for j = 1 : R %reps
        
        idx = j + R*(i-1);
        
        temp = data_tech{idx};
        
        bar_array_tech_avg(i,j) = mean(temp);
        bar_array_tech_std(i,j) = std(temp);
        bar_array_tech_SE(i,j) = std(temp)/sqrt(length(temp));
        
    end
end

%randomize technical replicates
rng('default');

%create random index array
for i = 1 : C
        
        idx_rand(i,:) = randperm(3);
        
end

bar_array_randtech_avg = zeros(C,R);
bar_array_randtech_std = zeros(C,R);

for i = 1 : C %cycles
    for j = 1 : R %reps
        
        idx = idx_rand(i,j);
        bar_array_randtech_avg(i,j) = bar_array_tech_avg(i,idx);
        bar_array_randtech_std(i,j) = bar_array_tech_std(i,idx);
        bar_array_randtech_SE(i,j) = bar_array_tech_SE(i,idx);
        
    end
end




%tech and exp rep data

array_LME_both = [];
array_string_tech_days = " ";
count = 1;

%save data for LME (tech and exp rep)
for i = 1 : y2
    
    cycle = file_tech_days(1,i);
    day = file_tech_days(2,i);
    tube = file_tech_days(3,i);
    temp1 = file_tech_days(4:end,i);
    temp2 = temp1(~isnan(temp1));
    L = length(temp2);
    
    if cycle == 19
        
        cycle = cycles(1);
        
    elseif cycle == 20
        
        cycle = cycles(2);
        
    elseif cycle == 21
        
        cycle = cycles(3);
        
    else
        
        cycle = cycles(4);
        
    end
    
    
    array_LME_tech_days(count:count+L-1,2) = cycle;
    array_LME_tech_days(count:count+L-1,3) = day;
    array_LME_tech_days(count:count+L-1,4) = tube;
    

    temp3 = temp2;
    
    array_LME_tech_days(count:count+L-1,1) = temp3;
    
    
    %change day into string
    if day == 1
        
        array_string_tech_days(count:count+L-1,1) = "A";
        
    elseif day == 2 
        
        array_string_tech_days(count:count+L-1,1) = "B";
        
    else
        
        array_string_tech_days(count:count+L-1,1) = "C";
        
    end
    
    
    %change tube into string
    if tube == 1
        
        array_string_tech_days(count:count+L-1,2) = "A";
        
    elseif tube == 2 
        
        array_string_tech_days(count:count+L-1,2) = "B";
        
    else
        
        array_string_tech_days(count:count+L-1,2) = "C";
        
    end
    
    
    
        
    
    
    count = count + L;
    
    %data_tech_days{i} = temp1(~isnan(temp1)); %data of technical replicates
    data_tech_days{i} = temp3; %data of technical replicates
   
    
end

%array_string_tech_days = array_string_tech_days'; %for LME

%save data for bar graph
C = length(cycles);
R = reps;
D = days;


for i = 1 : y2
    
    temp = data_tech_days{i};
    
    bar_tech_days_avg(i) = mean(temp);
    bar_tech_days_std(i) = std(temp);
    bar_tech_days_SE(i) = std(temp)/sqrt(length(temp));
    
end



%pull out day 1 and 2 (day 3 already calculated above)
for i = 1 : C
    
    
    bar_day1_avg(i) = bar_tech_days_avg(i + 4*(i-1));
    bar_day1_std(i) = bar_tech_days_std(i + 4*(i-1));
    bar_day1_SE(i) = bar_tech_days_SE(i + 4*(i-1));
    
    bar_day2_avg(i) = bar_tech_days_avg(i + 4*(i-1) + 1);
    bar_day2_std(i) = bar_tech_days_std(i + 4*(i-1) + 1);
    bar_day2_SE(i) = bar_tech_days_SE(i + 4*(i-1) + 1);
    
end

%adjust for removing day1 cycle22
bar_day2_avg(4) = bar_day1_avg(4);
bar_day2_std(4) = bar_day1_std(4);
bar_day2_SE(4) = bar_day1_SE(4);

bar_day1_avg(4) = [];
bar_day1_std(4) = [];
bar_day1_SE(4) = [];
    


%% LME

%tech reps only

%save data in table
tbl_tech = table(array_LME_tech(:,1),array_LME_tech(:,2),array_string_tech,'VariableNames',{'CPD','cycle','tube'});

%run LME
lme_tech = fitlme(tbl_tech,'CPD~cycle+(1|tube)');

%grab coefficients
lme_tech_coeff = lme_tech.Coefficients;
intercept_tech_mixed = double(lme_tech_coeff(1,2));
slope_tech_mixed = double(lme_tech_coeff(2,2));

intercept_tech_mixed_upper = double(lme_tech_coeff(1,8));
slope_tech_mixed_upper = double(lme_tech_coeff(2,8));

intercept_tech_mixed_lower = double(lme_tech_coeff(1,7));
slope_tech_mixed_lower = double(lme_tech_coeff(2,7));




%tech + exp reps

%save data in table
tbl_tech_days = table(array_LME_tech_days(:,1),array_LME_tech_days(:,2),array_string_tech_days(:,1),...
    array_string_tech_days(:,2),'VariableNames',{'CPD','cycle','day','tube'});

%run LME
lme_tech_days = fitlme(tbl_tech_days,'CPD~cycle+(1|day)+(1|tube)');

%grab coefficients
lme_tech_days_coeff = lme_tech_days.Coefficients;
intercept_tech_days_mixed = double(lme_tech_days_coeff(1,2));
slope_tech_days_mixed = double(lme_tech_days_coeff(2,2));

intercept_tech_days_mixed_upper = double(lme_tech_days_coeff(1,8));
slope_tech_days_mixed_upper = double(lme_tech_days_coeff(2,8));

intercept_tech_days_mixed_lower = double(lme_tech_days_coeff(1,7));
slope_tech_days_mixed_lower = double(lme_tech_days_coeff(2,7));



%% Figures

%colors
gray = linspecer('gray');
green = linspecer('green');
blue = linspecer('blue');
red = linspecer('red');

%single M gene population (tech reps, no random)
figure(1); clf(1)

h = barwitherr(bar_array_tech_std,bar_array_tech_avg);
set(h,'FaceColor',gray(64,:),'FaceAlpha',0.33,'linewidth',0.5);

yline(cpd_exp,':k','linewidth',1);
yline(cpd_exp*1.5,':k','linewidth',1);
yline(cpd_exp*2/3,':k','linewidth',1);



xlabel('Cycle')
ylabel('CPD')
ylim([0 350]);

box on
set(gca,'XTickLabel',{'19','20','21','22'})
set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','off',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 420,250])





%single M gene population (tech reps, random)
figure(2); clf(2)


h = barwitherr(bar_array_randtech_SE,bar_array_randtech_avg);
set(h,'FaceColor',gray(64,:),'FaceAlpha',0.33,'linewidth',0.5);


yline(cpd_exp,':k','linewidth',0.5);
yline(cpd_exp*1.5,':k','linewidth',1);
yline(cpd_exp*2/3,':k','linewidth',1);


xlabel('Cycle')
ylabel('CPD')
ylim([0 350]);

box on
set(gca,'XTickLabel',{'19','20','21','22'})
set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','off',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 420,250])





%single M gene population w/ LME model (tech reps, random)
figure(3); clf(3)

hold on

cycle_graph = [cycles; cycles; cycles];

offset = [-0.1 0 0.1];

%offset data points
for i = 1 : R
    
    cycles_offset = cycles + offset(i);
    
    errorbar(cycles_offset,bar_array_randtech_avg(:,i),bar_array_randtech_SE(:,i),'.k',...
        'markersize',5,'linewidth',0.5)
    
end

%LME model
x = 18:23;
y = slope_tech_mixed*x + intercept_tech_mixed;
plot(x,y,'--k','linewidth',0.5)

y_lower = slope_tech_mixed_lower*x + intercept_tech_mixed_lower;
plot(x,y_lower,'--k','linewidth',0.5)

y_upper = slope_tech_mixed_upper*x + intercept_tech_mixed_upper;
plot(x,y_upper,'--k','linewidth',0.5)

yline(cpd_exp,':k','linewidth',0.5);
yline(cpd_exp*1.5,':k','linewidth',1);
yline(cpd_exp*2/3,':k','linewidth',1);

hold off


xlabel('Cycle')
ylabel('CPD')
title('Tech Reps')
axis([18 23 0 300]);

box on

set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])
set(gcf, 'Position',  [100, 100, 420,250])






%single M gene population w/ LME model (tech and exp reps, random day 3)
figure(4); clf(4)

hold on

cycle_graph = [cycles; cycles; cycles];

offset = [-0.05 0 0.05];
  
plot(cycles(1:3),bar_day1_avg,'sk','markersize',10,'linewidth',0.5,'color',red(64,:))
    
plot(cycles,bar_day2_avg,'^k','markersize',8,'linewidth',0.5,'color',green(64,:))
    
%offset data points (day 3)
for i = 1 : R
    
    cycles_offset = cycles + offset(i);

    plot(cycles_offset,bar_array_randtech_avg(:,i),'ok','markersize',7,'linewidth',0.5,...
        'color',blue(64,:))
    
end

%LME model
x = cycle_start:23;
y = slope_tech_days_mixed*x + intercept_tech_days_mixed;
h2 = plot(x,y,'--k','linewidth',0.5);


y_lower = slope_tech_days_mixed_lower*x + intercept_tech_days_mixed_lower;
plot(x,y_lower,'--k','linewidth',0.5)

y_upper = slope_tech_days_mixed_upper*x + intercept_tech_days_mixed_upper;
plot(x,y_upper,'--k','linewidth',0.5)

h1 = yline(cpd_exp,':k','linewidth',1);
%yline(cpd_exp*1.5,':k','linewidth',1);
%yline(cpd_exp*2/3,':k','linewidth',1);
yline(cpd_exp*2,':k','linewidth',1);
yline(cpd_exp*0.5,':k','linewidth',1);

hold off




xlabel('Cycle')
ylabel('cpd')

% legend([h1 h2],{'171 cpd \pm 1.5-fold change','Fixed effects model w/ 95% CI'},...
%     'fontsize',6,'location','se')
legend([h1 h2],{'171 cpd \pm 2-fold change','Fixed effects model w/ 95% CI'},...
    'fontsize',6,'location','se')

%title('Tech + Exp Reps')
axis([cycle_start 23 0 370])

box on

set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on',...
    'yminortick','on','layer','top','ticklength',[0.015 1])

set(gcf, 'Position',  [100, 100, 405,230])

%print -painters -depsc FigS6C.eps


