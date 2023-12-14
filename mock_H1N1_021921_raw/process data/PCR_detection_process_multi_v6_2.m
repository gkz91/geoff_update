%Process multiplexed PCR detection data for MOI 0.1 burst size
%1-22-21 v6.2_multi added Cy5 data + cleaned code
%Geoff Zath

%delRn normalization

%mock drops 2-19-21

%%add split ratio to estimated positive drops


clear; clc

%% Inputs
dataname = 'detection_data_mock_021921.mat'; %rename me
nbins_curve = 1001; %number of bins for histogram curves (N+1)
nbins = 50; %number of bins for histograms

edges_FAM = linspace(0,1,nbins);
edges_CY5 = linspace(0.2,1.4,nbins);

%time inside thresholds
T_time = [4 7; 4 7; 4 8; 5 10; 5 10; 5 10; 5 10; 5 10];
%T_time = [0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10]; %no filter
edges_time = linspace(0,10,nbins);

%ROX thresholds
T_ROX = [0.3 0.5; 0.3 0.5; 0.4 0.7; 0.4 0.7; 0.4 0.7; 0.4 0.7; 0.4 0.7; 0.4 0.7];
%T_ROX = [0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10]; %no filter
edges_ROX = linspace(0.1,2,nbins);

%delRn (FAM) thresholds (final cleanup)
T_delRn_FAM = [-0.075 0.1; 0.1 3; 0.1 2; 0.1 2; 0.1 2; 0.1 2; 0.6 2; 0.6 2]; 
%T_delRn = [-2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10]; %no filter
edges_delRn_FAM = linspace(-0.1,1,nbins);

%delRn (CY5) thresholds (final cleanup)
T_delRn_CY5 = [-0.25 0.1; 0.2 1.5; 0.1 2; 0.1 2; 0.1 2; 0.1 2; 0.6 2; 0.6 2]; 
%T_delRn = [-2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10]; %no filter
edges_delRn_CY5 = linspace(-0.5,1,nbins);




%PMT gain correction
gain_correct = 'N'; %Y/N
m = 4.8572; %PMT gain slope (log10 scale)
b = 1.4367; %PMT gain intercept (log10 scale)
G_v_FAM_measured = 0.66; %V, PMT GAIN for measurements
G_v_FAM_base = 0.55; %V, PMT GAIN from standard curve
G_FAM_measured = 10^(m*G_v_FAM_measured + b); %convert PMT gain voltage to true gain (measured gain)
G_FAM_base = 10^(m*G_v_FAM_base + b); %convert PMT gain voltage to true gain (standard curve gain)

G_v_ROX_measured = 0.54; %V, PMT GAIN for measurements
G_v_ROX_base = 0.6; %V, PMT GAIN from standard curve
G_ROX_measured = 10^(m*G_v_ROX_measured + b); %convert PMT gain voltage to true gain (measured gain)
G_ROX_base = 10^(m*G_v_ROX_base + b); %convert PMT gain voltage to true gain (standard curve gain)

%% Load Data

%detection data
DETECTION = load(dataname);
FAM = DETECTION.FAM;
ROX = DETECTION.ROX;
CY5 = DETECTION.CY5;
TIME = DETECTION.Time_Inside;
cycle = DETECTION.Cycles;
C = length(cycle); %number of cycles

%% Post Processing

if gain_correct == 'Y'
    
    ratio_gain_FAM(1:3) = G_FAM_base/G_FAM_measured;
    ratio_gain_FAM(4:9) = 1;
    
    ratio_gain_ROX(1:3) = G_ROX_base/G_ROX_measured; %only 1 through 3 adjusted
    ratio_gain_ROX(4:9) = 1; 
    
else
    
    ratio_gain_FAM(1:9) = 1;
    ratio_gain_ROX(1:9) = 1;
    
end


%filter data
for i = 1 : C
    
    %load temp array
    FAM_temp_1 = FAM{i};
    FAM_temp_1 = FAM_temp_1*ratio_gain_FAM(i); %correct for change in PMT gain
    FAM_PMT_correct{i} = FAM_temp_1;
    CY5_temp_1 = CY5{i};
    ROX_temp_1 = ROX{i};
    ROX_temp_1 = ROX_temp_1*ratio_gain_ROX(i);
    TIME_temp_1 = TIME{i};
    drop_count_unfiltered(i) = length(FAM_temp_1);
    
%sort drops based on time inside
    FAM_temp_2 = FAM_temp_1(TIME_temp_1 >= T_time(i,1));
    CY5_temp_2 = CY5_temp_1(TIME_temp_1 >= T_time(i,1));
    ROX_temp_2 = ROX_temp_1(TIME_temp_1 >= T_time(i,1));
    TIME_temp_2 = TIME_temp_1(TIME_temp_1 >= T_time(i,1));
    
    FAM_temp_3 = FAM_temp_2(TIME_temp_2 < T_time(i,2));
    CY5_temp_3 = CY5_temp_2(TIME_temp_2 < T_time(i,2));
    ROX_temp_3 = ROX_temp_2(TIME_temp_2 < T_time(i,2));
    TIME_temp_3 = TIME_temp_2(TIME_temp_2 < T_time(i,2));
    
    %sort drops based on ROX
	FAM_temp_4 = FAM_temp_3(ROX_temp_3 >= T_ROX(i,1));
    CY5_temp_4 = CY5_temp_3(ROX_temp_3 >= T_ROX(i,1));
    TIME_temp_4 = TIME_temp_3(ROX_temp_3 >= T_ROX(i,1));
    ROX_temp_4 = ROX_temp_3(ROX_temp_3 >= T_ROX(i,1));
    
	FAM_temp_5 = FAM_temp_4(ROX_temp_4 < T_ROX(i,2));
    CY5_temp_5 = CY5_temp_4(ROX_temp_4 < T_ROX(i,2));
    TIME_temp_5 = TIME_temp_4(ROX_temp_4 < T_ROX(i,2));
    ROX_temp_5 = ROX_temp_4(ROX_temp_4 < T_ROX(i,2));
      
    FAM_corr{i} = FAM_temp_5;
    FAM_avg(i) = mean(FAM_temp_5);
    FAM_std(i) = std(FAM_temp_5);
    
    CY5_corr{i} = CY5_temp_5;
    CY5_avg(i) = mean(CY5_temp_5);
    CY5_std(i) = std(CY5_temp_5);
    
    ROX_corr{i} = ROX_temp_5;
    ROX_avg(i) = mean(ROX_temp_5);
    ROX_std(i) = std(ROX_temp_5);
    
    
    %Rn (FAM)
    Rn_FAM{i} = FAM_temp_5./ROX_temp_5;
    Rn_FAM_avg(i) = mean(FAM_temp_5./ROX_temp_5);
    Rn_FAM_std(i) = std(FAM_temp_5./ROX_temp_5);
    
    
    %delRn (FAM)
    delRn_FAM{i} = (FAM_temp_5./ROX_temp_5) - Rn_FAM_avg(1);
    delRn_FAM_avg(i) = mean((FAM_temp_5./ROX_temp_5) - Rn_FAM_avg(1));
    delRn_FAM_std(i) = std((FAM_temp_5./ROX_temp_5) - Rn_FAM_avg(1));
    
    %Rn (Cy5)
    Rn_CY5{i} = CY5_temp_5./ROX_temp_5;
    Rn_CY5_avg(i) = mean(CY5_temp_5./ROX_temp_5);
    Rn_CY5_std(i) = std(CY5_temp_5./ROX_temp_5);
    
    
    %delRn (Cy5)
    delRn_CY5{i} = (CY5_temp_5./ROX_temp_5) - Rn_CY5_avg(1);
    delRn_CY5_avg(i) = mean((CY5_temp_5./ROX_temp_5) - Rn_CY5_avg(1));
    delRn_CY5_std(i) = std((CY5_temp_5./ROX_temp_5) - Rn_CY5_avg(1));
    
    
    %filter delRn (FAM)
    delRn_FAM_temp_1 = delRn_FAM{i};
    
    delRn_FAM_temp_2 = delRn_FAM_temp_1(delRn_FAM_temp_1 > T_delRn_FAM(i,1));
    delRn_FAM_temp_3 = delRn_FAM_temp_2(delRn_FAM_temp_2 < T_delRn_FAM(i,2));
    
    delRn_FAM_filter{i} = delRn_FAM_temp_3;
    delRn_FAM_filter_avg(i) = mean(delRn_FAM_temp_3);
    delRn_FAM_filter_std(i) = std(delRn_FAM_temp_3);
    
    
    %filter delRn (Cy5)
    delRn_CY5_temp_1 = delRn_CY5{i};
    
    delRn_CY5_temp_2 = delRn_CY5_temp_1(delRn_CY5_temp_1 > T_delRn_CY5(i,1));
    delRn_CY5_temp_3 = delRn_CY5_temp_2(delRn_CY5_temp_2 < T_delRn_CY5(i,2));
    
    delRn_CY5_filter{i} = delRn_CY5_temp_3;
    delRn_CY5_filter_avg(i) = mean(delRn_CY5_temp_3);
    delRn_CY5_filter_std(i) = std(delRn_CY5_temp_3);
    
end

%percent amplified
for i = 1 : C;
    
    P_amp_FAM(i) = length(delRn_FAM_filter{i})/length(delRn_FAM{i})*100;
    P_amp_CY5(i) = length(delRn_CY5_filter{i})/length(delRn_CY5{i})*100;
    
end




%find expected percent positive drops for all cycles
for i = 1 : C
    
        
    pos_FAM_temp = delRn_FAM_filter{i};
    pos_CY5_temp = delRn_CY5_filter{i};
        
    %end
    
    %save positive drop threshold
    delRn_FAM_FINAL{i} = pos_FAM_temp;
    delRn_FAM_avg_FINAL(i) = mean(delRn_FAM_FINAL{i});
    delRn_FAM_std_FINAL(i) = std(delRn_FAM_FINAL{i});
    
    delRn_CY5_FINAL{i} = pos_CY5_temp;
    delRn_CY5_avg_FINAL(i) = mean(delRn_CY5_FINAL{i});
    delRn_CY5_std_FINAL(i) = std(delRn_CY5_FINAL{i});

end

savename = ['processed_delRn_' dataname];

save(savename,'delRn_FAM_avg_FINAL','delRn_FAM_std_FINAL','cycle','delRn_FAM_FINAL',...
    'delRn_CY5_avg_FINAL','delRn_CY5_std_FINAL','delRn_CY5_FINAL','FAM_PMT_correct',...
    'ROX','TIME','FAM_corr','ROX_corr','delRn_FAM','delRn_CY5')



 
%% Figures

%colors
red = linspecer('red');
blue = linspecer('blue');
red = linspecer('red');

%subplot of detection delRn (FAM) histograms
figure(1); clf(1)

hold on

for i = 1 : C
    
    subplot(C,1,i)
    histogram(delRn_FAM{i},edges_delRn_FAM)
    xline(T_delRn_FAM(i,1),'-k','linewidth',1);
    xline(T_delRn_FAM(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_delRn_FAM) max(edges_delRn_FAM) 0 inf])
    
end

xlabel('FAM \DeltaRn (a.u.)')
%xlabel('Intensity (a.u.)')
suptitle('FAM delRn Detection Data')

hold off





%subplot of detection time inside histograms
figure(2); clf(2)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    histogram(TIME{i},edges_time)
    xline(T_time(i,1),'-k','linewidth',1);
    xline(T_time(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_time) max(edges_time) 0 inf])
    
end

xlabel('Time Inside (ms)')
suptitle('Time Inside Detection Data')

hold off





%subplot of detection ROX histograms
figure(3); clf(3)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    histogram(ROX{i}*ratio_gain_ROX(i),edges_ROX)
    xline(T_ROX(i,1),'-k','linewidth',1);
    xline(T_ROX(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_ROX) max(edges_ROX) 0 inf])
    
end

xlabel('Fluoresence Intensity (a.u.)')
suptitle('ROX Detection Data')

hold off





%subplot of detection FAM histograms
figure(4); clf(4)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    histogram(FAM_PMT_correct{i},edges_FAM)
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_FAM) max(edges_FAM) 0 inf])
    
end

xlabel('Fluoresence Intensity (a.u.)')
suptitle('FAM Detection Data')

hold off





%subplot of detection CY5 histograms
figure(5); clf(5)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    histogram(CY5{i},edges_CY5)
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_CY5) max(edges_CY5) 0 inf])
    
end

xlabel('Fluoresence Intensity (a.u.)')
suptitle('CY5 Detection Data')

hold off





%subplot of detection delRn (Cy5) histograms
figure(6); clf(6)

hold on

for i = 1 : C
    
    subplot(C,1,i)
    histogram(delRn_CY5{i},edges_delRn_CY5)
    xline(T_delRn_CY5(i,1),'-k','linewidth',1);
    xline(T_delRn_CY5(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_delRn_CY5) max(edges_delRn_CY5) 0 inf])
    
end

xlabel('CY5 \DeltaRn (a.u.)')
%xlabel('Intensity (a.u.)')
suptitle('CY5 delRn Detection Data')

hold off













