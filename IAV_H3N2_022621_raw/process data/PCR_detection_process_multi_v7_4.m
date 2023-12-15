%Process PCR detection data for MOI 0.1 (H3N2) burst size IAV drops
%3-7-21 v7.4 update ROX and TIME threshold to work on post FAM filtered
%data + more conservative delRn C1 (or by %) and C40 avg +/- 3std
%delRn data so data is linked
%Geoff Zath

%delRn normalization

%IAV drops 2-26-21

%%to do list
%Otso threshold for std (done)  :)
%mean/mode of ROX/TIME data +/- X% (done)
%option for each data for manual or auto (done)


clear; clc

%% Inputs
dataname = 'detection_data_IAV_022621.mat'; %rename me
mgene_data = load('processed_delRn_detection_data_std_022621.mat'); %rename me

nbins = 100; %number of bins for histograms (100)
% Z = 3.392; %Z value for CI calc at cycle 1 (95% = 1.960, 99% = 2.576. 99.5% = 2.807, 99.9% = 3.392)



%CY5 theshold
edges_CY5 = linspace(0,1,nbins);



%FAM thresholds by drop #
FAM_auto = 'Y'; %auto filter = Y, manual filter = N
CV_auto = 'Y'; %auto CV filter = Y, manual filter = N
delX = 10; %delta X for drop# filter (10)
nbins_Otsu = 50; %for Otsu thresholding (50)
X_Otsu = [0.5 1 1 1 1 1]; %change Otsu threshold at cycle i (counter bias from noisy drops, 0.5 at cycle 1)
T_CV = [1, 0.3, 0.3, 0.3, 0.3, 0.3]; %manual CV threshold for automatic filtering
T_FAM = [1 0.5e4; 1 0.5e4; 1 0.5e4; 1 0.5e4; 1 0.5e4; 1 0.5e4]; %manual FAM threshold
edges_FAM = linspace(0,0.3,nbins);



%time inside thresholds
TIME_auto = 'Y'; %auto filter = Y, manual filter = N
P_TIME = 0.2; %+/- median range filter on TIME (0.2)
T_time = [4 6; 4.5 6.5; 4.5 6.5; 4.5 6.5; 4.5 6.5; 4.5 6.5]; %manual TIME threshold
%T_time = [0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10]; %no filter
edges_time = linspace(0,10,nbins);



%ROX thresholds
ROX_auto = 'Y'; %auto filter = Y, manual filter = N
P_ROX = 0.2; %+/- median range filter on ROX (0.2)
T_ROX = [0.6 1.2; 0.8 1.4; 0.8 1.4; 0.8 1.4; 0.8 1.4; 0.8 1.4; 0.8 1.4; 0.8 1.4]; %manual ROX threshold
%T_ROX = [0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10]; %no filter
edges_ROX = linspace(0,2,nbins);



%delRn (FAM) thresholds (final cleanup) (min set to ~1% positive in cycle
%1) (max set to C40 avg of mgene drops - one std)
delRn_auto = 'Y'; %auto filter = Y, manual filter = N
P_delRn_C1 = 0.005; %fraction of positive drops at cycle 1 used for min thresholed (0.005)
N_std_C1 = 3; %number of std above mean

%calculate max threshold
N_std_C40 = 3; %number of std below mean
mgene_delRn_avg = mgene_data.delRn_FAM_avg_FINAL;
mgene_delRn_std = mgene_data.delRn_FAM_std_FINAL;
T_C40 = mgene_delRn_avg(end) - N_std_C40*mgene_delRn_std(end);
T_delRn_FAM = [0.025 T_C40; 0.025 T_C40; 0.025 T_C40; 0.025 T_C40; 0.025 T_C40; 0.025 0.5]; %manual delRn threshold
%T_delRn = [-2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10]; %no filter
edges_delRn_FAM = linspace(-0.05,0.25,nbins);



%delRn (CY5) thresholds (final cleanup)
T_delRn_CY5 = [-0.25 0.25; -0.25 0.5; -0.25 0.5; -0.25 1; -0.25 1; -0.25 1.5; -0.25 1.5; -0.25 2]; %manual delRn threshold
%T_delRn = [-2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10; -2 10]; %no filter
edges_delRn_CY5 = linspace(-0.05,0.5,nbins);



%PMT gain correction (only use if gain difference between cycles)
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




%automatic drop# filter
for i = 1 : C

    raw_temp = FAM{i};
    
    L_raw = length(raw_temp);
    
    N_delx = floor(L_raw/delX); %number of sections
    
    count = 1;
    
    %calculate CV within each range
    for j = 1 : N_delx
        
        %section range
        start = 1 + delX*(j-1);
        stop = start + delX;
        
        std_section(i,j) = std(raw_temp(start:stop));
        avg_section(i,j) = mean(raw_temp(start:stop));
        CV_section(i,j) = std_section(i,j)/avg_section(i,j);
        
        start_array(i,j) = start;
        stop_array(i,j) = stop;
        
%         %save sections above threshold
%         if CV_section(i,j) >= T_CV(i)
%             
%             start_keep(count) = start;
%             stop_keep(count) = stop;
%             
%             count = count + 1;
%         
%         end
        
    end
    
    %find Otsu threshold
    temp_CV = CV_section(i,:);
    temp_CV = temp_CV(temp_CV>0);
    %N_CV = histcounts(temp_CV,nbins_Otsu); %bin data
    %T_CV_Otsu(i) = otsuthresh(N_CV); %calculate bin# threshold
    
    N_CV = histcounts(log10(temp_CV),nbins_Otsu); %bin data
    T_CV_Otsu(i) = otsuthresh(N_CV); %calculate bin# threshold
    
    %T_CV_Otsu(i) = graythresh(temp_CV);
    %T_CV_Otsu(i) = graythresh(log10(temp_CV));
    T2_CV_Otsu(i) = 10^(min(log10(temp_CV)) + T_CV_Otsu(i)*range(log10(temp_CV)));
    %T2_CV_Otsu(i) = min(temp_CV) + T_CV_Otsu(i)*range(temp_CV);
    
    %adjust cycle i threshold
    T2_CV_Otsu(i) = T2_CV_Otsu(i)*X_Otsu(i);
    
    if CV_auto == 'Y'
            
        CV_thresh(i) = T2_CV_Otsu(i);
        
    else
        
        CV_thresh(i) = T_CV(i);
        
    end

    
    %save section above Otsu or manual threshold
    for j = 1 : N_delx
        
        if CV_section(i,j) >= CV_thresh(i)
            
            start_keep(count) = start_array(i,j);
            stop_keep(count) = stop_array(i,j);
            
            count = count + 1;
        
        end
        
    end

    
    
    start_cell{i} = start_keep;
    stop_cell{i} = stop_keep;
    
    clear count start_keep stop_keep
     
    
    
end




% %automatic ROX and TIME filter
% for i = 1 : C
%     
%     %ROX
%     if ROX_auto == 'Y'
%         
%         ROX_median(i) = median(ROX{i});
%         T_ROX(i,1) = ROX_median(i)*(1 - P_ROX);
%         T_ROX(i,2) = ROX_median(i)*(1 + P_ROX);
%     
%     end
%     
%     %time inside
%     if TIME_auto == 'Y'
%         
%         TIME_median(i) = median(TIME{i});
%         T_time(i,1) = TIME_median(i)*(1 - P_TIME);
%         T_time(i,2) = TIME_median(i)*(1 + P_TIME);
%     
%     end
%     
%     
% end
        


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
    
    if FAM_auto == 'Y'
    
        %sort drops based on drop# (auto filter)
        %make shaded regions of filtered sections
        L = length(start_cell{i});

        start_temp = start_cell{i};
        stop_temp = stop_cell{i};

        %cell_temp = [];

        %filter by sections above CV threshold
        for j = 1 : L

            FAM_cell_temp{j} = FAM_temp_1(start_temp(j):stop_temp(j));
            CY5_cell_temp{j} = CY5_temp_1(start_temp(j):stop_temp(j));
            ROX_cell_temp{j} = ROX_temp_1(start_temp(j):stop_temp(j));
            TIME_cell_temp{j} = TIME_temp_1(start_temp(j):stop_temp(j));

        end

        %rearrange data
        FAM_temp_1 = [FAM_cell_temp{:}];
        FAM_temp_1 = FAM_temp_1(:);

        CY5_temp_1 = [CY5_cell_temp{:}];
        CY5_temp_1 = CY5_temp_1(:);

        ROX_temp_1 = [ROX_cell_temp{:}];
        ROX_temp_1 = ROX_temp_1(:);

        TIME_temp_1 = [TIME_cell_temp{:}];
        TIME_temp_1 = TIME_temp_1(:);

        clear FAM_cell_temp CY5_cell_temp ROX_cell_temp TIME_cell_temp
        
    else

        %sort drops based on drop# (manual filter)
        FAM_temp_1 = FAM_temp_1(T_FAM(i,1):T_FAM(i,2));
        CY5_temp_1 = CY5_temp_1(T_FAM(i,1):T_FAM(i,2));
        ROX_temp_1 = ROX_temp_1(T_FAM(i,1):T_FAM(i,2));
        TIME_temp_1 = TIME_temp_1(T_FAM(i,1):T_FAM(i,2));
    
    end
    
    %save FAM filtered data for graphing
    TIME_FAMfilter{i} = TIME_temp_1;
    ROX_FAMfilter{i} = ROX_temp_1;
    CY5_FAMfilter{i} = CY5_temp_1;
    FAM_FAMfilter{i} = FAM_temp_1;
    
    
    
    %automatic ROX and TIME filter
    %ROX
    if ROX_auto == 'Y'
        
        %ROX_median(i) = median(ROX{i});
        ROX_median(i) = median(ROX_FAMfilter{i});
        T_ROX(i,1) = ROX_median(i)*(1 - P_ROX);
        T_ROX(i,2) = ROX_median(i)*(1 + P_ROX);
    
    end
    
    %time inside
    if TIME_auto == 'Y'
        
        %TIME_median(i) = median(TIME{i});
        TIME_median(i) = median(TIME_FAMfilter{i});
        T_time(i,1) = TIME_median(i)*(1 - P_TIME);
        T_time(i,2) = TIME_median(i)*(1 + P_TIME);
    
    end
    
    
    
    
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
    
    %filter delRn by manual or auto
    if delRn_auto == 'Y'
        
        %set min to ~1% positive in cycle 1
        L_C1 = length(delRn_FAM{1});
        L_C1_T = floor(L_C1*P_delRn_C1);
        delRN_sort = flip(sort(delRn_FAM{1}));
        T_delRn_FAM(i,1) = delRN_sort(L_C1_T);

%         %set min to avg + N*std
%         T_delRn_FAM(i,1) = delRn_FAM_avg(1) + N_std_C1*delRn_FAM_std(1);
        
        %set max to mgene C40 avg - std
        T_delRn_FAM(i,2) = T_C40;

        %set C40 to manual value (all data)
        if i == C

            T_delRn_FAM(i,2) = T_delRn_FAM(end,end);

        end
        
    end
    
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
    
    
    %filter delRn (Cy5) by FAM
    delRn_CY5_FAM_temp_1 = delRn_CY5{i};
    
    delRn_CY5_FAM_temp_2 = delRn_CY5_FAM_temp_1(delRn_FAM_temp_1 > T_delRn_FAM(i,1));
    delRn_CY5_FAM_temp_3 = delRn_CY5_FAM_temp_2(delRn_FAM_temp_2 < T_delRn_FAM(i,2));
    
    delRn_CY5_FAM_filter{i} = delRn_CY5_FAM_temp_3;
    delRn_CY5_FAM_filter_avg(i) = mean(delRn_CY5_temp_3);
    delRn_CY5_FAM_filter_std(i) = std(delRn_CY5_temp_3);
    
    
end





%percent amplified
for i = 1 : C
    
    P_amp_FAM(i) = length(delRn_FAM_filter{i})/length(delRn_FAM{i})*100;
    P_amp_CY5(i) = length(delRn_CY5_filter{i})/length(delRn_CY5{i})*100;
    N_T(i) = length(delRn_FAM{i});
    N_P(i) = length(delRn_FAM_filter{i});
    
end



%save positive drops for all cycles
for i = 1 : C
    
        
    pos_FAM_temp = delRn_FAM_filter{i};
    pos_CY5_temp = delRn_CY5_filter{i};
    pos_CY5_FAM_temp = delRn_CY5_FAM_filter{i};
        
    %end
    
    %save positive drop threshold
    delRn_FAM_FINAL{i} = pos_FAM_temp;
    delRn_FAM_avg_FINAL(i) = mean(delRn_FAM_FINAL{i});
    delRn_FAM_std_FINAL(i) = std(delRn_FAM_FINAL{i});
    
    delRn_CY5_FINAL{i} = pos_CY5_temp;
    delRn_CY5_avg_FINAL(i) = mean(delRn_CY5_FINAL{i});
    delRn_CY5_std_FINAL(i) = std(delRn_CY5_FINAL{i});
    
    delRn_CY5_FAM_FINAL{i} = pos_CY5_FAM_temp;

end

savename = ['auto_processed_delRn_' dataname];

save(savename,'delRn_FAM_avg_FINAL','delRn_FAM_std_FINAL','cycle','delRn_FAM_FINAL',...
    'delRn_CY5_avg_FINAL','delRn_CY5_std_FINAL','delRn_CY5_FINAL','FAM_PMT_correct',...
    'ROX','TIME','FAM_corr','ROX_corr','delRn_CY5_FAM_FINAL','delRn_FAM','delRn_CY5')




 
%% Figures

%colors
red = linspecer('red');
blue = linspecer('blue');
green = linspecer('green');

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

xlabel('ROX Intensity (a.u.)')
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

xlabel('FAM Intensity (a.u.)')
suptitle('FAM Detection Data')

hold off






%subplot of filtered detection delRn (FAM) histograms
figure(5); clf(5)

hold on

for i = 1 : C
    
    subplot(C,1,i)
    histogram(delRn_FAM_filter{i},edges_delRn_FAM)
    %xline(T_delRn(i,1),'-k','linewidth',1);
    %xline(T_delRn(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_delRn_FAM) max(edges_delRn_FAM) 0 inf])
    
end

xlabel('FAM \DeltaRn (a.u.)')
suptitle('filtered FAM delRn Detection Data')

hold off






%subplot of detection CY5 histograms
figure(6); clf(6)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    histogram(CY5{i},edges_CY5)
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_CY5) max(edges_CY5) 0 inf])
    
end

xlabel('Cy5 Intensity (a.u.)')
suptitle('CY5 Detection Data')

hold off





%subplot of detection delRn (Cy5) histograms
figure(7); clf(7)

hold on

for i = 1 : C
    
    subplot(C,1,i)
    histogram(delRn_CY5{i},edges_delRn_CY5)
    %xline(T_delRn_CY5(i,1),'-k','linewidth',1);
    %xline(T_delRn_CY5(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_delRn_CY5) max(edges_delRn_CY5) 0 inf])
    
end

xlabel('CY5 \DeltaRn (a.u.)')
%xlabel('Intensity (a.u.)')
suptitle('CY5 delRn Detection Data')

hold off





%subplot of filtered detection delRn (Cy5) histograms
figure(8); clf(8)

hold on

for i = 1 : C
    
    subplot(C,1,i)
    histogram(delRn_CY5_FAM_filter{i},edges_delRn_CY5)
    %xline(T_delRn(i,1),'-k','linewidth',1);
    %xline(T_delRn(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Drop Count')
    axis([min(edges_delRn_CY5) max(edges_delRn_CY5) 0 inf])
    
end

xlabel('CY5 \DeltaRn (a.u.)')
suptitle('filtered CY5 delRn Detection Data')

hold off





%subplot of detection time inside scatter plot
figure(9); clf(9)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    plot(TIME_FAMfilter{i},'.k','markersize',1)
    %plot(TIME{i},'.k','markersize',1)
    %yline(T_time(i,1),'-r','linewidth',1);
    %yline(T_time(i,2),'-r','linewidth',1);
    %xline(T_FAM(i,1),'-r','linewidth',1);
    %xline(T_FAM(i,2),'-r','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('Time Inside (ms)')
    axis([0 inf min(edges_time) max(edges_time)])   
    
    %draw box
    xmin = 1;
    %xmax = length(TIME{i});
    xmax = length(TIME_FAMfilter{i});    
    
    x = [xmin xmax xmax xmin];
        
    ymin = min(T_time(i,1));
    ymax = max(T_time(i,2));
        
    y = [ymin ymin ymax ymax];
        
    patch(x,y,blue(96,:),'facealpha',0.25,'edgecolor','none')
    
    
end

xlabel('Drop #')
suptitle('Time Inside Detection Data')

hold off





%subplot of detection ROX scatter plot
figure(10); clf(10)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    %plot(ROX{i}*ratio_gain_ROX(i),'.k','markersize',1)
    plot(ROX_FAMfilter{i}*ratio_gain_ROX(i),'.k','markersize',1)
    %yline(T_ROX(i,1),'-r','linewidth',1);
    %yline(T_ROX(i,2),'-r','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('ROX Intensity (a.u.)')
    axis([0 inf min(edges_ROX) max(edges_ROX)])
    
    %draw box
    xmin = 1;
    %xmax = length(ROX{i});
    xmax = length(ROX_FAMfilter{i});
        
    x = [xmin xmax xmax xmin];
        
    ymin = min(T_ROX(i,1));
    ymax = max(T_ROX(i,2));
        
    y = [ymin ymin ymax ymax];
        
    patch(x,y,red(96,:),'facealpha',0.25,'edgecolor','none')
    
end

xlabel('Drop #')
suptitle('ROX Detection Data')

hold off





%subplot of detection FAM scatter plot
figure(11); clf(11)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    plot(FAM_PMT_correct{i},'.k','markersize',1)
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('FAM Intensity (a.u.)')
    axis([0 inf min(edges_FAM) max(edges_FAM)])
    
    
    if FAM_auto == 'Y'
        
        %make shaded regions of filtered sections
        L = length(start_cell{i});

        start_temp = start_cell{i};
        stop_temp = stop_cell{i};

        %draw boxes
        for j = 1 : L

            xmin = start_temp(j);
            xmax = stop_temp(j);

            x = [xmin xmax xmax xmin];

            ymin = min(edges_FAM);
            ymax = max(edges_FAM);

            y = [ymin ymin ymax ymax];

            patch(x,y,green(96,:),'facealpha',0.25,'edgecolor','none')

        end
        
    else
        
        xline(T_FAM(i,1),'--g','linewidth',1);
        xline(T_FAM(i,2),'--g','linewidth',1);
        
    end
    
    
  
    
    
end

xlabel('Drop #')
suptitle('FAM Detection Data')

hold off







%subplot of detection CY5 scatter plot
figure(12); clf(12)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    plot(CY5{i},'.k','markersize',1)
    %xline(T_FAM(i,1),'-k','linewidth',1);
    %xline(T_FAM(i,2),'-k','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('CY5 Intensity (a.u.)')
    axis([0 inf min(edges_CY5) max(edges_CY5)])
    
end

xlabel('Drop #')
suptitle('CY5 Detection Data')

hold off






%subplot of FAM std from sectioning
figure(13); clf(13)

hold on


for i = 1 : C
    
    subplot(C,1,i)
    plot(CV_section(i,:),'.k','markersize',3,'linewidth',1)
    %xline(T_FAM(i,1),'-k','linewidth',1);
    %xline(T_FAM(i,2),'-k','linewidth',1);
    %yline(T_CV(i),'-r','linewidth',1);
    yline(CV_thresh(i),'--r','linewidth',1);
    title(['Cycle ',num2str(cycle(i),'%2.0f')])
    ylabel('FAM CV (a.u.)')
    axis([0 inf 0 inf])
    
    set(gca,'yscale','log')
    
end

xlabel('Section #')
suptitle('FAM std sectioning')

hold off












