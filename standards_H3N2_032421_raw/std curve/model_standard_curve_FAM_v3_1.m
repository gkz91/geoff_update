%qPCR model standard curve based on efficiency model fit to drop detection
%data and efficiency model (no dilution series)
%1-22-21 v3.1_multi updated for FAM channel + cleaned code
%Geoff Zath

%3-24-21
%delRn normalization (FAM)
%M gene 10^2 cpd in 100 um drops (10^6 copies/uL pre-split)



clear; clc

%% Inputs

filename = 'processed_delRn_detection_data_std_032421.mat';

N_cycles = 40; %total number of thermal cycles
N_span = 1000; %number of model curves (increases resolution)
E_PCR = 0.951; %PCR efficiency fraction 
conc_fit = 171; %100 cpd fit post split (171 cpd) 
Rn_range = [0.5 2]; %range of linear region
cycle_graph = 28; %single cycle for figure 3
split_ratio = 1; %split ratio of EVO chip


%efficiency fitting
%N_fit = 31; %number of data points measured for fit
F_0 = 1; %initial fluorescence at cycle 1 (arbitrary, can't be zero)
a_eff_range = [0.00 0.1]; %initial guess of min efficiency
b_eff_range = [1 2]; %intitial guess of slope
c_eff_range = [21 23]; %initial guess of C_0.5 range
d_eff_range = [.9 1]; %initial guess of max efficiency
cycle_thresh = 1; %set cycle to start fitting curve
N_model_eff = 30; %mesh size (max at 30, or else PC will hurt)
cycle = 1:40; %PCR thermal cycles

%apply std curve
N_cycles_span = 100; %cycle range to calculate for model span
use_E_fit = 'Y'; %use estimated efficiency from curve fitting (Y/N), if no uses efficiency from dilution series
E_fit_offset = -3; %offset from C_0.5 to calculate PCR efficiency for model span
E_PCR_scale = 1; %scale PCR effiency calculated by dilution series


%% Load Data

%fitting data
PCR_curve = load(filename);
D_avg_curve = PCR_curve.delRn_FAM_avg_FINAL;
D_std_curve = PCR_curve.delRn_FAM_std_FINAL;
cycle_curve = PCR_curve.cycle;
C_curve = length(cycle_curve);




%% Best Fit PCR Curve (efficiency based PCR curve)
%sigmoidal model from "Sigmoidal curve-fitting redefines quantitative
%real-time PCR with the prospective of developing
%automated high-throughput applications"

%QS3 data model initial guess

a_eff = linspace(a_eff_range(1),a_eff_range(2),N_model_eff); %min efficiency
%a_eff = 0; %min efficiency
b_eff = linspace(b_eff_range(1),b_eff_range(2),N_model_eff); %slope
c_eff = linspace(c_eff_range(1),c_eff_range(2),N_model_eff); %C(1/2)
d_eff = linspace(d_eff_range(1),d_eff_range(2),N_model_eff); %max efficiency
%d_eff = 1; %max efficiency

cycle40 = 1:N_cycles;

y_bar = mean(D_avg_curve);
SS_total = sum((D_avg_curve - y_bar).^2);

%find best fit with minimum R^2

model_iter = zeros(N_cycles,1);
model_R = zeros(C_curve,1);
SS_res_eff = zeros(N_model_eff,N_model_eff);
R_sq_eff = SS_res_eff;

for m = 1 : N_model_eff %a
    for k = 1 : N_model_eff %b
        for j = 1 : N_model_eff %c
            for n = 1 : N_model_eff %d

                F_temp = zeros(N_cycles,1);
                F_temp(1) = F_0;

                for i = 1 : N_cycles-1 %all cycles modeled
                    d_n = d_eff(n);
                    E_rev_temp = d_eff(n)./(1 + exp((cycle40(i) - c_eff(j))/b_eff(k))) + a_eff(m); %efficiency decreasing with time

                    %model qPCR amp curve based on PCR efficiency
                    F_temp(i+1) = F_temp(i)*(1 + E_rev_temp);

                    model_iter(i+1) = F_temp(i+1);

                end

        
                for h = 1 : C_curve %cycles measured

                    model_R(h) = model_iter(cycle_curve(h));

                end
        

                F_model_scale = D_avg_curve(end)/model_iter(end); %scale fluorescence
                model_iter_scaled = model_iter*F_model_scale; %scale fluorescence
                model_R_scaled = model_R*F_model_scale;
                %model_thresh_scaled = model_thresh*F_model_scale; %scale fluorescence
                %model_iter_baseline = model_iter_scaled - F_0; %baseline

                SS_res(m,k,j,n) = sum((D_avg_curve' - model_R_scaled).^2);
                R_sq_eff(m,k,j,n) = 1 - SS_res(m,k,j,n)/SS_total;
              
            end
        end
    end
end

%find best fit model parameter location

%find max R_sq varying a
[M_a I_a] = max(R_sq_eff,[],1);
M_a = squeeze(M_a);
I_a = squeeze(I_a);

%find max R_sq varying b
[M_b I_b] = max(M_a,[],1);
M_b = squeeze(M_b);
I_b = squeeze(I_b);

%find max R_sq varying c
[M_c I_c] = max(M_b);

%fine max R_sq varying d
[M_d I_d] = max(M_c);

R_sq_eff_best = M_d;

%find location of best parameters
idx_d = I_d;
idx_c = I_c(idx_d);
idx_b = I_b(idx_c,idx_d);
idx_a = I_a(idx_b,idx_c,idx_d);

%Best fit model of bulk standard curve (max conc)
a_eff_best = a_eff(idx_a); %min efficiency
b_eff_best = b_eff(idx_b); %
c_eff_best = c_eff(idx_c); %
d_eff_best = d_eff(idx_d); %max efficiency


%best fit PCR efficiency
E_best_rev = d_eff_best./(1 + exp((cycle - c_eff_best)/b_eff_best)) + a_eff_best; %efficiency decreasing with time 

%model qPCR amp curve based on best fit PCR efficiency
F_best = zeros(N_cycles,1);
F_best(1) = F_0;

for i = 1 : N_cycles-1
    
    F_best(i+1) = F_best(i)*(1 + E_best_rev(i));
    
end

F_best_scale = D_avg_curve(end)/F_best(end);
F_best_scaled = F_best*F_best_scale;
model_best = F_best_scaled;

%span curve across all cycles
a_span_eff = a_eff_best; %min efficiency
b_span_eff = b_eff_best; 
c_span_eff = linspace(1,40,N_span);
d_span_eff = d_eff_best; %max efficiecny

%k_eff = dsearchn(c_span_eff,c_eff_best); %find location of C_0.5 in span range
cycle_span = 1:N_cycles_span; %cycles to use for spanning model
model_eff_span = zeros(N_cycles_span,N_span);
for j = 1 : N_span %curve
    
    F_temp = zeros(N_cycles_span,1);
    F_temp(1) = F_0;
    model_eff_span(1,j) = F_temp(1);
    
    for i = 1 : N_cycles_span-1 %cycles

        %model_span(i,j) = a_span + d_span/(1 + exp(-((cycle40(i) - c_span(j))/b_span)));
        
        %E_rev_temp = d_span_eff./(1 + exp((cycle40(i) - c_span_eff(j))/b_span_eff)) + a_span_eff; %efficiency decreasing with time
        E_rev_temp = d_span_eff./(1 + exp((cycle_span(i) - c_span_eff(j))/b_span_eff)) + a_span_eff; %efficiency decreasing with time

        %model qPCR amp curve based on PCR efficiency
        F_temp(i+1) = F_temp(i)*(1 + E_rev_temp);
        
        model_eff_span(i+1,j) = F_temp(i+1);

    end
    
    %not easy!!!
    %calculate C' (virtual cycle that intercects with C=40 of fitted line)
    C_prime(j) = 40 - (c_eff_best - c_span_eff(j));
    
    %linear interpolation to find F' (virtual fluoresecnce that intercects
    %with F(C=40) of fitted line
    C_i = floor(C_prime(j));
    C_ii = ceil(C_prime(j));
    
    %fix if C_prime is an integer
    if C_i == C_ii
        
        C_ii = C_i + 1;
        
    end
    
    F_i = model_eff_span(C_i,j);
    F_ii = model_eff_span(C_ii,j);
    
    F_prime(j) = F_i + (C_prime(j) - C_i)*((F_ii - F_i)/(C_ii - C_i)); %linear interpolation
    
    %scale data at virtual cycle
    F_prime_scale(j) = D_avg_curve(end)/F_prime(j);
    model_eff_scaled_span(:,j) = F_prime_scale(j)*model_eff_span(:,j);
    

    F_span_scale(j) = D_avg_curve(end)/model_eff_span(end,j);
    %model_eff_scaled_span(:,j) = F_span_scale(j)*model_eff_span(:,j);
    %model_eff_scaled_span(:,j) = F_best_scale*model_eff_span(:,j);
  
end


%standard curve calculation

%choose which PCR efficiency to use
if use_E_fit == 'Y'
    
    idx_offset = floor(c_eff_best + E_fit_offset);
    E_stdc = E_best_rev(idx_offset) 
    
else
    E_stdc = E_PCR*E_PCR_scale
end


slope_stdc = -1/log10(E_stdc + 1);
conc_stdc_eff = log10(conc_fit) - (c_eff_best - c_span_eff)/slope_stdc;
conc_stdc_eff = 10.^conc_stdc_eff;

conc_stdc = conc_stdc_eff;
model_scaled = model_eff_scaled_span;

%save data
stats_savename = ['eff_FAM_stdcurve' filename(10:end)];
save(stats_savename,'model_scaled','conc_stdc') 

%% Figures

blue = linspecer('blue');
red = linspecer('red');
gray = linspecer('gray');
color_span = linspecer(N_span);

%best fit model
figure(1); clf(1)



hold on

yyaxis right
plot(cycle40,model_best,'--','color',red(64,:),'linewidth',1.5)
errorbar(cycle_curve,D_avg_curve,D_std_curve,'.','color',gray(128,:),'linewidth',1,...
    'markersize',10)
ylabel('\DeltaR_N')
ylim([-0.05 inf])


yyaxis left
plot(cycle,E_best_rev,'--','linewidth',1.5,'color',blue(64,:))
ylabel('PCR Efficiency')
ylim([0 1.05])

hold off

box on
%title('Best Fit Sigmoidal Curve')
xlabel('Cycle')
%ylabel('\DeltaRn (a.u.)')
%ylabel('Normalized Fluoresence (a.u.)')
legend('PCR efficiency',['Amp curve fit (R^2 = ',num2str(R_sq_eff_best,'%1.4f'),')'],'fontsize',10,'location','se')
%set(gca,'fontsize',14,'linewidth',2,'yscale','log')
set(gca,'fontsize',14,'linewidth',1,'yscale','lin','layer','top','xminortick','on','yminortick','on')
xlim([0 41])
%axis([0 41 -0.2 inf])
%axis([0 41 1e-1 5e1])





%model span
figure(2); clf(2)

hold on

for i = 1 : N_span
    
    plot(cycle40,model_scaled(1:N_cycles,i),'-','linewidth',3,'color',color_span(i,:))
    
end

hold off

box on
title(['Curve Extended (N = ',num2str(N_span,'%3.0f'),')'])
xlabel('Cycle')
ylabel('\DeltaRn (a.u.)')
set(gca,'fontsize',14,'linewidth',1,'yscale','lin')
axis([0 41 -0.2 inf])





%cycle i standard curve (post split)
figure(3); clf(3)

semilogy(model_scaled(cycle_graph,:),conc_stdc,'-k','linewidth',1) %only need C_1/2 values (absolute)

box on

%xlabel('Fluorescence Intensity (a.u.)')
xlabel('\DeltaRn (a.u.)')
%xlabel('Normalized Fluorescence (a.u.)')
%ylabel('IAV RNA Copies per Drop (cpd)')
ylabel('C_{RNA} (cpd)')
title(['Standard Curve (Cycle ',num2str(cycle_graph,'%2.0f'),')'])
set(gca,'fontsize',14,'linewidth',1)
axis([0 2.5 1e-2 1e6])






%subplot of standard curves
figure(4); clf(4)
L = length(cycle_curve);

for i = 2 : L
    
    subplot(3,3,i)
    semilogy(model_scaled(cycle_curve(i),:),conc_stdc,'linewidth',1)
    xlabel('Fluorescence Intensity (a.u.)')
    ylabel('IAV RNA Copies per Drop (cpd)')
    title(['Standard Curve (Cycle ',num2str(cycle_curve(i),'%2.0f'),')'])
    set(gca,'fontsize',10,'linewidth',1)
    axis([0 2 1e-3 1e7])
    
end





























