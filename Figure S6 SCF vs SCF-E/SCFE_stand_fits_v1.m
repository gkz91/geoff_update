%plot SCF-E fits to all M gene RNA standards
%2-28-21 v1.1 add SCF fits
%Geoff Zath

clear; clc

%% Inputs




%load data

%SCFPE
A{1} = load('SCFPE_fit_std_1e9.mat');
A{2} = load('SCFPE_fit_std_1e8.mat');
A{3} = load('SCFPE_fit_std_1e7.mat');
A{4} = load('SCFPE_fit_std_1e6.mat');
A{5} = load('SCFPE_fit_std_1e5.mat');
A{6} = load('SCFPE_fit_std_1e4.mat');

B = load('Mgene_stand_conc.mat');
stand_conc = B.conc_stand;

C = load('Mgene_stand_avg.mat');
stand_avg = C.stand_avg;

D = load('Mgene_stand_std.mat');
stand_std = D.stand_std;



%SCF
E{1} = load('SCF_fit_std_1e9.mat');
E{2} = load('SCF_fit_std_1e8.mat');
E{3} = load('SCF_fit_std_1e7.mat');
E{4} = load('SCF_fit_std_1e6.mat');
E{5} = load('SCF_fit_std_1e5.mat');
E{6} = load('SCF_fit_std_1e4.mat');



%% Process Data

cycle = A{1}.cycle;


L = length(A);

for i = 1 : L
    
    SCFPE_fit(:,i) = A{i}.F_best_scaled;
    R_sq_SCFPE(i) = A{i}.R_sq_eff_best;
    
    SCF_fit(:,i) = E{i}.model_best;
    R_sq_SCF(:,i) = E{i}.maximum;
  
end

R_sq_SCFPE_avg = mean(R_sq_SCFPE);
R_sq_SCF_avg = mean(R_sq_SCF);

%% Figures

%setup colors
blue = linspecer('blue');
red = linspecer('red');
green = linspecer('green');
gray = linspecer('gray');
%color_span = linspecer(N_span);

%colorblind color palette from https://jfly.uni-koeln.de/color/
colorblind(1,:) = [204,121,167]/255; %reddish purple
colorblind(2,:) = [230,159,0]/255; %orange
colorblind(3,:) = [86,180,233]/255; %sky blue
colorblind(4,:) = [0,158,115]/255; %blueish green
colorblind(5,:) = [213,94,0]/255; %vermillion
colorblind(6,:) = [0,114,178]/255; %blue
colorblind(7,:) = [0,0,0]; %black
colorblind(8,:) = [240,228,66]/255; %yellow




%SCFPE fits to M gene dilutions
figure(1); clf(1)

hold on

for i = 1 : L
    
    plot(cycle,SCFPE_fit(:,i),':','linewidth',1.5,'color',gray(64,:))
    
    h(i) = errorbar(cycle,stand_avg(i,:),stand_std(i,:),'.','color',colorblind(i,:),'linewidth',0.5,...
        'markersize',7,'displayname','asdf');
    
    

end

hold off


ylabel('\DeltaR_N')
xlabel('Cycle')
title('SCF-E')



legend(h,[num2str(stand_conc(1),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(2),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(3),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(4),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(5),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(6),'%1.2E'),' c/uL'],...
    'fontsize',6,'location','nw')

xlim([0 41])
ylim([0 6])



box on

set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on','layer','top','ticklength',[0.015 1])

set(gcf, 'Position',  [100, 100, 420,220])


%print -painters -depsc SCF-E_bulk_standards.eps


%SCF fits to M gene dilutions
figure(2); clf(2)

hold on

for i = 1 : L
    
    plot(cycle,SCF_fit(:,i),':','linewidth',1.5,'color',gray(64,:))
    
    h(i) = errorbar(cycle,stand_avg(i,:),stand_std(i,:),'.','color',colorblind(i,:),'linewidth',0.5,...
        'markersize',7,'displayname','asdf');
    
    

end

hold off


ylabel('\DeltaR_N')
xlabel('Cycle')
title('SCF')




legend(h,[num2str(stand_conc(1),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(2),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(3),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(4),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(5),'%1.2E'),' c/uL'],...
    [num2str(stand_conc(6),'%1.2E'),' c/uL'],...
    'fontsize',6,'location','nw')

xlim([0 41])
ylim([0 6])



box on

set(gca,'fontsize',8,'linewidth',0.5,'yscale','lin','xminortick','on','yminortick','on','layer','top','ticklength',[0.015 1])

set(gcf, 'Position',  [100, 100, 420,220])

%print -painters -depsc SCF_bulk_standards.eps



