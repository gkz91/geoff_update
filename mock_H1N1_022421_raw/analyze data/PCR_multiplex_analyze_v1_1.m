%graph multiplex PCR data for FAM and CY5 probes in drops
%2-14-21 v1.1 include FAM filtered CY5 data and non-FAM filtered CY5 data
%Geoff Zath

%2-5-21 burst size experiment

clear; clc

%% Inputs

cycle = 28; %cycle data to graph
%nbins = 50;

%heatmap
Nx = 50; %x bin
Ny = 50; %y bin

xmin = -0.1;
xmax = 0.5;

ymin = -0.1;
ymax = 0.5;


%% Load data

PCR_data = load('processed_delRn_detection_data_mock_022421.mat');

PCR_cycles = PCR_data.cycle;
%FAM_filter = PCR_data.delRn_FAM_FINAL;
%CY5_filter = PCR_data.delRn_CY5_FAM_FINAL;
FAM = PCR_data.delRn_FAM;
CY5 = PCR_data.delRn_CY5;

cycle_loc = find(PCR_cycles == cycle);

%FAM_filter_graph = FAM_filter{cycle_loc};
%CY5_filter_graph = CY5_filter{cycle_loc};
FAM_graph = FAM{cycle_loc};
CY5_graph = CY5{cycle_loc};


%% Process data





%% Figures

%colors
map = colormap(linspecer);
%map(1,:) = 1; %background white = 1 or black = 0


% %Filtered Cy5 vs FAM
% figure(1); clf(1)
% 
% plot(FAM_filter_graph,CY5_filter_graph,'.k','markersize',5)
% 
% xlabel('FAM \DeltaRn (a.u.)')
% ylabel('CY5 \DeltaRn (a.u.)')
% 
% title(['Filtered data cycle ',num2str(cycle,'%2.0f')])
% 
% box on
% 
% axis([xmin xmax ymin ymax])
% 
% set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])
% 
% 
% 
% 
% 
% 
% %heatmap of filtered Cy5 vs FAM
% figure(2); clf(2)
% 
% %Make 2D histogram
% x = linspace(xmin, xmax, Nx);
% y = linspace(ymin, ymax, Ny);
% [N_2D Xedge Yedge] = histcounts2(FAM_filter_graph,CY5_filter_graph,x,y);
% imagesc(x,y,N_2D')
% set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');
% 
% box on
% colormap(map)
% c = colorbar('eastoutside','fontsize',14);
% c.Label.String = 'Number of Drops';
% 
% title(['Filtered data, Cycle ',num2str(cycle,'%2.0f')])
% xlabel('FAM \DeltaRn (a.u.)')
% ylabel('CY5 \DeltaRn (a.u.)')
% axis([xmin xmax ymin ymax])
% set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])






%All Cy5 vs FAM
figure(1); clf(1)

plot(FAM_graph,CY5_graph,'.k','markersize',5)

xlabel('FAM \DeltaRn (a.u.)')
ylabel('CY5 \DeltaRn (a.u.)')

title(['All data, Cycle ',num2str(cycle,'%2.0f')])

box on

axis([xmin xmax ymin ymax])

set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])






%heatmap of all Cy5 vs FAM
figure(2); clf(2)

%Make 2D histogram
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);
[N_2D Xedge Yedge] = histcounts2(FAM_graph,CY5_graph,x,y);
imagesc(x,y,N_2D')
set(gca, 'XLim', x([1 end]), 'YLim', y([1 end]), 'YDir', 'normal');

box on
colormap(map)
c = colorbar('eastoutside','fontsize',14);
c.Label.String = 'Number of Drops';

title(['All data, cycle ',num2str(cycle,'%2.0f')])
xlabel('FAM \DeltaRn (a.u.)')
ylabel('CY5 \DeltaRn (a.u.)')
axis([xmin xmax ymin ymax])
set(gca,'fontsize',14,'linewidth',1,'xminortick','on','yminortick','on','ticklength',[0.015 1])
