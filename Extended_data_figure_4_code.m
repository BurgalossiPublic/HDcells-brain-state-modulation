%% Analysis script to reproduce key findings of extended figure 4
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_4.mat');

% define cell types
hd_cells=strcmp(data.extData_figure_4.cell_type,'hd');
non_hd_cells=strcmp(data.extData_figure_4.cell_type,'non hd');


%% get pupil and whisking for hd and non hd

hd.pupil=data.extData_figure_4.pupil(hd_cells,:);
hd.whisking=data.extData_figure_4.whisker(hd_cells,:);

non_hd.pupil=data.extData_figure_4.pupil(non_hd_cells,:);
non_hd.whisking=data.extData_figure_4.whisker(non_hd_cells,:);

trigger_time=data.extData_figure_4.time(1,:);
%% plot 

xlim_val=[-3,4];

figure
subplot(1,2,1)
yyaxis right
plot(trigger_time,mean(hd.pupil),'g','LineWidth',1.5)
ylim([.35,.52])
ylabel 'pupil area (norm)'
yyaxis left
plot(trigger_time,mean(hd.whisking),'m','LineWidth',1.5)
xlim(xlim_val)
ylim([.1,.5])
ylabel 'whisking pad motion (norm)'
line([0,0],[0,.55],'color','r')
xlabel 'time after sound onset'

subplot(1,2,2)
yyaxis right
plot(trigger_time,mean(non_hd.pupil),'g','LineWidth',1.5)
ylim([.345,.45])
xlim(xlim_val)

ylabel 'pupil area (norm)'

yyaxis left
plot(trigger_time,mean(non_hd.whisking),'m','LineWidth',1.5)
xlim(xlim_val)
ylim([.11,.5])
ylabel 'whisking pad motion (norm)'
line([0,0],[0,.5],'color','r')
xlabel 'time after sound onset'
set(gcf,'renderer','Painters')
