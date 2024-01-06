%% Analysis script to reproduce key findings of extended figure 6
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_6.mat');

% define cell types
hd_cells=strcmp(data.extData_figure_6.cell_type,'hd');
non_hd_cells=strcmp(data.extData_figure_6.cell_type,'non hd');


%% get firing rate for hd and non hd

hd_FR=data.extData_figure_6.psth(hd_cells,:);
non_hd_FR=data.extData_figure_6.psth(non_hd_cells,:);

time_bins=data.extData_figure_6.psth_times(1,:);
%% plot

figure
subplot(1,2,1)
plot(time_bins,mean(hd_FR),'b','LineWidth',1.5)
hold on
plot(time_bins,mean(hd_FR)+std(hd_FR,[],1)/sqrt(sum(hd_cells)),':b')
plot(time_bins,mean(hd_FR)-std(hd_FR,[],1)/sqrt(sum(hd_cells)),':b')
ylim([0,1])
xlabel 'time (sec)'
ylabel 'norm rate'

subplot(1,2,2)
plot(time_bins,mean(non_hd_FR),'k','LineWidth',1.5)
hold on
plot(time_bins,mean(non_hd_FR)+std(non_hd_FR,[],1)/sqrt(sum(non_hd_cells)),':k')
plot(time_bins,mean(non_hd_FR)-std(non_hd_FR,[],1)/sqrt(sum(non_hd_cells)),':k')
ylim([0,1])
xlabel 'time (sec)'
ylabel 'norm rate'
