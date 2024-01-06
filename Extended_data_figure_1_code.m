%% Analysis script to reproduce key findings of extended figure 1
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_1.mat');

% define cell types
hd_cells=strcmp(data.tunning_curves.cell_type,'hd');
non_hd_cells=strcmp(data.tunning_curves.cell_type,'non hd');

%% plot head direction index hd vs non hd cells

index_hd=data.tunning_curves.hd_index(hd_cells);
index_non_hd=data.tunning_curves.hd_index(non_hd_cells);

hd_bin=linspace(0,1,30);

hd_bin_plot=medfilt1(hd_bin,2);
hd_bin_plot=hd_bin_plot(2:end);

[counts_hd,~]=histcounts(index_hd,hd_bin);
[counts_n_hd,~]=histcounts(index_non_hd,hd_bin);

figure
plot(hd_bin_plot,counts_hd,'b','LineWidth',1.5)
hold on
plot(hd_bin_plot,counts_n_hd,'k','LineWidth',1.5)
line([.4,.4],[0,90],'Color','r','LineWidth',2)
xlabel 'hd index'
ylabel 'cell count'
legend({'hd','non hd'})


%% tunning curves comparison

% hd cells
hd_rate_align=data.tunning_curves.tunning_curve(hd_cells,:);
% non hd cells
n_hd_rate_align=data.tunning_curves.tunning_curve(non_hd_cells,:);
% bins
tunning_bins=data.tunning_curves.tunning_bins(1,:);


figure
plot(tunning_bins,mean(hd_rate_align),'b')
hold on
plot(tunning_bins,mean(hd_rate_align)+std(hd_rate_align),':k')
plot(tunning_bins,mean(hd_rate_align)-std(hd_rate_align),':k')
xlim([-180 180])
ylim([-5,100])

plot(tunning_bins,mean(n_hd_rate_align),'k')
plot(tunning_bins,mean(n_hd_rate_align)+std(n_hd_rate_align),':k')
plot(tunning_bins,mean(n_hd_rate_align)-std(n_hd_rate_align),':k')
xlim([-180 180])
ylim([-5,100])
xlabel 'HD (degrees)'
ylabel 'FR (Hz)'


