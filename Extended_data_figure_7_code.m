%% Analysis script to reproduce key findings of extended figure 7
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_7.mat');

% define cell types
hd_cells_motion=strcmp(data.motion_sound_response.cell_type,'hd');
non_hd_cells_motion=strcmp(data.motion_sound_response.cell_type,'non hd');

%% get pupil and whisker motion for hd and non hd cells


hd.pupil=data.motion_sound_response.pupil(hd_cells_motion,:);
hd.whisking=data.motion_sound_response.whisking(hd_cells_motion,:);

non_hd.pupil=data.motion_sound_response.pupil(non_hd_cells_motion,:);
non_hd.whisking=data.motion_sound_response.whisking(non_hd_cells_motion,:);

pupil_time=data.motion_sound_response.times_pupil(1,:);
whisk_time=data.motion_sound_response.times_whisking(1,:);

%% plot event average responses


figure
subplot(2,2,1)
plot(pupil_time,mean(hd.pupil),'g')
title 'HD pupil'
ylabel 'norm amplitude'
xlabel 'time after onset(sec)'
ylim ([0,.55])
xlim ([-2,11])

subplot(2,2,2)
plot(whisk_time,smooth(mean(hd.whisking),3),'m')
title 'HD whisker'
xlabel 'time after onset(sec)'
ylim ([0,.45])
xlim ([-3,4])

subplot(2,2,3)
plot(pupil_time,mean(non_hd.pupil),'g')
title 'non HD pupil'
ylabel 'norm amplitude'
xlabel 'time after onset(sec)'
ylim ([0,.55])
xlim ([-2,11])

subplot(2,2,4)
plot(whisk_time,smooth(mean(non_hd.whisking),3),'m')
title 'non HD whisker'
xlabel 'time after onset(sec)'
ylim ([0,.45])
xlim ([-3,4])

%% define cell types for pupil cycle

hd_cells_cycle=strcmp(data.pupil_cycle_lfp_pow.cell_type,'hd');
non_hd_cells_cycle=strcmp(data.pupil_cycle_lfp_pow.cell_type,'non hd');


%% get hd and non hd lfp pow on pupil cycle

hd_lfp_pow=data.pupil_cycle_lfp_pow.lfp_pow(hd_cells_cycle,:);
non_hd_lfp_pow=data.pupil_cycle_lfp_pow.lfp_pow(non_hd_cells_cycle,:);

phase_bins=data.pupil_cycle_lfp_pow.pupil_phase(1,:);

%% plot 1-10 pow across pupil cycle


figure

plot(phase_bins,mean(hd_lfp_pow,'omitnan'),'b','LineWidth',1.5)
hold on
plot(phase_bins,mean(hd_lfp_pow,'omitnan')+std(hd_lfp_pow,'omitnan'),':b')
plot(phase_bins,mean(hd_lfp_pow,'omitnan')-std(hd_lfp_pow,'omitnan'),':b')

plot(phase_bins,mean(non_hd_lfp_pow,'omitnan'),'g','LineWidth',1.5)
plot(phase_bins,mean(non_hd_lfp_pow,'omitnan')+std(non_hd_lfp_pow,'omitnan'),':g')
plot(phase_bins,mean(non_hd_lfp_pow,'omitnan')-std(non_hd_lfp_pow,'omitnan'),':g')
set(gca,'TickDir','out');

legend({'HD','','','nonHD','',''},'Location','best')

title '1-10 Hz lfp pow'
axis tight
ylabel 'zscore'
xlabel 'pupil phase'
set(gcf,'renderer','Painters')


%% box plot pupil cycle phases hd versus non hd

% define windows in phase bins
left_win=[];
left_win=phase_bins>=-2.5 & phase_bins<=-1.5;

right_win=[];
right_win=phase_bins>=1.5 & phase_bins<=2.5;

% select side of the pupil cycle
left_comp=[mean(hd_lfp_pow(:,left_win),2); mean(non_hd_lfp_pow(:,left_win),2)];
right_comp=[mean(hd_lfp_pow(:,right_win),2); mean(non_hd_lfp_pow(:,right_win),2)];

% group data
var_plot=[left_comp;right_comp];
var_grp=[ones(sum(hd_cells_cycle),1);ones(sum(non_hd_cells_cycle),1)*2 ;ones(sum(hd_cells_cycle),1)*3 ;ones(sum(non_hd_cells_cycle),1)*4];

% boxplot
figure
boxplot(var_plot,var_grp)
xticklabels({'hd dila','non hd dila','hd const','non hd const'})
ylabel '1-10 Hz lfp pow (zscore)'
set(gca,'TickDir','out');
set(gcf,'renderer','Painters')

% get stats for multicomparison
[p_ks,tbl_ks,stats_ks]=kruskalwallis(var_plot,var_grp,'off');
[multi_comp,multi_m]=multcompare(stats_ks,'Display','off');

title (['p KS= ',num2str(p_ks)])

% plot the results from multicomparisons
figure
rectangle('Position',[0 0 10 10])
axis([0 10 0 10])
text(1.5,5,num2str(multi_comp),'FontSize',9)
title(['ks p= ',num2str(p_ks)])









