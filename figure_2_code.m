%% Analysis script to reproduce key findings of figure 2
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b
%% loading data

clear all
data=load('Figure_2.mat');

% define cell types
hd_cells=strcmp(data.psth_motion.cell_type,'hd');
non_hd_cells=strcmp(data.psth_motion.cell_type,'non hd');

%% get PSTH for pupil and whisking trigger motion

% define the binning of the raster window
% binning was set to 250ms for pupil and 50ms for whisking

halftime_mot=[15,4];%half time psth in sec
mov_binning=[halftime_mot(1)*4,halftime_mot(2)*20];

for mm=1:length(halftime_mot)
    
    % define parameters
    half_window=halftime_mot(mm)*1000;% in ms
    n_bins=mov_binning(mm);
    time_bin=(half_window/n_bins)/1000;% in sec
    
    % define raster binning
    rasterEdges=linspace(-half_window,half_window,n_bins*2);
    
    % define rasters to use
    if mm==1
        rasters_time=data.psth_motion.spk_raster_time_pupil;
        rasters_stim=data.psth_motion.spk_raster_trigg_pupil;
    else
        rasters_time=data.psth_motion.spk_raster_time_whisker;
        rasters_stim=data.psth_motion.spk_raster_trigg_whisker;
    end
    
    % alocate raster variable
    rate_mov=zeros(size(rasters_time{1},1),length(rasterEdges)-1);
    
    for cc=1:size(rasters_time,1)
        
        if any(rasters_time{cc})
            
            [count_mov,~]=histcounts(rasters_time{cc},rasterEdges);
        else
            continue
        end
        % get firing rate
        if sum(count_mov)~=0
            fr_mov=count_mov/(max(rasters_stim{cc})*time_bin);
            rate_mov(cc,:)=fr_mov;
        end
        
    end
    
    if mm==1
        rate_pupil=rate_mov;
        pupil_bin_plot=medfilt1(rasterEdges,2);
        pupil_bin_plot=pupil_bin_plot(2:end)/1000;
        
    else
        rate_whisk=rate_mov;
        whisk_bin_plot=medfilt1(rasterEdges,2);
        whisk_bin_plot=whisk_bin_plot(2:end)/1000;
    end
    
end

%% normalize the firing rate to pre stimulus period

% pupil
t_pre_pupil=[];
t_pre_pupil=pupil_bin_plot<=0;
zs_pupil=(rate_pupil - mean(rate_pupil(:,t_pre_pupil),2))./std(rate_pupil(:,t_pre_pupil),[],2);

% whisking
t_pre_whisk=[];
t_pre_whisk=whisk_bin_plot<=0;
zs_whisk=(rate_whisk - mean(rate_whisk(:,t_pre_whisk),2))./std(rate_whisk(:,t_pre_whisk),[],2);

%% compute modulation index

% define windows

% pupil
t_post_pupil=[];
t_post_pupil=pupil_bin_plot>=2 & pupil_bin_plot<=8;


%whisking
t_post_whisk=[];
t_post_whisk=whisk_bin_plot>=0 & whisk_bin_plot<=1;


% get pre and post means
mean_pre_pupil=mean(rate_pupil(:,t_pre_pupil),2);
mean_post_pupil=mean(rate_pupil(:,t_post_pupil),2);

mean_pre_whisk=mean(rate_whisk(:,t_pre_whisk),2);
mean_post_whisk=mean(rate_whisk(:,t_post_whisk),2);


% compute modulation index
mi_pupil=(mean_post_pupil-mean_pre_pupil)./(mean_post_pupil+mean_pre_pupil);
mi_whisk=(mean_post_whisk-mean_pre_whisk)./(mean_post_whisk+mean_pre_whisk);

% get modulation index hd vs non hd

mi_hd_pupil=mi_pupil(hd_cells);
mi_n_hd_pupil=mi_pupil(non_hd_cells);

mi_hd_whisk=mi_whisk(hd_cells);
mi_n_hd_whisk=mi_whisk(non_hd_cells);

%% get the normalize rate for hd cells and non hd cells

% hd cells
hd_zs_pupil=zs_pupil(hd_cells,:);
hd_zs_whisk=zs_whisk(hd_cells,:);

% non hd cells
n_hd_zs_pupil=zs_pupil(non_hd_cells,:);
n_hd_zs_whisk=zs_whisk(non_hd_cells,:);

%% sort by max

% sorting the pupil rate
t_post_pupil=[];
t_post_pupil=pupil_bin_plot>=2 & pupil_bin_plot<=8;

% take firing post to sort
hd_post_pupil=mean(hd_zs_pupil(:,t_post_pupil),2);
n_hd_post_pupil=mean(n_hd_zs_pupil(:,t_post_pupil),2);

% indexes to sort
[~,s_hd_pupil]=sort(hd_post_pupil,'descend');
[~,s_n_hd_pupil]=sort(n_hd_post_pupil,'descend');


% sorting the whisking rate
t_post_whisk=[];
t_post_whisk=whisk_bin_plot>=0 & whisk_bin_plot<=1;

% take firing post to sort
hd_post_whisk=mean(hd_zs_whisk(:,t_post_whisk),2);
n_hd_post_whisk=mean(n_hd_zs_whisk(:,t_post_whisk),2);

% indexes to sort
[~,s_hd_whisk]=sort(hd_post_whisk,'descend');
[~,s_n_hd_whisk]=sort(n_hd_post_whisk,'descend');

%% plot psth

% define limits for the plot
x_lim_pupil=[-4,12];
x_lim_whisk=[-1,2];

y_lim_pupil=[-2,6];
y_lim_whisk=[-2,11];

c_lim_pupil=[-2,5];
c_lim_whisk=[-2,8];

% head direction cells

figure
subplot(2,2,1)
imagesc(pupil_bin_plot,[],hd_zs_pupil(s_hd_pupil,:))
caxis(c_lim_pupil)
xlim(x_lim_pupil)
title('pupil trigg')
ylabel 'cells'

subplot(2,2,2)
imagesc(whisk_bin_plot,[],hd_zs_whisk(s_hd_whisk,:))
caxis(c_lim_whisk)
xlim(x_lim_whisk)
title('whisk trigg')
ylabel 'cells'

subplot(2,2,3)
plot(pupil_bin_plot,mean(hd_zs_pupil),'b','LineWidth',1.5)
hold on
plot(pupil_bin_plot,mean(hd_zs_pupil)+std(hd_zs_pupil,[],1),':b')
plot(pupil_bin_plot,mean(hd_zs_pupil)-std(hd_zs_pupil,[],1),':b')
ylim(y_lim_pupil)
xlim(x_lim_pupil)
line([0,0],y_lim_pupil,'color','r')
ylabel 'zscore'
xlabel 'time after onset (sec)'

subplot(2,2,4)
plot(whisk_bin_plot,mean(hd_zs_whisk),'b','LineWidth',1.5)
hold on
plot(whisk_bin_plot,mean(hd_zs_whisk)+std(hd_zs_whisk,[],1),':b')
plot(whisk_bin_plot,mean(hd_zs_whisk)-std(hd_zs_whisk,[],1),':b')
ylim(y_lim_whisk)
xlim(x_lim_whisk)
line([0,0],y_lim_whisk,'color','r')
ylabel 'zscore'
xlabel 'time after onset (sec)'
set(gcf,'renderer','Painters')

% non direction cells
figure
subplot(2,2,1)
imagesc(pupil_bin_plot,[],n_hd_zs_pupil(s_n_hd_pupil,:))
caxis(c_lim_pupil)
xlim(x_lim_pupil)
title('pupil trigg')
ylabel 'cells'

subplot(2,2,2)
imagesc(whisk_bin_plot,[],n_hd_zs_whisk(s_n_hd_whisk,:))
caxis(c_lim_whisk)
xlim(x_lim_whisk)
title('whisk trigg')
ylabel 'cells'

subplot(2,2,3)
plot(pupil_bin_plot,mean(n_hd_zs_pupil),'k','LineWidth',1.5)
hold on
plot(pupil_bin_plot,mean(n_hd_zs_pupil)+std(n_hd_zs_pupil,[],1),':k')
plot(pupil_bin_plot,mean(n_hd_zs_pupil)-std(n_hd_zs_pupil,[],1),':k')
ylim(y_lim_pupil)
xlim(x_lim_pupil)
line([0,0],y_lim_pupil,'color','r')
ylabel 'zscore'
xlabel 'time after onset (sec)'

subplot(2,2,4)
plot(whisk_bin_plot,mean(n_hd_zs_whisk),'k','LineWidth',1.5)
hold on
plot(whisk_bin_plot,mean(n_hd_zs_whisk)+std(n_hd_zs_whisk,[],1),':k')
plot(whisk_bin_plot,mean(n_hd_zs_whisk)-std(n_hd_zs_whisk,[],1),':k')
ylim(y_lim_whisk)
xlim(x_lim_whisk)
line([0,0],y_lim_whisk,'color','r')
ylabel 'zscore'
xlabel 'time after onset (sec)'
set(gcf,'renderer','Painters')


%% boxplot modulation index
% note that p values are different from those reported in the paper, as
% the latter were optained with a linear mix model.

figure
subplot(1,2,1)
mi_whisk_val=[mi_hd_whisk;mi_n_hd_whisk];
mi_whisk_grp=[ones(size(mi_hd_whisk)); ones(size(mi_n_hd_whisk))*2];
boxplot(mi_whisk_val,mi_whisk_grp)
p=ranksum(mi_hd_whisk,mi_n_hd_whisk);
xticklabels({'hd','non hd '})
ylabel 'modulation index'
ylim([-.4,.4])
set(gca,'TickDir','out');
title(['whisk ',num2str(p)])


subplot(1,2,2)
mi_pupil_val=[mi_hd_pupil;mi_n_hd_pupil];
mi_pupil_grp=[ones(size(mi_hd_pupil)); ones(size(mi_n_hd_pupil))*2];
boxplot(mi_pupil_val,mi_pupil_grp)
p=ranksum(mi_hd_pupil,mi_n_hd_pupil);
xticklabels({'hd','non hd '})
ylabel 'modulation index'
ylim([-.4,.8])
set(gca,'TickDir','out');
title(['pupil ',num2str(p)])

%% Pupil cycle analysis

% define cell types
hd_cells=strcmp(data.pupil_cycle.cell_type,'hd');
non_hd_cells=strcmp(data.pupil_cycle.cell_type,'non hd');


%% plot pupil cycle

% define hd and non hd matrix
cycle_FR_hd=data.pupil_cycle.FR_zscore(hd_cells,:);
cycle_FR_non_hd=data.pupil_cycle.FR_zscore(non_hd_cells,:);

% bins
phase_bins=data.pupil_cycle.cycle_bins(1,:);

figure

plot(phase_bins,mean(cycle_FR_hd,'omitnan'),'b','LineWidth',1.5)
hold on
plot(phase_bins,mean(cycle_FR_hd,'omitnan')+std(cycle_FR_hd,'omitnan'),':b')
plot(phase_bins,mean(cycle_FR_hd,'omitnan')-std(cycle_FR_hd,'omitnan'),':b')

plot(phase_bins,mean(cycle_FR_non_hd,'omitnan'),'g','LineWidth',1.5)
plot(phase_bins,mean(cycle_FR_non_hd,'omitnan')+std(cycle_FR_non_hd,'omitnan'),':g')
plot(phase_bins,mean(cycle_FR_non_hd,'omitnan')-std(cycle_FR_non_hd,'omitnan'),':g')
set(gca,'TickDir','out');

legend({'HD','','','nonHD','',''},'Location','best')

title 'firing rate'
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
left_comp=[mean(cycle_FR_hd(:,left_win),2); mean(cycle_FR_non_hd(:,left_win),2)];
right_comp=[mean(cycle_FR_hd(:,right_win),2); mean(cycle_FR_non_hd(:,right_win),2)];

% group data
var_plot=[left_comp;right_comp];
var_grp=[ones(sum(hd_cells),1);ones(sum(non_hd_cells),1)*2 ;ones(sum(hd_cells),1)*3 ;ones(sum(non_hd_cells),1)*4];

% boxplot
figure
boxplot(var_plot,var_grp)
xticklabels({'hd dila','non hd dila','hd const','non hd const'})
ylabel 'firing rate (zscore)'
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

%% Analysis pupil amplitude vs firing rate during prefer direction window

% define cell types
hd_cells=strcmp(data.per_pass_data.cell_type,'hd');
non_hd_cells=strcmp(data.per_pass_data.cell_type,'non hd');

% separate hd non hd pases mean FR
hd_table=data.per_pass_data(hd_cells,:);
non_hd_table=data.per_pass_data(non_hd_cells,:);

%% compute LMM 

% LMM 
lme_intercept_slope_hd = fitlme(hd_table,...
    'norm_FR ~ 1 + norm_Pupil + (1+norm_Pupil|subjects)','FitMethod','REML');

lme_intercept_slope_non_hd = fitlme(non_hd_table,...
    'norm_FR ~ 1 + norm_Pupil + (1+norm_Pupil|subjects)','FitMethod','REML');

% get fixed effects p-vals residuals 
[~,~,stats_hd] = fixedEffects(lme_intercept_slope_hd,'alpha',0.05);
[~,~,stats_non_hd] = fixedEffects(lme_intercept_slope_non_hd,'alpha',0.05);



%% ploting

% get coefficients from LMM models 
coefs_hd=lme_intercept_slope_hd.Coefficients.Estimate;
coefs_non_hd=lme_intercept_slope_non_hd.Coefficients.Estimate;


% ploting the LMM fit per subject
figure
subplot(1,2,1)
scatter(hd_table.norm_Pupil,hd_table.norm_FR,30,'.k');
hold on
lfit=refline(coefs_hd(2),coefs_hd(1));
lfit.Color='r';
lfit.LineWidth=2;
ylim([0,1])
xlabel('norm Pupil'), ylabel('norm FR')

title(['HD',' LMM pval=',num2str(stats_hd.pValue(2)),' R2=',num2str(lme_intercept_slope_hd.Rsquared.Adjusted)])


subplot(1,2,2)
scatter(non_hd_table.norm_Pupil,non_hd_table.norm_FR,30,'.k');
hold on
lfit=refline(coefs_non_hd(2),coefs_non_hd(1));
lfit.Color='r';
lfit.LineWidth=2;
ylim([0,1])
xlabel('norm Pupil'), ylabel('norm FR')

title(['nonHD ','LMM pval=',num2str(stats_non_hd.pValue(2)),' R2=',num2str(lme_intercept_slope_non_hd.Rsquared.Adjusted)])



