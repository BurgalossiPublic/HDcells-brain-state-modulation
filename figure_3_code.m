%% Analysis script to reproduce key findings of figure 3
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Figure_3.mat');

% define cell types
hd_cells=strcmp(data.social_data.cell_type,'hd');
non_hd_cells=strcmp(data.social_data.cell_type,'non hd');

%% compute psth

% define the binning of the raster window 
% the raw rasters window is defined by 'halfwindow' above (5sec original)

% define parameters 
half_window=5*1000;% in ms
n_bins=5*10;

time_bin=(half_window/n_bins)/1000;% in sec

% define raster binning
rasterEdges=linspace(-half_window,half_window,n_bins*2);

% alocate raster variable
rate_social=zeros(size(data.social_data,1),length(rasterEdges)-1);

for cc=1:size(data.social_data,1)
    
    if any(data.social_data.raster_spike_time{cc})
        
        [count,~]=histcounts(data.social_data.raster_spike_time{cc},rasterEdges);
        
        % get firing rate
        if sum(count)~=0
            fr_social=count/(max(data.social_data.raster_spike_interaction{cc})*time_bin);
            rate_social(cc,:)=fr_social;
        end
        
    end
end

% bins to plot
raster_bin=medfilt1(rasterEdges,2);
raster_bin=raster_bin(2:end);


%% get psth head direction and non direction

time_pre=[];
time_pre=raster_bin<=0;

time_post=[];
time_post=raster_bin>=0 & raster_bin<=2000;

% get psth
psth_hd=rate_social(hd_cells,:);
psth_non_hd=rate_social(non_hd_cells,:);

%% normalize the data to pre interaction

% normalize psth 
zs_hd=(psth_hd - mean(psth_hd(:,time_pre),2,'omitnan'))./std(psth_hd(:,time_pre),[],2,'omitnan');
zs_n_hd=(psth_non_hd - mean(psth_non_hd(:,time_pre),2,'omitnan'))./std(psth_non_hd(:,time_pre),[],2,'omitnan');

% replace inf if any
zs_hd(isinf(zs_hd))=nan;
zs_n_hd(isinf(zs_n_hd))=nan;


% sorting the post stim period
mean_post_hd=mean(zs_hd(:,time_post),2);
mean_post_n_hd=mean(zs_n_hd(:,time_post),2);

[~,s_hd]=sort(mean_post_hd,'descend');
[~,s_n_hd]=sort(mean_post_n_hd,'descend');

%% compute modulation index

% get mean rate pre and post stimuli
mean_pre=mean(rate_social(:,time_pre),2);
mean_post=mean(rate_social(:,time_post),2);

% compute modulation index
mi=(mean_post-mean_pre)./(mean_post+mean_pre);

% modulation index head direction non head direction
mi_hd=mi(hd_cells);
mi_n_hd=mi(non_hd_cells);


%% plot psth

bins2plot_sec=raster_bin/1000;

figure

subplot(2,2,1)
imagesc(bins2plot_sec,[],zs_hd(s_hd,:))
caxis([-2,4])
title 'HD'
ylabel 'cells'

subplot(2,2,2)
imagesc(bins2plot_sec,[],zs_n_hd(s_n_hd,:))
caxis([-2,4])
title 'non HD'
ylabel 'cells'

subplot(2,2,3)
plot (bins2plot_sec,mean(zs_hd,'omitnan'),'b','LineWidth',1.5)
hold on
plot (bins2plot_sec,mean(zs_hd,'omitnan')+std(zs_hd,'omitnan'),':b')
plot (bins2plot_sec,mean(zs_hd,'omitnan')-std(zs_hd,'omitnan'),':b')
ylim([-2,6])
line([0,0],[-2,6],'color','r')
xlabel 'time after interaction (sec)'
set(gca,'TickDir','out');

subplot(2,2,4)
plot (bins2plot_sec,mean(zs_n_hd,'omitnan'),'k','LineWidth',1.5)
hold on
plot (bins2plot_sec,mean(zs_n_hd,'omitnan')+std(zs_n_hd,'omitnan'),':k')
plot (bins2plot_sec,mean(zs_n_hd,'omitnan')-std(zs_n_hd,'omitnan'),':k')
ylim([-2,6])
line([0,0],[-2,6],'color','r')
xlabel 'time after interaction (sec)'
set(gca,'TickDir','out');

set(gcf,'renderer','Painters')

%% plot modulation index
% note that p values are different from those reported in the paper, as
% the latter were optained with a linear mix model.

% group data
mi_comp=[mi_hd;mi_n_hd];
mi_grp=[ones(size(mi_hd));ones(size(mi_n_hd))*2];

% test
[p,~]=ranksum(mi_hd,mi_n_hd);

figure
boxplot(mi_comp,mi_grp)
xticklabels ({'hd','no hd'})
title(['pval=', num2str(p)])
ylim([-.4,.7])
ylabel 'modulation index'
set(gcf,'renderer','Painters')
set(gca,'TickDir','out');

