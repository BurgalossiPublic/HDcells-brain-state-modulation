%% Analysis script to reproduce key findings of extended figure 3
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_3.mat');

% define cell types
in_PD=strcmp(data.extData_figure_3.cell_type,'inPD');
out_PD=strcmp(data.extData_figure_3.cell_type,'outPD');


%% get raster histograms

% define the binning of the raster window
half_window=1000;% in ms
n_bins=400;
time_bin=(half_window/n_bins)/1000;% in sec

% define raster binning
rasterEdges=linspace(-half_window,half_window,n_bins*2);


% alocate raster variables
rate_sound=zeros(length(data.extData_figure_3.cell_index),length(rasterEdges)-1);

for cc=1:length(data.extData_figure_3.cell_index)
    
    % get firing rate
    if ~isempty(data.extData_figure_3.spk_raster_time{cc})
        [count_on_sound,~]=histcounts(data.extData_figure_3.spk_raster_time{cc},rasterEdges);
        rate_sound(cc,:)=count_on_sound/(length(data.extData_figure_3.spk_raster_stim{cc})*time_bin);
        
    end
    
end

% bins to plot
raster_bin=medfilt1(rasterEdges,2);
raster_bin=raster_bin(2:end);

%% normalize rate
% zscore from pre stimulus period
time_pre=raster_bin<=0;

zs_sound=(rate_sound - mean(rate_sound(:,time_pre),2,'omitnan'))./std(rate_sound(:,time_pre),[],2,'omitnan');

%% get rasters for in and out side field stimulation

zs_sound_in=zs_sound(in_PD,:);
zs_sound_out=zs_sound(out_PD,:);

%% sort

% take firing post to sort
time_post=[];
time_post=raster_bin>=5 & raster_bin<=20;

% indexes to sort head direction
[~,sort_in]=sort(max(zs_sound_in(:,time_post),[],2),'descend');
[~,sort_out]=sort(max(zs_sound_out(:,time_post),[],2),'descend');


%% plotting head direction cells sound stimuli

left_lim=-20;
right_lim=100;

figure
subplot(2,2,1)
imagesc(raster_bin,[],zs_sound_in(sort_in,:))
caxis([-2,10])
title 'sound in PD hd'
ylabel 'cells'
xlim([left_lim,right_lim])

subplot(2,2,2)
imagesc(raster_bin,[],zs_sound_out(sort_out,:))
caxis([-2,35])
title 'sound off PD hd'
ylabel 'cells'
xlim([left_lim,right_lim])

subplot(2,2,3)
plot(raster_bin,mean(zs_sound_in,'omitnan'),'r','LineWidth',1.5)
hold on
plot(raster_bin,mean(zs_sound_in,'omitnan')+std(zs_sound_in,'omitnan'),':k')
plot(raster_bin,mean(zs_sound_in,'omitnan')-std(zs_sound_in,'omitnan'),':k')
ylim([-2,12])
xlim([left_lim,right_lim])
xlabel 'time after stim (ms)'
ylabel 'zscore (from pre stim)'

subplot(2,2,4)
plot(raster_bin,mean(zs_sound_out,'omitnan'),'r','LineWidth',1.5)
hold on
plot(raster_bin,mean(zs_sound_out,'omitnan')+std(zs_sound_out,'omitnan'),':k')
plot(raster_bin,mean(zs_sound_out,'omitnan')-std(zs_sound_out,'omitnan'),':k')
ylim([-2,35])
xlim([left_lim,right_lim])
xlabel 'time after stim (ms)'
ylabel 'zscore (from pre stim)'
set(gcf,'renderer','Painters')

%% compute modulation index  on/off PD


% get mean rate pre and post stimuli
mean_pre=mean(rate_sound(:,time_pre),2,'omitnan');
mean_post=mean(rate_sound(:,time_post),2,'omitnan');

% compute modulation index
mi_all=(mean_post-mean_pre)./(mean_post+mean_pre);

% get modulation index hd non hd
mi_in=mi_all(in_PD);
mi_out=mi_all(out_PD);


% plot figure
figure
scatter(mean_pre(in_PD),mean_post(in_PD),'g')
hold on
scatter(mean_pre(out_PD),mean_post(out_PD),'k')
ylim([0,max(mean_post(in_PD))])
xlim([0,max(mean_post(in_PD))])
line([0,max(mean_post(in_PD))],[0,max(mean_post(in_PD))],'Color','r')
xlabel 'firing rate pre'
ylabel 'firing post'
legend ({'in PD','out PD'})
%% modulation index comparison on PD vs off PD
% note that p values are different from those reported in the paper, as
% those were optained with a linear mix model.

[p,~]=ranksum(mi_in,mi_out);

figure
boxplot(mi_all,data.extData_figure_3.cell_type)
set(gca,'XTickLabel',{'In PD','Out PD'});
title(['pval=',num2str(p)])
ylabel 'modulation index'
ylim([-1,1])






