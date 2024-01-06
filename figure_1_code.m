%% Analysis script to reproduce key findings of figure 1
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Figure_1.mat');

% define cell types
sound_hd=strcmp(data.psth_sound.cell_type,'hd');
sound_non_hd=strcmp(data.psth_sound.cell_type,'non hd');

%% compute psth for sound stimuli

% define parameters 
half_window=1000;% in ms
n_bins=400;
time_bin=(half_window/n_bins)/1000;% in sec

% define raster binning
rasterEdges=linspace(-half_window,half_window,n_bins*2);

% alocate raster variable
rate_sound=zeros(size(data.psth_sound,1),length(rasterEdges)-1);

for cc=1:size(data.psth_sound,1)
    
    [count_sound,~]=histcounts(data.psth_sound.spk_raster_time{cc},rasterEdges);
        
    % get firing rate
    if sum(count_sound)~=0
        fr_sound=count_sound/(max(data.psth_sound.spk_raster_stim{cc})*time_bin);
        rate_sound(cc,:)=fr_sound;
    end
    
end

% define raster bins
raster_bin_sound=medfilt1(rasterEdges,2);
raster_bin_sound=raster_bin_sound(2:end);

%% compute psth for whisker stimuli

% define parameters 
half_window=1000;% in ms
n_bins=400;
time_bin=(half_window/n_bins)/1000;% in sec

% define raster binning
rasterEdges=linspace(-half_window,half_window,n_bins*2);

% alocate raster variable
rate_whisker=zeros(size(data.psth_whisker,1),length(rasterEdges)-1);

for cc=1:size(data.psth_whisker,1)
    
    [count_sound,~]=histcounts(data.psth_whisker.spk_raster_time{cc},rasterEdges);
        
    % get firing rate
    if sum(count_sound)~=0
        fr_whisker=count_sound/(max(data.psth_whisker.spk_raster_stim{cc})*time_bin);
        rate_whisker(cc,:)=fr_whisker;
    end
    
end

% remove empty bins
zero_bins=sum(rate_whisker)==0;
rate_whisker(:,zero_bins)=[];

% define raster bins
raster_bin_whisker=medfilt1(rasterEdges,2);
raster_bin_whisker=raster_bin_whisker(2:end);
raster_bin_whisker(zero_bins)=[];

%% zscore normalization

% select window pre stimulus
pre_sound=raster_bin_sound<=0;
pre_whisker=raster_bin_whisker<=0;

% compute zscore
zs_sound=(rate_sound - mean(rate_sound(:,pre_sound),2))./std(rate_sound(:,pre_sound),[],2);
zs_whisker=(rate_whisker - mean(rate_whisker(:,pre_whisker),2))./std(rate_whisker(:,pre_whisker),[],2);

% get hd and non hd cells (only sound)
zs_sound_hd=zs_sound(sound_hd,:);
zs_sound_non_hd=zs_sound(sound_non_hd,:);

%% compute modulation index (for sound only)

% define window post stimulus
time_post=[];
time_post=raster_bin_sound>=5 & raster_bin_sound<=20;

% get mean rate pre and post stimuli
mean_pre=mean(rate_sound(:,pre_sound),2);
mean_post=mean(rate_sound(:,time_post),2);

% compute modulation index
mi=(mean_post-mean_pre)./(mean_post+mean_pre);

%% sorting psth

% sound

% take firing post to sort
mean_post=max(zs_sound_hd(:,time_post),[],2);
mean_post_non_hd=max(zs_sound_non_hd(:,time_post),[],2);

% indexes to sort
[~,sort_hd]=sort(mean_post,'descend');
[~,sort_non_hd]=sort(mean_post_non_hd,'descend');

% whisker

% take firing post to sort
time_post_whisker=raster_bin_whisker>=5 & raster_bin_whisker<=20;
mean_post_whisker=max(zs_whisker(:,time_post_whisker),[],2);

% indexes to sort
[~,sort_whisker]=sort(mean_post_whisker,'descend');

%% ploting psth for sound
x_axis_lim=[-20,100];
y_axis_lim=[-1,20];

c_lim_left=-2;
c_lim_right=20;

figure
subplot(2,2,1)
imagesc(raster_bin_sound,[],zs_sound_hd(sort_hd,:))
xlim(x_axis_lim)
ylabel 'cells'
caxis([c_lim_left,c_lim_right])

subplot(2,2,2)
imagesc(raster_bin_sound,[],zs_sound_non_hd(sort_non_hd,:))
xlim(x_axis_lim)
ylabel 'cells'
caxis([c_lim_left,c_lim_right])

subplot(2,2,3)
plot(raster_bin_sound,mean(zs_sound_hd),'r','LineWidth',1.5)
hold on
plot(raster_bin_sound,mean(zs_sound_hd)+std(zs_sound_hd),':k')
plot(raster_bin_sound,mean(zs_sound_hd)-std(zs_sound_hd),':k')
xlim(x_axis_lim)
ylim(y_axis_lim)
ylabel 'FR zscore'
xlabel 'time after onset (ms)'

subplot(2,2,4)
plot(raster_bin_sound,mean(zs_sound_non_hd),'r','LineWidth',1.5)
hold on
plot(raster_bin_sound,mean(zs_sound_non_hd)+std(zs_sound_non_hd),':k')
plot(raster_bin_sound,mean(zs_sound_non_hd)-std(zs_sound_non_hd),':k')
xlim(x_axis_lim)
ylim(y_axis_lim)
ylabel 'FR zscore'
xlabel 'time after onset (ms)'

%% ploting psth for whisker
% note that in figure 1 Firing rate at bins -+2ms around stimuli were omitted due the
% artifact produced by the piezo electric device.

% two stimuli are visble one at zero and another at 500ms, last one indicates
% the returning of the piezo to original position.

x_axis_lim=[-100,650];
y_axis_lim=[-1,15];

figure
subplot(2,1,1)
imagesc(raster_bin_whisker,[],zs_whisker(sort_whisker,:))
xlim(x_axis_lim)
ylabel 'cells'


subplot(2,1,2)
plot(raster_bin_whisker,mean(zs_whisker),'r','LineWidth',1.5)
hold on
plot(raster_bin_whisker,mean(zs_whisker)+std(zs_whisker),':k')
plot(raster_bin_whisker,mean(zs_whisker)-std(zs_whisker),':k')
xlim(x_axis_lim)
ylim(y_axis_lim)
ylabel 'FR zscore'
xlabel 'time after onset (ms)'
%% ploting modulation index (for sound only)

% get modulation index hd non hd
hd_mi=mi(sound_hd);
non_hd_mi=mi(sound_non_hd);

% define bins
mi_bins=linspace(-1,1,20);

% get histogram
[hd_mi_hist,~]=histcounts(hd_mi,mi_bins,'Normalization','Probability');
[non_hd_mi_hist,~]=histcounts(non_hd_mi,mi_bins,'Normalization','Probability');

% correct bins
mi_bins=medfilt1(mi_bins,2);
mi_bins=mi_bins(2:end);

figure
plot(mi_bins,hd_mi_hist,'LineWidth',1.5)
hold on
plot(mi_bins,non_hd_mi_hist,'LineWidth',1.5)
xlabel 'modulation index'
ylabel 'probability'
%% overlay sound and whisker response

figure

plot(raster_bin_sound,mean(zs_sound_hd),'b','LineWidth',1.5)
hold on
plot(raster_bin_whisker,mean(zs_whisker),'g','LineWidth',1.5)
xlim([-10,50])
ylabel 'FR zscore'
xlabel 'time after onset (ms)'
legend({'sound','whisker'})


