%% Analysis script to reproduce key findings of extended figure 5
% Sensory and behavioral modulation of the internal compass 
% Eduardo Blanco-Hernández, Giuseppe Balsamo, Patricia Preston-Ferrer and Andrea Burgalossi

% note: the script is meant to be run per sections in the order here provided.
% MATLAB version:2018b

%% loading data

clear all
data=load('Extended_Data_figure_5.mat');

% define cell types
hd_cells=strcmp(data.extData_figure_5.cell_type,'hd');
non_hd_cells=strcmp(data.extData_figure_5.cell_type,'non hd');

%% compute psth

% define parameters
half_window=1*1000;% in ms
n_bins=400;
time_bin=(half_window/n_bins)/1000;% in sec

% define raster binning
rasterEdges=linspace(-half_window,half_window,n_bins*2);

% alocate raster variable
rate_sound=zeros(size(data.extData_figure_5,1),length(rasterEdges)-1);

for cc=1:size(data.extData_figure_5,1)
    
    if ~isempty(data.extData_figure_5.raster_spike_time{cc})
        
        [count_sound,~]=histcounts(data.extData_figure_5.raster_spike_time{cc},rasterEdges);
        
        
        % get firing rate
        
        fr_sound=count_sound/(max(data.extData_figure_5.raster_spike_stim{cc})*time_bin);
        rate_sound(cc,:)=fr_sound;
        
    end
end

% get bins for plotting
bins_plot=medfilt1(rasterEdges,2);
bins_plot=bins_plot(2:end);

%% normalize rate
% zscore from pre stimulus period
time_pre=bins_plot<=0;

zs_sound=(rate_sound - mean(rate_sound(:,time_pre),2))./std(rate_sound(:,time_pre),[],2);

%% sorted by max response

bin_window=bins_plot>=5 & bins_plot<=30;
tmp=[];
for ss=1:size(zs_sound,1)
    
   [~,idx]=max(zs_sound(ss,bin_window));
   tmp(ss)=idx;
end

[~,sort_idx]=sort(tmp);

%% plotting 

figure
subplot(2,1,1)
imagesc(bins_plot,[],zs_sound(sort_idx,:))
ylabel 'cells'
xlim ([-20,500])
caxis ([0,30])

subplot(2,1,2)
plot(bins_plot,mean(zs_sound,'omitnan'),'r')
hold on
plot(bins_plot,mean(zs_sound,'omitnan')+std(zs_sound,'omitnan'),':k')
plot(bins_plot,mean(zs_sound,'omitnan')-std(zs_sound,'omitnan'),':k')
line([0,0],[-20,60],'color','r')
xlim ([-20,500])
ylim ([-2,30])
xlabel 'time after onset (ms)'
ylabel 'zscore'




