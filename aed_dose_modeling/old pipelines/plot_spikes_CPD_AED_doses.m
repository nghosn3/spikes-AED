%% visualize the change point detection on DOWNSAMPLED SPIKE RATE with the seizure times, AED dosage, and alcohol
close all;clear;
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/')
addpath( genpath( 'spikes-AED'))

%load spike rate 
load('all_spike_rate_ds.mat')

% normalize AED doses for plotting/stats?
aed_norm =1;

for x = 1:length(ptIDs)
    ptID = ptIDs{x};
    
    %get the medication data
    if aed_norm
        [med_times,med_names,days,each_drug_dose,total_dose,alc_times] = get_med_info_norm(ptID);
    else
        [med_times,med_names,days,each_drug_dose,total_dose,alc_times] = get_med_info(ptID);
    end
    
    % calculated spike rate fpr patient in list
    spike_rate=all_spike_rate{x};
    
    %line up time axes for plotting:
    timelim = [0 max(med_times)];
    %% plot the spike rate
    figure();
   
    subplot(2,1,1)
    time = (1:length(spike_rate))*.5;
    %plot(time,spike_rate,'linewidth',1.2); xlim(timelim);
    plot(time,spike_rate,'linewidth',1.2); xlim(timelim);
    ylabel('spikes/30min'); xlabel('time (hrs)')
    title([ptID ': Spike rate and medication dose over EMU stay']);
    
    % plot the AEDS
    for i =1:length(med_times), xline(med_times(i),'b'); end
        hold on;
        
    % plot the seizures
    seizure_times = get_seizure_times(ptID);
    for j =1:length(seizure_times)
        xline(seizure_times(j,1),'--r','linewidth',1);hold on;
    end
    
    %plot alcohol times
    for j =1:height(alc_times)
        xline(alc_times.start(1)./3600,'--g','linewidth',1);hold on; %alc times are in seconds, convert to hours
    end
    
    % shade the plot based on the detected breakpoints (from python)-- is it possibe to run a python script from matlab?
    load(['[''' ptIDs{x} ''']' '_bkpts.mat'])
    bkpts = [1 bkpts];
    blue1 = [0 0.4470 0.7410];
    blue2 = [0.3010 0.7450 0.9330 ];
    
    grey1  = [127 127 127]./255;
    grey2 = [100 100 100]./255;
    for j=1:length(bkpts)-1
        hold on;
        if mod(j,2)==0 
            color = blue1;
        else
            color =blue2;
        end
        area(time(bkpts(j):bkpts(j+1)),spike_rate(bkpts(j):bkpts(j+1)),'basevalue',0,'FaceColor',color);
    end
     
    %% plot doses of medication over EMU stay
    subplot(2,1,2)
    % convert time axis to hours
    plot(days*24,each_drug_dose,'.-'); hold on;
    plot(days*24,total_dose,'k','LineWidth',1.5);xlim(timelim)
    legend([med_names; {'total AED dosage'}]);
    xlabel('time (hrs) from iEEG record start'); ylabel('dose (mg)');
    legend_names = [med_names; {'total AED dosage'}];
    %plot alcohol times
    for j =1:height(alc_times)
        xline(alc_times.start(1)./3600,'--g','linewidth',1);hold on; %alc times are in seconds, convert to hours
        legend_names = [legend_names; alc_times.medication(j)];
    end
   
    l = legend(legend_names,'Box','off','Location','West');
    % ,'Orientation','horizontal','Location','SouthOutside',
    
    %plot seizure times
    for j =1:length(seizure_times)
        xline(seizure_times(j,1),'--r','linewidth',1);hold on;
    end
    
    
   
end

% addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/tools')
% cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/results/AED_spikes_CPD_plots')
% print_all_figures_to_eps()
