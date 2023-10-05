close all;clear;
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')

% Get which patients
cohort_info = readtable('AED-connectivity-cohort-list.xlsx');
ptIDs = cohort_info.PtIDs;

%load spike rate
load('all_spike_rate_ds.mat');

%fully sampled spikes
%load('fs_spike_rate21-Oct-2021.mat')

[all_dose_curves,all_Hr,ptIDs,~] = get_aed_curve(ptIDs);

for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID =  ptIDs{ipt};
    [med_times,med_names,days,each_drug_dose,total_dose,alc_times,meds] = get_all_med_info(ptID);
    figure(ipt);
    %% plot the spikes
    spike_rate=all_spike_rate{ipt}; % calculated spike rate fpr patient in list
    %line up time axes for plotting:
    timelim = [0 max(med_times)];
    subplot(2,1,1)
    time = (1:length(spike_rate))*.5;
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
    %addpath([pwd '/spikes-AED/results/change_point/full_samp_spikes'])
    %load(['[''' ptID ''']' '_fbkpts.mat'])
    load(['[''' ptID ''']' '_bkpts.mat'])
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
    
    for x = 1:length(all_dose_curves{ipt})
        c_drug = all_dose_curves{ipt}{x};
        
        
        %% plot the AED curve
        %make tHr the length of the Spike data, the truncate the blood curves
        tHr = all_Hr{ipt}{x};
        subplot(2,1,1)
        
        %% plot the AED curves
        
        subplot(2,1,2)
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'--r','linewidth',1);hold on;
        end
        hold on;
        c_drug_norm = c_drug ./ max(c_drug);
        plot(tHr,c_drug_norm,'linewidth',1.4); hold on;
        xlim(timelim);%xlim([0 tHr(end)]); % this should be the range of the spike data, when that is incorperated
        title(['Normalized approx. plasma concentration of each AED: ' ptIDs{ipt}]);xlabel('hours since EMU admission');
        legend(med_names);
        
        
    end
    
    
end

% to save figs:
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/tools');