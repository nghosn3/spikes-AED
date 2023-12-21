%% Get the results spike detection from 5 min block every 30min (10-22-21: every 10min)
% we want to get the spike rate relative to the AED administration and
% seizure onset times

% get and save spike rate. make into function?
close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])



% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
soz_info = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve','Sheet','SOZ');

% save file?
save_file = 1;
use_all_elecs = 0;

%% 
all_spike_rate = cell(1,length(ptIDs));
all_spike_times = cell(1,length(ptIDs));
file_inds = cell(1,length(ptIDs));
spike_file_dur = zeros(3,length(ptIDs));

for x = 26%1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(x))];
    %load the spike data
    addpath('/Volumes/USERS/erinconr/projects/fc_toolbox/results/analysis/intermediate/');
    fname = [ptID '.mat'];

    soz_chans = soz_info.("SOZ electrode")(contains(soz_info.name,num2str(ptIDs(x))));
    soz_chans = strsplit(soz_chans{:},', ');
    %try
        load(fname);
        % later can remove spikes that occur during seizures
        soz_inds = zeros(length(soz_chans),1);
        bad_inds =[];
        for i =1:length(soz_chans)
            this_ind = find(strcmp(summ.labels,soz_chans{i}));
            if ~isempty(this_ind) 
                soz_inds(i) = this_ind;
            else 
                bad_inds = [bad_inds i];
            end 
        end
        soz_inds(bad_inds) =[];

        if ~isempty(soz_inds) && use_all_elecs ==0
            spike_rate = nansum(summ.spikes(soz_inds,:),1);
            
        else

            spike_rate = nansum(summ.spikes(:,:),1); %spikes across SOZ channels!!
            disp(['used all elecs for ' ptID])
        end
        %spike_rate =  nansum(summ.spikes); spikes across all channels
        all_spike_rate(x) = {spike_rate};
        all_spike_times(x) = {summ.times};
        file_inds(x) = {summ.file_index};
        
        file_ind = find(diff(summ.file_times) <0);
        if ~isempty(file_ind)
            spike_file_dur(1:length(file_ind),x) =file_ind; %in blocks
        end
%     catch
%         ptIDs(x)=NaN;
%     end
end

%% save file in Pioneer/DATA
if save_file
    dirr = '/Volumes/users/nghosn3/Pioneer/DATA/spikes_rates_SOZ_102723.mat';
    save(dirr,'all_spike_rate','ptIDs','spike_file_dur','all_spike_times','file_inds');
    disp 'saved downsampled spikes: 10min blocks'
end
