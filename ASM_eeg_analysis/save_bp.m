%% Get the results spike detection from 5 min block every 30min (10-22-21: every 10min)
% we want to get the spike rate relative to the AED administration and
% seizure onset times

% get and save spike rate. make into function?
close all;clear;
curr_path = pwd;
addpath(curr_path)
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED'])


% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;

% save file?
save_file = 1;

all_bp = cell(1,length(ptIDs));
all_bp_times = cell(1,length(ptIDs));
file_inds = cell(1,length(ptIDs));
bp_file_dur = zeros(3,length(ptIDs));
freqs = [0.5 4;...
    4, 8;
    8, 12;
    12, 30;
    30, 80];

for x = 1%:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(x))];
    %load the spike data
    addpath('/Volumes/USERS/erinconr/projects/fc_toolbox/results/analysis/intermediate/');
    fname = [ptID '.mat'];
    try
        load(fname);
          %average across all channels!!
        delta =  (squeeze(summ.bp(:,1,:)));
        theta =  nanmean(squeeze(summ.bp(:,2,:)));
        alpha = nanmean(squeeze(summ.bp(:,3,:)));
        beta =  nanmean(squeeze(summ.bp(:,4,:)));
        gamma =  nanmean(squeeze(summ.bp(:,5,:)));

        
        bps = [delta; theta; alpha; beta; gamma];
        all_bp(x) = {bps};
        all_bp_times(x) = {summ.times};
        file_inds(x) = {summ.file_index};
        
        file_ind = find(diff(summ.file_times) <0);
        if ~isempty(file_ind)
            spike_file_dur(:,x) =[file_ind]; %in blocks
        end
    catch
        ptIDs(x)=NaN;
    end
end

%% save file in Pioneer/DATA
% if save_file
%     dirr = '/Volumes/USERS/nghosn3/Pioneer/DATA/spikes_rates_021323.mat';
%     save(dirr,'all_spike_rate','ptIDs','spike_file_dur','all_spike_times','file_inds');
%     disp 'saved downsampled spikes: 10min blocks'
% end
