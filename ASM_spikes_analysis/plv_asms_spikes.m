%% calculate phase locking value over time between ASM load and Spikes, and see if it changes as ASM levels decrease

close all;clear;
cd('/Volumes/USERS/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/ASM-spikes-analysis'])
addpath('/Volumes/users/nghosn3/tools/wavelet-coherence')

%load spike rate and med data - new from 2/13/23 (samp/10min)
spikes_fname = 'spikes_rates_021323.mat';
load(spikes_fname);
% get the seizure and SOZ localization information
soz_info = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve','Sheet','SOZ');

cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

meds_fname = 'MAR_032122.mat';
load(meds_fname);

[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes(meds_fname,spikes_fname,all_ieeg_offset, all_dose_curves, all_tHr,ptIDs);

%%

all_plv = cell(1,length(ptIDs));
all_asm_plv = cell(1,length(ptIDs));
for ipt = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];

    % get spikes and ASM signals
    spikes = all_spike_rate{ipt};
    asm_load = sum(all_pts_drug_samp{ipt});
   

    %only use data up to first sz -  check if enough data before first seizure, otherwise exclude.
    sz_inds = get_spike_seizure_inds(ptID,file_inds{ipt});
    time_thresh = 6*24*2; % two days of data
    if sz_inds(1) > time_thresh
        asm_load = asm_load(1:sz_inds(1));
        spikes = spikes(1:sz_inds(1));
        signal1 = spikes;
        signal2 = asm_load;

        % Define analysis parameters - can adjust
        segment_length = 12;  % segment can be 1 hour - 6, 10minute samples
        overlap = 6;         % Overlap between segments can be 50% (adjust as needed)

        % Calculate the number of segments
        num_segments = floor((length(signal1) - overlap) / (segment_length - overlap));

        % Initialize an array to store PLV values over time
        plv_over_time = zeros(1, num_segments);
        asm_over_time = zeros(1, num_segments);
        spikes_over_time = zeros(1, num_segments);

        % Iterate through the segments
        for i = 1:num_segments
            % Define the start and end indices for the current segment
            start_idx = (i - 1) * (segment_length - overlap) + 1;
            end_idx = start_idx + segment_length - 1;

            % Extract the data for the current segment
            segment_signal1 = signal1(start_idx:end_idx);
            segment_signal2 = signal2(start_idx:end_idx);

            % Calculate PLV for the segment (as described in the previous response)
            signal1_analytic = hilbert(segment_signal1);
            signal2_analytic = hilbert(segment_signal2);
            phase_difference = angle(signal1_analytic ./ signal2_analytic);
            plv_over_time(i) = abs(mean(exp(1i * phase_difference)));
            asm_over_time(i) = mean(signal2(start_idx:end_idx));
            spikes_over_time(i) = mean(signal1(start_idx:end_idx));
        end

        % Create a time vector for the PLV values
        time_vector = (segment_length - overlap) * (0:num_segments - 1);

        % Plot the PLV values over time
        plot(time_vector, plv_over_time);
        xlabel('Time (ms)');
        ylabel('PLV');
        title('Phase Locking Value Over Time');

        all_plv{ipt} = plv_over_time;
        all_asm_plv{ipt} = asm_over_time;
        all_spikes_plv{ipt} = spikes_over_time;
    end


end

%%
plvs = [all_plv{:}];
asms = [all_asm_plv{:}];

figure;
tiledlayout('flow');
for ipt = 1:length(all_plv)
    if ~isempty(all_plv{ipt})
        nexttile;
        plot(all_plv{ipt},all_asm_plv{ipt},'.','markersize',10);
        title(['HUP' num2str(ptIDs(ipt))])
        [r,p] = corr(all_plv{ipt}',all_asm_plv{ipt}');
        legend(['R = ' num2str(r) ', p = ' num2str(p)],'Location', 'southoutside')
        
    end
end 

%% plot spikes and ASM load per patient

figure;
tiledlayout('flow');
for ipt = 1:length(all_plv)
    if ~isempty(all_plv{ipt})
        nexttile;
        plot(all_asm_plv{ipt},all_spikes_plv{ipt},'.','markersize',10);
        title(['HUP' num2str(ptIDs(ipt))])
        [r,p] = corr(all_asm_plv{ipt}',all_spikes_plv{ipt}');
        legend(['R = ' num2str(r) ', p = ' num2str(p)],'Location', 'southoutside')
        
    end 
end 