%% look at spike sequences
close all;clear;

cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
%addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])

seq_dir = '/Volumes/users/aguilac/Interictal_Spike_Analysis/HUMAN/working_feat_extract_code/working features/clean_spikeleads/clean_atlas_leaders.csv';
spike_leaders_table = readtable(seq_dir);

cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;
meds_fname = 'MAR_032122.mat';
spikes_fname = 'spikes_rates_SOZ_102723.mat';

load(meds_fname)
[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes(meds_fname,spikes_fname,all_ieeg_offset, all_dose_curves, all_tHr,ptIDs);


%% curate table with ptIDs and sort
close all;

spikes_seqs = table();
seq_ptids = [];
tiledlayout('flow');
LCFs = nan(length(ptIDs),2);
for i = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    pt_inds = contains(spike_leaders_table.pt_id,ptID);
    if sum(pt_inds)>0

        pt_info = spike_leaders_table(pt_inds,:);
        asm_curves = all_pts_drug_samp{i};
        asm_tvec = all_pts_tvec{i}{1}*3600; % convert to seconds
        [med_names,~,~,starts_eeg,starts_emu] = parse_MAR(ptID,all_meds);


        % find day of lowest dose and day of highest dose - first time align spike times to match ASM CURVES
        offsets = all_ieeg_offset{2,i};
     

        if length(offsets)==1
            
            [~,sorted_inds] = sort(pt_info.spike_sequence);
            pt_info = pt_info(sorted_inds,:);
            
            peak_times = pt_info.peak_time_usec./(1e6); % convert to seconds to match offsets
            peak_sizes = pt_info.seq_total_dur;

            spikes_seqs = [spikes_seqs; pt_info];
            seq_ptids = [seq_ptids ptIDs(i)];

            % align ieeg times for each file with emu medication times
            x=1; % change to reflect all files used
            peak_times = peak_times - starts_eeg(x) + starts_emu(x); % only useful if more than one file
            peak_times = (peak_times + offsets(1)); %shift for t=0 to be start of emu stay, not start of ieeg recording. convert to hours


            % find the mean load every 24hours
            averageEveryN = @(v, n) arrayfun(@(i) mean(v(i:min(i + n - 1, end))), 1:n:length(v)-n+1);
            deltaT = 24*6; % in samples
            avg_load = averageEveryN(sum(asm_curves), deltaT);
            [~,max_day] = max(avg_load);
            [~,min_day] = min(avg_load);

            high_load_bounds = [1+(max_day-1)*deltaT , (deltaT+ ((max_day-1) *deltaT))]; %* (10*60); % convert to seconds, 10 = samp rate
            high_load_bounds = asm_tvec(high_load_bounds);
            low_load_bounds = [1+(min_day-1)*deltaT , (deltaT+ ((min_day-1) *deltaT))]; 
            low_load_bounds = asm_tvec(low_load_bounds);

            % generate distributions for spike sequences in high and low ASM load days

            % make probability plot for high ASM load day
            nexttile();
            high_asm_inds = peak_times>=high_load_bounds(1) & peak_times<=high_load_bounds(2);
            sizes = peak_sizes(high_asm_inds);
%             h = histogram((sizes),'normalization','probability');
%             prob_vals = h.Values;
%             cas_size = h.BinEdges;
            [prob_vals,cas_size] = histcounts(sizes,'normalization','probability');
            cas_size = cas_size(1:end-1) + diff(cas_size) / 2;

            loglog(cas_size,prob_vals,'linewidth',2,'Color','blue')
            title(ptID);
            xlabel('sequence size (samples)');
            ylabel('P(size)');
            set(gca,'fontsize',14)

            hold on;
            low_asm_inds = peak_times>=low_load_bounds(1) & peak_times<=low_load_bounds(2);
            sizes = peak_sizes(low_asm_inds);
%             h = histogram(sizes,'normalization','probability');
%             prob_vals = h.Values;
%             cas_size = h.BinEdges;
            [prob_vals,cas_size] = histcounts(sizes,'normalization','probability');
            cas_size = cas_size(1:end-1) + diff(cas_size) / 2;
            loglog(cas_size,prob_vals,'linewidth',2,'Color','red')

            %legend('highest ASM load day','lowest ASM load day')

            %  calculate LCF for high and low loads
            sys_size = max(pt_info.channel_index);
            LCFs(i,1)=[]; % high ASM load
            LCFs(i,2)=[]; % low ASM load






        end
    end
end


% get average ASM level during each spike sequence (first peak index to last peak index)




%%

% files_to_run = readtable('filenames_w_ids.csv');
% 
% for i = 1:length(ptIDs)
%     ptID = ['HUP' num2str(ptIDs(i))];
%     pt_inds = contains(files_to_run.filename,ptID);
%     files_to_run.to_use_nina(pt_inds) = 1;
% end
% 
% 
% writetable(files_to_run,'filenames_w_ids_nina.csv')