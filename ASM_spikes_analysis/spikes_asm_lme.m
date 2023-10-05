%% sample the medication curve at spike times to run time bin and correlation analyses

close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])

tic


%load spike rate - new from 2/13/23 (samp/10min)
load('spikes_rates_021323.mat');
spikes_fname = 'spikes_rates_021323.mat';
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;
load('MAR_032122.mat')
meds_fname = 'MAR_032122.mat';

[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes(meds_fname,spikes_fname,all_ieeg_offset, all_dose_curves, all_tHr,ptIDs);


%% time since last seizure

all_t_last_sz = cell(1,length(ptIDs));
all_t_implant = cell(1,length(ptIDs));

for ipt = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt};
    % Get time since last seizure
    t_last_seizure = nan(1,length(spikes));

    [sz_inds] = get_spike_seizure_inds(ptID,file_inds{ipt});
    seizure_inds = sz_inds;
    start = seizure_inds(1);
    ind=1;
    for j=start:seizure_inds(end)
        if ind<length(seizure_inds)
            if j>=seizure_inds(ind) && j < seizure_inds(ind+1)
                t_last_seizure(j)=j-seizure_inds(ind);
                
            else
                ind=ind+1;
                t_last_seizure(j)=j-seizure_inds(ind);
            end
        end
    end
    t_last_seizure(j+1:end) = (j+1:length(t_last_seizure))-seizure_inds(end);
    
    % assume last seizure happened 'first_sz' hours before emu stay and add that in
    % beginning
    nan_inds = find(isnan(t_last_seizure));
    first_sz = -24*2; 
    first_inds = nan_inds - first_sz;
    t_last_seizure(nan_inds)=first_inds;
    all_t_last_sz(ipt) = {t_last_seizure};
    all_t_implant(ipt) = {1:length(spikes)};
    
    
    
end

%% create feature matrix - model w/ response of spike rate to see coefficients - lme?
feat_names = [{'ptID'}, {'asm_load'}, {'t_last_sz'} {'t_implant'} {'spike_rate'}];
features = cell(4,length(ptIDs));
for i = 1:length(ptIDs)
    % add ptID
    features{1,i} = i*ones(1,length(all_pts_drug_samp{i}));
    % add asm load 
    features{2,i} = sum(all_pts_drug_samp{i}); % add all drugs together - averaging with zero decreases level
    % add t_last_sz
    features{3,i} = all_t_last_sz{i};
    % add time since implant
    features{4,i} = all_t_implant{i};
    % add spike rate
    features{5,i} = all_spike_rate{i};%./max(all_spike_rate{i}); % normalize spike rate within each patient
end 

drug_table = table();
for i =1:length(feat_names)
    feat_name = feat_names{i};
    this_feat = features(i,:);
    this_feat = horzcat(this_feat{:});
    drug_table.(feat_name)=this_feat';
end
%drug_table.ptID = categorical(drug_table.ptID);

%% fit linear mixed effects regression model 
modelspec = 'spike_rate~asm_load+t_last_sz+(1|ptID)';
mdl =fitglme(drug_table,modelspec) % 'Distribution','Poisson'

%% compare spike rate in low and high ASM states before seizures 

% ASM load vs spikes before first seizure
spikes_before_sz = cell(1,length(ptIDs));
for ipt = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt};%./max(all_spike_rate{ipt});
    % Get time since last seizure
    t_last_seizure = nan(1,length(spikes));

    [sz_inds] = get_spike_seizure_inds(ptID,file_inds{ipt});
    spikes_before_sz(ipt) = {spikes(1:sz_inds(1))};
    
end 

%% compare time to first seizure and average spikes 
t_first_sz = cellfun(@length,spikes_before_sz); % in hours
asm_before_sz = cell(1,length(ptIDs));
asm_state = cell(1,length(ptIDs)); % 0 for low load, 1 for high load

% for each patient, catagorize ASM loads to low and high 
for i = 1:length(ptIDs)
    asm = nanmean(all_pts_drug_samp{i});
    asm_before_sz{i} = asm(1:t_first_sz(i));
    states = nan(1,t_first_sz(i));
    states(asm_before_sz{i} > prctile(asm_before_sz{i},75) )= 1;
    states(asm_before_sz{i} < prctile(asm_before_sz{i},25) )= 0;

    asm_state{i} = states;

end

% patients with early vs late seizure - seizure within the first two days 
thresh = 2 * 24 * 6; % two days * 24 hours * 60min/hour * 1 sample/10min;
early_sz = t_first_sz <= thresh;

subplot(1,2,1)
scatter([asm_before_sz{:}],[spikes_before_sz{:}]);
xlabel('asm load','fontsize',18);
ylabel('spike rate','fontsize',18);
title('Before first seizure: all pts','fontsize',18)
axis square

subplot(1,2,2)
plot([spikes_before_sz{:}],'o');
xlabel('time (10min)','fontsize',18);
ylabel('spike rate','fontsize',18);
title('Before first seizure: all pts','fontsize',18)
axis square

all_asm_states = [asm_state{:}];
all_spikes_preictal = [spikes_before_sz{:}];

figure;
max_len = max([length(all_spikes_preictal(all_asm_states==1)'),length(all_spikes_preictal(all_asm_states==0)')]);
data = nan(max_len,2);
data(1:length(all_spikes_preictal(all_asm_states==1)),1)= all_spikes_preictal(all_asm_states==1)';
data(1:length(all_spikes_preictal(all_asm_states==0)),2)= all_spikes_preictal(all_asm_states==0)';

boxplot(data,'Labels',{'high asm','low asm'});
title('before first seizure','fontsize',18);

%%
figure;
subplot(1,2,1)
early_asm_states = [asm_state{early_sz}];
early_spikes_preictal = [spikes_before_sz{early_sz}];

late_asm_states = [asm_state{~early_sz}];
late_spikes_preictal= [spikes_before_sz{~early_sz}];

max_len = max([length(early_spikes_preictal(early_asm_states==1)'),length(early_spikes_preictal(early_asm_states==0)')]);
data = nan(max_len,2);
data(1:length(early_spikes_preictal(early_asm_states==1)),1)= early_spikes_preictal(early_asm_states==1)';
data(1:length(early_spikes_preictal(early_asm_states==0)),2)= early_spikes_preictal(early_asm_states==0)';

boxplot(data,'Labels',{'high asm','low asm'});
title('early seizures: before first seizure','fontsize',18);
axis square

%%
subplot(1,2,2)
max_len = max([length(late_spikes_preictal(late_asm_states==1)'),length(late_spikes_preictal(late_asm_states==0)')]);
data = nan(max_len,2);
data(1:length(late_spikes_preictal(late_asm_states==1)),1)= late_spikes_preictal(late_asm_states==1)';
data(1:length(late_spikes_preictal(late_asm_states==0)),2)= late_spikes_preictal(late_asm_states==0)';
boxplot(data,'Labels',{'high asm','low asm'});
title('late seizures: before first seizure','fontsize',18);
axis square




