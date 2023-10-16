%% sample the medication curve at spike times to run time bin and correlation analyses

close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
%addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])



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

% get AD sleep labels
[sleep_labels,sleep_times,ptIDs] = ad_label_sleep(ptIDs);

%%
% spikes before and after first seizure

sz_info = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve');

ind_first_seizure = zeros(length(ptIDs),1);
avg_rate_before_sz = zeros(length(ptIDs),1);
avg_rate_after_sz = zeros(length(ptIDs),1);
pre_ictal_rate = zeros(length(ptIDs),1);
post_ictal_rate = zeros(length(ptIDs),1);
interictal_rate = zeros(length(ptIDs),1);

for ipt =1:length(ptIDs)
    % get seizure times in ieeg and convert to idx of spike rate
    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt}';
    spikes = spikes./max(spikes); %normalize to compare across patients

    [med_names,meds,explant_date,starts_eeg,starts_emu] = parse_MAR(ptID,all_meds);
    [taper_info,any_taper] = get_taper_info(meds, med_names);


    [sz_inds] = get_spike_seizure_inds(ptID,file_inds{ipt});
    % index of seizure occurance, use data before
    first_ind = sz_inds(1);
    ind_first_seizure(ipt) = first_ind; % seconds / 60sec/min /10min/point

    % split the spikes, and remove data around the seizures - 1hr before
    % and 1 hour after
    two_hr = 6*2; % 10 min * 12 = 2hours of time
    before_spikes = spikes(1:(first_ind-two_hr));

    rm_inds = [];
    preictal_inds = [];
    postictal_inds =[];
    seizure_inds = sz_inds;
    for n = 1:length(seizure_inds)
        if seizure_inds(n)>two_hr
            rm_inds = [rm_inds seizure_inds(n)-two_hr:seizure_inds(n)+two_hr];
            preictal_inds = [preictal_inds seizure_inds(n)-two_hr:seizure_inds(n)];
            postictal_inds = [postictal_inds seizure_inds(n):seizure_inds(n)+two_hr];
        end
    end

    after_spikes = spikes;
    after_spikes(sz_inds) = [];
    after_spikes = after_spikes(first_ind+two_hr:end);
    
    avg_rate_before_sz(ipt) = mean(before_spikes);
    avg_rate_after_sz(ipt) =  mean(after_spikes);

    pre_ictal_rate(ipt) = mean(spikes(preictal_inds));
    postictal_inds(postictal_inds >length(spikes))=[];
    post_ictal_rate(ipt) = mean(spikes(postictal_inds));
    %     interictal_rate(ipt) = ;


end

% may want to exclude patients that haqd very early seizures
ex_pts = ~(ind_first_seizure < 6*6); % first seizure occured in first 6/12 hours


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
feat_names = [{'ptID'}, {'asm_load'}, {'t_last_sz'} {'t_implant'} {'awake'} {'spike_rate'}];
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
    % add sleep labels
     features{5,i} = sleep_labels{i};
    % add spike rate
    features{6,i} = all_spike_rate{i};%./max(all_spike_rate{i}); % normalize spike rate within each patient
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
modelspec = 'spike_rate~asm_load+t_last_sz+ awake + (1|ptID)';
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


%% spikes beginning and end of emu stay
early_late_emu_rate = zeros(length(ptIDs),2);
for ipt = 1:length(all_spike_rate)
    spikes = all_spike_rate{ipt}; spikes = spikes./max(spikes);
    early_late_emu_rate(ipt,:) = [mean(spikes(1:round(length(spikes)./5))) mean(spikes(round(length(spikes)./5)+1:end))];
end

% is the spike rate higher in the second half of the emu stay?
[h,p] = ttest(early_late_emu_rate(:,1),early_late_emu_rate(:,2)) % H = 0; no, the spike rate is not greater at the end


%% PLOTTING
% spike rate in low vs high ASM load states - paired analysis
high_asm_spikes = zeros(length(ptIDs),1);
low_asm_spikes = zeros(length(ptIDs),1);
med_asm_spikes = zeros(length(ptIDs),1);

for ipt = 1:length(ptIDs)

    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt}';
    spikes = (spikes)./max(spikes); %normalize to compare across patients
    asm_load = sum(all_pts_drug_samp{ipt});

    high_asm_inds = asm_load >= prctile(asm_load,75);
    low_asm_inds = asm_load <= prctile(asm_load,25);
    %med_asm_inds = asm_load > prctile(asm_load,25) &  asm_load < prctile(asm_load,75);


    high_asm_spikes(ipt) = mean(spikes(high_asm_inds));
    low_asm_spikes(ipt) = mean(spikes(low_asm_inds));
    %med_asm_spikes(ipt) = mean(spikes(med_asm_inds));

end


tiledlayout(2,4)
[h,p,ci,stats] = ttest(low_asm_spikes, high_asm_spikes)

nexttile([1, 2]);
plot([1 2],[low_asm_spikes(low_asm_spikes>high_asm_spikes) high_asm_spikes(low_asm_spikes>high_asm_spikes)],'b.-','linewidth',1)
hold on;
plot([1 2],[low_asm_spikes(low_asm_spikes<high_asm_spikes) high_asm_spikes(low_asm_spikes<high_asm_spikes)],'r.-', 'linewidth',1)
legend('N.S')
xlim([0.5 2.5]);
ylabel('norm spike rate','fontsize',14)
xticks([1,2])
xticklabels({'high ASM load','low ASM load'})
set(gca, 'FontSize', 12);
title('Spikes by ASM load state','fontsize',14)

%axis square;


nexttile([1, 2]);
plot(drug_table.asm_load,drug_table.spike_rate,'ko','markersize',5);
set(gca, 'FontSize', 14);
ylabel('spike rate');
xlabel('ASM load');
%axis square;
set(gca, 'FontSize', 12);
title('Spikes vs ASM load','fontsize',14)


% spikes in sleep and wake within each patient
sleep_spikes_asm = zeros(length(ptIDs),2);
awake_spikes_asm = zeros(length(ptIDs),2);
med_asm_spikes = zeros(length(ptIDs),1);

for ipt = 1:length(ptIDs)

    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt}';
    spikes = (spikes)./max(spikes); %normalize to compare across patients
    asm_load = sum(all_pts_drug_samp{ipt});
    labels = sleep_labels{ipt}';

    sleep_spikes_asm(ipt,1) = median(spikes(labels==0));
    sleep_spikes_asm(ipt,2) = median(asm_load(labels==0));
    
    awake_spikes_asm(ipt,1) = median(spikes(labels==1));
    awake_spikes_asm(ipt,2) = median(asm_load(labels==1));

end

nexttile();
higher_inds = awake_spikes_asm(:,1)>sleep_spikes_asm(:,1);
plot(awake_spikes_asm(higher_inds,1), sleep_spikes_asm(higher_inds,1),'b.','markersize',15);
hold on;
plot(awake_spikes_asm(~higher_inds,1), sleep_spikes_asm(~higher_inds,1),'rv','markersize',5,'MarkerFaceColor','k');
plot([0 0.4],[0 0.4],'r--','linewidth',1.5)
xlabel('Spikes in wake','fontsize',14)
ylabel('Spikes in sleep','fontsize',14)
%axis square;
set(gca, 'FontSize', 12);
title('Spikes by wakefulness','fontsize',14)
[h,p] = ttest(awake_spikes_asm(:,1),sleep_spikes_asm(:,1))
legend(['p = ' num2str(p)])


nexttile();
higher_inds = awake_spikes_asm(:,2)>sleep_spikes_asm(:,2);
plot(awake_spikes_asm(higher_inds,2), sleep_spikes_asm(higher_inds,2),'b.','markersize',15);
hold on;
plot(awake_spikes_asm(~higher_inds,2), sleep_spikes_asm(~higher_inds,2),'rv','markersize',5,'MarkerFaceColor','k');
plot([1 4],[1 4],'r--','linewidth',1.5)
xlabel('ASM load in wake','fontsize',14)
ylabel('ASM load in sleep','fontsize',14)
%axis square;
set(gca, 'FontSize', 12);
title('ASM load by wakefulness','fontsize',14)
[h,p] = ttest(awake_spikes_asm(:,2),sleep_spikes_asm(:,2))
legend(['p = ' num2str(p)])



% plot spikes before and after first seizure
nexttile([1, 2]);
plot(avg_rate_before_sz(ex_pts),avg_rate_after_sz(ex_pts),'.b','markersize',20);
hold on;
plot(linspace(0,0.3,100),linspace(0,0.3,100),'--k','linewidth',3)
ylabel('norm spike rate after first sz','fontsize',14);
xlabel('norm spike rate before first sz','fontsize',14);
limits = [0 0.3];
ylim(limits); xlim(limits); axis square;
title('inter-ictal spike rate before vs after first seizure','fontsize',16)

% compare pre and post ictal spike rate across patients (aggregated across seizures?)
[h,p] = ttest(pre_ictal_rate,post_ictal_rate); % post ictal spike rate is higher??









%%
% figure;
% subplot(1,2,1)
% early_asm_states = [asm_state{early_sz}];
% early_spikes_preictal = [spikes_before_sz{early_sz}];
% 
% late_asm_states = [asm_state{~early_sz}];
% late_spikes_preictal= [spikes_before_sz{~early_sz}];
% 
% max_len = max([length(early_spikes_preictal(early_asm_states==1)'),length(early_spikes_preictal(early_asm_states==0)')]);
% data = nan(max_len,2);
% data(1:length(early_spikes_preictal(early_asm_states==1)),1)= early_spikes_preictal(early_asm_states==1)';
% data(1:length(early_spikes_preictal(early_asm_states==0)),2)= early_spikes_preictal(early_asm_states==0)';
% 
% boxplot(data,'Labels',{'high asm','low asm'});
% title('early seizures: before first seizure','fontsize',18);
% axis square
% 
% %%
% subplot(1,2,2)
% max_len = max([length(late_spikes_preictal(late_asm_states==1)'),length(late_spikes_preictal(late_asm_states==0)')]);
% data = nan(max_len,2);
% data(1:length(late_spikes_preictal(late_asm_states==1)),1)= late_spikes_preictal(late_asm_states==1)';
% data(1:length(late_spikes_preictal(late_asm_states==0)),2)= late_spikes_preictal(late_asm_states==0)';
% boxplot(data,'Labels',{'high asm','low asm'});
% title('late seizures: before first seizure','fontsize',18);
% axis square





