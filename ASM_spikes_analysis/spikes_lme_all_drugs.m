%% sample the medication curve at spike times to run time bin and correlation analyses

close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
%addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])


%load spike rate - new from 2/13/23 (samp/10min)
%load('spikes_rates_021323.mat');
%spikes_fname = 'spikes_rates_021323.mat';
spikes_fname = 'spikes_rates_SOZ_102723.mat';
load(spikes_fname);

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
all_t_last_sz = cell(1,length(ptIDs));
all_t_implant = cell(1,length(ptIDs));


for ipt =1:length(ptIDs)
    % get seizure times in ieeg and convert to idx of spike rate
    ptID = ['HUP' num2str(ptIDs(ipt))];
    spikes = all_spike_rate{ipt}';
    spikes = spikes./max(spikes); %normalize to compare across patients
    t_last_seizure = nan(1,length(spikes));

    [med_names,meds,explant_date,starts_eeg,starts_emu] = parse_MAR(ptID,all_meds);
    [taper_info,any_taper] = get_taper_info(meds, med_names);


    [sz_inds] = get_spike_seizure_inds(ptID,file_inds{ipt});
    % index of seizure occurance, use data before
    first_ind = sz_inds(1);
    ind_first_seizure(ipt) = first_ind; % seconds / 60sec/min /10min/point

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

% find meds that are given to at least 10 patients 
[GC,GR]=groupcounts(vertcat(all_med_names{:}));
med_feats = GR(GC>5);


feat_names = [{'ptID'}, {'t_last_sz'} {'t_implant'} {'awake'} med_feats' {'spike_rate'}];
features = cell(length(feat_names),length(ptIDs));
for i = 1:length(ptIDs)
    % add ptID
    features{1,i} = i*ones(1,length(all_pts_drug_samp{i}));
    % add asm loads 
    % add t_last_sz
    features{2,i} = all_t_last_sz{i};
    % add time since implant
    features{3,i} = all_t_implant{i};
    % add sleep labels
    features{4,i} = sleep_labels{i};
    % add spike rate
    features{end,i} = all_spike_rate{i}./max(all_spike_rate{i}); % normalize spike rate within each patient 
    asm_load = all_pts_drug_samp{i};
    pt_med_names = all_med_names{i};

        for m = 1:length(med_feats)
            ind = contains(pt_med_names,med_feats(m));
            feat_ind = contains(feat_names,med_feats(m));
            if sum(ind) > 0  
                features{feat_ind,i} = asm_load(ind,:);
            else 
                features{feat_ind,i} = zeros(1,width(asm_load));
            end 

        end 

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
med_char=[];
for i = 1:length(med_feats)
    med_char = [med_char, med_feats{i} '+'];
end 


modelspec =[ 'spike_rate~' med_char 't_last_sz+ awake + (1|ptID)'];
mdl =fitglme(drug_table,modelspec); % 'Distribution','Poisson'

%% plotting
% Assuming 'coefs' is a table or array with coefficient estimates and CIs
close all;
coefs =table;
coefs.Estimate = double(mdl.Coefficients(:,2));
coefs.Upper95CI = double(mdl.Coefficients(:,end));
coefs.Lower95CI = double(mdl.Coefficients(:,end-1));
coefs.Variable = mdl.CoefficientNames';  

OddsRatio = exp(coefs.Estimate);
LowerCI = exp(coefs.Lower95CI);
UpperCI = exp(coefs.Upper95CI);

% Plotting
figure;
errorbar(OddsRatio, 1:height(coefs), OddsRatio - LowerCI, UpperCI - OddsRatio, 'ko','horizontal', 'LineWidth', 1.5);
yticks(1:height(coefs));
yticklabels(coefs.Variable);
xlabel('Odds Ratio');
title('Odds Ratios with 95% Confidence Intervals');
xline(1,'--','linewidth',1.5,'color',[.5 .5 .5 .5])
set(gca,'fontsize',14,'Box','off')
axis square
grid on;
%xlim([.9 1.2])


