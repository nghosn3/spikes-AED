
close all;clear;
cd('/Volumes/USERS/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])


%load spike rate - new from 2/13/23 (samp/10min)
load('spikes_rates_021323.mat');
inds = cellfun('isempty',all_spike_rate);

cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;
load('MAR_032122.mat')

[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes('MAR_032122.mat','spikes_rates_021323.mat',all_ieeg_offset, all_dose_curves, all_tHr,ptIDs)

%% calculate time since last seizure
all_t_last_sz = cell(1,length(ptIDs));

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
    
    
end


%% time series analysis - VAR model for granger causality

all_gc_results = nan(1,length(ptIDs));
all_best_lag = nan(1,length(ptIDs));

for ipt = 1:length(ptIDs)
    % test for stationarity and make data stationary - spikes are generally
    % stationary but asm load is not
    asm_load = mean(all_pts_drug_samp{ipt},1)'; % averaged ASM load across medications
    spikes = all_spike_rate{ipt}';

    % make the asm load stationary by takeing the differende between points
    % (derivative)
    asm_load_diff = asm_load;%[diff(asm_load); 0];% add a point at the beginning since it will be one less

    % need to normalize the spike rates
    spikes_norm = spikes ./max(spikes);

    % get a time series of time since last seizure
    t_last_seizure = all_t_last_sz{ipt}./max(all_t_last_sz{ipt});

    % store pre-processed data into a table for model fitting and selection -
    % do this for each patient - start with example 1

    tbl = table(asm_load_diff,spikes_norm);
    T = size(tbl,1);

    % now need to find the optimal number of lags
    numseries = 2;
    numlags = (1:3:15)'; %(lags are in 10min blocks, so try 1:5 hours in increments of 1 hr)
    nummdls = numel(numlags);

    % Partition time base.
    maxp = max(numlags) * 5; % Maximum number of required presample responses - try making the pre estimate time base longer
    idxpre = 1:maxp;
    idxest = (maxp + 1):T;

    % Preallocation
    EstMdl(nummdls) = varm(numseries,0);
    aic = zeros(nummdls,1);

    % Fit VAR models to data.
    Y0 = tbl{idxpre,:}; % Presample
    Y = tbl{idxest,:};  % Estimation sample
    for j = 1:numel(numlags)
        Mdl = varm(numseries,numlags(j));
        Mdl.SeriesNames = string(tbl.Properties.VariableNames);
        EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
        results = summarize(EstMdl(j));
        aic(j) = results.AIC;
    end

    [~,bestidx] = min(aic);
    p = numlags(bestidx); %if p=6 then the best number of lags is 6 which is the same as one hour

    % now we take the best model and do a granger causality test

    BestMdl = EstMdl(bestidx);
    [h,Summary] = gctest(BestMdl,'Display',false); % gc test shows that the asm load granger causes the spike rate???

    all_gc_results(ipt) = h(1); % we want the result of testing the hyp for the affect of asm_load on spikes rate 
    all_best_lag(ipt) = p;

end


% can we conduct a test across all patients to see if this is generally
% true?

%% SARIMAX - seasonal arima model with exogenous variable

