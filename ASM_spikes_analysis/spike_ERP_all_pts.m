%% perform an ERP style analysis within and across patients at time of ASM dosing


close all;clear;
cd('/Volumes/USERS/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/ASM-spikes-analysis'])
addpath('/Volumes/users/nghosn3/tools/wavelet-coherence')

% % add fooof stuff to path
addpath(genpath('/Volumes/users/nghosn3/Pioneer/fooof_mat'))
% addpath('/Volumes/users/nghosn3/Pioneer/fooof_mat/fooof/')
% addpath('/Volumes/users/nghosn3/Pioneer/fooof_mat/fooof_mat/fooof/fooof/')

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

%%  find the dosing times and match across all drugs - only use patients that had q12 dosing

spike_trials = cell(1,length(ptIDs));
drug_trials = cell(1,length(ptIDs));
trial_labels = cell(1,length(ptIDs));

for ipt =1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    asm_curves = all_pts_drug_samp{ipt};
    drug_sum = sum(asm_curves);
    med_names = all_med_names{ipt};
    tvec = all_pts_tvec{ipt}{:};
    spikes = (all_spike_rate{ipt}+1);
    [~,meds,~,~,~] = parse_MAR(ptID,all_meds);

    %find when drugs are administered...
    [pks,locs] = findpeaks(-drug_sum);
    admin_times(1:length(locs)) = locs;
    max_admins = length(locs);

    % create day/night labels for the trials
    startDate =meds.date(1);
    duration = meds.date(end) - meds.date(1);
    numDates = length(drug_sum);  % Change this to the desired number of dates
    dateVector = linspace(startDate, startDate+duration, numDates);

    all_tod_labels = dateVector.Hour > 5 & dateVector.Hour <19; %label of 1 means dayime


    admin_times = admin_times(:,1:max_admins);
   
    
    [sz_inds] = get_spike_seizure_inds(ptID,file_inds{ipt});
    sz_inds = sort(sz_inds);

%     % find the trials that include seizures and remove
%     for s = 1:length(sz_inds)
%         admin_times(admin_times > (sz_inds(s)-6) & admin_times < (sz_inds(s)+6)) = [];
% 
%     end

    % get rid of all trials after the first seizure
    admin_times(admin_times >= sz_inds(1)) =[];


    % find the admin times to be used as trials

    trial_labels{ipt} = all_tod_labels(admin_times);


    all_drugs_admined = admin_times;
    pt_trials = zeros(length(all_drugs_admined),(12*6)+1); % 12hours by number of trials
    pt_drug_trials = zeros(length(all_drugs_admined),(12*6)+1); % 12hours by number of trials

    avg_asm_curves = sum(asm_curves);
    offset = 12*6; %xhours
    for t = 1:length(all_drugs_admined)
        % get 6hrs before and 6hrs after admin for each trial
        t_zero = all_drugs_admined(t); %
        t_start = t_zero - offset;
        t_end = t_zero +offset;

        if t_end < length(spikes) && t_start>0
            trial = spikes(t_start:t_end);
            pt_trials(t,1:length(trial)) = trial;
            pt_drug_trials(t,1:length(trial)) = avg_asm_curves(t_start:t_end);
        end

    end

    spike_trials(ipt) = {pt_trials};
    drug_trials(ipt) = {pt_drug_trials};


end

%% plot all trials across patients
all_trials = [];
all_trials_day = [];
all_trials_night = [];
all_asms_day = [];
all_asms_night =[];
all_ptids_night =[];
all_ptids_day =[];
plot_trials=1;

inc=0;
for i = 74%1:length(spike_trials)
    pt_trials = spike_trials{i};
    asm_trials = drug_trials{i};
    labels = trial_labels{i};
    zero_inds = (all(pt_trials == 0, 2));
    pt_trials(zero_inds,:)=[];
    labels(zero_inds)=[];

    if plot_trials
        if ~isempty(pt_trials)
            inc = inc+1;
            
            % plot the day trials
            %subplot(3,1,1)
            trials = pt_trials(labels,:);

            if ~isempty(trials)
                x = linspace(-6,6,width(trials));
                nexttile
                plot(x',trials);
                set(gca, 'box', 'off')
                hold on;
                plot(x',mean(trials),'k','linewidth',1)
                %title(['Day spike rate: example HUP' num2str(ptIDs(i))])
                title('Spike rate','fontsize',14)
                xline(x(x==0),'--r','linewidth',1)

                nexttile;
                plot(x',asm_trials); hold on;
                set(gca, 'box', 'off')
                xline(x(x==0),'--r','linewidth',1)
                xlabel('time (hrs)','fontsize',14)
                title('ASM administrations','fontsize',14)


            end

            % plot the evening trials
%             subplot(3,1,2)
%             trials = pt_trials(~labels,:);
%             if ~isempty(trials)
%                 x = linspace(-6,6,width(trials));
%                 plot(x',trials+(inc*.25));
%                 hold on;
%                 plot(x',mean(trials)+(inc*.25),'k','linewidth',1)
%                 title(['Night spike rate: HUP' num2str(ptIDs(i))])
%             end
% 
%             % plot the avg asm load over each other to make sure its indexing properly?
%             subplot(3,1,3)
%             plot(x',asm_trials); hold on;
%             plot(x',mean(asm_trials));
%             xline(0,'r--')
%             title(['ASM load: HUP' num2str(ptIDs(i))])
        end
    end
    max_val = max(max(pt_trials));
    z_scored_trials = zscore(pt_trials);
%     all_trials_day = [all_trials_day; mean(pt_trials(labels,:)./max_val,1)]; %the average normalized spike rate
%     all_trials_night = [all_trials_night; mean(pt_trials(~labels,:)./max_val,1)]; %the average normalized spike rate
    all_trials_day = [all_trials_day; z_scored_trials(labels,:)]; %the average normalized spike rate
    all_trials_night = [all_trials_night; z_scored_trials(~labels,:)]; %the average normalized spike rate
    
    all_asms_day = [all_asms_day; asm_trials(labels,:)];
    all_asms_night = [all_asms_night; asm_trials(~labels,:)];
    
    all_ptids_night = [all_ptids_night; ones(height(asm_trials(~labels,:)),1)*ptIDs(i)];
    all_ptids_day =   [all_ptids_day;   ones(height(asm_trials(labels,:)),1)*ptIDs(i)];
    
 
end

%%
% figure;
% plot(x',(all_trials),'b','linewidth',.1); hold on;
% plot(x',mean(all_trials),'k','linewidth',3); hold on;
% title('all trials across patients')
% ylabel('normalized spikes')
% xlabel('time hrs')


% Calculate the mean signal and standard deviation
meanSignal = nanmean(all_trials_day, 1);
stdSignal = nanstd(all_trials_day, 1);

% Calculate the confidence interval bounds
alpha = 0.05;  % Significance level (1 - confidence levelc)
n = size(all_trials_day, 1);
t = tinv(1 - alpha/2, n - 1);  % t-value for two-tailed t-distribution
ciLower = meanSignal - t * stdSignal / sqrt(n);
ciUpper = meanSignal + t * stdSignal / sqrt(n);

% Plot the mean day spike rate and confidence interval
figure;
subplot(2,1,1)
x = linspace(-6,6,length(meanSignal));
time = x;  % Time axis
plot(time, meanSignal, 'b', 'LineWidth', 2);
hold on;
fill([time, fliplr(time)], [ciLower, fliplr(ciUpper)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
[bkpts,res] = findchangepts(meanSignal,MaxNumChanges=1,statistic = 'mean');
xline(time(bkpts),'r','LineWidth',2)
xline(0,'k--')
legend('mean spike rate','95% CI','change in mean','','medication admin')
% Customize the plot
xlabel('Time relative to ASM admin');
ylabel('normalized spikes');
title('Average spike rate after ASM administration: day');
grid on;

%now plot the night trials
% Calculate the mean signal and standard deviation
meanSignal = nanmean(all_trials_night, 1);
stdSignal = nanstd(all_trials_night, 1);

% Calculate the confidence interval bounds
alpha = 0.05;  % Significance level (1 - confidence levelc)
n = size(all_trials_day, 1);
t = tinv(1 - alpha/2, n - 1);  % t-value for two-tailed t-distribution
ciLower = meanSignal - t * stdSignal / sqrt(n);
ciUpper = meanSignal + t * stdSignal / sqrt(n);

% Plot the mean day spike rate and confidence interval
subplot(2,1,2)
time = x;  % Time axis
plot(time, meanSignal, 'b', 'LineWidth', 2);
hold on;
fill([time, fliplr(time)], [ciLower, fliplr(ciUpper)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
[bkpts,res] = findchangepts(meanSignal,MaxNumChanges=1);
xline(time(bkpts),'r','LineWidth',2)
xline(0,'k--')
legend('mean spike rate','95% CI','change in mean','','medication admin')
% Customize the plot
xlabel('Time relative to ASM admin');
ylabel('normalized spikes');
title('Average spike rate after ASM administration: Night');
grid on;

%% make a plot comparing the pre and post admin spike rate for day and night. color code by mean asm load

% figure;
tbl_day = table;
tbl_day.pre_asm_spikes = mean(all_trials_day(:,1:end/2),2); % mean of the pre period for each trial - day
tbl_day.post_asm_spikes = mean(all_trials_day(:,end/2+1:end),2); % mean of the post period for each trial - day 
tbl_day.asm_load = mean(all_asms_day(:,end/2+1:end),2);
tbl_day.ptID = all_ptids_day;
tbl_day.change_spike_rate = tbl_day.post_asm_spikes - tbl_day.pre_asm_spikes;% should be a negative number if higher before admin

tbl_night = table;
tbl_night.pre_asm_spikes = mean(all_trials_night(:,1:end/2),2);
tbl_night.post_asm_spikes = mean(all_trials_night(:,end/2+1:end),2);
tbl_night.asm_load = mean(all_asms_night(:,end/2+1:end),2);
tbl_night.ptID = all_ptids_night;
tbl_night.change_spike_rate = tbl_night.post_asm_spikes - tbl_night.pre_asm_spikes;% should be a negative number if higher before admin

% do stats across patients 
plot([1 2], [tbl_day.pre_asm_spikes tbl_day.post_asm_spikes],'.k','markersize',10); hold on;
plot([1 2], [tbl_day.pre_asm_spikes tbl_day.post_asm_spikes],'-k');


[h,p,ci,stats] = ttest(tbl_day.pre_asm_spikes, tbl_day.post_asm_spikes);
[h,p,ci,stats] = ttest(tbl_night.pre_asm_spikes, tbl_night.post_asm_spikes);
% result: no significant change in spike rate before and after (??)

% do stats within patients across trials - 
% ASM load does not affect the change in spike rate around administration
modelspec = 'change_spike_rate ~asm_load+(1|ptID)';
mdl_day =fitglme(tbl_day,modelspec) % 'Distribution','Poisson';
mdl_night =fitglme(tbl_night,modelspec) % 'Distribution','Poisson';

% which patients show decrease in spikes after admin
patients_decrease = unique(tbl_day.ptID(tbl_day.change_spike_rate <  0));

%% save the results
% dirr = '/Volumes/users/nghosn3/Pioneer/spikes-AED/ASM-spikes-analysis/';
% save([dirr 'asm_erp_analysis_results.mat'],'patients_decrease','tbl_day','tbl_night','all_trials_day','all_trials_night','all_ptids_night','all_ptids_day','spike_trials','drug_trials','trial_labels')

%% subplot(2,1,2)

tiledlayout('flow');
meanSignal = nanmean(all_trials_day, 1);
stdSignal = nanstd(all_trials_day, 1);

% Calculate the confidence interval bounds
alpha = 0.05;  % Significance level (1 - confidence levelc)
n = size(all_trials_day, 1);
t = tinv(1 - alpha/2, n - 1);  % t-value for two-tailed t-distribution
ciLower = meanSignal - t * stdSignal / sqrt(n);
ciUpper = meanSignal + t * stdSignal / sqrt(n);


time = x;  % Time axis
nexttile([1 2])
plot(time, meanSignal, 'b', 'LineWidth', 2);
hold on;
fill([time, fliplr(time)], [ciLower, fliplr(ciUpper)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xline(0,'r--','linewidth',2)

legend('mean spike rate','95% CI','medication admin')
% Customize the plot
xlabel('Time relative to ASM admin');
ylabel('normalized spikes');
title('Average spike rate after ASM administration: day');
grid on;

nexttile
histogram(tbl_day.change_spike_rate);
title('change after ASM administration','fontsize',14);
xlabel('change in spike rate (z-score)','fontsize',12);
ylabel('density','fontsize',12)




