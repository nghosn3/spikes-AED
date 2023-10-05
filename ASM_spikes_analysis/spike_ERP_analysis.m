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
q12_only =  [225	224	219	217	216	215	214	213	211	210	209	208	207	205	196	195	194	193	191	190	189 186	185	182	179	178	177	174	173	171	170	168	165	164	163	159	156	155	148	147	145	144	143	142	141	140	138];


spike_trials = cell(1,length(q12_only));
drug_trials = cell(1,length(q12_only));

for ipt =1:length(q12_only)
    ptID = ['HUP' num2str(q12_only(ipt))];
    asm_curves = all_pts_drug_samp{ipt};
    med_names = all_med_names{ipt};
    tvec = all_pts_tvec{ipt}{:};
    spikes = (all_spike_rate{ipt}+1);
    spikes = spikes./max(spikes);
    %spikes = smoothdata(spikes);

    %find when each drug was stopped and find period of halt
    admin_times = zeros(length(med_names),50); %pre-allocate to store admin locs
    max_admins =0;
    for n=1:length(med_names)
        if ~(contains(med_names{n},'lorazepam'))
            med_curve = asm_curves(n,:);
            [pks,locs] = findpeaks(-med_curve);
            admin_times(n,1:length(locs)) = locs;
            max_admins = max([max_admins length(locs)]);

        end
    end

    % find the admin times to be used as trials
    admin_times = admin_times(:,1:max_admins);
    % remove all zero rows
    admin_times((all(admin_times == 0, 2)),:)=[];
    % remove rows for drugs that were started late in the emu stay

    admins = zeros(height(admin_times),length(spikes));
    for r = 1:height(admins)
        inds = admin_times(r,:); inds(inds==0)=[];
        admins(r,inds) = 1;

    end
    % find times where they are administed at the same time
    all_drugs_admined  = find(sum(admins)==r);
    % find times where all drugs are administered within an hour of each
    % other
    sum1 = (sum(admins));
    inds = find(sum1~=0 & sum1~=r);

    add_inds = find(diff(inds)<=6);

    for id = add_inds'
        if sum(sum1(inds(id):inds(id)+6)) == r
            last_id = find(sum1(inds(id):inds(id)+6) ~= 0);
            all_drugs_admined = [all_drugs_admined inds(id)+(last_id(end)-1)];
        end
    end

    all_drugs_admined = sort(all_drugs_admined);
    pt_trials = zeros(length(all_drugs_admined),(12*6)+1); % 12hours by number of trials
    pt_drug_trials = zeros(length(all_drugs_admined),(12*6)+1); % 12hours by number of trials

    avg_asm_curves = mean(asm_curves);
    offset = 6*6; %6hours
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
close all;
all_trials = [];
figure;
inc=0;
for i = 1:length(spike_trials)
    pt_trials = spike_trials{i};
    asm_trials = drug_trials{i};
    pt_trials((all(pt_trials == 0, 2)),:)=[];


    all_trials = [all_trials; pt_trials];
    if ~isempty(pt_trials)
        inc = inc+1;
        figure;
        subplot(2,1,1)
        x = linspace(-6,6,width(pt_trials));
        plot(x',pt_trials+(inc*.25));
        hold on;
        plot(x',mean(pt_trials)+(inc*.25),'k','linewidth',1)
        title(['spike rate: HUP' num2str(ptIDs(i))])

                % plot the avg asm load over each other to make sure its indexing properly?
                subplot(2,1,2)
                plot(x',asm_trials); hold on;
                plot(x',mean(asm_trials));
                xline(0,'r--')
                title(['ASM load: HUP' num2str(ptIDs(i))])

    end
end
xline(0,'r--')
ylabel('patient#');
xlabel('time hrs')
title('spike rate')
%%
figure;
plot(x',(all_trials),'b','linewidth',.1); hold on;
plot(x',mean(all_trials),'k','linewidth',3); hold on;
title('all trials across patients')
ylabel('normalized spikes')
xlabel('time hrs')


% Calculate the mean signal and standard deviation
meanSignal = mean(all_trials, 1);
stdSignal = std(all_trials, 1);

% Calculate the confidence interval bounds
alpha = 0.05;  % Significance level (1 - confidence level)
n = size(all_trials, 1);
t = tinv(1 - alpha/2, n - 1);  % t-value for two-tailed t-distribution
ciLower = meanSignal - t * stdSignal / sqrt(n);
ciUpper = meanSignal + t * stdSignal / sqrt(n);

% Plot the mean signal and confidence interval
figure;
time = x;  % Time axis
plot(time, meanSignal, 'b', 'LineWidth', 2);
hold on;
fill([time, fliplr(time)], [ciLower, fliplr(ciUpper)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
[bkpts,res] = findchangepts(meanSignal,MaxNumChanges=2);
xline(time(bkpts),'r','LineWidth',2)
xline(0,'k--')
legend('mean spike rate','95% CI','change in mean','','medication admin')
% Customize the plot
xlabel('Time relative to ASM admin');
ylabel('normalized spikes');
title('Average spike rate after ASM administration');
grid on;





% need to split between day and night dosing, maybe something happens
