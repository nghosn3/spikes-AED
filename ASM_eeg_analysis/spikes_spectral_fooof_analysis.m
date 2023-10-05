%% spectral analysis of spike rate using wavelet transform and FOOOF

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

%% find patients for spectral analysis: 1) taper and 'off' of at least one to all medications for >2days
pt_data_clips = table();
pt_data_clips.ptID = (ptIDs);
pt_data_clips.med_names = all_med_names';
pt_data_clips.inds = cell(length(ptIDs),1);
fooof_results = struct();

days=2; % threshold in days
taper_thresh = days*(24*60)./10; % 10min segments in 2 day period

for ipt =1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))]
    asm_curves = all_pts_drug_samp{ipt};
    med_names = all_med_names{ipt};
    offsets = all_ieeg_offset{2,ipt};
    ieeg_offset_datasets = all_ieeg_offset{1,ipt};

    %find when each drug was stopped and find period of halt
    inds = nan(length(med_names),3);
    for n=1:length(med_names)
        if ~(contains(med_names{n},'lorazepam'))
            med_curve = asm_curves(n,:);
            [pks,locs] = findpeaks(med_curve);
            diffs = diff(locs);
            if length(diffs)>1
                locs(end+1) = length(med_curve);
                diffs = diff(locs);
            end
            last_admin = find(diffs >= taper_thresh);
            if ~isempty(last_admin)
                stop_med =locs(last_admin(1))+6; %ind+1 of last administration
                restart_med =locs(last_admin(1)+1)-1; % ind -1 of next peak after stopping

                inds(n,1) = 1;
                inds(n,2) = stop_med;
                inds(n,3) = restart_med;

                %              %manual testing: HUP 137
                %                 inds(n,1) = 1;
                %                 inds(n,2) = 106*6;
                %                 inds(n,3) = 146*6;


            end
        end
    end
    pt_data_clips.inds(ipt) = {inds};

    %     % add seizure times
    %     [seizure_times] = get_seizure_times_from_sheet(ptID);
    %
    %     sz_inds = zeros(length(seizure_times),1);
    %     for sz = 1:length(seizure_times)
    %         [~,sz_ind] = min(abs(all_spike_times{ipt}-seizure_times(sz,1)));
    %         sz_inds(sz) = sz_ind;
    %     end
    %     pt_data_clips.sz_inds(ipt) = {sz_inds};

end


% get soz localization (interested in temporal patients)
for i=1:(height(pt_data_clips))
    soz_loc = soz_info.region(contains(soz_info.name,num2str(pt_data_clips.ptID(i))));
    pt_data_clips.soz_loc(i) = soz_loc;
    pt_data_clips.temporal(i) = contains(lower(soz_loc),'temporal');
end

%%

% set constants
plot_spec = 0;
plot_fooofs = 0;
do_fooof = 1;
fs = 6 ; % 6 sample every hour - units are in hours
prom_thresh = 0.2; % for grabbing peaks on periodogram
err = 0.3; % bounds around frequency of interest

% For inverse Wavelet Transform: Setting Constants from Torrence and Compo (TBL 2)
dt = 1/fs;
Cd = 0.776;
Psi0 = pi^(-.25);
dj = 1/12;
modifier = (dj*dt^(1/2))/(Cd*Psi0);
% Defining inverse wavelet function
invcwt = @ (wave,scale,pk_locs) modifier*sum(real(wave(pk_locs,:))./(scale(pk_locs)'.^(1/2)),1);


close all;

spec_data = table();
excluded = false(length(ptIDs),1);
tic
for ipt =1:length(ptIDs)

    ptID = ['HUP' num2str(ptIDs(ipt))];

    %if ~(strcmp(ptID,'HUP177') || strcmp(ptID,'HUP166') || strcmp(ptID,'HUP191') )
        med_names = all_med_names{ipt};
        asm_load = all_pts_drug_samp{ipt};
        spikes = (all_spike_rate{ipt}+1);
        spikes = spikes./max(spikes);

        spec_data.ptID(ipt) = {ptID};

        signal = zscore(spikes);
        exact_time = all_spike_times{ipt};
        time =linspace(0,max(exact_time),length(signal))./3600;


        % at least one medication tapered for >2days - pre/post condition breakpoints
        all_asm_inds = pt_data_clips.inds(ipt); all_asm_inds = all_asm_inds{:};
        dosing_periods = zeros(length(med_names),1);
        dosing_periods_amp = zeros(length(med_names),1);
        before_after_peak_periods = zeros(length(med_names),2);
        spike_peaks_before = [];
        spike_peaks_after = [];
        pre_post_med_spike_peak = zeros(length(med_names),2); % get the freq that the peak occur in, pre/post

        if ~all(all(isnan(all_asm_inds)))
            ind1 = max(all_asm_inds(:,2));
            ind2 = min(all_asm_inds(:,3));

            if (ind2-ind1) <= taper_thresh
                excluded(i) = true;

            else
                exact_time = all_spike_times{ipt};
                time =linspace(0,max(exact_time),length(spikes))./3600;

                asm_before = asm_load(:,1:ind1);
                asm_after = asm_load(:,ind1+1:ind2);

                spikes_before = signal(:,1:ind1);
                spikes_after = signal(:,ind1+1:ind2);

                % determine dosing freqs in before period using wavelet transform

                for m =  1:length(med_names)
                    sig = asm_before(m,:)';
                    [pks,locs] = findpeaks(sig);
                    if length(diff(locs)) > 1 % at least 2 administrations in the dosing period
                        target_dose_int = mean(diff(time(locs))); % in hours from med curve
                        if ~(contains(med_names{m},'lorazepam')) && ~isempty(pks) && ~isnan(target_dose_int)
                            [wave,periods,~,~] = wt([time(1:ind1)',sig],hours(1/fs));
                            periodogram = mean(abs(wave).^2,2);
                            [pks,locs] = findpeaks(periodogram./max(periodogram)); %,'minpeakprominence',.3
                            if ~isempty(locs) %&& round(max(periods))>=20
                                [per_diff,loc_int] = min(abs(periods(locs)-target_dose_int));
                                if (per_diff < target_dose_int*.2) % verify that it is a peak period in the signal - 33% error
                                    pk_periods = (periods(locs));
                                    dosing_freqs = (pk_periods(loc_int));
                                    dosing_periods(m) = round(dosing_freqs*10)./10;
                                    [~,ind]=min(abs(periods-dosing_freqs));
                                    dosing_periods_amp(m) = periodogram(ind);
                                else
                                    dosing_periods(m) = target_dose_int;
                                    [~,ind]=min(abs(periods-target_dose_int));
                                    dosing_periods_amp(m) = periodogram(ind);

                                end
                            end
                        end
                    end

                end

                if plot_spec % plot the taper period, the asms, & spikes
                    figure;clf; hold on %#ok<*UNRCH>
                    subplot(2,1,1)
                    plot(time,signal)
                    xline(time(ind1),'r-','linewidth',3); xline(time(ind2),'r-','linewidth',3)
                    xlabel('Time(hours)')
                    subplot(2,1,2)
                    plot(time, asm_load'); hold on;
                    xline(time(ind1),'r-','linewidth',3); xline(time(ind2),'r-','linewidth',3)
                    legend(med_names);title(ptID)
                end

                % amplitude of periodogram peak pre/post - normalized
                % across subjects
                [wave,period,scale,coi] = wt([time',signal'],hours(1/fs));
                pre_periodogram = mean(abs(wave(:,1:ind1)).^2,2);
                post_periodogram = mean(abs(wave(:,ind1:ind2)).^2,2);

                % FOOOF settings
                freqs = period(period <=48)';
                settings = struct();  % Use defaults
%                 settings.min_peak_height = .2;
%                 settings.peak_threshold = 1;
%                 settings.aperiodic_mode = 'fixed';


                f_range = [min(freqs), 40]; %40hrs

                % Run FOOOF, also returning the model
                if do_fooof
                    try
                        fooof_results_pre = fooof(freqs, pre_periodogram(period' <=48)', f_range, settings, true);
                        fooof_results_post = fooof(freqs, post_periodogram( period' <=48)', f_range, settings, true);


                        % Plot the resulting models
                        if plot_fooofs
                            fooof_plot(fooof_results_pre);
                            title([ptID ' spikes during dosing'])
                            fooof_plot(fooof_results_post);
                            title([ptID ' spikes after dosing'])
                        end

                        %store the fooof results - the peak frequencies in the
                        %periodic fit, and the height (logscale) from the baseline
                        %aperiodic 1/f fit
                        [pre_peaks,pre_locs,width_pre,~] = findpeaks(fooof_results_pre.fooofed_spectrum);
                        model_fit_pre = fooof_results_pre.fooofed_spectrum-fooof_results_pre.ap_fit;
                        pre_spike_peaks = model_fit_pre(pre_locs);

                        [post_peaks,post_locs,width_post,~] = findpeaks(fooof_results_post.fooofed_spectrum);
                        model_fit_post = fooof_results_post.fooofed_spectrum-fooof_results_post.ap_fit;
                        post_spike_peaks = model_fit_post(post_locs);

                        spike_peaks_before = [freqs(pre_locs) pre_spike_peaks' width_pre'];
                        spike_peaks_after = [freqs(post_locs) post_spike_peaks' width_post'];

                    catch
                        disp([ptID ' fooof error'])

                    end
                end

            end

            pt_data_clips.dosing_periods(ipt) = {dosing_periods};
            pt_data_clips.dosing_periods_amp(ipt) = {dosing_periods_amp};
            pt_data_clips.spike_fooof_peakfreqs_pre(ipt) = {spike_peaks_before};
            pt_data_clips.spike_fooof_peakfreqs_post(ipt) = {spike_peaks_after};
        end
    %end
end
toc

%% test all spike amp in  pre/post dosing period:
close all;
% temporal patients only
temporal_pts = pt_data_clips(pt_data_clips.temporal==1,:);
med_label =  temporal_pts.med_names(:);
med_label =  vertcat(med_label{:});

err = 0.1; %small error for identifying 12hr peak to differentiate dosing.. ex dosing at 13hr or 15 isnot circadian

% is there a peak in the dosing period? binomial test (y/n) - does ASM
% dosing influence spikes in a periodic fashion?
has_dosing_peak = nan(1,length(ptIDs));
for i = 1:length(ptIDs)
    dosing_periods = pt_data_clips.dosing_periods{i};
    spike_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i};
    if any(dosing_periods~= 0) && ~isempty(spike_peaks)
        has_peak = 0;
        spike_peaks = spike_peaks(:,1);
        for p = 1:length(dosing_periods)
            has_peak = has_peak | any(abs(dosing_periods(p)-spike_peaks) <= err*spike_peaks);
        end
        has_dosing_peak(i)=has_peak;
    end

end

dosed_q12 = nan(1,length(ptIDs));
has_12_peak =nan(1,length(ptIDs));
for i = 1:length(ptIDs)
    dosing_periods = pt_data_clips.dosing_periods{i};
    spike_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i};
    if any(dosing_periods~= 0) && ~isempty(spike_peaks)
        dosed_q12(i) = any(abs(dosing_periods-12) <= err*12);
        spike_peaks = spike_peaks(:,1);
        has_12_peak(i) = any(abs(spike_peaks-12) <= err*12);
    end

end


% conduct binomial test
n = sum(~isnan(has_dosing_peak)); % Total number of trials (patients)
k = nansum(has_dosing_peak);  % Number of successes (peaks at dosing frequency)

% Define the null hypothesis probability
p0 = 0.5; % Expected probability under the null hypothesis

% Perform the binomial test
p_spike_peak = 2 * min(binocdf(k, n, p0), 1 - binocdf(k - 1, n, p0))

% conduct a fishers test: dosed every 12 and has 12hr peak
table_12hr = table();
table_12hr.dosed_q12 = dosed_q12';
table_12hr.has_12_peak = has_12_peak';
[tbl,~,~,labels] = crosstab(table_12hr.dosed_q12,table_12hr.has_12_peak)
[p_12hr, ~] = fishertest(tbl, 'Tail', 'both')
heatmap(table_12hr,'dosed_q12','has_12_peak')



%% all patients

all_dosing_periods = pt_data_clips.dosing_periods(:); 
spike_fooof_peaks_pre =  vertcat(pt_data_clips.spike_fooof_peakfreqs_pre{:});
spike_fooof_peaks_post =  vertcat(pt_data_clips.spike_fooof_peakfreqs_post{:});

% for each patient, see what happens to the width and amp - normalize  columns
amp_max = max([spike_fooof_peaks_pre(:,2); spike_fooof_peaks_post(:,2) ]);
width_max = max([spike_fooof_peaks_pre(:,3); spike_fooof_peaks_post(:,3) ])

spike_fooof_peaks_pre(:,2) = spike_fooof_peaks_pre(:,2)./amp_max;
spike_fooof_peaks_pre(:,3) = spike_fooof_peaks_pre(:,3)./width_max;

spike_fooof_peaks_post(:,2) = spike_fooof_peaks_post(:,2)./amp_max;
spike_fooof_peaks_post(:,3) = spike_fooof_peaks_post(:,3)./width_max;

edges = 0:3:ceil(max([spike_fooof_peaks_pre(:,1); spike_fooof_peaks_post(:,1) ]));

freq_bins_pre = discretize(spike_fooof_peaks_pre(:,1),edges);
freq_bins_post = discretize(spike_fooof_peaks_post(:,1),edges);

half_max_pre = spike_fooof_peaks_pre(:,2)./spike_fooof_peaks_pre(:,3);
half_max_post = spike_fooof_peaks_post(:,2)./spike_fooof_peaks_post(:,3);

%% compare the peakiness of the spike peak at dosing period to distance from
% 12hr dosing during pre period
tbl=table();
tbl.ptid=zeros(0);
tbl.dosing_period=zeros(0);

tbl.fooof_period_pre=zeros(0);
tbl.fooof_width_pre =zeros(0);
tbl.fooof_amp_pre =zeros(0);

tbl.fooof_period_post=zeros(0);
tbl.fooof_width_post=zeros(0);
tbl.fooof_amp_post=zeros(0);
err = .1;
next=1;
for i = 1:height(pt_data_clips)
    dosing_periods = pt_data_clips.dosing_periods{i};
    fooof_pre_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i};
    fooof_post_peaks = pt_data_clips.spike_fooof_peakfreqs_post{i};
    if any(dosing_periods~=0) && ~isempty(fooof_pre_peaks)
        for d = 1:length(dosing_periods)
            tbl.dosing_period(next) = dosing_periods(d);
            tbl.ptid(next) = pt_data_clips.ptID(i);
            ind = find(abs(dosing_periods(d) - fooof_pre_peaks(:,1)) <= err*fooof_pre_peaks(:,1)); % find the dosing peak on the pre-spec

            if ~isempty(ind)
                tbl.fooof_period_pre(next) = fooof_pre_peaks(ind,1);
                tbl.fooof_width_pre(next) = fooof_pre_peaks(ind,3);
                tbl.fooof_amp_pre(next) = fooof_pre_peaks(ind,2);
            end 
             ind_post = find(abs(dosing_periods(d) - fooof_post_peaks(:,1)) <= err*fooof_post_peaks(:,1)); % find the dosing peak on the post-spec if it exists

                if isempty(ind_post) && ~isempty(ind)
                    tbl.fooof_amp_post(next) = 0;
                    tbl.fooof_width_post(next) = 0;
                    tbl.fooof_period_post(next)=0;
                elseif  ~isempty(ind_post)
                    tbl.fooof_period_post(next) = fooof_post_peaks(ind_post,1);
                    tbl.fooof_width_post(next)=fooof_post_peaks(ind_post,3);
                    tbl.fooof_amp_post(next)=fooof_post_peaks(ind_post,2);
                end 
next=next+1;
        end
        
    end

end

% look at the peaki-ness of the peaks -p, prominance, w, width
% p/w ratio: the bigger it is, the peakier the peak
peri_err = 0.3;
peri_12_dosed  = abs(tbl.dosing_period-12) <= 12*peri_err;

x = (tbl.dosing_period(peri_12_dosed)-12);
y = tbl.fooof_width_pre(peri_12_dosed);
[h,p,ci,stats] = ttest(x,y)
subplot(3,1,1)
plot(x,y,'.','markersize',15);
title('peak width vs dosing period-12hrs')
ylabel('peak width')
xlabel('dosing period - 12');
[r,p]=corr(x,y,'type','pearson')
legend(['r = ' num2str(r) ', p = ' num2str(p)])

subplot(3,1,2)
x = (tbl.dosing_period(peri_12_dosed)-12);
y = tbl.fooof_amp_pre(peri_12_dosed);
plot(x,y,'.','markersize',15);
title('peak amplitude vs dosing period-12hrs')
ylabel('peak amplitude')
xlabel('dosing period - 12');
[r,p]=corr(x,y,'type','pearson')
legend(['r = ' num2str(r) ', p = ' num2str(p)])


subplot(3,1,3)
x = (tbl.dosing_period(peri_12_dosed)-12);
y = tbl.fooof_amp_pre(peri_12_dosed)./tbl.fooof_width_pre(peri_12_dosed);
plot(x,y,'.','markersize',15);
title('peak amplitude/width vs dosing period-12hrs')
ylabel('peak amp/width')
xlabel('dosing period - 12');
[r,p]=corr(x,y,'type','pearson')
legend(['r = ' num2str(r) ', p = ' num2str(p)])


figure;
scatter(tbl.dosing_period(tbl.dosing_period~= 0),tbl.fooof_period_pre(tbl.dosing_period~= 0))
hold on; plot([0 30],[0 30]);
title('dosing period vs fooof peak period');
xlabel('dosing period'); ylabel('peak in fooofed model')
ylim([0 30])
xlim([0 30])
% try a linear mixed model to tease out the effect of medication taper on
% spikes at dosing frequency


%% try a monte carlo simulation to see if 12hr peak is associated with 12hour dosing
% dosingPeak= dosed_q12(~isnan(dosed_q12));
% spectrogramPeak =has_12_peak(~isnan(dosed_q12));
% 
% % Observed test statistic
% observedStatistic = sum(dosingPeak & spectrogramPeak) / sum(spectrogramPeak);
% 
% % Set the number of simulations
% numSimulations = 10000;
% 
% % Perform the Monte Carlo simulation
% nullDistribution = zeros(1, numSimulations);
% for i = 1:numSimulations
%     % Permute the dosingPeak variable
%     permutedDosingPeak = dosingPeak(randperm(length(dosingPeak)));
% 
%     % Calculate the test statistic for the permuted data
%     nullDistribution(i) = sum(permutedDosingPeak & spectrogramPeak) / sum(spectrogramPeak);
% end
% 
% % Calculate the p-value
% pValue = sum(nullDistribution >= observedStatistic) / numSimulations
% 
% 



