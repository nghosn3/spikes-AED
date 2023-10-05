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

%%  include all patients for spectral analysis: dosed with medication before taper
pt_data_clips = table();
pt_data_clips.ptID = (ptIDs);
pt_data_clips.med_names = all_med_names';
pt_data_clips.inds = cell(length(ptIDs),1);
fooof_results = struct();

days=2; % threshold in days
taper_thresh = days*(24*60)./10; % 10min segments in 2 day period

for ipt =1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    asm_curves = all_pts_drug_samp{ipt};
    med_names = all_med_names{ipt};
    offsets = all_ieeg_offset{2,ipt};
    ieeg_offset_datasets = all_ieeg_offset{1,ipt};

    %find when each drug was stopped and find period of halt
    inds = nan(length(med_names),2);
    for n=1:length(med_names)
        if ~(contains(med_names{n},'lorazepam'))
            med_curve = asm_curves(n,:);
            [pks,locs] = findpeaks(med_curve);
            if length(locs)>2 && (locs(1) <= taper_thresh)
                diffs = diff(locs);
%                 if length(diffs)>1
%                     locs(end+1) = length(med_curve);
%                     diffs = diff(locs);
%                 end
                last_admin = find(diffs >= taper_thresh);
                if ~isempty(last_admin)
                    stop_med =locs(last_admin(1))+12; %ind+1 of last administration
                else
                    stop_med = locs(end)+12;
                    stop_med = min([stop_med length(med_curve)]);
                end

                inds(n,1) = 1;
                inds(n,2) = stop_med;
            end
        end
    end
    pt_data_clips.inds(ipt) = {inds};

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

excluded = false(length(ptIDs),1);
tic
for ipt =1:length(ptIDs)%[67,79]

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
    dosing_periods_wavelet = zeros(length(med_names),1);
    dosing_periods = zeros(length(med_names),1);
    
    dosing_periods_var = zeros(length(med_names),1);
    dosing_periods_amp = zeros(length(med_names),1);
    
    before_after_peak_periods = zeros(length(med_names),2);
    spike_peaks_before = [];
    spike_peaks_after = [];
    pre_post_med_spike_peak = zeros(length(med_names),2); % get the freq that the peak occur in, pre/post

    if ~all(all(isnan(all_asm_inds)))
        ind1 = max(all_asm_inds(:,2)); % min ; when last med is stopped


        exact_time = all_spike_times{ipt};
        time =linspace(0,max(exact_time),length(spikes))./3600;

        asm_before = asm_load(:,1:ind1);
        spikes_before = signal(:,1:ind1);

        % determine dosing freqs in before period using wavelet transform

        for m =  1:length(med_names)
            sig = asm_before(m,:)';
            [pks,locs] = findpeaks(sig);
            if length(diff(locs)) >= 2 % at least 3 administrations in the dosing period
                target_dose_int = median(diff(time(locs))); % in hours from med curve
                dosing_periods(m) = target_dose_int;
                dosing_periods_var(m) = std(diff(time(locs))); % standard deviation 
                if ~(contains(med_names{m},'lorazepam')) && ~isempty(pks) && ~isnan(target_dose_int)
                    [wave,periods,~,~] = wt([time(1:ind1)',sig],hours(1/fs));
                    periodogram = mean(abs(wave).^2,2);
                    [pks,locs] = findpeaks(periodogram./max(periodogram)); %,'minpeakprominence',.3
                    if ~isempty(locs) %&& round(max(periods))>=20
                        [per_diff,loc_int] = min(abs(periods(locs)-target_dose_int));
                        if (per_diff < target_dose_int*.2) % verify that it is a peak period in the signal - 33% error
                            pk_periods = (periods(locs));
                            dosing_freqs = (pk_periods(loc_int));
                            dosing_periods_wavelet(m) = round(dosing_freqs*10)./10;
                            [~,ind]=min(abs(periods-dosing_freqs));
                            dosing_periods_amp(m) = periodogram(ind);
                        else
                            dosing_periods_wavelet(m) = NaN;
                        end
                    end
                end
            end
        end

        if plot_spec % plot the taper period, the asms, & spikes
            figure;clf; hold on %#ok<*UNRCH>
            subplot(2,1,1)
            plot(time,signal)
            xline(time(ind1),'r-','linewidth',3); %xline(time(ind2),'r-','linewidth',3)
            xlabel('Time(hours)')
            subplot(2,1,2)
            plot(time, asm_load'); hold on;
            xline(time(ind1),'r-','linewidth',3); %xline(time(ind2),'r-','linewidth',3)
            legend(med_names);title(ptID)
        end

        % amplitude of periodogram peak pre/post - normalized
        % across subjects
        [wave,period,scale,coi] = wt([time',signal'],hours(1/fs));
        pre_periodogram = mean(abs(wave(:,1:ind1)).^2,2);

        % FOOOF settings
        freqs = period(period <=48)';
        settings = struct();  % Use defaults
        settings.min_peak_height = .2;
        settings.peak_threshold = 1;
        settings.aperiodic_mode = 'fixed';


        f_range = [min(freqs), 40]; %40hrs

        % Run FOOOF, also returning the model
        if do_fooof
            try
                try % try with pre-set settings that are more generous 
                    try
                    fooof_results_pre = fooof(freqs, pre_periodogram(period' <=48)', f_range, settings, true);
                    catch
                        settings.aperiodic_mode = 'knee';
                         fooof_results_pre = fooof(freqs, pre_periodogram(period' <=48)', f_range, settings, true);
                    end 

                    % Plot the resulting models
                    if plot_fooofs
                        fooof_plot(fooof_results_pre);
                        title([ptID ' spikes during dosing'])
                    end

                    [pre_peaks,pre_locs,width_pre,~] = findpeaks(fooof_results_pre.fooofed_spectrum);
                    model_fit_pre = fooof_results_pre.fooofed_spectrum-fooof_results_pre.ap_fit;
                    pre_spike_peaks = model_fit_pre(pre_locs);

                    spike_peaks_before = [freqs(pre_locs) pre_spike_peaks' width_pre'];

                catch % of that doesnt work then use default settings 
                    settings = struct();
                    fooof_results_pre = fooof(freqs, pre_periodogram(period' <=48)', f_range, settings, true);

                    % Plot the resulting models
                    if plot_fooofs
                        fooof_plot(fooof_results_pre);
                        title([ptID ' spikes during dosing'])
                    end

                    %store the fooof results 
                    [pre_peaks,pre_locs,width_pre,~] = findpeaks(fooof_results_pre.fooofed_spectrum);
                    model_fit_pre = fooof_results_pre.fooofed_spectrum-fooof_results_pre.ap_fit;
                    pre_spike_peaks = model_fit_pre(pre_locs);

                    spike_peaks_before = [freqs(pre_locs) pre_spike_peaks' width_pre'];
                end
            catch % if that doesnt work, quit
                disp([ptID ' fooof error'])

            end
        end

    end

    pt_data_clips.dosing_periods(ipt) = {dosing_periods};
    pt_data_clips.dosing_periods_wavelet(ipt) = {dosing_periods_wavelet};
    pt_data_clips.dosing_periods_amp(ipt) = {dosing_periods_amp};
    pt_data_clips.dosing_periods_var(ipt) = {dosing_periods_var};
    pt_data_clips.spike_fooof_peakfreqs_pre(ipt) = {spike_peaks_before};
end


toc

%% compare the peakiness of the spike peak at dosing period to distance from
% 12hr dosing during pre period
tbl=table();
tbl.ptid=zeros(0);
tbl.dosing_period=zeros(0);
tbl.dosing_period_var=zeros(0);
tbl.med_name = cell(0);


tbl.fooof_period_pre=zeros(0);
tbl.fooof_width_pre =zeros(0);
tbl.fooof_amp_pre =zeros(0);

tbl.fooof_period_post=zeros(0);
tbl.fooof_width_post=zeros(0);
%tbl.fooof_amp_post=zeros(0);
err = .3;
next=1;
for i = 1:height(pt_data_clips)
   
    % find dosing periods that were not detected and add them manually
    zeros_manual = pt_data_clips.dosing_periods{i} == 0 |  isnan(pt_data_clips.dosing_periods{i});
    zeros_wavelet = pt_data_clips.dosing_periods_wavelet{i} == 0 |  isnan(pt_data_clips.dosing_periods_wavelet{i});
    use_manual = find(zeros_manual ~= 1 & zeros_wavelet == 1);% when dosing is not zero but wavlet is zero
    
    temp = pt_data_clips.dosing_periods_wavelet{i};%% CHOOSE TO USE WAVELET DOSING PERIOD OR MANUAL
    temp_manual  = pt_data_clips.dosing_periods{i};

    dosing_periods = temp;
    dosing_periods(use_manual) = temp_manual(use_manual);
    
    dosing_periods_var = pt_data_clips.dosing_periods_var{i};
    fooof_pre_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i};
    med_names = pt_data_clips.med_names{i};
    %fooof_post_peaks = pt_data_clips.spike_fooof_peakfreqs_post{i};
    if any(dosing_periods~=0) && ~isempty(fooof_pre_peaks)
        for d = 1:length(dosing_periods)
            tbl.dosing_period(next) = dosing_periods(d);
            tbl.dosing_period_var(next) = dosing_periods_var(d);
            tbl.ptid(next) = pt_data_clips.ptID(i);
            tbl.med_name(next) = med_names(d);
            ind = find(abs(dosing_periods(d) - fooof_pre_peaks(:,1)) <= err*fooof_pre_peaks(:,1)); % find the dosing peak on the pre-spec

            if ~isempty(ind)
                tbl.fooof_period_pre(next) = fooof_pre_peaks(ind(1),1);
                tbl.fooof_width_pre(next) = fooof_pre_peaks(ind(1),3);
                tbl.fooof_amp_pre(next) = fooof_pre_peaks(ind(1),2);
            end
            %ind_post = find(abs(dosing_periods(d) - fooof_post_peaks(:,1)) <= err*fooof_post_peaks(:,1)); % find the dosing peak on the post-spec if it exists

            %                 if isempty(ind_post) && ~isempty(ind)
            %                     tbl.fooof_amp_post(next) = 0;
            %                     tbl.fooof_width_post(next) = 0;
            %                     tbl.fooof_period_post(next)=0;
            %                 elseif  ~isempty(ind_post)
            %                     tbl.fooof_period_post(next) = fooof_post_peaks(ind_post,1);
            %                     tbl.fooof_width_post(next)=fooof_post_peaks(ind_post,3);
            %                     tbl.fooof_amp_post(next)=fooof_post_peaks(ind_post,2);
            %                 end
            next=next+1;
        end

    end

end

% remove entries that correspons to ativan 
ativan_inds = contains(tbl.med_name,'lorazepam');
tbl(ativan_inds,:) = [];

%% look at the peaki-ness of the peaks -p, prominance, w, width
% p/w ratio: the bigger it is, the peakier the peak

close all;
peri_err = 0.3;

has_peak = tbl.fooof_period_pre ~= 0;
ptIDs_has_peak = unique(tbl.ptid(tbl.fooof_period_pre ~= 0));
no_peak = ~has_peak;


%x = abs(tbl.dosing_period(has_peak)-tbl.fooof_period_pre(has_peak));
x = tbl.dosing_period_var(has_peak);
y = tbl.fooof_width_pre(has_peak);
[r,p]=corr(x,y,'type','pearson');

figure;
subplot(1,3,1)
plot(x,y,'.','markersize',15);
axis square

title('width')
ylabel('peak width')
xlabel('std(dosing intervals)');

legend(['r = ' num2str(r) ', p = ' num2str(p)])

subplot(1,3,2)
y = tbl.fooof_amp_pre(has_peak);
plot(x,y,'.','markersize',15);
axis square

title('amplitude')
ylabel('peak amplitude')
xlabel('std(dosing intervals)');
[r,p]=corr(x,y,'type','pearson');
legend(['r = ' num2str(r) ', p = ' num2str(p)])


subplot(1,3,3)
y = tbl.fooof_period_pre(has_peak)-tbl.dosing_period(has_peak)
%y = (tbl.fooof_width_pre(has_peak)./max(tbl.fooof_width_pre(has_peak)) ./ (tbl.fooof_amp_pre(has_peak)./max(tbl.fooof_amp_pre(has_peak))));
%[f,gof] = fit(x,y,'exp2')
plot(x,y,'.','markersize',15)%'.','markersize',15);
axis square

title('fatness of peak')
ylabel('peak width/amp')
xlabel('std(dosing intervals)');
[r,p]=corr(x,y,'type','pearson');
legend(['r = ' num2str(r) ', p = ' num2str(p)])

% peak width and difference in peak period from dosing period - new plot

x = abs(tbl.dosing_period(has_peak)-tbl.fooof_period_pre(has_peak));
y = tbl.fooof_width_pre(has_peak);

figure;
subplot(1,3,1)
plot(x,y,'.','markersize',15);
axis square

title('width')
ylabel('peak width')
xlabel('dosing period - peak period');
[r,p]=corr(x,y,'type','pearson');
legend(['r = ' num2str(r) ', p = ' num2str(p)])

subplot(1,3,2)
y = tbl.fooof_amp_pre(has_peak);
plot(x,y,'.','markersize',15);
axis square

title('amplitude')
ylabel('peak amplitude')
xlabel('dosing period - peak period');
[r,p]=corr(x,y,'type','pearson');
legend(['r = ' num2str(r) ', p = ' num2str(p)])


subplot(1,3,3)
y = tbl.fooof_width_pre(has_peak) ./ tbl.fooof_amp_pre(has_peak);
%[f,gof] = fit(x,y,'exp2')
plot(x,y,'.','markersize',15)%'.','markersize',15);
axis square

title('fatness of peak')
ylabel('peak width/amp')
xlabel('dosing period - peak period');
[r,p]=corr(x,y,'type','pearson');
legend(['r = ' num2str(r) ', p = ' num2str(p)])



% difference in peak period & dosing period vs variation in dosing 
figure;
x = tbl.dosing_period_var(has_peak);
y = abs(tbl.dosing_period(has_peak)-tbl.fooof_period_pre(has_peak));
plot(x,y,'.','markersize',15); hold on;
xlabel('std of dosing intervals','fontsize',16);
ylabel('dosing period - spike spectral peak','fontsize',16)
[rho,p]=corr(x,y);

p1 = polyfit(x,y,1);
x1 = linspace(min(x),max(x));
y1 = polyval(p1,x1);
plot(x1,y1);
legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit')
hold off



% dosing period vs peak frequency detected by fooof
figure;
subplot(1,3,1)

peri_12_inds = tbl.dosing_period >= 0 & tbl.dosing_period <= 16;
non_zero_inds = tbl.fooof_period_pre ~=0;
x = tbl.dosing_period(peri_12_inds & non_zero_inds);
y = tbl.fooof_period_pre(peri_12_inds & non_zero_inds);
scatter(x,y)
xlim([0 20])
ylim([0 20])
axis square;
hold on;
title('>0hrs - 16hrs');
xlabel('dosing period'); ylabel('peak in fooofed model')
% ylim([0 40])
% xlim([0 40])

p1 = polyfit(x,y,1);
x1 = linspace(min(x),max(x));
y1 = polyval(p1,x1);
plot(x1,y1,'linewidth',2);
%plot([0 30],[0 30]);
[rho,p]=corr(x,y);
legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit','y=x')
hold off



subplot(1,3,2)
peri_24_inds = tbl.dosing_period > 16 & tbl.dosing_period <= 35;
x = tbl.dosing_period(peri_24_inds & non_zero_inds); %& non_zero_inds
y = tbl.fooof_period_pre(peri_24_inds& non_zero_inds);

scatter(x,y)
xlim([16 35])
ylim([16 35])
axis square;
hold on; 
title('16hr - 35hrs');
xlabel('dosing period'); ylabel('peak in fooofed model')

p1 = polyfit(x,y,1);
x1 = linspace(min(x),max(x));
y1 = polyval(p1,x1);
plot(x1,y1,'linewidth',2);
%plot([0 30],[0 30]);
[rho,p]=corr(x,y);
legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit','y=x')
hold off


subplot(1,3,3)

% aggregate data at a patient level- range of dosing periods
ptids= unique(tbl.ptid);
dosing_range = zeros(1,length(ptids));
has_any_peak = zeros(1,length(ptids));
has_any_dosing_period =  zeros(1,length(ptids));
has_no_dosing_period = zeros(1,length(ptids));
for i = 1:length(ptids)
    pt_inds = tbl.ptid == ptids(i);
    dosing_range(i) = range(tbl.dosing_period(pt_inds));
    has_any_peak(i) =( (sum(tbl.fooof_period_pre(pt_inds) ~= 0)./ height(tbl.fooof_period_pre(pt_inds))))*100;
    has_no_dosing_period(i) = all(tbl.dosing_period(pt_inds) == 0 | isnan(tbl.dosing_period(pt_inds)) );
end

% 
% x = dosing_range';
% y = has_any_peak';
% scatter(x,y); hold on;
% xlabel('range of dosing periods (hrs)');
% ylabel('% of dosing periods w/ peaks')
% title('interaction between dosing periods')
% axis square;
% 
% p1 = polyfit(x,y,1);
% x1 = linspace(0,max(x));
% y1 = polyval(p1,x1);
% plot(x1,y1);
% [rho,p]=corr(x,y);
% legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit')
% hold off

x = tbl.dosing_period; %& non_zero_inds
y = tbl.fooof_period_pre;

scatter(x,y)
axis square;
hold on; 
title('all dosing periods');
xlabel('dosing period'); ylabel('peak in fooofed model')

p1 = polyfit(x,y,1);
x1 = linspace(min(x),max(x));
y1 = polyval(p1,x1);
plot(x1,y1,'linewidth',2);
%plot([0 30],[0 30]);
[rho,p]=corr(x,y);
legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit','y=x')
hold off



%% look at patients that have a 24hr peak

% how many patients with a 24hr peak also have a 12 hr and which of those
% were also dosed q12
err = 0.3; % setting error threshold low 

has_24hr_peak = false(1,height(pt_data_clips));
has_12hr_peak = false(1,height(pt_data_clips));

for i = 1:height(pt_data_clips)
    try
    spec_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i}(:,1);
    has_12hr_peak(i) = any(abs(spec_peaks -12) <= err*12);
    has_24hr_peak(i) = any(abs(spec_peaks -24) <= err*24);
    catch
    end
end

dosed_q12 = false(1,height(pt_data_clips));
only_dosed_q12 = false(1,height(pt_data_clips));
dosed_q24= false(1,height(pt_data_clips));

for i = 1:length(ptIDs)
    pt_inds = tbl.ptid == ptIDs(i);
    dosing_periods = tbl.dosing_period(pt_inds);

    dosed_q24(i) = any(abs(dosing_periods -24) <= err*24);
    dosed_q12(i) = any(abs(dosing_periods -12) <= err*12);
    only_dosed_q12(i) = all(abs(dosing_periods -12) <= err*12);

end

%% perform fishers tests 
[conttbl,chi2,p] = crosstab(has_24hr_peak,has_12hr_peak);
peak_table = table;
peak_table.has_24hr_peak = has_24hr_peak'; peak_table.has_12hr_peak = has_12hr_peak';

[h,p,stats] = fishertest(conttbl);
heatmap(peak_table,'has_24hr_peak','has_12hr_peak');

disp(['results for association between a 24hr peak and a 12hr peak, p = ' num2str(p)])

[conttbl,chi2,p] = crosstab(dosed_q12,has_12hr_peak);
[h,p,stats] = fishertest(conttbl);
disp(['results for association between a 12hr peak and 12hr dosing, p = ' num2str(p)])

[conttbl,chi2,p] = crosstab(only_dosed_q12,has_12hr_peak);
[h,p,stats] = fishertest(conttbl);
disp(['results for association between a 12hr peak and *only* 12hr dosing, p = ' num2str(p)])

[conttbl,chi2,p] = crosstab(has_24hr_peak&has_12hr_peak,~dosed_q12);
[h,p,stats] = fishertest(conttbl);
disp(['results for association between presence of both 24hr peak and a 12hr peak with  q12 dosing, p = ' num2str(p)])

[conttbl,chi2,p] = crosstab(has_24hr_peak,dosed_q24);
[h,p,stats] = fishertest(conttbl);
disp(['results for association between 24hr peak and a 24hr dosing, p = ' num2str(p)])

%% some post hoc testing 
pts_12hr_dosing_peak = ptIDs(has_12hr_peak & dosed_q12);
pts_12hr_dosing_nopeak = ptIDs(~has_12hr_peak & dosed_q12);

std_peak = zeros(1,length(pts_12hr_dosing_peak));
for i = 1:length(pts_12hr_dosing_peak)
    pt = pts_12hr_dosing_peak(i);
    std_peak(i) = mean(tbl.dosing_period_var(tbl.ptid == pt & (abs(tbl.dosing_period-12) <=err*12)));

end 

std_nopeak = zeros(1,length(pts_12hr_dosing_nopeak));
for i = 1:length(pts_12hr_dosing_nopeak)
   pt = pts_12hr_dosing_nopeak(i);
   std_nopeak(i) = sum(tbl.dosing_period_var(tbl.ptid == pt & (abs(tbl.dosing_period-12) <=err*12))); % sum of the std for drugs dosed at 12hrs
end 

figure;
subplot(1,2,1)
[p,h,stats]=ranksum(std_nopeak,std_peak)
data = [mean(std_peak) mean(std_nopeak)];
bar(data);

errhigh =mean(std_peak) + std(mean(std_peak)/length(mean(std_peak)));
errlow  = mean(std_peak) - std(mean(std_peak)/length(mean(std_peak)));
hold on

er = errorbar([1,2],data,errlow,errhigh);    

hold off


ylabel('total std of 12hr dosing');
xticklabels({'dosed q12 & has peak','dosed q12 & no peak'})
title('12hr dosing')
axis square;

% do same thing for 24hr peak

pts_24hr_dosing_peak = ptIDs(has_24hr_peak & dosed_q24);
pts_24hr_dosing_nopeak = ptIDs(~has_24hr_peak & dosed_q24);

std_peak = zeros(1,length(pts_24hr_dosing_peak));
for i = 1:length(pts_24hr_dosing_peak)
    pt = pts_24hr_dosing_peak(i);
    std_peak(i) = mean(tbl.dosing_period_var(tbl.ptid == pt & (abs(tbl.dosing_period-24) <=err*24)));

end 

std_nopeak = zeros(1,length(pts_24hr_dosing_nopeak));
for i = 1:length(pts_24hr_dosing_nopeak)
   pt = pts_24hr_dosing_nopeak(i);
   std_nopeak(i) = sum(tbl.dosing_period_var(tbl.ptid == pt & (abs(tbl.dosing_period-24) <=err*24))); % sum of the std for drugs dosed at 12hrs
end 


subplot(1,2,2)
[p,h,stats]=ranksum(std_nopeak,std_peak)
data = [mean(std_peak) mean(std_nopeak)];
bar(data);
title('24hr dosing')
axis square;

errhigh =mean(std_peak) + std(mean(std_peak)/length(mean(std_peak)));
errlow  = mean(std_peak) - std(mean(std_peak)/length(mean(std_peak)));
hold on

er = errorbar([1,2],data,errlow,errhigh);    

hold off


ylabel('total std of 24hr dosing');
xticklabels({'dosed q24 & has peak','dosed q24 & no peak'})









% figure;
% ptids =unique(pt_data_clips.ptID);
% colors = turbo(length(ptids));
% for i = 1:length(ptids)
%     inds = find(pt_data_clips.ptID==ptids(i));
%     jitter = rand(length(inds),1)./5;
%     plot(zeros(length(inds),1)+(i)+jitter-mean(jitter),pt_data_clips.dosing_periods{inds},'.','Color',colors(i,:),'Markersize',25); hold on;
% end






