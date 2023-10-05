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

    % use drug sum
    asm_curves = all_pts_drug_samp{ipt};
    med_curve = sum(asm_curves);
    med_names = all_med_names{ipt};

    %find when/if all drugs are stopped
    [pks,locs] = findpeaks(med_curve);
    if length(locs)>2 && (locs(1) <= taper_thresh)
        diffs = diff(locs);

        last_admin = find(diffs >= taper_thresh);
        if ~isempty(last_admin)
            stop_med =locs(last_admin(1)); %ind+1 of last administration
        else
            stop_med = locs(end)+12;
            stop_med = min([stop_med length(med_curve)]);
        end

    end

    pt_data_clips.inds(ipt) = {[locs(1) stop_med]};

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

tic
for ipt =1:length(ptIDs)

    ptID = ['HUP' num2str(ptIDs(ipt))];

    med_names = all_med_names{ipt};
    asm_load = all_pts_drug_samp{ipt};
    asm_load = sum(asm_load);
    spikes = (all_spike_rate{ipt}+1);
    spikes = spikes./max(spikes);

    exact_time = all_spike_times{ipt};
    time =linspace(0,max(exact_time),length(spikes))./3600;


    % all medication tapered for >2days - pre/post condition breakpoints
    all_asm_inds = pt_data_clips.inds{ipt};
    dosing_periods_wavelet = [];
    dosing_periods =[];

    dosing_periods_var = [];
    dosing_periods_amp =[];

    ind1 = all_asm_inds(1); % min ; when first med is started
    ind2 = all_asm_inds(2);

    exact_time = all_spike_times{ipt};
    time =linspace(0,max(exact_time),length(spikes))./3600;

    asm_before = asm_load(:,ind1:ind2);

    spikes_before = spikes(:,ind1:ind2);

    % determine dosing freqs in before period using wavelet transform


    sig = asm_before;
    [pks,locs] = findpeaks(sig);
    if length(diff(locs)) >= 2 % at least 3 administrations in the dosing period
        target_dose_int = mean(diff(time(locs))); % in hours from med curve
        dosing_periods = target_dose_int;
        dosing_periods_var = std(diff(time(locs))); % standard deviation
            [wave,periods,~,~] = wt([time(ind1:ind2)',sig'],hours(1/fs));
            periodogram = mean(abs(wave).^2,2);
            [pks,locs] = findpeaks(periodogram./max(periodogram)); %,'minpeakprominence',.3
            if ~isempty(locs) %&& round(max(periods))>=20
                [per_diff,loc_int] = min(abs(periods(locs)-target_dose_int));
                if (per_diff < target_dose_int*.3) % verify that it is a peak period in the signal - 33% error
                    pk_periods = (periods(locs));
                    dosing_freqs = (pk_periods(loc_int));
                    dosing_periods_wavelet = round(dosing_freqs*10)./10;
                    [~,ind]=min(abs(periods-dosing_freqs));
                    dosing_periods_amp = periodogram(ind);
                else
                    dosing_periods_wavelet = NaN;
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
        xline(time(ind2),'r-','linewidth',3); %xline(time(ind2),'r-','linewidth',3)
        xline(time(ind1),'r-','linewidth',3); %xline(time(ind2),'r-','linewidth',3)
        legend(med_names);title(ptID)
    end

    % amplitude of periodogram peak pre/post - normalized
    % across subjects
    [wave,period,scale,coi] = wt([time',spikes'],hours(1/fs));
    pre_periodogram = mean(abs(wave(:,ind1:ind2)).^2,2);

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
    dosing_periods = pt_data_clips.dosing_periods{i};
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

has_peak = tbl.fooof_period_pre == tbl.fooof_period_pre;
no_peak = ~has_peak;


%x = abs(tbl.dosing_period(has_peak)-tbl.fooof_period_pre(has_peak));
x = tbl.dosing_period_var(has_peak );
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
xlabel('std of dosing intervals');
ylabel('dosing period - spike spectral peak')
[rho,p]=corr(x,y);

p1 = polyfit(x,y,1);
x1 = linspace(min(x),max(x));
y1 = polyval(p1,x1);
plot(x1,y1);
legend(['R = ' num2str(rho) ', p = ' num2str(p)],'line of best fit')
hold off



% dosing period vs peak frequency detected by fooof
figure;

peri_12_inds = tbl.dosing_period >= 8 & tbl.dosing_period <= 16;
non_zero_inds = tbl.fooof_period_pre ~=0;
x = tbl.dosing_period(peri_12_inds & non_zero_inds);
y = tbl.fooof_period_pre(peri_12_inds & non_zero_inds);
scatter(x,y)
axis square;
hold on;
title('dosing period vs fooof peak period');
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


%% look at patients that have a 24hr peak

% how many patients with a 24hr peak also have a 12 hr and which of those
% were also dosed q12
err = 0.3; % setting error threshold low

has_24hr_peak = false(1,height(pt_data_clips));
has_12hr_peak = false(1,height(pt_data_clips));

dosed_q12 = false(1,height(pt_data_clips));
dosed_q24= false(1,height(pt_data_clips));

for i = 1:height(pt_data_clips)
    try
        spec_peaks = pt_data_clips.spike_fooof_peakfreqs_pre{i}(:,1);
        dosing_periods = pt_data_clips.dosing_periods{i};

        dosed_q12(i) = any(abs(dosing_periods -12) <= err*12);
        has_12hr_peak(i) = any(abs(spec_peaks -12) <= err*12);

        dosed_q24(i) = any(abs(dosing_periods -24) <= err*24);
        has_24hr_peak(i) = any(abs(spec_peaks -24) <= err*24);
    catch
    end

end

%% perform fishers tests
[conttbl,chi2,p] = crosstab(has_24hr_peak,has_12hr_peak);
peak_table = table;
peak_table.has_24hr_peak = has_24hr_peak'; peak_table.has_12hr_peak = has_12hr_peak';

[h,p,stats] = fishertest(conttbl);
figure;
heatmap(peak_table,'has_24hr_peak','has_12hr_peak');

disp(['results for association between a 24hr peak and a 12hr peak, p = ' num2str(p)])

[conttbl,chi2,p] = crosstab(dosed_q12,has_12hr_peak);
[h,p,stats] = fishertest(conttbl);
disp(['results for association between a 12hr peak and 12hr dosing, p = ' num2str(p)])
figure;
peak_table = table;
peak_table.dosed_q12 = dosed_q12'; peak_table.has_12hr_peak = has_12hr_peak';
heatmap(peak_table,'dosed_q12','has_12hr_peak');


%% some post hoc testing
pts_12hr_dosing_peak = (has_12hr_peak & dosed_q12);
pts_12hr_dosing_nopeak = (~has_12hr_peak & dosed_q12);

std_peak = pt_data_clips.dosing_periods_var(pts_12hr_dosing_peak); std_peak = [std_peak{:}];
std_nopeak = pt_data_clips.dosing_periods_var(pts_12hr_dosing_nopeak);std_nopeak = [std_nopeak{:}]


figure;
[p,h,stats]=ranksum(std_nopeak,std_peak)
data = [mean(std_peak) mean(std_nopeak)];
bar(data);

errhigh =mean(std_peak) + std(mean(std_peak)/length(mean(std_peak)));
errlow  = mean(std_peak) - std(mean(std_peak)/length(mean(std_peak)));
hold on

er = errorbar([1,2],data,errlow,errhigh);

hold off


ylabel('std of 12hr dosing');
xticklabels({'dosed q12 & has peak','dosed q12 & no peak'})
title('12hr dosing')
axis square;



% figure;
% ptids =unique(pt_data_clips.ptID);
% colors = turbo(length(ptids));
% for i = 1:length(ptids)
%     inds = find(pt_data_clips.ptID==ptids(i));
%     jitter = rand(length(inds),1)./5;
%     plot(zeros(length(inds),1)+(i)+jitter-mean(jitter),pt_data_clips.dosing_periods{inds},'.','Color',colors(i,:),'Markersize',25); hold on;
% end






