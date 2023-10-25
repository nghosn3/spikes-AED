%% spectral analysis of spike rate using wavelet transform and FOOOF

close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
%addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])
addpath('/Volumes/users/nghosn3/Pioneer/DATA/mat files')
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


%%

% set constants
fs = 6 ; % 6 sample every hour - units are in hours
prom_thresh = 0.2; % for grabbing peaks on periodogram

% For inverse Wavelet Transform: Setting Constants from Torrence and Compo (TBL 2)
dt = 1/fs;
Cd = 0.776;
Psi0 = pi^(-.25);
dj = 1/12;
modifier = (dj*dt^(1/2))/(Cd*Psi0);
% Defining inverse wavelet function
invcwt = @ (wave,scale,pk_locs) modifier*sum(real(wave(pk_locs,:))./(scale(pk_locs)'.^(1/2)),1);


tiledlayout('flow');
for ipt =79%1:length(ptIDs)%[67,79]

    ptID = ['HUP' num2str(ptIDs(ipt))];

    med_names = all_med_names{ipt};
    med_plot_inds = ~(contains(med_names,'lorazepam'));

    asm_load = all_pts_drug_samp{ipt};
    spikes = (all_spike_rate{ipt}+1);
    spikes = spikes./max(spikes);

    signal = zscore(spikes);
    exact_time = all_spike_times{ipt};
    time =linspace(0,max(exact_time),length(signal))./3600;

    asm_before = asm_load;
    spikes_before = signal;

%     % plot spikes
%     nexttile;
%     plot(time,smoothdata(signal),'linewidth',1.5)
%     ylabel('Spike rate (#/10min)')
%     title('Spike rate','fontsize',14)
%     xlabel('Time(hours)','fontsize',14)
%     set(gca,'fontsize',14)
%     set(gca, 'box', 'off')
% 
%     % plot ASM load over time 
%     nexttile();
%     plot(time, asm_load(med_plot_inds,:)','linewidth',1.5); hold on;
%     legend(med_names(med_plot_inds));title(ptID,'fontsize',14)
%     title('ASM load','fontsize',14);
%     ylabel('normalized BPL','fontsize',14);
%     xlabel('time (hrs)','fontsize',14);
%     set(gca,'fontsize',14)
%     set(gca, 'box', 'off')
%     

      % plot spike spectrum
    nexttile();
    [wave,period,scale,coi] = wt([time',signal'],hours(1/fs));
    pre_periodogram = mean(abs(wave(:,:)).^2,2);
    plot(period,pre_periodogram./max(pre_periodogram),'linewidth',2,'color','black');
    title('Periodogram: Spike rate');
    set(gca,'fontsize',14)
    ylabel('prominance','fontsize',14);
    xlabel('period (hrs)','fontsize',14)
    set(gca, 'box', 'off')

  
    % plot ASM spectrum : calc and plot periodograms of ASMs
    hold on;
    for m =  1:length(med_names)
        sig = asm_before(m,:)';
        if  ~contains(med_names(m),'lorazepam') % at least 3 administrations in the dosing period
            [wave,periods,~,~] = wt([time',sig],hours(1/fs));
            periodogram = mean(abs(wave).^2,2);
            plot(periods,periodogram./max(periodogram),'linewidth',1.5);hold on;
        end
    end

    
    title('ASM load and spike rate','fontsize',14);
    ylabel('Normalized prominance','fontsize',14);
    xlabel('period (hrs)','fontsize',14)
    set(gca,'fontsize',14)
    set(gca, 'box', 'off')
    xline(12,'--','linewidth',2,'Color',[.5 .5 .5])
    xline(24,'--','linewidth',2,'Color',[.7 .7 .7])
    legend([{'spike rate'} ; med_names(med_plot_inds)]); hold off;
    xticks([0 12 24 48]);



end


%% compare the peakiness of the spike peak at dosing period to distance from
% 12hr dosing during pre period

% load the fooof results and
load('fooof_results_091223.mat');

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
         
            next=next+1;
        end

    end

end

% remove entries that correspons to ativan
ativan_inds = contains(tbl.med_name,'lorazepam');
tbl(ativan_inds,:) = [];


%% look at patients that have a 24hr peak

err = 0.3; % setting error threshold lows

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

%% figure showing these peaks are likely not related ...
peak_table = table();
peak_table.has_24hr_peak = has_24hr_peak'; 
peak_table.has_12hr_peak = has_12hr_peak';
peak_table.dosed_q12 = dosed_q12';

% nexttile()
% [conttbl,~,~] = crosstab(has_24hr_peak,has_12hr_peak);
% [~,p_24,~] = fishertest(conttbl);
% heatmap(peak_table,'has_24hr_peak','has_12hr_peak');
% title('')
% ylabel('12hr spike peak')
% xlabel('24hr spike peak')
% set(gca,'fontsize',14)

nexttile();
[conttbl,~,~] = crosstab(dosed_q12,has_12hr_peak);
[~,p_12,~] = fishertest(conttbl);
heatmap(peak_table,'has_12hr_peak','dosed_q12');
ylabel('12hr ASM dosing')
xlabel('12hr spike peak')
set(gca,'fontsize',14)
title('12hr ASM and spike rate peaks')





%% now look at plv of spikes and ASMs
[all_plv,all_asm_plv,all_spikes_plv] = calc_plv_spikes_asm(all_spike_rate,all_pts_drug_samp,file_inds,ptIDs);


%%
r_plv = nan(length(ptIDs),1);

%tiledlayout('flow')
for i = 1:length(ptIDs)
    plv = all_plv{i};
    nan_inds =isnan(plv);
    plv(nan_inds)=[];
    asm_load = all_asm_plv{i};
    asm_load(nan_inds)=[];

    r = corrcoef(asm_load,plv);
    r_plv(i) = r(1,2);
    %     nexttile;
    %     plot(plv,asm_load,'.')


end

% figure;
% tiledlayout('flow');
nexttile;
histogram(r_plv, 'FaceAlpha', 0.5,'Facecolor','black','EdgeColor', 'black', 'linewidth', 2);
set(gca,'fontsize',14);
ylabel('# patients')
xlabel('correlation');
title('PLV and ASM load')

[~,p] = ttest(r_plv);




