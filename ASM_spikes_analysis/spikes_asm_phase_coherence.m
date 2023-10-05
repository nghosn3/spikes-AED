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
load('pt_data_clips.mat')

[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes('MAR_032122.mat','spikes_rates_021323.mat',all_ieeg_offset, all_dose_curves, all_tHr,ptIDs);

%% find overall dosing_period

for ipt =  1:length(ptIDs)

    ptID = ['HUP' num2str(ptIDs(ipt))];
    asm_load = all_pts_drug_samp{ipt};
    fs = 6 ; % 6 sample every hour - units are in hours


    % at least one medication tapered for >2days - pre/post condition breakpoints
    all_asm_inds = pt_data_clips.inds(ipt); all_asm_inds = all_asm_inds{:};

    dosing_period_overall = [];
    dosing_periods_var_overall = [];


    if ~all(all(isnan(all_asm_inds)))
        ind1 = length(asm_load);%max(all_asm_inds(:,2)); % min ; when last med is stopped


        exact_time = all_spike_times{ipt};
        time =linspace(0,max(exact_time),length(asm_load))./3600;

        asm_before = asm_load(:,1:ind1);
        sig = asm_before(1,:)';

        [pks,locs] = findpeaks(sig);
        if length(diff(locs)) >= 2 % at least 3 administrations in the dosing period
            target_dose_int = median(diff(time(locs))); % in hours from med curve
            dosing_periods_var_overall = nanstd(diff(time(locs))); % standard deviation
            if  ~isempty(pks) && ~isnan(target_dose_int)
                [wave,periods,~,~] = wt([time(1:ind1)',sig],hours(1/fs));
                periodogram = mean(abs(wave).^2,2);
                [pks,locs] = findpeaks(periodogram./max(periodogram)); %,'minpeakprominence',.3
                if ~isempty(locs) %&& round(max(periods))>=20
                    [per_diff,loc_int] = min(abs(periods(locs)-target_dose_int));
                    if (per_diff < target_dose_int*.2) % verify that it is a peak period in the signal - 33% error
                        pk_periods = (periods(locs));
                        dosing_freqs = (pk_periods(loc_int));
                        dosing_periods_overall = round(dosing_freqs*10)./10;

                    else
                        dosing_periods_overall= NaN;
                    end
                end
            end
            pt_data_clips.dosing_period_overall(ipt) = dosing_periods_overall;
            pt_data_clips.dosing_period_overall_var(ipt) = dosing_periods_var_overall;
        end
    end



end

%%
close all;
plot_figs = 1;

xcorr_results = nan(length(ptIDs),2);
xcorr_results_12_24hr = nan(length(ptIDs),2); 
for ipt =  1:length(ptIDs)
    fs = 6 ; % 6 sample every hour - units are in hours
    exact_time = all_spike_times{ipt};
    signal = zscore(all_spike_rate{ipt});
    asm_load = sum(all_pts_drug_samp{ipt});
    time =linspace(0,max(exact_time),length(signal))./3600;

    % pre and post condition breakpoints
    all_asm_inds = pt_data_clips.inds(ipt); all_asm_inds = all_asm_inds{:};
    ind1 = max(all_asm_inds(:,2));
    ind1 = length(asm_load);

    ind2 = min(all_asm_inds(:,3));

    if ~isnan(ind1)
        time = time(1:ind1);
        signal = signal(1:ind1);

        % Calculate wavelet energy
        [wave,period,scale,coi] = wt([time',signal'],hours(1/fs));
        prom_thresh = 0.5;
        all_periodogram = mean(abs(wave).^2,2);

        % Inverse Wavelet Transform
        % Setting Constants from Torrence and Compo (TBL 2)
        dt = 1/fs;
        Cd = 0.776;
        Psi0 = pi^(-.25);
        dj = 1/12;
        modifier = (dj*dt^(1/2))/(Cd*Psi0);

        % Defining inverse wavelet function
        invcwt = @ (wave,scale,pk_locs) modifier*sum(real(wave(pk_locs,:))./(scale(pk_locs)'.^(1/2)),1);

        %Finding 5 sec peak period (with +/- 33% per Baud et. al., 2018)
        pk_per = pt_data_clips.dosing_period_overall(ipt);% pt_data_clips.dosing_period(ipt); % average dosing
        if isnan(pk_per) | pk_per == 0
            pk_per = 12;
        end 

        pk_per_12 = 12;
        pk_per_24 = 24;

        % You can also use the pre_pks found above by filtering for the specific peak of
        % interest and instead of 5, pk_per = period(pk_loc) where pk_loc is the idx of peak of interest
        err = 0.05; % bounds around frequency of interest
        per_bounds = [pk_per-err*pk_per,pk_per+err*pk_per];
        per_mask = period > per_bounds(1) & period < per_bounds(2);

        % Plotting inverse wavelet and signal (scales changed for comparison)
        inverse_wavelet = invcwt(wave,scale,per_mask);
        invw_phase = angle(hilbert(inverse_wavelet));


        % phase of ASM load
        asm_signal = asm_load(1:ind1);

        [wave,period,scale,coi] = wt([time(1:ind1)',asm_signal'],hours(1/fs));
        inverse_wavelet_asm = invcwt(wave,scale,per_mask);
        invw_phase_asm = angle(hilbert(inverse_wavelet_asm));

        max_lags = 24*6; % 12hours - max lag is one cycle - 6hrs is half a cycle
        [xr,lags] = xcorr(invw_phase,invw_phase_asm,max_lags,'normalized');
        [max_xr,max_ind] = max(abs(xr));
        max_xr_lags  = [(lags(max_ind)./6) max_xr]; % lags in hours

        xcorr_results(ipt,:) = max_xr_lags;

        % 12hr wavelet and 24hr wavelet
        per_bounds = [pk_per_12-err*pk_per_12,pk_per_12+err*pk_per_12];
        per_mask = period > per_bounds(1) & period < per_bounds(2);

        % Plotting inverse wavelet and signal (scales changed for comparison)
        inverse_wavelet_12hr = invcwt(wave,scale,per_mask);
        
        %24
        per_bounds = [pk_per_24-err*pk_per_24,pk_per_24+err*pk_per_24];
        per_mask = period > per_bounds(1) & period < per_bounds(2);
        inverse_wavelet_24hr = invcwt(wave,scale,per_mask);
        
        [xr,lags] = xcorr(inverse_wavelet_12hr,inverse_wavelet_24hr,max_lags,'normalized');
        [max_xr,max_ind] = max(abs(xr));
        max_xr_lags  = [(lags(max_ind)./6) max_xr]; % lags in hours

        xcorr_results_12_24hr(ipt,:) = max_xr_lags;


%         figure;
%         plot(time(1:ind1),inverse_wavelet_12hr,'LineWidth',1.5,'LineStyle','-');hold on;
%         plot(time(1:ind1),inverse_wavelet_24hr,'LineWidth',1)
%         legend(['12hr ICWT spikes'],[ '24hr ICWT of spikes'])



        if plot_figs
            figure;
            subplot(2,2,1)
            plot(period,all_periodogram); hold on;
            xlabel('Period (hrs)')
            title('periodogram: pre-taper')
            % Find the peaks (filter noise out with prominence)
            [pks,locs] = findpeaks(all_periodogram);%,"MinPeakProminence",prom_thresh);
            % if using different signal, will have to tune this prominence val
            stem(period(locs),pks)
            xticks(round(period(locs)*10)./10)

            subplot(2,2,2)

            plot(time(1:ind1),asm_signal); hold on;
            plot(time(1:ind1),(inverse_wavelet_asm+1) *2,'LineWidth',1.5,'LineStyle','-');hold on;
            legend('asm load','wavelet transform of asm load')
            title(['ASM load and ' num2str(pk_per) 'hr wavelet'])

            subplot(2,2,3)
            plot(time,signal,'LineWidth',.5);hold on;
            plot(time,inverse_wavelet*5,'LineWidth',2,'LineStyle','-')
            xlabel('Time(hours)')
            title(['HUP ' num2str(ptIDs(ipt)) ': wavelet transform of spikes'])
            ylabel('normalized amplitude')
            legend('Raw Signal', ['' num2str(pk_per) ' ICWT of spikes'])


            subplot(2,2,4)
            plot(time(1:ind1),inverse_wavelet*50,'LineWidth',1.5,'LineStyle','-');hold on;
            plot(time(1:ind1),invw_phase,'LineWidth',1)
            legend(['' num2str(pk_per) ' ICWT of spikes'],[ num2str(pk_per) ' phase of spikes'])
            xlabel('Time (hrs)')
            ylabel('normalized amplitude')
            title('wavelet transform of spikes and phase')



            figure;
            subplot(3,1,1)

            plot(time(1:ind1),inverse_wavelet_asm,'LineWidth',1.5,'LineStyle','-');hold on;
            plot(time(1:ind1),inverse_wavelet,'LineWidth',1)
            legend(['' num2str(pk_per) ' ICWT of asm load'],[ num2str(pk_per) 'ICWT of spikes'])
            
            xlabel('Time (hours)')
            ylabel('normalized amplitude')
            title(['HUP ' num2str(ptIDs(ipt)) ': wavelet transform of spikes and ASMs'])

            subplot(3,1,2)
            plot(time(1:ind1),invw_phase_asm,'LineWidth',1.5,'LineStyle','-');hold on;
            plot(time(1:ind1),invw_phase,'LineWidth',1.5,'LineStyle','-');hold on;
            legend(["phase of asm load","phase of spikes"])
            xlabel('Time (hours)')
            title(['HUP ' num2str(ptIDs(ipt)) ':phase of spikes and ASMs'])


            subplot(3,1,3)
            stem(lags./6,xr);
            xlabel('lags (hours)')
            title(['HUP ' num2str(ptIDs(ipt)) ': lagged cross correlation'])


        end


    end

end
%%
figure;
histogram(xcorr_results(:,1))
xlabel('lag (hrs)','fontsize',14)
ylabel('# patients','fontsize',14)
title('lag of maximal correlation','fontsize',14)

figure;
histogram(xcorr_results_12_24hr(:,1),'numbins',15)
xlabel('lag (hrs)','fontsize',14)
ylabel('# patients','fontsize',14)
title('lag of maximal correlation between 12hr and 24hr wavelets','fontsize',14)
