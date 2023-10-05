%% dose as a function of pt home dose and minimum therapuetic dose at seizure onset
close all;clear;
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling/figures_code')

% load medication data
load('MAR_032122.mat')
load('home_meds.mat');

% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;

% Get the drug curve
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_min_dose(ptIDs);

% load AED meta data for min effective dose
% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

%%
min_dose_sz = cell(1,length(ptIDs));
home_dose_sz = cell(1,length(ptIDs));

total_home_dose = cell(1,length(ptIDs));
total_min_dose = cell(1,length(ptIDs));
for i = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    [med_names,meds,explant_date] = parse_MAR(ptID,all_meds);
    
    % get dose per day in the EMU 
    [out] = get_dose_per_day(meds, med_names);

    % grab home meds for pt
    pt_home_meds = home_meds(home_meds.HUP_ID==str2double(ptID(4:end)),:);
    % only for meds within 3days of explant date
    pt_home_meds= pt_home_meds(pt_home_meds.EVENT_TIME>= explant_date & pt_home_meds.EVENT_TIME<=explant_date+days(4),:);
    [~,ia] = unique(pt_home_meds.medication);
    pt_home_meds = pt_home_meds(ia,:);
    
    % find the dose of the home meds from the EMU meds if NaN
    %add meds given in first day ofstay
    missed_meds=cell(1,length(med_names));
    ind=1;
    if ~isempty(pt_home_meds)
        for j =1:length(med_names)
            if ~contains(pt_home_meds.medication,med_names(j))
                missed_meds{ind} =med_names(j);
            end
        end
    else
        pt_home_meds(1:length(med_names),'medication') = med_names;
    end
    
    % if dose for med not included in data
    home_med_names = unique(pt_home_meds.medication);
    for m=1:length(home_med_names)
        if isnan(pt_home_meds.DOSE(m)) || pt_home_meds.DOSE(m)==0
            temp_med_doses = [meds.dose(contains(meds.medication,home_med_names(m))); 0];
            pt_home_meds.dose(m)= temp_med_doses(1);
        else 
            pt_home_meds.dose(contains(pt_home_meds.medication,home_med_names(m)))= pt_home_meds.DOSE(m);
        end 
        % get medication frequency
        if ~isempty(pt_home_meds.SIG{m})
            [~, spd] = parse_med_frequencies(pt_home_meds.SIG(m));
            pt_home_meds.frequency(m) = spd;
        else 
            pt_home_meds.frequency(m)=NaN;
        end
    end 
    
    % calculate the dose per day as a function of min dose and home med dose
    emu_days = max(cellfun(@length,out.days));
    pt_total_min_dose =  zeros(1,emu_days);
    pt_total_home_dose = zeros(1,emu_days);
    
    for n = 1:length(out.med_name)
        if ~isempty(out.med_name{n})
            med_ind = contains(aed_params.medication,out.med_name{n});
            min_dose = aed_params.min_dose_single_mg(med_ind);
            pt_total_min_dose = nansum([pt_total_min_dose;(out.dose_schedule{n}./min_dose)],1);
             c
            med_ind = contains(pt_home_meds.medication,out.med_name{n});
            if sum(med_ind)>0
                home_dose = pt_home_meds.frequency(med_ind)*pt_home_meds.dose(med_ind);
                if ~isnan(home_dose)
                    pt_total_home_dose = pt_total_home_dose + (out.dose_schedule{n}./home_dose);
                end
            end
        end
    end
    
    total_home_dose{i} =pt_total_home_dose;
    total_min_dose{i} =pt_total_min_dose;
    
    % get the dose load at seizure times for that patient
    
    % get seizure times in EMU time 
    offsets = ieeg_offset{2,i};
    ieeg_offset_datasets = ieeg_offset{1,i};
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    sz_inds = floor(seizure_times./24); sz_inds(sz_inds(:,1)==0)=1;
    
    min_dose_sz{i} = pt_total_min_dose(sz_inds(:,1));
    home_dose_sz{i} = pt_total_home_dose(sz_inds(:,1));
    
end 

%% make plots
subplot(1,2,1)
histogram(horzcat(min_dose_sz{:})); title('fraction of minimum therapuetic dose at seizure'); axis square
hold on; histogram(horzcat(total_min_dose{:})); legend({'daily dose at seizure','all daily doses'})

subplot(1,2,2)
% remove outliers from med parsing
home_doses = horzcat(home_dose_sz{:});
home_doses(home_doses>100)=NaN;
histogram(home_doses); title('fraction of prescribed home dose at seizure'); axis square;
total_home_doses = horzcat(total_home_dose{:}); 
total_home_doses(total_home_doses>100)=NaN;
hold on; histogram(total_home_doses); legend({'daily dose at seizure','all daily doses'})

figure;
[out,pval_binom_alt,successes,successes_alt] = plot_orders(total_min_dose,min_dose_sz);


