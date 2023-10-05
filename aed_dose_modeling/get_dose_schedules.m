% inputs: meds, med_names from parse_MAR.m
function [all_dose_schedules,all_dose_times,all_med_names,all_ieeg_offset,emu_dur] = get_dose_schedules(ptIDs)
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% load all med data
load('MAR_032122.mat');
all_med_names = cell(1,length(ptIDs));
all_dose_times = cell(1,length(ptIDs));
all_dose_schedules = cell(1,length(ptIDs));

all_ieeg_offset=cell(3,length(ptIDs));
emu_dur = zeros(1,length(ptIDs));
for ipt = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    disp(ptID);
    [med_names,meds,~] = parse_MAR(ptID,all_meds);
    all_med_names(ipt) = {med_names};
    emu_dur(ipt) = hours(meds.date(end)-meds.date(1)) +24;
    
    pt_dose_curves = cell(1,length(med_names));
    pt_time = cell(1,length(med_names));
    
    
    %get offset in ieeg
    eeg_diff = (meds.admin_time*3600 - meds.OffsetSecondsInIeeg);
    eeg_round = unique(round(eeg_diff*10^5)/10^5);
    eeg_offsets = eeg_round(~isnan(eeg_round));
    
    eeg_datasets = unique(meds.dataset);
    eeg_datasets = eeg_datasets(contains(eeg_datasets,'D'));
    all_ieeg_offset(1,ipt) = {eeg_datasets};
    all_ieeg_offset(2,ipt) = {eeg_offsets};
    
    datasetStarts = unique(meds.DatasetStartIeeg);
    datasetStarts = datasetStarts(~isnat(unique(meds.DatasetStartIeeg)));
    all_ieeg_offset(3,ipt) = {datasetStarts};
    for n = 1:length(med_names)
        
        %find the specific medication
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        
        sched = hours(diff(this_med_info.date));
        drug_doses = this_med_info.dose;
        drug_times = this_med_info.admin_time;
        
        skipped_days = find(sched > 24);
        for i =1:length(skipped_days)
            drug_doses = [drug_doses(1:skipped_days(i)); 0 ; drug_doses(skipped_days(i)+1:end)];
            drug_times = [drug_times(1:skipped_days(i)); mean([drug_times(skipped_days(i)) drug_times(skipped_days(i)+1)]) ; drug_times(skipped_days(i)+1:end)];
        end
        
        % interpolate drug_doses and drug_times to have one minute bins for the whole EMU stay
        time = linspace(1,emu_dur(ipt),emu_dur(ipt)*60);
        interp_doses = nan(1,length(time));
        
        dose_inds = zeros(1,length(drug_doses));
        for i = 1:length(drug_times)
            [~,dose_inds(i)] = min(abs(time-drug_times(i)));
        end
        interp_doses(dose_inds) = drug_doses;
        
        % interpolate at the zeros, fill in with previous nonzero value
        for i = 1:length(dose_inds)-1
            val = drug_doses(i);
            ind1 = dose_inds(i);
            ind2 = dose_inds(i+1);
            interp_doses(ind1:ind2) = val;
        end
        
        pt_dose_curves(n) = {interp_doses};
        pt_time(n) = {time};
        
    end
    all_dose_schedules{ipt}=pt_dose_curves;
    all_dose_times{ipt}=pt_time;
    
    
end
end



