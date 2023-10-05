close all;clear;
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% Get which patients
cohort_info = readtable('AED-connectivity-cohort-list.xlsx');
ptIDs = cohort_info.PtIDs;

% Get their medication BPL curves 
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve(ptIDs);

%% loop through patients: get info on medication taper
all_taper_info=cell(1,length(ptIDs));
for i=1:length(ptIDs)
    ptID = ptIDs{i};
    [~,med_names,~,meds] = get_med_info_master(ptID);
    all_taper_info(i) = {get_taper_info(meds, med_names)};
end 

