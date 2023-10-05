close all;clear;
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');

% Get which patients
cohort_info = readtable('AED-connectivity-cohort-list.xlsx');
ptIDs = cohort_info.PtIDs;
% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

all_dose_curves =cell(1,length(ptIDs));
all_thr_curves =cell(1,length(ptIDs));

for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID =  ptIDs{ipt};
    [med_times,med_names,days,each_drug_dose,total_dose,alc_times,meds] = get_all_med_info(ptID);
    
    % will want to do this for each drug when modeling all of them... wont need to select which med, just loop through all
    c_all=cell(1,length(med_names));
    figure(ipt);
    for n = 1:length(med_names)
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        %placeholder for adding the last timepoint in the emu stay, to get
        %the concentration to the end.
        this_med_info.start(end+1) = 400*3600;
        this_med_schedule = diff(this_med_info.start)./3600; %in hours
        
        % parameters that will be drawn from the spreadsheet: example parameters of levetiracetam
        med_ind =contains(aed_params.medication,med_names(n));
        tmax = aed_params(med_ind,:).t_max_c;
        tHalf = aed_params(med_ind,:).t_half_e;
        F = aed_params(med_ind,:).F;
        
        %set initial dose to zero, or to presurgical levels (plus decay since measured)
        %c_total = zeros(1,60*300);
        c_drug=[];
        dose = this_med_info.dose_mg(1);
        for i = 2:height(this_med_info)
            tInt = this_med_schedule(i-1);
            
            [c_t,t]= get_single_dose_curve(dose,tHalf,tInt,F);
            dose = this_med_info.dose_mg(i) + c_t(end);
            c_drug=[c_drug c_t];
        end
        
        c_all(ipt) = {c_drug};
        all_dose_curves(ipt) ={c_all};
        %% get the curve- preset the total curve by 60*tHr(of emu stay)
        
        tHr = linspace(this_med_info.start(1)./3600,this_med_info.start(end)./3600,length(c_drug));
        all_tHr(ipt) = {tHr};
        subplot(2,1,1)
        plot(tHr,c_drug,'linewidth',1.4); hold on;
        xlim([0 tHr(end)]); % this should be the range of the spike data, when that is incorperated
        title(['approx. plasma concentration of each AED: ' ptIDs{ipt}]);xlabel('hours since EMU admission')
        legend(med_names);
        
        subplot(2,1,2)
        c_drug_norm = c_drug ./ max(c_drug);
        plot(tHr,c_drug_norm,'linewidth',1.4); hold on;
        xlim([0 tHr(end)]); % this should be the range of the spike data, when that is incorperated
        title(['Normalized approx. plasma concentration of each AED: ' ptIDs{ipt}]);xlabel('hours since EMU admission')
        
    end
    %get seizures and plot them on normalized curve?
    seizure_times = get_seizure_times(ptID);
    for j =1:length(seizure_times)
        xline(seizure_times(j,1),'--r','linewidth',1);hold on;
    end
end
