function [all_dose_curves,all_tHr,ptIDs,all_med_names] = get_aed_curve(ptIDs)

%turn warning off
id = 'MATLAB:table:RowsAddedExistingVars';
warning('off',id)

addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

all_dose_curves =cell(1,length(ptIDs));
all_tHr =cell(1,length(ptIDs));
all_med_names = cell(1,length(ptIDs));

for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID =  ptIDs{ipt};
    [~,med_names,~,~,~,~,meds] = get_all_med_info(ptID);
    all_med_names{ipt} = med_names;
    % will want to do this for each drug when modeling all of them... wont need to select which med, just loop through all
    c_all=cell(1,length(med_names));
    t_all=cell(1,length(med_names));
    for n = 1:length(med_names)
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        %placeholder for adding the last timepoint in the emu stay, to get
        %the concentration to the end.
        this_med_info.start(end+1) = 400*3600;
        this_med_schedule = diff(this_med_info.start)./3600; %in hours
        
        % grab drug parameters for model:
        med_ind =contains(aed_params.medication,med_names(n));
    
        F = aed_params(med_ind,:).F;
        vd=aed_params(med_ind,:).vd;
        ka=aed_params(med_ind,:).ka;
        tmax = aed_params(med_ind,:).t_max;
        % take mean of range for current model:
        tHalf = [aed_params(med_ind,:).t_half_e];
        tHalf=strsplit(tHalf{1},',');tHalf = cellfun(@str2double,tHalf);
        tHalf = mean(tHalf);

        %set initial dose to zero, or to presurgical levels (plus decay since measured)
        c_drug=[];
        dose = this_med_info.dose_mg(1); 
        c0=0;
        
        for i = 2:height(this_med_info)
            tInt = this_med_schedule(i-1);
            [c_t,~]= get_single_dose_curve(c0,dose,tHalf,tInt,F,vd,tmax,ka);
            c_drug=[c_drug c_t];
            c0 = c_drug(end);
            
            dose = this_med_info.dose_mg(i);
        end
        
        c_all(n) = {c_drug};
        all_dose_curves(ipt) ={c_all};
        
        tHr = linspace(this_med_info.start(1)./3600,this_med_info.start(end)./3600,length(c_drug));
        t_all(n) = {tHr};
        all_tHr(ipt) = {t_all};
        
        
        
    end
end


end