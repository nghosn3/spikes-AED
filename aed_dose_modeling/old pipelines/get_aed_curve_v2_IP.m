function [all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_v2_IP(ptIDs)

%turn warning off
id = 'MATLAB:table:RowsAddedExistingVars';
warning('off',id)

addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);
%aed_params.medication=convert_AED_to_generic(aed_params.medication);

% load lab and home med data- data with med names to convert - update spreadsheets when full patient cohort arrives
labs=readtable('outpatient_meds_labs.xlsx','Sheet','labs');
home_meds=readtable('outpatient_meds_labs.xlsx','Sheet','home_meds');
implant_info=readtable('outpatient_meds_labs.xlsx','Sheet','implant_time');

% clean up medication name field in each data sheet:
labs.medication=cellfun(@lower,labs.medication,'UniformOutput',false);
temp_list = cellfun(@strsplit,home_meds.medication,'UniformOutput',false);
for i =1:length(temp_list)
    home_meds.medication(i) = temp_list{i}(1);
end
home_meds.medication = cellfun(@lower,home_meds.medication,'UniformOutput',false);

%convert names to generic
home_meds.medication=convert_AED_to_generic(home_meds.medication);
labs.medication = convert_AED_to_generic(labs.medication);

%initialize dose curves for EMU stay
all_dose_curves =cell(1,length(ptIDs));
all_tHr =cell(1,length(ptIDs));
all_med_names=cell(1,length(ptIDs));

% get the ieeg_offset for plotting and the dataset the offset is for
all_ieeg_offset=cell(3,length(ptIDs));

emu_dur =zeros(1,length(ptIDs));
for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID =  ptIDs{ipt};
    disp(ptID);
    [~,med_names,~,meds] = get_med_info_master(ptID);% will want to do this for each drug when modeling all of them... wont need to select which med, just loop through all
    all_med_names(ipt) = {med_names};
    emu_dur(ipt) = hours(meds.Date(end)-meds.Date(1)) +24; %approx
    
    %get offset in ieeg
    eeg_diff = (meds.admin_time*3600 - meds.OffsetSecondsInIeeg);
    eeg_round = unique(round(eeg_diff*10^5)/10^5);
    eeg_offsets = eeg_round(~isnan(eeg_round));
    
    eeg_datasets = unique(meds.Dataset);
    eeg_datasets = eeg_datasets(contains(eeg_datasets,'D'));
    all_ieeg_offset(1,ipt) = {eeg_datasets};
    all_ieeg_offset(2,ipt) = {eeg_offsets};
    
    datasetStarts = unique(meds.DatasetStart);
    datasetStarts = datasetStarts(~isnat(unique(meds.DatasetStart)));
    all_ieeg_offset(3,ipt) = {datasetStarts};
    
    % initialize BPL curves
    c_all=cell(1,length(med_names));
    t_all=cell(1,length(med_names));
    
    % grab home meds and labs for pt
    pt_home_meds = home_meds(contains(home_meds.HUP,ptID),:);
    pt_labs = labs(contains(labs.HUP,ptID),:);
    
    % get implant date and time to use for medication times
    implant_time = implant_info.time(implant_info.HUP_ID==str2double(ptID(4:end)));
    implant_date = implant_info.date(implant_info.HUP_ID==str2double(ptID(4:end)));
    
    for n = 1:length(med_names)
        disp(med_names{n});
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        %labs or home dosing for this med:
        med_labs =pt_labs(contains(pt_labs.medication,med_names(med_name_ind)),:);
        med_routine = pt_home_meds(contains(pt_home_meds.medication,med_names(med_name_ind)),:);
        
        %placeholder for adding the last timepoint in the emu stay, to get
        %the concentration to the end....using 450 can change later
        max_dur =450;
        this_med_info.admin_time(end+1) = max_dur; %hrs
        this_med_schedule = diff(this_med_info.admin_time); %in hours
        
        % grab drug parameters for model:
        med_ind =contains(aed_params.medication,med_names(n));
        
        % if med_ind is zero, then drug not exist or its wine
        if sum(med_ind)==0
            disp(['warning: drug ' med_names{n} ' not found...']);
        else
            F = aed_params(med_ind,:).F;
            vd=aed_params(med_ind,:).vd;
            ka=aed_params(med_ind,:).ka;
            tmax = aed_params(med_ind,:).t_max;
            % take mean of range for current model:
            tHalf = [aed_params(med_ind,:).t_half_e]; tHalf=strsplit(tHalf{1},',');tHalf = cellfun(@str2double,tHalf);
            tHalf = mean(tHalf);
            
            %% titrate the half life of that drug if a lab result of that drug exists
            
            
            % if exists,
            % tHalf_pt = titrate_tHalf()
            % then set tHalf = tHalf_pt
            %store the patient tHalf to compare tHalfs of patients to
            %the published ranges
            
            
            % calculate c0, the initial concentration at the time of their EMU
            % stay
            
            %% get the initial concentration at EMU admission following surgery
            %get home meds
            % load spreadsheet of outpatient labs and meds early, and find and
            % parse entries for that patient [drug dose frequency]
            % find med_ind of this drug in their home meds
            ind=contains(pt_home_meds.medication,med_names(n));
            if sum(ind)~=0
                dose = pt_home_meds.dose1(ind); % later take into account nonuniform dose schedule
                frequency = home_meds.frequency(ind);
                
                %% need to use implant date and time and date/time of first implant to pick c0 from c_ss
                [~,css]=get_c0(dose,frequency,tHalf,F,vd,tmax,ka); %this is curve of one day
                c0=css(round((mod(this_med_schedule(1),24./frequency))*60)); %this is the c0 at the first admin,indexed from t=0 (midnight) so need to adjust based on sleep/wake estimate
            else
                c0=0;
            end
            
            %%
            c_drug=[];
            dose = this_med_info.Dose(1);
            for i = 2:height(this_med_info)
                tInt = this_med_schedule(i-1);
                [c_t,~]= get_single_dose_curve(c0,dose,tHalf,tInt,F,vd,tmax,ka);
                c_drug=[c_drug c_t];
                c0 = c_drug(end);
                
                dose = this_med_info.Dose(i);
            end
            
            c_all(n) = {c_drug};
            all_dose_curves(ipt) ={c_all};
            
            tHr = linspace(this_med_info.admin_time(1),this_med_info.admin_time(end),length(c_drug));
            t_all(n) = {tHr};
            all_tHr(ipt) = {t_all};
            
            
        end
    end
end


end