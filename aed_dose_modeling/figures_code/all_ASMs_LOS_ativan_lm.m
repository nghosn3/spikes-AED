close all;clear;
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling/figures_code')
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');


% load the aed metadata
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

% Get which patients have AED data, load the data
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% load home medications
load('home_meds.mat');
load('MAR_032122.mat')

% get overall decrease and decrease in first 24hrs
[tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds);
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,0,0); 


%%
% convert doses to proportion of minimum effective dose
taper_info =all_taper_info; %with adjusted for minimum effective dose

medications = vertcat(all_med_names{:});
all_medications = unique(medications);

all_asm= zeros(length(ptIDs),length(all_medications));

aed_decrease = zeros(1,length(taper_info));
aed_day1_decrease = zeros(1,length(taper_info));
for i = 1:length(taper_info)
    
    func1 = @(x) x(2);
    min_func = @(x) min(x(2:end));
    norm_func = @(x) x./max(x);
    curr_pt = taper_info{i};
    pt_meds = {curr_pt.med_name};
    
    for n=1:length(pt_meds)
        med_ind = strcmp(pt_meds{n},all_medications);
        if ~isempty(curr_pt(n).dose_schedule)
        curr_pt_doses = [curr_pt(n).dose_schedule{:}];
        
        day1_doses = curr_pt_doses(2);
        day2_doses = func2(curr_pt_doses);
        
        day1_decrease = (day1_doses - day2_doses) ./ day1_doses;
        inf_inds = day1_decrease == -Inf; day1_decrease(inf_inds)=0;
        all_asm(i,med_ind) = day1_decrease;
        end 
        
    end
end


%% find the seizures followed by ativan administration
all_seizures =table(); % add ptID, seizure ID, preictal AED load, and binary yes/no for ativan admin w/in 1hr

all_times = [];
seizure_offsets = [];
start_ind =1;
for ipt=1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(ipt))];
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    
    %get drug curve to grab pre-ictal levels
    [med_names,~,~] = parse_MAR(ptID,all_meds);
    
    % get the total AED dose over time
    drugs =zeros(length(med_names),max_dur*60); %450 hours of EMU stay in minutes
    for i =1:length(med_names)
        drug=all_dose_curves{ipt}{i};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(all_tHr{ipt}{i}(1)*60)-1;
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum=drug_sum./max(drug_sum); %normalize for number of drugs
    drug_sum(drug_sum==0) =NaN; %not include all zeros in average, and in histogram
    
    % get seizure times
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    % remove clustered seizures <2hrs apart
    cluster_diff = diff(seizure_times(:,1));
    noncluster_inds = [1; cluster_diff>1];
    seizure_times=seizure_times(logical(noncluster_inds),:);
    
    end_ind = start_ind + length(seizure_times);
    
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    %remove seizure inds that are greater than emu recording
    %time =emu_dur(ipt)*60;%number of minutes of emu stay
    %seizure_inds(seizure_inds>time)=[];
    
    
    % med curves are sampled at one point per minute
    pt_preictal_1hr = zeros(1,length(seizure_inds));
    for i=1:length(seizure_inds)
        pt_preictal_1hr(i) = nanmean(drug_sum(seizure_inds(i)-60:seizure_inds(i)));        
    end
    
    %Get ativan times
    [~,meds,~] = parse_MAR(ptID,all_meds);
    ativan_inds = strcmp(meds.medication,'lorazepam');
    times = meds.admin_time(ativan_inds);
    
    
    time_to_closest_sz = nan(length(seizure_times(:,1)),2);
    for n=1:length(seizure_times(:,1))
        sz_diffs = seizure_times(n,1)-times;
        before_ativan = sz_diffs< 0;
        if ~isempty(times)
            if ~isempty(before_ativan) && ~(sum(before_ativan)==0)
                [mval,~]=min(abs(sz_diffs(before_ativan)));
                ind = find(abs(sz_diffs)==mval);
                time_to_closest_sz(n,:)=[times(ind(1)) seizure_times(n,1)];
            end
        end
    end
    all_seizures.ptID(start_ind:end_ind-1)=ptIDs(ipt)*ones(1,end_ind-start_ind);
    all_seizures.seizureEEC(start_ind:end_ind-1) = seizure_times(:,1);
    all_seizures.preictal_aed_load(start_ind:end_ind-1) = pt_preictal_1hr;
    all_seizures.t_closest_ativan(start_ind:end_ind-1) =  time_to_closest_sz(:,1);
    
     start_ind = end_ind;
end

% remove the columns for drugs with less than 20 patients, and remove ativan
med_counts=cell(1,length(all_medications));
for n =1:length(all_medications)
    med_counts{n}=sum(contains(medications,all_medications{n}));
end
med_info = [all_medications med_counts'];
exclude_meds = [med_counts{:}]<20;
all_medications(contains(all_medications,'valproic acid'))={'valproic_acid'};

top_meds=all_medications(~exclude_meds);


all_seizures.ativan_sz = double(all_seizures.t_closest_ativan - all_seizures.seizureEEC <= 1 &  all_seizures.t_closest_ativan - all_seizures.seizureEEC >=0);

% get binary of that patient had a convulsion or not 
has_conv = false(length(ptIDs),1);
for i=1:length(ptIDs)
   pt_inds = all_seizures.ptID == ptIDs(i);
   has_conv(i) = any(all_seizures.ativan_sz(pt_inds));
end

% get baseline seizure frequencies 
baseline_sz_freqs = readtable('no_phi_baseline_sz_freq.xlsx');
sz_freqs = zeros(length(ptIDs),1);
for i = 1:length(ptIDs)
    ind = baseline_sz_freqs.HUP_ID == ptIDs(i);
    sz_freqs(i) = baseline_sz_freqs.sz_per_month(ind);
end

%% run model

% fill table
tbl = table();
for n=1:length(all_medications)
    name = all_medications(n);
    tbl.(name{:})=all_asm(:,n);
end 
tbl(:,exclude_meds)=[];

tbl.baseline_sz_freq = sz_freqs;
tbl1 = tbl;
tbl1.has_conv = double(has_conv);

mdl_severity = fitglm(tbl1,'distribution','binomial')
coef_names=mdl_severity.CoefficientNames(2:end);
coefs = mdl_severity.Coefficients.Estimate(2:end);
SEs = mdl_severity.Coefficients.SE(2:end);
ORs = (coefs);
CIs = ([coefs+(1.96*SEs) coefs-(1.96*SEs)]) %SE*1.96 for 95% CI
inds = length(ORs):-1:1;

figure;
neg = ORs-CIs(:,2);
pos = CIs(:,1)-ORs;
errorbar(ORs,inds,neg,pos,'.k','horizontal','linewidth',2); hold on;
plot((ORs),inds,'.r','markersize',20);

xline(0,'--r','linewidth',1.5)

yticks(1:length(ORs))
yticklabels(coef_names(inds))
xlabel('model coefficient'); 
%xlim([-2 20])

sig_coeffs = mdl_severity.Coefficients.pValue(2:end)   <0.05;
plot(-1,inds(sig_coeffs),'*k');


%%
tbl2 = tbl;
tbl2.length_stay = hours(cohort_info.Explant_Date - cohort_info.Implant_Date);

mdl_LOS = fitlm(tbl2)

% can also try time to first seizure 

function out = func2(x)
if length(x)>2
    out = x(3);
else
    out =NaN;
end
end