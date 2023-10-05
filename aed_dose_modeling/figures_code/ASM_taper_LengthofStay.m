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
%[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,0,0); 


%%
% convert doses to proportion of minimum effective dose
taper_info = cell(1,length(all_taper_info)); %with adjusted for minimum effective dose
for i=1:length(all_taper_info)
    pt_info=all_taper_info{i};
    for n=1:width(pt_info)
        if ~isempty(pt_info(n).med_name)
            med_ind = contains(aed_params.medication,pt_info(n).med_name);
            min_dose = aed_params(med_ind,:).min_dose_single_mg * aed_params(med_ind,:).min_dose_freq ;
            doses = pt_info(n).dose_schedule;
            if iscell(doses)
                doses = doses{:}; %something weird with cell array
            end
            pt_info(n).min_dose_schedule =  {doses ./ min_dose};            
        end
    end
    taper_info{i} = pt_info;
end


aed_decrease = zeros(1,length(taper_info));
aed_day1_decrease = zeros(1,length(taper_info));
for i = 1:length(taper_info)
    
    func1 = @(x) x(2);
    min_func = @(x) min(x(2:end));
    norm_func = @(x) x./max(x);
    curr_pt = taper_info{i};
    
    curr_pt_doses = [curr_pt(:).dose_schedule];
    curr_pt_doses_norm = cellfun(norm_func,curr_pt_doses,'UniformOutput',false);
    len = min(cellfun(@length,curr_pt_doses_norm));
    curr_pt_doses_norm = cellfun(@(x)x(1:len),curr_pt_doses_norm,'UniformOutput',false);%match lengths
    [min_daily_dose,day_ind] = min_func(sum(vertcat(curr_pt_doses_norm{:}),1)); day_ind = day_ind+1; % since excluding first day in indexing    
    
    day1_doses = cellfun(func1,curr_pt_doses,'UniformOutput',false);
    day2_doses = cellfun(@func2,curr_pt_doses,'UniformOutput',false);
    
    get_min_dose = @(x) x(day_ind);
    min_doses = cellfun(get_min_dose,curr_pt_doses,'UniformOutput',false); %first half day excluded from analysis
    
    all_aed_decreases = ([day1_doses{:}]- [min_doses{:}]) ./ [day1_doses{:}];
    inf_inds = all_aed_decreases ==-Inf; all_aed_decreases(inf_inds) = 0;
    aed_decrease(i) = nanmean(all_aed_decreases);
    day1_decrease = ([day1_doses{:}] - [day2_doses{:}]) ./ [day1_doses{:}];
    aed_day1_decrease(i) = nanmean(day1_decrease);
    
end

%% get which seizures are followed by ativan
all_seizures =table(); % add ptID, seizure ID, preictal AED load, and binary yes/no for ativan admin w/in 1hr
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);

%% find the seizures followed by ativan administration
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


% get time to first seizure (in EMU time for all)
% time_to_first_seizure = zeros(length(ptIDs),1);
% for i = 1:length(ptIDs)
%     ptID = ['HUP' num2str(ptIDs(i))];
%     offsets = ieeg_offset{2,i};
%     ieeg_offset_datasets = ieeg_offset{1,i};
%     [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID);
%     time_to_first_seizure(i)=seizure_times(1); %the first seizure
% end

%% run model for severity 
tbl = table();
tbl.asm_decrease = aed_day1_decrease'; %aed_decrease';
tbl.baseline_sz_freq = sz_freqs;%./max(sz_freqs); % normalize like other feature
tbl1 = tbl;
tbl1.has_conv = double(has_conv);
%tbl.time_to_first_seizure = time_to_first_seizure;

% find -inf and negative values and make them zero for no decrease
zero_inds = tbl1.asm_decrease == -inf | tbl1.asm_decrease <0;
tbl1.asm_decrease(zero_inds)=0;


mdl_severity = fitglm(tbl1,'distribution','binomial')

%% run model for length of stay
tbl2 = tbl;
tbl2.length_stay = hours(cohort_info.Explant_Date - cohort_info.Implant_Date);
zero_inds = tbl2.asm_decrease == -Inf | tbl2.asm_decrease <0;
tbl2(zero_inds,:)=[];
outliers = find(isoutlier(tbl2.baseline_sz_freq,'quartiles'));
tbl2(39,:)=[];
mdl_LOS = fitlm(tbl2)

%% plot stuff 
figure;
subplot(1,2,1);
x1=tbl2.asm_decrease;
x2 = tbl2.baseline_sz_freq;
y= tbl2.length_stay;
plot3(x1,x2,y,'.k','markersize',15); axis square; hold on;
 
ylim([0 300])
title('linear model for length of stay')
xlabel('ASM decrease');ylabel('baseline seizure frequency');zlabel('length of stay (hrs)');

% plot mesh on top
b = [mdl_LOS.Coefficients.Estimate];
x1fit = linspace(min(x1),max(x1),50);
x2fit = linspace(min(x2),max(x2),50);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)


subplot(1,2,2)
decrease_had_conv =aed_day1_decrease(has_conv);
decrease_no_conv =aed_day1_decrease(~has_conv);
len = max([length(decrease_had_conv) length(decrease_no_conv)]);
data = nan(2,len);
data(1,1:length( decrease_had_conv))= decrease_had_conv;
data(2,1:length( decrease_no_conv))= decrease_no_conv; 

boxplot(data',[{'no ativan sz'},{'had ativan sz'}]); axis square;
p=ranksum(decrease_had_conv,decrease_no_conv)
save_path='/Users/ninaghosn/Documents/Litt_Lab/projects/Pioneer/AED-taper-networks/results-figures/';
%print([save_path 'fig04_los_conv.eps'],'-depsc2','-painters', '-tiff', '-r300', '-f')


% coefs = mdl_severity.Coefficients.Estimate(2:end);
% SEs = mdl_severity.Coefficients.SE(2:end);
% ORs = exp(coefs);
% CIs = exp([coefs+(1.96*SEs) coefs-(1.96*SEs)]);
% coef_names=mdl_severity.CoefficientNames(2:end);
% inds = length(ORs):-1:1;
% 
% neg = ORs-CIs(:,2);
% pos = CIs(:,1)-ORs;
% errorbar(ORs,inds,neg,pos,'.k','horizontal','linewidth',2); hold on;
% plot((ORs),inds,'.r','markersize',20); axis square;
% xline(1,'--r','linewidth',1.5)
% yticks(1:length(ORs))
% yticklabels(coef_names(inds))
% xlabel('Odds Ratio'); 

function out = func2(x)
if length(x)>2
    out = x(3);
else
    out =NaN;
end
end