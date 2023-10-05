%% run linear mixed model for curves and drugs

close all;clear;
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% get AED curves
load('MAR_032122.mat')

% using the AED BPL model
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);

% using only interpolation dosage 
%all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);

%%
exclude_ativan=0;

medications = vertcat(all_med_names{:});
all_medications = unique(medications);
if exclude_ativan
    ativan_ind = contains(all_medications,'lorazepam');
    all_medications(ativan_ind)=[];
end 

% make matrix with all variables (levels for each drug)
num_feats=length(all_medications)+2;
sz=num_feats;
features = cell(num_feats,length(ptIDs));
winLen =120;

for i=1:length(ptIDs) % exclude HUP136
    ptID = ['HUP' num2str(ptIDs(i))];
    % for grabbing seizure times 
    offsets = ieeg_offset{2,i};
    ieeg_offset_datasets = ieeg_offset{1,i};
    
    med_names = all_med_names{i};
    pt_curves = all_dose_curves{i};
    tHr=all_tHr{i};
    
    % exclude ativan in drug sum
    if exclude_ativan
        ativan_ind = contains(med_names,'lorazepam'); %#ok<*UNRCH>
        med_names(ativan_ind)=[];
        pt_curves(ativan_ind)=[];
        tHr(ativan_ind)=[];
    end
    
    % align drugs on same time scale
    end_times = cellfun(@(x)x(end),tHr);
    curve_lens = cellfun(@length,pt_curves);
    max_time = max(end_times);
    drugs =zeros(length(med_names),ceil(max_time*60)); 
    for n =1:length(med_names)
        drug=pt_curves{n};
        %drug=drug./nanmax(drug); %normalize drug curve again?
        if ~isempty(drug)
            dStart = round(tHr{n}(1)*60)-1;
            drugs(n,dStart+1:dStart+length(drug))=drug;
        end
    end
    
    % get binned drug values
    nbins=ceil(length(drugs)./winLen);
    drug_wins = zeros(1,nbins);
    
    for r = 1:length(med_names)
        ind =1;
        for c = 1:winLen:length(drugs)-winLen
            drug_wins(r,ind)=nanmean(drugs(r,c:c+winLen));
            ind=ind+1;
        end
    end
    
    for n=1:length(med_names)
        ind = contains(all_medications,med_names{n});
        features(ind,i)={drug_wins(n,:)};
    end
    
    empty_inds =cellfun(@isempty,features(:,i));
    features(empty_inds,i)={zeros(1,nbins)};
    features(end-1,i)={i*ones(1,nbins)}; %ptID
    
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    % adjust seizure times to winLen minute bins
    seizure_times = seizure_times*60./winLen;
    seizure_times(seizure_times(:,1)>(max_time*60./winLen),:)=[]; %this is for you, HUP175!
    
    sz_bins = zeros(1,nbins);
    sz_bins(round(seizure_times(:,1))) =1;
    
    features(sz,i)={sz_bins};
    
    
end

%%
all_medications(contains(all_medications,'valproic acid'))={'valproic_acid'};

% get the medication counts, exclude the meds with less than 5 subjects - exclude subjects with none of the used meds
med_counts=cell(1,length(all_medications));
for n =1:length(all_medications)
    med_counts{n}=sum(contains(medications,all_medications{n}));
end
med_info = [all_medications med_counts'];
exclude_meds = [med_counts{:}]<5;
top_meds=all_medications(~exclude_meds);
med_feats=features;
med_feats(exclude_meds,:)=[];

feat_names = [top_meds ;{'ptID'}; {'seizure'}];

drug_table = table();
for i =1:length(feat_names)
    feat_name = feat_names{i};
    this_feat = med_feats(i,:);
    this_feat = horzcat(this_feat{:})';
    drug_table.(feat_name)=this_feat;
end 

% Use Lasso regression to choose features for mixed model 
X = table2array(drug_table(:,top_meds));
y = table2array(drug_table(:,'seizure'));
[B,FitInfo] = lasso(X,y,'CV',10);

lam=FitInfo.Index1SE;   
B(:,lam);
[row, col] = find(B==0);
[~, ia, ~] = unique(row, 'rows');
zeroed = horzcat(row, col);
zeroed = sortrows(zeroed(ia, :), 2,'descend');

new_feat_names = top_meds(zeroed(:,1));
new_feat_names = [new_feat_names; {'ptID'}; {'seizure'}]; % chose top 9 features since they died >90 lambda
% run linear regression on best features
fixed_effects=[];
%random_effects=[];
for n=1:length(new_feat_names)-2 % do not include seizure or ptID, the response variable
    name =new_feat_names{n};
    fixed_effects = [fixed_effects ' + ' name];
    %random_effects = [random_effects ' + (' name ' | ptID)'];
end 

random_effects = ' + (1 | ptID)';
% change ptID feature to nominal 
drug_table.ptID=categorical(drug_table.ptID);
%%
lme = fitlme(drug_table(:,new_feat_names),['seizure ~ 1'  fixed_effects random_effects]);

