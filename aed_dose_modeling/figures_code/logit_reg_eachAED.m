%% run regression for all patients and each drug seperatly
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
%[all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);

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
    drugs =nan(length(med_names),ceil(max_time*60));
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
    drug_wins = nan(1,nbins);
    
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
    
    % fill drugs with no levels 
    empty_inds =cellfun(@isempty,features(:,i));
    features(empty_inds,i)={zeros(1,nbins)};
    
    % add ptID feature
    features(end-1,i)={i*ones(1,nbins)}; %ptID
    
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    % adjust seizure times to winLen minute bins
    seizure_times = seizure_times*60./winLen;
    seizure_times(seizure_times(:,1)>(max_time*60./winLen),:)=[]; %this is for you, HUP175!
    
    % shift seizure bins left so 1 = seizure in next time bin
    sz_bins = zeros(1,nbins+1);
    sz_bins(round(seizure_times(:,1))) =1;
    sz_bins(1)=[];
    
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



%% run logistic regression
features = med_feats; %[all_feats; all_pt_bins; all_sz_bins];
feat_names = [top_meds; {'ptID'}; {'seizure'}];
drug_table = table();
for i =1:length(feat_names)
    feat_name = feat_names{i};
    this_feat = features(i,:);
    this_feat = horzcat(this_feat{:})';
    
    % replace NaN with zeros
    %this_feat(isnan(this_feat))=0;
    
    drug_table.(feat_name)=this_feat;
end
%drug_table.seizure = categorical(drug_table.seizure);


figure;
for n = 1:length(top_meds)
    
    pts = 1:length(ptIDs);
    feats = features;
    
    % remove the pts that didnt take this aed
    has_aed = cellfun(@nansum,feats) > 0;
    has_this_aed = has_aed(n,:);
    this_pts = pts(has_this_aed);
    
    pt_inds =[];
    for ipt = 1:length(this_pts)
        pt_inds = [pt_inds; find(drug_table.ptID == this_pts(ipt))]; % exclude patients without that AED
            
    end
    
    iter=10;
    pts = length(this_pts);
    
    AUCs = zeros(1,iter);
    xROCs = cell(iter,1);
    yROCs = cell(iter,1);
    accuracies = zeros(1,iter);
    
    this_feat_names = [{top_meds{n}} {'ptID'},{'seizure'}];
    this_table = drug_table(pt_inds,[{top_meds{n}} {'ptID'},{'seizure'}]); % last two rows are ptID and seizure;
    this_table.ptID = categorical(this_table.ptID);

    for i = 1:iter
        pt_inds = randperm(pts);
        split=floor(.7*pts);
        train_inds = this_pts(pt_inds(1:split));
        test_inds = this_pts(split+1:end);
        
        test_pts = cellstr(strsplit(num2str(test_inds)));
        train_pts =  cellstr(strsplit(num2str(train_inds)));
        
        all_feats_test = [];
        for j = 1:length(test_pts)
            all_feats_test = [all_feats_test; this_table(contains(cellstr(this_table.ptID),test_pts(j)),:)];
        end
        feats_test = all_feats_test(:,1:end-1);
        labels_test = all_feats_test.seizure;
        
        feats_train = [];
        for j = 1:length(train_pts)
            feats_train = [feats_train; this_table(contains(cellstr(this_table.ptID),train_pts(j)),:)];
        end
        
        %modelspec = 'seizure ~ 1 + aed_load + ptID';
        fixed_effects=[];
        %random_effects=[];
        for m=1:length(this_feat_names)-2 % do not include seizure or ptID, the response variable
            name =this_feat_names{m};
            fixed_effects = [fixed_effects ' + ' name];
            %random_effects = [random_effects ' + (' name ' | ptID)'];
        end
        
        random_effects = ' + (1 | ptID)';
        
        
        modelspec = ['seizure ~ 1 ' fixed_effects random_effects];
        tbl = feats_train(:,:);
        
        mdl = fitglme(tbl,modelspec,'Distribution','binomial');
        
        tbl = feats_test(:,:);
        sz_pred = predict(mdl,feats_test);
        
        intervals= linspace(0, 1, 1000);
        
        subplot(3,4,n)
        [Xroc, Yroc, ~, AUCs(i)]= perfcurve(labels_test,sz_pred,1,'Xvals',intervals);
    
        plot(Xroc, Yroc, 'LineWidth', 1); hold on;
        
        %for getting an average of the ROC curves from each fold
        x_adj= adjust_unique_points(Xroc); %interp1 requires unique points
        if i==1 %if is the first fold
            mean_curve= (interp1(x_adj, Yroc, intervals))/iter;
        else
            mean_curve= mean_curve+ (interp1(x_adj, Yroc, intervals))/iter;
        end

    end
    
    
end


plot(intervals, mean_curve, 'Color', 'Black', 'LineWidth', 3.0);
axis square; title([feat_names{n} ': AUC = ' num2str(mean(AUCs))])
ylabel('TPR'); xlabel('FPR'); xlim([0 1]);






function x= adjust_unique_points(Xroc)
x= zeros(1, length(Xroc));
aux= 0.0001;
for i=1: length(Xroc)
    if i~=1
        x(i)= Xroc(i)+aux;
        aux= aux+0.000001;
    end
    
end
end

