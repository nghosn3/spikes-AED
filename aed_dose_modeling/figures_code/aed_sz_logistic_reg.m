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

% using AED BPLs
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);

% using dosage schedule
%[all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);


%% calculate features 

all_feats = cell(1,length(ptIDs));
all_sz_bins = cell(1,length(ptIDs));

tapered = 0;
rm_ativan=0;
winLen=60; %minutes 
for i = 1:length(ptIDs)
        pt = ['HUP' num2str(ptIDs(i))];
        [med_names,meds,explant_date] = parse_MAR(pt,all_meds);
        pt_curves = all_dose_curves{i};
        tHr = all_tHr{i};
        
        %only use tapered drugs to get drug sum:
        if tapered
            tapered_drugs = get_taper_info(meds, med_names);
            if isfield(tapered_drugs,'med_name')
            med_names = {tapered_drugs.med_name};
            else 
                med_names = [];
            end 
        end
        
        if rm_ativan
            ativan =contains(med_names,'lorazepam');
            med_names(ativan)=[];
            pt_curves(ativan)=[];
            tHr(ativan)=[];
        end
        
        % get drug sum curve (AED load)
        drugs =zeros(length(med_names),max_dur*60); %450 hours of EMU stay in minutes
        for n =1:length(med_names)
            drug=pt_curves{n};
            drug=drug./nanmax(drug); %normalize each drug curve
            if ~isempty(drug)
                dStart = round(tHr{n}(1)*60)-1;
                drugs(n,dStart+1:dStart+length(drug))=drug;
            end
        end
        drug_sum=nansum(drugs,1);
        drug_sum=drug_sum./length(med_names);
        
        % cut off drug curve to only be length of emu stay
        time =emu_dur(i)*60; % minutes
        if time<length(drug_sum)
            drug_sum = drug_sum(1:time);
        end
        
        % F2: now average AED load in winLen minute, nonoverlapping windows
        nbins=ceil(length(drug_sum)./winLen);
        ind =1;
        drug_sum_wins = zeros(1,nbins);
        
        for j = 1:winLen:length(drug_sum)-winLen
            drug_sum_wins(ind)=mean(drug_sum(j:j+winLen));
            ind=ind+1;
        end
        
        %get seizure times
        offsets = ieeg_offset{2,i};
        ieeg_offset_datasets = ieeg_offset{1,i};
        [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,pt); % seizure times in hours of EMU stay
        
        seizure_times = seizure_times*60./(winLen); %convert to winLen minute blocks
        sz_inds = round(seizure_times(:,1));
        
        sz_bins = zeros(1,length(drug_sum_wins));
        sz_bins(sz_inds) =1;
        
        % gut check 1 for drug_sum_wins
%         tvec = linspace(1,time./60,length(drug_sum_wins));
%         plot(tvec,drug_sum_wins); hold on;
%         plot(tvec(sz_inds),drug_sum_wins(sz_inds),'o'); title(ptID); hold off;
        
        % construct feature matrix
        features = drug_sum_wins;
        
        zero_inds = drug_sum_wins==0;
        features(zero_inds)=NaN;
        %sz_bins(zero_inds)=[];
        
        % gut check 2 for drug_sum_wins
%         tvec = linspace(1,time./60,length(features));
%         plot(tvec,features); hold on;
%         plot(tvec(sz_inds),features(sz_inds),'o'); title(pt); hold off;
        
        
        all_feats{i}=features;
        all_sz_bins{i}= (sz_bins);
        
    
end

%% run logistic regression

pts = ptIDs;
med_feat = all_feats;
sz_bins = all_sz_bins;
iter=50;

AUCs = zeros(1,iter);
xROCs = cell(iter,1);
yROCs = cell(iter,1);
accuracies = zeros(1,iter);

figure(1)
for i = 1:iter
    pt_inds = randperm(length(pts));
    split=floor(.7*length(pts));
    train_inds = pt_inds(1:split);
    test_inds = pt_inds(split+1:end);
    
    feats_test = [];
    labels_test=horzcat(all_sz_bins{test_inds})+1;
    for l =1:height(med_feat)
        %feat_name = feat_names{i};
        this_feat = med_feat(l,test_inds);
        this_feat = horzcat(this_feat{:})';
        feats_test = [feats_test; this_feat'];
        
    end
    
    feats_train = [];
    labels_train=horzcat(all_sz_bins{train_inds})+1;
    for l =1:height(med_feat)
        %feat_name = feat_names{i};
        this_feat = med_feat(l,train_inds);
        this_feat = horzcat(this_feat{:})';
        feats_train = [feats_train; this_feat'];
    end
    
    
    
    B = mnrfit(feats_train',labels_train');
    sz_hat = mnrval(B,feats_test');
    
    
    sz_prob = sz_hat(:,1);
    
    sz_and_pred = [labels_test', sz_prob];
    sz_inds = labels_test == 2;
    % to get accuracy, need to choose a seizure prediction threshold - base it off the distribution?
    
    sz_preds = sz_and_pred(sz_inds,:);
    nonsz_preds = sz_and_pred(~sz_inds,:);
    
    intervals= linspace(0, 1, 3000);
    [Xroc, Yroc, ~, AUCs(i)]= perfcurve(labels_test,sz_prob,1,'Xvals',intervals);
    plot(Xroc, Yroc, 'LineWidth', 1); hold on;
    
    %for getting an average of the ROC curves from each fold
    x_adj= adjust_unique_points(Xroc); %interp1 requires unique points
    if i==1 %if is the first fold
        mean_curve= (interp1(x_adj, Yroc, intervals))/iter;
    else
        mean_curve= mean_curve+ (interp1(x_adj, Yroc, intervals))/iter;
    end
end

plot(intervals, mean_curve, 'Color', 'Black', 'LineWidth', 3.0);
axis square; title(['All pts, all AEDs: AUC= ' num2str(mean(AUCs))])
ylabel('TPR'); xlabel('FPR'); xlim([0 1]); hold off;

figure(2)
edges = .95:.0005:1;
histogram(sz_preds(:,2),'Normalization','probability','BinEdges',edges); hold on;
histogram(nonsz_preds(:,2),'Normalization','probability','BinEdges',edges); hold on;
legend('seizure','non-seizure')
ylim([0 .2])

function x= adjust_unique_points(Xroc)
x= zeros(1, length(Xroc));
aux= 0.0001;
for i=1: length(Xroc)
    if i~=1
        x(i)= Xroc(i)+aux;
        aux= aux+0.0001;
    end
    
end
end
