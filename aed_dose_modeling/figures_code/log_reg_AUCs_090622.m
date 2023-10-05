% run logistic regression to get AUCs and check what is broken 09/06/22

%% get AUC and operating point
% estimates 
load('coefficient_estimates_090622.mat')
null_names = [{'baseline_sz_freq'};{'t_last_seizure'}; {'ptID'}; {'seizure'}];

pts = 1:length(ptIDs);

pt_inds_master = cell(1,iter);
for i = 1:iter
    %random sample of 80% of patients
    pt_inds_master{i} =randperm(length(pts)); 
end

model = 1;
[OPs,AUCs]  = model_bootstrap_ROC(drug_table,feat_names,model,iter,pt_inds_master);
disp('done with model AUCs')
model =0;
[OPs_null,AUCs_null,] =model_bootstrap_ROC(drug_table,null_names,model,iter,pt_inds_master);
disp('done with null model AUCs')

AUCout = bootstrap_ci_and_p(AUCs,AUCs_null,1);
save('results_AUC_stats_090622.mat');
