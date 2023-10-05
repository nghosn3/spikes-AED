function [mdl,feats_train,feats_test] = run_NULL_logistic_on_AED(feats_train,feat_names)

fixed_effects=[];
for n=1:length(feat_names)-2 % do not include seizure or ptID, the response variable
    name =feat_names{n};
    fixed_effects = [fixed_effects ' + ' name];
end

random_effects = ' + (1 | ptID)';
modelspec = ['seizure ~ 1 ' fixed_effects random_effects];
tbl = feats_train;
mdl = fitglme(tbl,modelspec,'Distribution','binomial');

end