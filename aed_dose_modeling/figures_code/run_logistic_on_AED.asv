function [mdl] = run_logistic_on_AED(feats_train,feat_names)

fixed_effects=[];
for n=1:length(feat_names)-2 % do not include seizure or ptID, the response variable
    name =feat_names{n};
    fixed_effects = [fixed_effects ' + ' name];
end

random_effects = ' + (1 | ptID)';
modelspec = ['seizure ~ 1 ' fixed_effects random_effects];
%mdl = fitglme(feats_train,modelspec,'Distribution','binomial');

try
    mdl = fitglme(feats_train,modelspec,'Distribution','binomial');
    
catch
    disp('there was an error')
    
    mdl = [];
end


end
