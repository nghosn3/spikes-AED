% class_sz and class_no_sz are for seizure labels and non seizure labels
% disc is the cutoff for classification of sz/nonsz

function out = make_conf_matrix(class_sz,class_no_sz,disc)

tp = nansum(class_sz>disc);
fn = nansum(class_sz<disc);
fp = nansum(class_no_sz>disc);
tn = nansum(class_no_sz<disc);

sensitivity = tp/(tp+fn); % what % of those truly positive test positive?
specificity = tn/(tn+fp); % what % of those truly negative test negative?
ppv = tp/(tp+fp); % what % of those predicted to be positive are positive?
npv = tn/(tn+fn); % what % of those predicted to be negative are negative?

mat = [tp fn;fp tn];

out.mat = mat;
out.sensitivity = sensitivity;
out.specificity = specificity;
out.ppv = ppv;
out.npv = npv;

end