%% machine learnng models to eluciate the relationship between ASMs and Spike rate and seizures 

close all;clear;
cd('/Volumes/USERS/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/spikes-AED/aed_dose_modeling/figures_code'])
addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/ASM-spikes-analysis'])
addpath('/Volumes/users/nghosn3/tools/wavelet-coherence')

% % add fooof stuff to path
addpath(genpath('/Volumes/users/nghosn3/Pioneer/fooof_mat'))
% addpath('/Volumes/users/nghosn3/Pioneer/fooof_mat/fooof/')
% addpath('/Volumes/users/nghosn3/Pioneer/fooof_mat/fooof_mat/fooof/fooof/')

%load spike rate and med data - new from 2/13/23 (samp/10min)
spikes_fname = 'spikes_rates_021323.mat';
load(spikes_fname);
% get the seizure and SOZ localization information
soz_info = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve','Sheet','SOZ');

cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

meds_fname = 'MAR_032122.mat';
load(meds_fname);

[all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
[all_pts_drug_samp, all_pts_tvec] = get_samp_asm_spikes(meds_fname,spikes_fname,all_ieeg_offset, all_dose_curves, all_tHr,ptIDs);

%% 

features=table();
features.ptid =zeros(0);
features.asm_load = zeros(0);
features.tbin = zeros(0);
features.sz =zeros(0);
features.spike_rate =zeros(0);

% later- patient selection. exclude patients that had too many seizures with no medication taper
start_ind = 1;
for i=1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    % get average asm_load an times 
    avg_asm_load = mean(all_pts_drug_samp{i})';
    tbin = (1:length(avg_asm_load))';
    ptids = ptIDs(i) * ones(size(tbin));
    spike_rate = (all_spike_rate{i})./max(all_spike_rate{i});

    sz = zeros(size(tbin));
    % add seizure times
    [seizure_times] = get_seizure_times_from_sheet(ptID);

    for j = 1:length(seizure_times)
        [~,sz_ind] = min(abs(all_spike_times{i}-seizure_times(j,1)));
        sz(sz_ind) = 1;
    end

    
    ind1 = start_ind;
    ind2= start_ind+length(ptids)-1;
    
    features.ptid(ind1:ind2) =ptids;
    features.asm_load(ind1:ind2) = avg_asm_load;
    features.tbin(ind1:ind2) = tbin;
    features.sz(ind1:ind2) =sz;
    features.spike_rate(ind1:ind2) =spike_rate;

    start_ind = ind2+1;

end 


%% try a random forest to predict spike rate 
data = features;


% Specify the predictor variables and response variable
predictorVars = {'ptid', 'asm_load', 'tbin','sz'};
%predictorVars = {'ptid', 'tbin','sz'};

responseVar = 'spike_rate';

% Convert categorical variables to dummy variables
%features = dummyvar(features, 'categoricalVars', {'medication_name'});

% down sample the data to a data point each hour
% Group data by patient ID (ptid) and downsample within each group
% Define the downsampling factor (number of observations to be condensed)
downsampleFactor = 1;

% Group data by patient ID (ptid) and downsample within each group
downsampledData = table();
groups = findgroups(data.ptid);
uniqueGroups = unique(groups);
for i = 1:numel(uniqueGroups)
    groupIdx = groups == uniqueGroups(i);
    groupData = data(groupIdx, :);

    numGroups = ceil(size(groupData, 1) / downsampleFactor);
    for j = 1:numGroups
        startIndex = (j - 1) * downsampleFactor + 1;
        endIndex = min(j * downsampleFactor, size(groupData, 1));

        subgroupData = groupData(startIndex:endIndex, :);

        ptid = unique(subgroupData.ptid);
        tbin = max(subgroupData.tbin);
        sz = max(subgroupData.sz);
        spike_rate = sum(subgroupData.spike_rate);
        asm_load = nanmean(subgroupData.asm_load);

        downsampledData = [downsampledData; table(ptid, asm_load, tbin, sz, spike_rate)];
    end
end
%% Split the data into training, validation, and testing sets (adjust the ratios as needed)

trainRatio = 0.6;
valRatio = 0.2;
testRatio = 0.2;
%[trainData, valData, testData] = splitData(data, trainRatio, valRatio, testRatio);
[trainData, valData, testData] = splitData_within(data, trainRatio, valRatio, testRatio);
% Extract predictor variables and response variable from the training data
trainPredictors = trainData(:, predictorVars);
trainResponse = trainData.(responseVar);

%% Create a random forest model (adjust the number of trees and other parameters as needed)
numTrees = 100;
randomForest = TreeBagger(numTrees, trainPredictors, trainResponse,'Method','Regression','OOBPrediction','ON','OOBPredictorImportance','ON');

% Extract predictor variables and response variable from the validation data
valPredictors = valData(:, predictorVars);
valResponse = valData.(responseVar);

%% Predict spike rates for the validation data and calculate some measures of fit

% Evaluate model performance on the validation data (using mean squared error as an example)
valPredictions = predict(randomForest, valPredictors);
mse = mean((valPredictions - valResponse).^2);
fprintf('Mean Squared Error (Validation): %.4f\n', mse);

% Extract predictor variables and response variable from the testing data
testPredictors = testData(:, predictorVars);
testResponse = testData.(responseVar);

% Predict spike rates for the testing data
testPredictions = predict(randomForest, testPredictors);

% Evaluate model performance on the testing data
mse = mean((testPredictions - testResponse).^2);
fprintf('Mean Squared Error (Testing): %.4f\n', mse);

predictedValues = testPredictions;
trueResponse = testResponse;

% Calculate Root Mean Squared Error (RMSE)
rmse = sqrt(mse);

% Calculate R-squared (R²)
ssTotal = sum((trueResponse - mean(trueResponse)).^2);
ssResidual = sum((trueResponse - predictedValues).^2);
rSquared = 1 - (ssResidual/ssTotal);

% Display the evaluation metrics
fprintf('Mean Squared Error (MSE): %.4f\n', mse);
fprintf('Root Mean Squared Error (RMSE): %.4f\n', rmse);
fprintf('R-squared (R²): %.4f\n', rSquared);

%% look at the feature importance

% Get the predictor names
predictorNames = randomForest.PredictorNames;

% Initialize an array to store feature importance values
featureImportance = zeros(1, numel(predictorNames));

% Loop over each tree in the random forest
for i = 1:numTrees
    % Get the feature importance values for the current tree
    tree = randomForest.Trees{i};
    imp = predictorImportance(tree);
    
    % Accumulate the feature importance values
    featureImportance = featureImportance + imp;
end

% Calculate the average feature importance across all trees
featureImportance = featureImportance / numTrees;

% Sort the feature importance in descending order
[sortedImportance, sortedIndices] = sort(featureImportance, 'descend');

% Plot the feature importance
figure;
bar(sortedImportance);
xlabel('Features');
ylabel('Importance');
title('Feature Importance');
xticklabels(predictorNames(sortedIndices));
xtickangle(45);
%%
% Example patient ID
examplePatientID = 171;

% Filter data for the example patient
examplePatientData = data(data.ptid == examplePatientID, :);

% Extract predictor variables and true response values for the example patient
examplePredictors = examplePatientData(:, {'ptid','tbin', 'asm_load', 'sz'});
exampleTrueResponse = examplePatientData.spike_rate;

% Predict the response values using the trained model
examplePredictedResponse = randomForest.predict(examplePredictors);
% Plot the true and predicted values
figure;
plot(exampleTrueResponse, 'b', 'LineWidth', 2);  % True response values
hold on;
plot(examplePredictedResponse, 'r--', 'LineWidth', 2);  % Predicted response values
hold off;
xlabel('Observation');
ylabel('Response');
title(['True vs Predicted Response for Example Patient HUP' num2str(examplePatientID)]);
legend('True', 'Predicted');


%% perform leave one out cross validation rather than test/train split
% Initialize an array to store R-squared values for each patient
numPatients = length(ptIDs);
rSquaredValues = zeros(numPatients, 1);

% Perform leave-one-out cross-validation
for i = 1:numPatients
    % Split the data into training and testing sets
    trainData = data(data.ptid ~= ptIDs(i), :);
    testData = data(data.ptid == ptIDs(i), :);

    % Extract predictor variables and response variable from the training data
    trainPredictors = trainData(:, predictorVars);
    trainResponse = trainData.(responseVar);

    % Create a random forest model
    randomForest = TreeBagger(numTrees, trainPredictors, trainResponse,'Method','Regression','OOBPrediction','ON','OOBPredictorImportance','ON');

    % Extract predictor variables and response variable from the testing data
    testPredictors = testData(:, predictorVars);
    testResponse = testData.(responseVar);

    % Predict response values for the testing data
    testPredictions = predict(randomForest, testPredictors);

    % Calculate R-squared for the current patient
    ssTotal = sum((testResponse - mean(testResponse)).^2);
    ssResidual = sum((testResponse - testPredictions).^2);
    rSquared = 1 - (ssResidual / ssTotal);

    % Store the R-squared value
    rSquaredValues(i) = rSquared;
end

% Compute the overall mean R-squared value
meanRSquared = mean(rSquaredValues);

% Display the R-squared values for each patient
for i = 1:numPatients
    fprintf('R-squared for Patient %d: %.4f\n', ptIDs(i), rSquaredValues(i));
end

% Display the overall mean R-squared value
fprintf('Mean R-squared: %.4f\n', meanRSquared);

figure;bar(rSquaredValues);
xticks(1:length(ptIDs))
xticklabels(ptIDs)



