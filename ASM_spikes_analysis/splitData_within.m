function [trainData, valData, testData] = splitData_within(data, trainRatio, valRatio, testRatio)
    % Unique patient IDs
    uniquePatients = unique(data.ptid);
    
    trainData = table();  % Initialize empty table for training data
    valData = table();    % Initialize empty table for validation data
    testData = table();   % Initialize empty table for testing data

    % Iterate over each unique patient
    for i = 1:numel(uniquePatients)
        patientData = data(data.ptid == uniquePatients(i), :);  % Subset data for current patient

        % Calculate the number of observations for each split within the patient
        numObservations = height(patientData);
        numTrainObservations = round(trainRatio * numObservations);
        numValObservations = round(valRatio * numObservations);
        numTestObservations = numObservations - numTrainObservations - numValObservations;

        % Split the data within the patient based on observation indices
        trainIndices = 1:numTrainObservations;
        valIndices = numTrainObservations+1:numTrainObservations+numValObservations;
        testIndices = numTrainObservations+numValObservations+1:numObservations;

        % Append the split data for the current patient to the respective subsets
        trainData = [trainData; patientData(trainIndices, :)];
        valData = [valData; patientData(valIndices, :)];
        testData = [testData; patientData(testIndices, :)];
    end
end
