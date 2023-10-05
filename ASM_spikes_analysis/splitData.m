function [trainData, valData, testData] = splitData(data, trainRatio, valRatio, testRatio)
    % Unique patient IDs
    uniquePatients = unique(data.ptid);

    % Calculate the number of patients for each split
    numPatients = numel(uniquePatients);
    numTrainPatients = round(trainRatio * numPatients);
    numValPatients = round(valRatio * numPatients);
    numTestPatients = numPatients - numTrainPatients - numValPatients;

    % Split the patient IDs into train, validation, and test sets
    trainPatients = uniquePatients(1:numTrainPatients);
    valPatients = uniquePatients(numTrainPatients+1:numTrainPatients+numValPatients);
    testPatients = uniquePatients(numTrainPatients+numValPatients+1:end);

    % Split the data based on patient IDs while preserving the order of observations
    trainData = data(ismember(data.ptid, trainPatients), :);
    valData = data(ismember(data.ptid, valPatients), :);
    testData = data(ismember(data.ptid, testPatients), :);
end