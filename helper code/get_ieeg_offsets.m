%% get ieeg offsets and dataset info for patients in ptIDs

function [all_ieeg_offset] = get_ieeg_offsets(ptIDs,all_meds)

all_ieeg_offset=cell(3,length(ptIDs));
for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID = ['HUP' num2str(ptIDs(ipt))];
    disp(ptID);
    [~,meds,~] = parse_MAR(ptID,all_meds);    
    
    %get offset in ieeg
    eeg_diff = (meds.admin_time*3600 - meds.OffsetSecondsInIeeg);
    eeg_round = unique(round(eeg_diff*10^5)/10^5);
    eeg_offsets = eeg_round(~isnan(eeg_round));
    
    eeg_datasets = unique(meds.dataset);
    eeg_datasets = eeg_datasets(contains(eeg_datasets,'D'));
    all_ieeg_offset(1,ipt) = {eeg_datasets};
    all_ieeg_offset(2,ipt) = {eeg_offsets};
    
    datasetStarts = unique(meds.DatasetStartIeeg);
    datasetStarts = datasetStarts(~isnat(unique(meds.DatasetStartIeeg)));
    all_ieeg_offset(3,ipt) = {datasetStarts};
    
end
end