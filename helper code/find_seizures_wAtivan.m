%% identify which seizures had ativan admistered after 

close all;clear;
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');

% Get which patients have AED data, load the data
load('MAR_032122.mat')
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;


% Get the patients that also have seizure severity scores
seizure_metadata = readtable('seizure_metadata_with_severity.xlsx',"ReadVariableNames", true);
sz_pts = unique(seizure_metadata.Patient);

% patients included in analysis
pt_inds=zeros(1,length(ptIDs));
for x = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(x))];
    pt_inds(x) = sum(contains(sz_pts,ptID));
end

ptIDs = ptIDs(logical(pt_inds));
[all_ieeg_offset] = get_ieeg_offsets(ptIDs,all_meds);
%% find the seizures followed by ativan administration
all_times = [];
seizure_offsets = [];
all_pt_inds =[];
for ipt=1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(ipt))];    
    
    pt_inds = find(contains(seizure_metadata.Patient,ptID));
    pt_szData = seizure_metadata(pt_inds,:);
    seizure_times = pt_szData.SeizureEEC; %in seconds
    
    % get the datset name for the offset
    dataset=cellfun(@(x) strsplit(x, '_'),pt_szData.iEEGFilename,'UniformOutput',false);
    dataset_name = cell(length(dataset),1);
    for i =1:length(dataset)
        dataset_name(i) = dataset{i}(end);
    end
    dataset_name(~contains(dataset_name,'D'))={'one file'};
    seizure_dataset = dataset_name;
    
    offsets = all_ieeg_offset{2,ipt};
    all_ieeg_offset_datasets=all_ieeg_offset{1,ipt};
    these_offsets =[];
    for j =1:length(seizure_times)
        % check which dataset the seizure is from, and add appropriate offset
        if isequal(seizure_dataset{j},'D01') || isequal(seizure_dataset{j},'one file')
            dataset_offset = offsets(1);
            seizure_times(j,1)= (seizure_times(j,1) + dataset_offset)./3600; %hours
        else
            ind = contains(all_ieeg_offset_datasets,['D0' seizure_dataset{j}(end)]);
            dataset_offset = offsets(ind);
            seizure_times(j,1)= (seizure_times(j,1) + dataset_offset)./3600; %hours
        end
        these_offsets = [these_offsets dataset_offset];
        
    end
    
    seizure_offsets = [seizure_offsets these_offsets];
    %Get ativan times
    [~,meds,~] = parse_MAR(ptID,all_meds);
    ativan_inds = strcmp(meds.medication,'lorazepam');
    times = meds.admin_time(ativan_inds);
    
    
    time_to_closest_sz = nan(length(seizure_times(:,1)),2);
    for n=1:length(seizure_times(:,1))
        sz_diffs = seizure_times(n,1)-times;
        before_ativan = sz_diffs< 0;
        if ~isempty(times)
            if ~isempty(before_ativan) && ~(length(before_ativan)==1 || sum(before_ativan)==0)
                [mval,~]=min(abs(sz_diffs(before_ativan)));
                ind = find(abs(sz_diffs)==mval);
                time_to_closest_sz(n,:)=[times(ind(1)) seizure_times(n,1)];
            end
        end
    end
    all_times=[all_times; time_to_closest_sz];
    all_pt_inds=[all_pt_inds pt_inds'];
end

all_sz_data = seizure_metadata(all_pt_inds,:);
all_sz_data.closest_ativan = all_times(:,1);
all_sz_data.seizure_time_emu = all_times(:,2);

%%
save('seizure_meta_data_with_ativan_times','all_sz_data')