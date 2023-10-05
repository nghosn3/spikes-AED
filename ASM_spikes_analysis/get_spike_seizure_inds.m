function [sz_inds] = get_spike_seizure_inds(ptID,file_inds)

[seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID);

%file inds are for each 10 minute block - convert seizure times to a
%10minute block
sz_inds = zeros(1,length(seizure_times));
for x = 1:length(seizure_times)
    sz = round((seizure_times(x,1) ./60)./10); % convert to units of 10 minutes 
    sz_dataset = str2double(seizure_dataset{x}(end));
    if contains(seizure_dataset{x},'one file')
        sz_inds(x)=sz;
    else 
    dataset_inds = find(file_inds == sz_dataset);
    sz_inds(x)  = (dataset_inds(1)-1) + sz;
    end 
end 


end 