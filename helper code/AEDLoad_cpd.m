
close all;clear;
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
% Get which patients
cohort_info = readtable('AED-connectivity-cohort-list.xlsx');
ptIDs = cohort_info.PtIDs;

[all_dose_curves,all_Hr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_v2_IP(ptIDs); %EDIT!! TESTING SCRIPT
%%
for ipt = 1:length(ptIDs)
    ptID = ptIDs{ipt};
    med_names = all_med_names{ipt};
    offsets = ieeg_offset{2,ipt};
    
    % get the total AED dose over time
    drugs =zeros(length(med_names),max_dur*60); %400 hours of EMU stay in minutes
    for i =1:length(med_names)
        drug=all_dose_curves{ipt}{i};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(all_Hr{ipt}{i}(1)*60)-1;
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1)./length(med_names);
    
    drug_sum=drug_sum(1:emu_dur(ipt)*60); %emu_dur is in hours, index in minutes
  % change point - using mean and slope of data
    n_pts = 5; %optimize parameter to using sse of segments
    findchangepts(drug_sum,'Statistic','linear','MaxNumChanges',n_pts)
    
    %plot the seizure times with it
    [seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID);
    for j =1:length(seizure_times)
        % check which dataset the seizure is from, and add appropriate offset
        if ~(isequal(seizure_dataset{j},'D01') || isequal(seizure_dataset{j},'one file'))
            ind = str2double(seizure_dataset{j}(end));
            dataset_offset = offsets(ind);
            %dataset_offset = seconds(diff(ieeg_offset{3}));
            seizure_times(j,1)= (seizure_times(j,1) + dataset_offset); 
        end 
        seizure_times(j,:) = seizure_times(j,:)./60;
        xline(seizure_times(j,1),'--r','linewidth',1);hold on;
    end
end

%plot sleep wake mask?





    
    
    
    
    
%%
% [sleep_labels,sleep_times,ptIDs] = ad_label_sleep(ptIDs);
% sleep_time=sleep_times{x};
% labels = sleep_labels{x};
%%
% ERROR
% x = sleep_time(labels);
% y = 1000;
% 
% x=0:0.1:10;
% y1=exp(-x/2);
% y2=exp(-x/3);
% figure
% hold all
% plot(x,y1)
% plot(x,y2)
% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')
% hold off