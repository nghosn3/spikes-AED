

% using the AED BPL model
[all_dose_curves_ex,all_tHr_ex,~,all_med_names_ex,ieeg_offset_ex,max_dur_ex,emu_dur_ex] = get_aed_curve_kg(example_ptIDs,example_weights);


figure('Position', [10 10 900 500])
for ipt = 1:length(example_ptIDs)
    % patient medication administration
    ptID = ['HUP' num2str(example_ptIDs(ipt))];
    offsets = ieeg_offset_ex{2,ipt};
    ieeg_offset_datasets = ieeg_offset_ex{1,ipt};
    
    % plot dose curve
    % subplot(m,n,ipt)
    
    subplot(3,1,ipt);
    h=zeros(1,length(all_med_names_ex{ipt}));
    med_colors = lines(length(all_med_names_ex{ipt}));
    for i=1:length(all_med_names_ex{ipt})
        curve=all_dose_curves_ex{ipt}{i};
        if ~isempty(curve)
            h(i)=plot(all_tHr_ex{ipt}{i},curve,'LineWidth',2,'Color',[med_colors(i,:) .5]);hold on;
            %h(i)=plot(all_Hr{ipt}{i},curve);hold on;
        else
            all_med_names_ex{ipt}(i)=[];
            h(i)=NaN;
        end
    end
    h(isnan(h))=[];
    % plot the seizures -- in ieeg times
    
    %get seizure times
    [seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID);
    if ~isempty(offsets)
        for j =1:height(seizure_times)
            % check which dataset the seizure is from, and add appropriate offset
            if isequal(seizure_dataset{j},'D01') || isequal(seizure_dataset{j},'one file')
                seizure_times(j,1)= (offsets(1)+(seizure_times(j,1)))./3600;
            else
                %ind = str2double(seizure_dataset{j}(end));
                ind = contains(ieeg_offset_datasets,['D0' seizure_dataset{j}(end)]);
                dataset_offset = offsets(ind);
                seizure_times(j,1)= (seizure_times(j,1) + dataset_offset)./3600; %convert to hours
            end
            xline(seizure_times(j,1),'--r','linewidth',2);hold on;
        end
    end
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    
    ylabel('Normalized AED BPL','FontSize',14);xlabel('time (Hr)','FontSize',14)
    %xlim([0 time(end)]);
    if ~isempty(seizure_times)
        legend(h(:),all_med_names_ex{ipt}');
    else
        legend(all_med_names_ex{ipt}');
    end
    title(ptID,'FontSize',16);
    
end


