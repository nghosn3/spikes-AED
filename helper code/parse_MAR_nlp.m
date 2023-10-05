%% parse medication info from MAR (2/24/22) sheet 
function [med_names,meds,explant_date,implant_date] = parse_MAR_nlp(ptID,all_meds)

%pt_inds = all_meds.ptID== ptID;
pt_inds = all_meds.Record_ID== ptID;
meds=all_meds(pt_inds,:);

% load implant dates for all patients and find info for ptID
% implant_info=readtable('HUP_implant_dates.xlsx','VariableNamingRule','preserve');
% implant_info.implant_time=implant_info.implant_time*24;

explant_date = unique(meds.HOSP_DISCH_DATE); explant_date = explant_date(~isnat(explant_date));
implant_date = unique(meds.HOSP_ADMSN_DATE); implant_date = implant_date(~isnat(implant_date));

pt_implant.Implant_date = implant_date;
pt_implant.Explant_date = explant_date;

%%
EMU_med_inds=false(height(meds),1);
for i = 1:length(implant_date) % if there is more than one hosp admission
    inds= ((meds.TAKEN_DATETIME-timeofday(meds.TAKEN_DATETIME)) <= pt_implant.Explant_date(i)) & (meds.TAKEN_DATETIME-timeofday(meds.TAKEN_DATETIME)) >= pt_implant.Implant_date(i);
    EMU_med_inds = EMU_med_inds | inds;
end
EMU_meds = meds(EMU_med_inds,:);

% filter AEDs only 
emu_med_names = cellfun(@lower,EMU_meds.DISPLAY_NAME,'UniformOutput',false);
emu_med_names = cellfun(@strsplit,emu_med_names,'UniformOutput',false);
emu_med_names = cellfun(@(x)x{1},emu_med_names,'UniformOutput',false);


[gen_names] = get_AEDs_from_list(emu_med_names);
del_inds = cell2mat(cellfun(@isempty,gen_names,'UniformOutput',false));
EMU_meds(del_inds,:)=[];
EMU_meds.generic_name = gen_names(~del_inds);

% filter medications that were not given 
given = contains(EMU_meds.MAR_ACTION,'Given') | contains(EMU_meds.MAR_ACTION,'New Bag') |  contains(EMU_meds.MAR_ACTION,'Administered by Provider') | contains(EMU_meds.MAR_ACTION,'Performed') ;
EMU_meds = EMU_meds(given,:);

% now get field values for meds table
meds=table();

%meds.HUP = EMU_meds.ptID;
meds.HUP = EMU_meds.Record_ID;
meds.medication = EMU_meds.generic_name;
meds.dose = EMU_meds.DOSE;
meds.dose_unit = cellfun(@lower,EMU_meds.UNITS,'UniformOutput',false);
meds.date = EMU_meds.TAKEN_DATETIME;
meds.ActionTime = timeofday(EMU_meds.TAKEN_DATETIME);
meds.hosp_admsn = EMU_meds.HOSP_ADMSN_DATE;

source = EMU_meds.DATA_SOURCE;

% get the hr of the total EMU stay for the admin
[~,sorted_inds] = sort(meds.date);
meds = meds(sorted_inds,:);
source = source(sorted_inds);

for i=1:length(implant_date)
    inds = meds.hosp_admsn == implant_date(i);
    admin_times = hours(meds.date(inds)-implant_date(i));
    meds.admin_time(inds) = admin_times;
end

%drop duplcate entries from the same data source and re-sort
if ~isempty(meds)
    meds = unique(meds);
    [~,sorted_inds] = sortrows(meds,'date');
    meds.source =source(sorted_inds);
    med_names = unique(meds.medication);
else
    med_names = [];
end

% sources = unique(meds.source);
% sorted_meds =table();
% if ~isempty(meds)
%     for i=1:length(sources)
%     source_inds = contains(meds.source,sources(i));
%     unique_meds = unique(meds(source_inds,:));
%     unique_meds =sortrows(unique_meds,'date');
%     sorted_meds = [sorted_meds; unique_meds];
%     end 
% end
% 
% meds = sorted_meds;


%% now get ieeg times for medication - not all pts with ieeg, later add for HUP patients
% ieeg_info = readtable('HUP_ieeg_conversion.xlsx','VariableNamingRule' , 'preserve');
% pt_inds = ieeg_info.HUP == str2double(ptID(4:end));
% pt_ieeg = ieeg_info(pt_inds,:); 
% 
% %temp dates of dataset start
% t = pt_ieeg.dataset_start_ieeg;
% dataset_starts_ieeg = t-timeofday(t); %ieeg time
% day_offsets = cumsum([0; days(diff(dataset_starts_ieeg))]);
% for i = 1:height(pt_ieeg)
%     if isnat(pt_ieeg.dataset_start_EMU(i))
%         pt_ieeg.dataset_start_EMU(i) = pt_ieeg.EMU_admission_date(i)+days(day_offsets(i)) + timeofday(pt_ieeg.dataset_start_ieeg(i));
%     end
% end
% 
% each_dataset_start =NaT(height(meds),1);
% each_dataset_start_ieeg=NaT(height(meds),1);
% each_dataset_name = cell(height(meds),1);
% 
% dataset_starts= pt_ieeg.dataset_start_EMU; %emu time
% dataset_names = unique(pt_ieeg.dataset);
% for n=1:length(dataset_starts)
%     each_dataset_start_ieeg(meds.date > dataset_starts(n))=dataset_starts_ieeg(n);
%     each_dataset_start(meds.date > dataset_starts(n))=dataset_starts(n);
%     each_dataset_name(meds.date > dataset_starts(n)) = dataset_names(n);
% end
% % for the last entries, assume it is part of the last dataset; later if error occurs then admin is after end of ieeg recording
% meds.OffsetSecondsInIeeg = seconds(meds.date-each_dataset_start);
% meds.dataset=each_dataset_name; 
% meds.dataset(cellfun(@isempty,meds.dataset))={'before start of ieeg recording'};
% meds.DatasetStartIeeg = each_dataset_start_ieeg;
% 
% % other outputs
% med_names = unique(meds.medication);
% explant_date = pt_implant.Explant_Date;
% end 
