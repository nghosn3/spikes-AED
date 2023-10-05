%% function that loads and calculates all the medication (AED) data for a specified patient. need to point to
% the directories in borel
% update 9/16/21 to also find the times of alcohol admin
% function that uses the master spreadsheet - which will be updated when
% more data is recieved.

function [med_times,med_names,alc_times,meds] = get_med_info_master(ptID)

% load the medication data
%cd( '/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer') % '/Volumes/g/public/USERS/nghosn3'
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data');

% get patient medication information from master spreadsheet
all_meds = readtable('Pioneer_InPatientMedications.xlsx');
pt_inds = all_meds.HUP== str2double(ptID(4:end));
meds=all_meds(pt_inds,:);

% convert med_names to generic name
meds.medication = arrayfun(@lower,meds.medication);
temp_list = cellfun(@strsplit,meds.medication,'UniformOutput',false);
temp_names = cell(1,length(meds.medication));
for i =1:length(temp_list)
    temp_names{i} = temp_list{i}{1};
end
meds.medication = convert_AED_to_generic(temp_names);

% dose in mg to dose as double
meds.Dose= replace(meds.Dose,'mg','');
meds.Dose = cellfun(@str2double,meds.Dose);

%% find the alcohol times, if any and remove from the meds table
alc_times =table; 
alc_inds = zeros(1,height(meds));
ind =1;
for i =1:height(meds)
    if contains(meds.medication{i},'wine') || contains(meds.medication{i},'vodka')
        alc_times(ind,:) = meds(i,:); 
        alc_inds(ind) = i;
        ind = ind +1;
    end
end
alc_inds(alc_inds==0)=[];
% remove the entries from the meds
meds(alc_inds,:) = [];

%% retrieve the rest of the med data
hours_admin = floor(meds.ActionTime./100) + mod(meds.ActionTime,100)./60;
meds.admin_time =hours(meds.Date-meds.Date(1)) + hours_admin;


med_names = unique(meds.medication);
med_times = meds.admin_time;


end