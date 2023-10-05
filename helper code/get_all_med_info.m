%% function that loads and calculates all the medication (AED) data for a specified patient. need to point to
% the directories in borel
% update 9/16/21 to also find the times of alcohl admin

function [med_times,med_names,days,each_drug_dose,total_dose,alc_times,meds] = get_all_med_info(ptID)

% load the medication data
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data');


meds = readtable([ptID '_meds.xlsx']);
meds.medication = arrayfun(@lower,meds.medication);


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

%% remove the med times from the other EMU records that are not being analyzed, if any
if sum(strcmp('dataset',meds.Properties.VariableNames)) >0 && any(strcmp('D02',meds.dataset))
    inds = ~strcmp('D02',meds.dataset);
    meds(inds,:) = [];   
end 

%% retrieve the rest of the med data
%days = meds.day;

med_names = unique(meds.medication);
med_times = [meds.start]./3600; %start times in seconds from ieeg, convert to hr for plotting
days=ceil(med_times./24); %get the day of ieeg recording that the meds were given on
meds.day = days;

nan_inds = isnan(med_times); med_times(nan_inds) = 0;
%% get total dose per day, including acute drugs for breakthrough seizures
each_drug_dose = zeros(length(days),length(med_names));
for i=1:length(days)
    day_inds = find(meds.day == days(i));
    med_doses = meds.dose(day_inds);
    for x=1:length(med_names)
        med_inds = find(contains(meds.medication(day_inds),med_names(x)));
        doses = cell2mat(med_doses(med_inds)');
        if ~isempty(doses)
            doses=replace(doses,'mg',',');doses=strsplit(doses,',');
            doses = cell2mat(cellfun(@str2num,doses,'uniform',0));
            each_drug_dose(i,x)=sum(doses);
        else
            each_drug_dose(i,x)=0;
        end
    end
end

%add the doses for each admin time to the table in double
all_doses= cell2mat(meds.dose(:)');
all_doses=replace(all_doses,'mg',',');all_doses=replace(all_doses,'hr',',');all_doses=strsplit(all_doses,',');
all_doses = cell2mat(cellfun(@str2num,all_doses,'uniform',0));
meds.dose_mg(:)=all_doses';

%days=ceil(med_times./24); %get the day of ieeg recording that the meds were given on
total_dose = sum(each_drug_dose,2); %including ativan

end