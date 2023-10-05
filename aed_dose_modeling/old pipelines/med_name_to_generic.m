%% load home medications and labs, and replace meication names with generic from container map of all AEDs
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data')


% load lab and home med data- data with med names to convert
labs=readtable('outpatient_meds_labs.xlsx','Sheet','labs');
home_meds=readtable('outpatient_meds_labs.xlsx','Sheet','home_meds');
implant_info=readtable('outpatient_meds_labs.xlsx','Sheet','implant_time');

% clean up medication name field in each data sheet:
labs.medication=cellfun(@lower,labs.medication,'UniformOutput',false);
temp_list = cellfun(@strsplit,home_meds.medication,'UniformOutput',false);
for i =1:length(temp_list)
    home_meds.medication(i) = temp_list{i}(1);
end
home_meds.medication = cellfun(@lower,home_meds.medication,'UniformOutput',false);

%convert names to generic
home_meds.medication=convert_AED_to_generic(home_meds.medication);
labs.medication = convert_AED_to_generic(labs.medication);

