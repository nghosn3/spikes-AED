% function that pulls all the seizure times for a patient, and returns the
% start and end times in HOURS
function [seizure_times] = get_seizure_times(ptID)

% load info from borel:
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA')
load('patient_localization_final.mat');

addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/ieeg-metadata')
fname = 'DATA_MASTER.json';
fid = fopen(fname); raw = fread(fid,inf); str = char(raw');
fclose(fid);
data_master = jsondecode(str);

%get the seizure times
seizures = data_master.PATIENTS.(ptID).Events.Ictal;
    seizure_times =[];
    vars = fields(seizures);
    for j =1:length(vars)
        start = seizures.(vars{j}).SeizureEEC;
        try stop = seizures.(vars{j}).SeizureEnd;
        catch; stop =NaN;
        end
        seizure_times = [seizure_times;start stop];
    end
    seizure_times = seizure_times ./3600; %convert to hours

end

