
% load the medication data
%cd( '/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer') % '/Volumes/g/public/USERS/nghosn3'
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA/med-data');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')

%%
all_labs = readtable('full_report_DI.xlsx','Sheet','Labs','VariableNamingRule', 'preserve');
save('all_labs_042622.mat','all_labs')