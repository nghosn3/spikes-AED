%% run Erin's spike detector 

addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED')

% List of patients to include in analysis
%whichPts = [140 144 146 157 160 164 166 177 181 187 190];
whichPts = [140 144]; %need to run with two patients at a time 

new_spikes(whichPts)
