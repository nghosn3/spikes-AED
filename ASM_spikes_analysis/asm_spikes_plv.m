%% create a figure that shows the 

close all;clear;
cd('/Volumes/users/nghosn3/Pioneer')
curr_path = pwd;
addpath([curr_path '/DATA']);
addpath([curr_path '/DATA/mat files']);

addpath([curr_path '/spikes-AED/aed_dose_modeling'])
addpath([curr_path '/spikes-AED/helper code'])
addpath([curr_path '/spikes-AED/ASM_spikes_analysis'])


%load spike rate - new from 2/13/23 (samp/10min)
load('spikes_rates_021323.mat');
spikes_fname = 'spikes_rates_021323.mat';
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;
load('MAR_032122.mat')
meds_fname = 'MAR_032122.mat';

% load the fooof results and 
load('fooof_results_091223.mat')

[all_plv,all_asm_plv,all_spikes_plv] = calc_plv_spikes_asm(all_spike_rate,all_pts_drug_samp,file_inds,ptIDs); % plv before first seizure

%%
r_plv = nan(length(ptIDs),1);

tiledlayout('flow')
for i = 75%1:length(ptIDs)
    plv = all_plv{i};
    nan_inds =isnan(plv);
    plv(nan_inds)=[];
    asm_load = all_asm_plv{i};
    asm_load(nan_inds)=[];

    r = corrcoef(asm_load,plv);
    r_plv(i) = r(1,2);
%     nexttile;
%     plot(plv,asm_load,'.')


end 


nexttile;
histogram(r_plv, 'FaceAlpha', 0.5,'Facecolor','black','EdgeColor', 'black', 'linewidth', 2);
set(gca,'fontsize',14);
ylabel('# patients')
xlabel('correlation');
title('PLV of ASM load and spike rate vs. ASM load')

[~,p] = ttest(r_plv)

tiledlayout('flow')
[sorted_plv,inds] = sort(r_plv);
bar(sorted_plv,'faceColor',[0.5 0.5 0.5],'facealpha',0.1, 'edgecolor',[0.5 0.5 0.5],'linewidth',2);
xticks(1:length(ptIDs))
xticklabels(ptIDs(inds))
set(gca,'box','off','fontsize',14)
ylabel('R')
xtickangle(60)
set(gca, 'FontSize', 14,'linewidth',2);
set(gca, 'box', 'off')
title('ASM load vs. PLV(asm load, spikes)','fontsize',14)




