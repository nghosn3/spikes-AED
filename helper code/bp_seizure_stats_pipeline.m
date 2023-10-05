%% bandpower over EMU stay relative to blood plasma model: phase plots and statistics
close all;clear;
base_path = '/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/';
cd(base_path)
addpath([pwd '/spikes-AED'])
addpath([pwd '/DATA']);

% load cohort list
load('ptIDs_10_13.mat');

%load patient localization file
load('patient_localization_final.mat');
fname = 'DATA_MASTER.json';
fid = fopen(fname); raw = fread(fid,inf); str = char(raw');
fclose(fid);
data_master = jsondecode(str);

bandpower_pts = cell(1,length(ptIDs));
all_t = cell(1,length(ptIDs));
for x =1:length(ptIDs)
    tic
    ptID = ptIDs{x};
    
    %load band power
    %get their non-soz/soz electrodes that were included in targetRegions
    cd('/Volumes/borel.seas.upenn.edu/public/USERS/pattnaik/hmm-emu-state-space/data')
    addpath([cd '/' ptID])
    load('bandpower-windows-5win.mat');
    load('target-electrodes-regions.mat')
    
    %find patient in struct
    for k = 1:length(patient_localization)
        if patient_localization(k).patient == ptIDs{x}
            pt_ind = k; end
    end
    
    elecs = zeros(1,length(patient_localization(pt_ind).soz));
    elecs(targetElectrodesRegionInds) =1; %for getting labels of used electrodes
    soz_elecs = patient_localization(pt_ind).soz;
    if sum(soz_elecs)<1
        soz_elecs = patient_localization(pt_ind).resect; %use resected electrodes when the soz elecs are not identified
    end
    
    soz_elecs_target = logical(soz_elecs(targetElectrodesRegionInds)); %included in calculations
    
    %band power in target regions
    numElecs =height(Regions);
    band_names = [{'delta'} {'theta'} {'alpha'} {'beta'} {'gamma'} {'high gamma'}];
    allFeats = real(10*log10(allFeats));
    deltas = mean(allFeats(:,1:numElecs),2);
    thetas = mean(allFeats(:,numElecs+1:numElecs*2),2);
    alphas = mean(allFeats(:,numElecs*2+1:numElecs*3),2);
    betas = mean(allFeats(:,numElecs*3+1:numElecs*4),2);
    gamma1 = mean(allFeats(:,numElecs*4+1:numElecs*5),2);
    gamma2 = mean(allFeats(:,numElecs*5+1:numElecs*6),2);
    bands = [deltas'; thetas'; alphas'; betas'; gamma1'; gamma2';];
    
    bandpower_pts{x}=bands;
    all_t{x} =  (entireT*1e-6)./3600; %time in hours
    toc
end


%% BAND POWER PHASE PLOTS?
% get medication curves
load('all_bp.mat')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')
[all_dose_curves,all_Hr,ptIDs,all_med_names] = get_aed_curve(ptIDs);
%%
band_names = [{'delta'} {'theta'} {'alpha'} {'beta'} {'gamma'} {'high gamma'}];

b = 1; %freq bands
figure;

for ipt = 1:length(ptIDs)
    bp = bandpower_pts{ipt}(b,:);
    time = all_t{ipt};
    
    % get total aed dose again
    pad =max(cellfun('length',all_dose_curves{ipt}));
    drug_sum =zeros(1,pad);
    for i =1:length(all_med_names{ipt})
        drug = [all_dose_curves{ipt}{i} zeros(1,pad-length(all_dose_curves{ipt}{i}))];
        drug = drug./max(drug); %normalize each drug curve
        drug_sum = drug_sum+drug;
    end
    
    %cut off drug curve to only be length of emu stay
    tmax =round(time(end)*60); %number of minutes of emu stay ( of bp data)
    drug_sum = drug_sum(1:tmax);
    drug_sum_plot = drug_sum./length(all_med_names{ipt});
    
    % resample data to match
    bInds = floor(linspace(1,length(bp),length(drug_sum)));
    bp=bp(bInds);
    time = time(bInds);
    
    %convert seizure times to inds to get bandpower value and drug value
    seizure_times = get_seizure_times(ptIDs{ipt});
    seizure_inds = zeros(length(seizure_times),2);
    for n = 1:length(seizure_times)
        [~,si]=min(abs(time-seizure_times(n,1)));
        seizure_inds(n,1)=mean(bp([si-2,si-1,si]));
        seizure_inds(n,2)=drug_sum(si)./length(all_med_names{ipt});
    end
    
    
    subplot(3,4,ipt)
    c = colormap(gray(length(drug_sum_plot)));
    for i = 3:length(drug_sum_plot)-2
        plot(nanmean(bp([i-2,i-1,i])),drug_sum_plot(i),'.','Color',c(i,:));hold on;
    end
    plot(seizure_inds(:,1),seizure_inds(:,2),'or','LineWidth',2); axis square;
    xlabel([band_names{b} ' power (db)']);ylabel('normalized AED BPL');
    title(ptIDs{ipt}); ylim([0 1]);
    
    
end


