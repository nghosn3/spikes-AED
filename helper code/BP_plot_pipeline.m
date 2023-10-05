%% pipeline for plotting/packaging bandpower over emu stay for comparison with AED data and for change point detection.
% averages bandpower in target regions in delta, theta, alpha, beta, high
% and low gamma. Bandpower is averaged in the soz and non-soz regions.
close all;clear;

% plot average time course or just spectrogram?
plot_avg_bp = 0;

% load updated data master file from ieeg-metadata
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/ieeg-metadata')
fname = 'DATA_MASTER.json';
fid = fopen(fname); raw = fread(fid,inf); str = char(raw');
fclose(fid);
data_master = jsondecode(str);

% load rest of data needed
cd('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA')

% load localization file for soz
load('patient_localization_final.mat');
cohort = readtable('AED-connectivity-cohort-list.xlsx');
ptIDs = cohort.PtIDs;

ptIDs ={'HUP190'};


for y = 1:length(ptIDs)
    
    %find patient in struct
    for k = 1:length(patient_localization)
        if patient_localization(k).patient == ptIDs{y}
            pt_ind = k; end
    end
    ptID = ptIDs{y};
    
    %get their seizure times
    seizure_times = get_seizure_times(ptID);
    %seizure_times_plot = seizure_times*24;%+ days(1)*24; %shift for plotting to align with EMU DAY#
    
    %get their non-soz/soz electrodes that were included in targetRegions
    cd('/Volumes/borel.seas.upenn.edu/public/USERS/pattnaik/hmm-emu-state-space/data')
    addpath([cd '/' ptID])
    load('bandpower-windows-5win.mat');
    load('target-electrodes-Regions.mat')
    
    elecs = zeros(1,length(patient_localization(pt_ind).soz));
    elecs(targetElectrodesRegionInds) =1; %for getting labels of used electrodes
    soz_elecs = patient_localization(pt_ind).soz;
    if sum(soz_elecs)<1
        soz_elecs = patient_localization(pt_ind).resect; %use resected electrodes when the soz elecs are not identified
    end
    
    soz_elecs_target = logical(soz_elecs(targetElectrodesRegionInds)); %included in calculations
    
    
    %% bandpower
    
    tSec = entireT*1e-6;
    %convert to decibels
    allFeats = 10*log10(allFeats);
    %band power in target regions
    numElecs =height(Regions);
    bands = [{'delta'} {'theta'} {'alpha'} {'beta'} {'gamma'} {'high gamma'}];
    deltas = allFeats(:,1:numElecs);
    thetas = allFeats(:,numElecs+1:numElecs*2);
    alphas = allFeats(:,numElecs*2+1:numElecs*3);
    betas = allFeats(:,numElecs*3+1:numElecs*4);
    gamma1 = allFeats(:,numElecs*4+1:numElecs*5);
    gamma2 = allFeats(:,numElecs*5+1:numElecs*6);
    
    %% average band power in nonsoz vs. soz elecs
    tHr = tSec./3600;
    
    if plot_avg_bp
        soz_deltas = mean(deltas(:,soz_elecs_target),2);
        nonsoz_deltas = mean(deltas(:,~soz_elecs_target),2);
        
        soz_thetas = mean(thetas(:,soz_elecs_target),2);
        nonsoz_thetas = mean(thetas(:,~soz_elecs_target),2);
        
        soz_alphas = mean(alphas(:,soz_elecs_target),2);
        nonsoz_alphas = mean(alphas(:,~soz_elecs_target),2);
        
        soz_betas = mean(betas(:,soz_elecs_target),2);
        nonsoz_betas = mean(betas(:,~soz_elecs_target),2);
        
        soz_gamma1 = mean(gamma1(:,soz_elecs_target),2);
        nonsoz_gamma1 = mean(gamma1(:,~soz_elecs_target),2);
        
        soz_gamma2 = mean(gamma2(:,soz_elecs_target),2);
        nonsoz_gamma2 = mean(gamma2(:,~soz_elecs_target),2);
        
        %% plot the bandpower
        figure();title(ptID)
        ylimits =[0 1e3];
        subplot(3,2,1)
        plot(tHr,(nonsoz_deltas));
        hold on;plot(tHr,(soz_deltas));
        title([ptID 'avg delta power']); legend([{'non-soz electrdoes'},{'soz electrodes'}])
        %ylim(ylimits)
        
        subplot(3,2,2)
        plot(tHr,(nonsoz_thetas));
        hold on;plot(tHr,(soz_thetas));
        title('avg theta power');
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
        %ylim(ylimits)
        
        subplot(3,2,3)
        plot(tHr,(nonsoz_alphas));
        hold on;plot(tHr,(soz_alphas));
        title('avg alpha power');
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
        %ylim(ylimits)
        
        subplot(3,2,4)
        plot(tHr,(nonsoz_betas));
        hold on;plot(tHr,(soz_betas));
        title('avg beta power');
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
        %ylim(ylimits)
        
        subplot(3,2,5)
        plot(tHr,(nonsoz_gamma1));
        hold on;plot(tHr,(soz_gamma1));
        title('avg gamma (30-50Hz) power');
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
        %ylim(ylimits)
        
        subplot(3,2,6)
        plot(tHr,(nonsoz_gamma2));
        hold on;plot(tHr,(soz_gamma2));
        title('avg gamma (50-80Hz?) power');
        
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
        %ylim(ylimits)
        
        %electrodes included in the averages
        soz_elec_labels = Regions(soz_elecs_target,:);
        nonsoz_elec_labels = Regions(~soz_elecs_target,:);
        
    else 
        figure;
        imagesc(allFeats');
        title([ptID ': bandpower in target regions'])
        %yticks(1:26);yticklabels(Regions);
        xticks(1:7200:length(tSec));
        xticklabels(round(tHr(1:7200:end)));
        c = colorbar(); c.Label.String = 'dB';
        
        hold on;
        %plot seizure times
        for j =1:length(seizure_times)
            xline(seizure_times(j,1),'r--');hold on;
        end
       
    end
end