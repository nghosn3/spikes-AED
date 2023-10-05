close all;clear;

addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer')
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA')

load('patient_localization_final.mat');

fname = 'DATA_MASTER.json';
fid = fopen(fname); raw = fread(fid,inf); str = char(raw');
fclose(fid);
data_master = jsondecode(str);

ptIDs = {'HUP173','HUP166','HUP157','HUP160','HUP181','HUP190','HUP144','HUP146','HUP140','HUP164','HUP177','HUP187'};
%show plots of band power for all regions, or for (non)soz, or none
plot_bp=0;
plot_bp_soz=1;


%%
for y = 1:length(ptIDs)
    
    %find patient in struct
    for k = 1:length(patient_localization)
        if patient_localization(k).patient == ptIDs{y}
            pt_ind = k; end
    end
    ptID = ptIDs{y};
    
    
    %% plot medication dosage over EMU stay to define wean curve
    % load the medication data
    [med_times,med_names,days,each_drug_dose,total_dose,alc_times] = get_med_info(ptID);
    
    
    %get their seizure times
    seizure_times = get_seizure_times(ptID);
    seizure_times_plot = seizure_times+ days(1)*24; %shift for plotting to align with EMU DAY#
    
    %get their non-soz/soz electrodes that were included in targetRegions
    cd('/Volumes/borel.seas.upenn.edu/public/USERS/pattnaik/hmm-emu-state-space/data')
    addpath([cd '/' ptID])
    load('bandpower-windows-5win.mat');
    load('target-electrodes-regions.mat')
    
    
    elecs = zeros(1,length(patient_localization(pt_ind).soz));
    elecs(targetElectrodesRegionInds) =1; %for getting labels of used electrodes
    soz_elecs = patient_localization(pt_ind).soz;
    if sum(soz_elecs)<1
        soz_elecs = patient_localization(pt_ind).resect; %use resected electrodes when the soz elecs are not identified
    end
    
    soz_elecs_target = logical(soz_elecs(targetElectrodesRegionInds)); %included in calculations
    
    
    %% bandpower and AED administration
    if plot_bp_soz
        tSec = entireT*1e-6;
        
        rem_times = []; %artifact times
        % Remove artifacts
        for i = 1:height(rem_times)
            t1 = rem_times(i,1); t2 = rem_times(i,2);
            [~,ind1] = min(abs(tSec-t1)); [~,ind2] = min(abs(tSec-t2));
            % replace with NaN
            allFeats(ind1:ind2,:)=NaN;
        end
        
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
        movAvg = dsp.MovingAverage(20); %for 5 second windows, this is a one minute moving average
        tHr = tSec./3600;
        
        soz_deltas = mean(deltas(:,soz_elecs_target),2);
        nonsoz_deltas = mean(deltas(:,~soz_elecs_target),2);
        %
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
        plot(tHr,movAvg(nonsoz_deltas));
        hold on;plot(tHr,movAvg(soz_deltas));
        title([ptID 'avg delta power']); legend([{'non-soz electrdoes'},{'soz electrodes'}])
        ylim(ylimits)
        
        subplot(3,2,2)
        plot(tHr,movAvg(nonsoz_thetas));
        hold on;plot(tHr,movAvg(soz_thetas));
        title('avg theta power');
        %plot seizure times
        for j =1:length(seizure_times_plot)
            xline(seizure_times_plot(j,1),'r');hold on;
        end
        ylim(ylimits)
        
        subplot(3,2,3)
        plot(tHr,movAvg(nonsoz_alphas));
        hold on;plot(tHr,movAvg(soz_alphas));
        title('avg alpha power');
        %plot seizure times
        for j =1:length(seizure_times_plot)
            xline(seizure_times_plot(j,1),'r');hold on;
        end
        ylim(ylimits)
        
        subplot(3,2,4)
        plot(tHr,movAvg(nonsoz_betas));
        hold on;plot(tHr,movAvg(soz_betas));
        title('avg beta power'); 
        %plot seizure times
        for j =1:length(seizure_times_plot)
            xline(seizure_times_plot(j,1),'r');hold on;
        end
        ylim(ylimits)
        
        subplot(3,2,5)
        plot(tHr,movAvg(nonsoz_gamma1));
        hold on;plot(tHr,movAvg(soz_gamma1));
        title('avg gamma (30-50Hz) power'); 
        %plot seizure times
        for j =1:length(seizure_times_plot)
            xline(seizure_times_plot(j,1),'r');hold on;
        end
        ylim(ylimits)
        
        subplot(3,2,6)
        plot(tHr,movAvg(nonsoz_gamma2));
        hold on;plot(tHr,movAvg(soz_gamma2));
        title('avg gamma (50-80Hz?) power'); 
        %plot seizure times
        for j =1:length(seizure_times_plot)
            xline(seizure_times_plot(j,1),'r');hold on;
        end
        ylim(ylimits)
        
    end
    
    %% plot band power for all channels
    
    if plot_bp == 1
        figure();
        title(ptID)
        subplot(3,2,1);
        tHr = tSec./3600;
        plot(tHr,movAvg(deltas)); hold on;
        for i =1:length(med_times), xline(med_times(i),'r'); end
        ylim(ylimits)
        ylabel('power (dB)');xlabel('time (hr)');title('delta power')
        
        subplot(3,2,2);
        tHr = tSec./3600;
        plot(tHr,thetas); hold on;
        for i =1:length(med_times), xline(med_times(i),'r'); end
        ylim(ylimits)
        ylabel('power (dB)');xlabel('time (hr)');title('theta power')
        
        subplot(3,2,3);
        tHr = tSec./3600;
        plot(tHr,alphas); hold on;
        for i =1:length(med_times), xline(med_times(i),'r'); end
        ylim(ylimits)
        ylabel('power (dB)');xlabel('time (hr)');title('Alpha power')
        
        subplot(3,2,4);
        tHr = tSec./3600;
        plot(tHr,betas); hold on;
        for i =1:length(med_times), xline(med_times(i),'r'); end
        ylim(ylimits)
        ylabel('power (dB)');xlabel('time (hr)');title('beta power')
        
        subplot(3,2,5);
        tHr = tSec./3600;
        plot(tHr,gamma1); hold on;
        for i =1:length(med_times), xline(med_times(i),'r'); end
        ylim(ylimits)
        ylabel('power (dB)');xlabel('time (hr)');title('gamma power')
        
        subplot(3,2,6);
        tHr = tSec./3600;
        plot(tHr,gamma2); hold on;
        for i =1:length(med_times), xline(med_times(i),'k'); end
        ylim(ylimits); legend(Regions);
        ylabel('power (dB)');xlabel('time (hr)');title('high gamma power')
        
    end
    
end
