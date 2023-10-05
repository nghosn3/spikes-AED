%% eeg data pipeline

addpath(genpath('/Volumes/users/nghosn3/tools/ieeg-matlab-1.14.49'));
addpath('/Volumes/users/nghosn3/Pioneer/spikes-AED/ASM_eeg_analysis');
addpath(genpath('/Volumes/users/nghosn3/Pioneer/DATA'));

for_borel=1;
if for_borel
    addpath(genpath( '/mnt/leif/littlab/users/nghosn3/tools/ieeg-matlab-1.14.49'));
    addpath( '/mnt/leif/littlab/users/nghosn3/Pioneer/spikes-AED/ASM_eeg_analysis');
    addpath(genpath( '/mnt/leif/littlab/users/nghosn3/Pioneer/DATA'));

end

HUP_eeg_info = readtable('HUP_ieeg_conversion.xlsx');
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;

login_name = 'nghosn3';
pwfile = 'ngh_ieeglogin.bin';
extras =1;

frequencyBands = [0.5, 4; 4, 8; 8, 12; 12, 30; 30, 100];

for ipt =2:length(ptIDs)
    % strart the session to see data duration
    pt_file_inds = find(HUP_eeg_info.HUP == ptIDs(ipt));
   
    pt_data = table();
    pt_data.theta = zeros(0);
    pt_data.delta = zeros(0);
    pt_data.alpha = zeros(0);
    pt_data.beta = zeros(0);
    pt_data.gamma = zeros(0);
    pt_data.times = zeros(0);
    pt_data.file = zeros(0);


    for f = 1:length(pt_file_inds)
        fname = HUP_eeg_info.dataset(pt_file_inds(f));
        f_ind = str2double(fname{:}(end));
        if ~isnumeric(f_ind) || isnan(f_ind)
            f_ind = 1;
        end 
        
        session = IEEGSession(fname, login_name, pwfile);
        duration = session.data.rawChannels(1).get_tsdetails.getDuration/1e6; %in seconds
        fs = session.data.sampleRate;
        tvec = linspace(1,(duration),duration*fs);

        % calculate the number of pulls - should grab data in chunks of 2hrs
        tstep = 3600*.5; %two hour pulls
        runs = 1:tstep:(duration);

        for t = 1:length(runs)-1

            run_times = [runs(t) runs(t+1)]; % want all of the eeg data, so can change function to download all
            data = download_ieeg_data(fname, login_name, pwfile, run_times, extras);
            fs = data.fs;
            %data_vals = data.values;
            data_vals = nanmean(data.values,2);% the mean across all channels

            % filter the data
            filterOrder = 4; % Filter order (adjust as needed)
            nyquist = fs / 2;
            freqRange = [0.5 250];
            [b, a] = butter(filterOrder, [freqRange(1)/nyquist, freqRange(2)/nyquist], 'bandpass');

            % take care of NaN and 0 values
            % Replace NaN values with mean of data
            data_vals(isnan(data_vals)) = nanmean(data_vals);

            % Replace empty values with mean of data
            data_vals(isempty(data_vals)) = nanmean(data_vals);

            % Apply forward and backward filtering
            filteredSignal = filtfilt(b, a, data_vals);

            % Sliding window parameters
            windowSize = 1; % Window size in seconds
            windowSamples = windowSize * fs;
            overlap = 0; % Overlap between consecutive windows (50%)

            % Calculate number of windows
            numWindows = floor((length(filteredSignal) - windowSamples) / (windowSamples * (1 - overlap))) + 1;

            % Initialize bandpower values matrix
            bandpowerValues = zeros(size(frequencyBands, 1), numWindows);

            % Perform sliding window analysis
            for i = 1:numWindows
                startIdx = 1 + floor((i - 1) * windowSamples * (1 - overlap));
                endIdx = startIdx + windowSamples - 1;

                % Extract windowed signal
                windowedSignal = filteredSignal(startIdx:endIdx);

                % Calculate bandpower for each frequency band
                for j = 1:size(frequencyBands, 1)

                    bandpowerValues(j, i) = bandpower(windowedSignal, fs,[frequencyBands(j,:)]);
                end

            end


        
    
            % append the bandpower values for that run to the patient data
            % table.
            tnew= table();
            tnew.theta = bandpowerValues(1,:)';
            tnew.delta = bandpowerValues(2,:)';
            tnew.alpha = bandpowerValues(3,:)';
            tnew.beta = bandpowerValues(4,:)';
            tnew.gamma =  bandpowerValues(5,:)';
            tnew.times = linspace(run_times(1),run_times(2),numWindows)'; % an estimate of the time vector 
            tnew.file = ones(numWindows,1) * f_ind;
            
            pt_data = [pt_data; tnew];

            disp(['run ' num2str(t) ' out of ' num2str(length(runs)) ' completed for bandpower calculation'])
        end

    end
    
    pt_data_file_name = ['HUP' num2str(ptIDs(ipt)) '_bandpower.mat'];
    save_path = '/mnt/leif/littlab/users/nghosn3/Pioneer/DATA/bandpower_060623/';
    save([save_path pt_data_file_name],'pt_data')

end



% %% plot some data - NOTE: data.duration is duration of entire file, in seconds
% plot_chans = length(data.chLabels)-11:length(data.chLabels);
% plot_len = 60*data.fs;
% plot_start = 120*data.fs;
% data_plot = data.values(plot_start:plot_start+plot_len,plot_chans);
% tvec = linspace(plot_start,plot_start+plot_len,length(data_plot));
% spacing = ones(length(plot_chans),1) * 1:length(plot_chans);
% plot(tvec,zscore(data_plot')+spacing')