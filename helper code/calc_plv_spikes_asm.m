
%% calculate phase locking value over time between ASM load and Spikes before the first seizure
function [all_plv,all_asm_plv,all_spikes_plv] = calc_plv_spikes_asm(all_spike_rate,all_pts_drug_samp,file_inds,ptIDs)


all_plv = cell(1,length(ptIDs));
all_asm_plv = cell(1,length(ptIDs));
all_spikes_plv = cell(1,length(ptIDs));

for ipt = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];

    % get spikes and ASM signals
    spikes = all_spike_rate{ipt};
    asm_load = sum(all_pts_drug_samp{ipt});
   

    %only use data up to first sz -  check if enough data before first seizure, otherwise exclude.
    sz_inds = get_spike_seizure_inds(ptID,file_inds{ipt});
    time_thresh = 6*24*1; % two days of data
    sz_inds = length(asm_load);
    if sz_inds(1) > time_thresh
        asm_load = asm_load(1:sz_inds(1));
        spikes = spikes(1:sz_inds(1));
        signal1 = spikes;
        signal2 = asm_load;

        % Define analysis parameters - can adjust
        segment_length = 12;  % segment can be 1 hour - 6, 10minute samples
        overlap = 6;         % Overlap between segments can be 50% (adjust as needed)

        % Calculate the number of segments
        num_segments = floor((length(signal1) - overlap) / (segment_length - overlap));

        % Initialize an array to store PLV values over time
        plv_over_time = zeros(1, num_segments);
        asm_over_time = zeros(1, num_segments);
        spikes_over_time = zeros(1, num_segments);

        % Iterate through the segments
        for i = 1:num_segments
            % Define the start and end indices for the current segment
            start_idx = (i - 1) * (segment_length - overlap) + 1;
            end_idx = start_idx + segment_length - 1;

            % Extract the data for the current segment
            segment_signal1 = signal1(start_idx:end_idx);
            segment_signal2 = signal2(start_idx:end_idx);

            % Calculate PLV for the segment (as described in the previous response)
            signal1_analytic = hilbert(segment_signal1);
            signal2_analytic = hilbert(segment_signal2);
            phase_difference = angle(signal1_analytic ./ signal2_analytic);
            plv_over_time(i) = abs(mean(exp(1i * phase_difference)));
            asm_over_time(i) = mean(signal2(start_idx:end_idx));
            spikes_over_time(i) = mean(signal1(start_idx:end_idx));
        end

        % Create a time vector for the PLV values
        time_vector = (segment_length - overlap) * (0:num_segments - 1);

%         % Plot the PLV values over time
%         plot(time_vector, plv_over_time);
%         hold on;
%         xlabel('Time (ms)');
%         ylabel('PLV');
%         title('Phase Locking Value Over Time');

        all_plv{ipt} = plv_over_time;
        all_asm_plv{ipt} = asm_over_time;
        all_spikes_plv{ipt} = spikes_over_time;
    end


end

end 