%% calculate phase locking value over time between two signals 
function [plv_over_time,time_vector] = calc_plv(signal1,signal2)


% Define analysis parameters - can adjust
segment_length = 12;  % segment can be 1 hour - 6, 10minute samples
overlap = 6;         % Overlap between segments can be 50% (adjust as needed)

% Calculate the number of segments
num_segments = floor((length(signal1) - overlap) / (segment_length - overlap));

% Initialize an array to store PLV values over time
plv_over_time = zeros(1, num_segments);

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

end

% Create a time vector for the PLV values
time_vector = (segment_length - overlap) * (0:num_segments - 1);


end






