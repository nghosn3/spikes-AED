function [pps, spd] = parse_med_frequencies_v1(med_string)

pps = NaN;
spd = NaN;
is_daily = false;

%maps text numbers like one to numbers like 1
num_str = containers.Map({'one' 'two' 'three' 'four' 'five' 'six' 'seven' ...
    'eight' 'nine' 'ten', 'once', 'twice'}, {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '1', '2'});

%split the string into words
med_split = split(med_string);

%keep track of where the times and numbers are
num_inds = zeros(1,length(med_split));

for i = 1:length(med_split)
    
    %check if the number for freq is represented in text form and replace as needed
    if isKey(num_str, lower(med_split{i}))
        med_split{i} = num_str(lower(med_split{i}));
    end
    
    if ~isnan(str2double(med_split{i}))
        num_inds(i) = 1;
    elseif contains(lower(med_split{i}), 'hour') || contains(lower(med_split{i}), 'day') || contains(lower(med_split{i}), 'dai')
        num_inds(i) = 2;
        %         if contains(lower(med_split{i}), 'day') || contains(lower(med_split{i}), 'dai')
        %             is_daily = true;
    elseif contains(lower(med_split{i}), 'morning') || contains(lower(med_split{i}), 'am') || contains(lower(med_split{i}), 'pm') || contains(lower(med_split{i}), 'mid') || contains(lower(med_split{i}), 'evening') || contains(lower(med_split{i}), 'bed') || contains(lower(med_split{i}), 'night')
         num_inds(i) = 2;
    end
end



%try to parse the frequency spd
%find where times can occur
time_idx = find(num_inds==2);


if contains(lower(med_string), 'need') || contains(lower(med_string), 'seiz') || contains(lower(med_string), 'sz')
    spd = NaN;
    fprintf(strcat("As needed found: ", med_string))
else
    
    %if a time wasn't found, then assume tablets in a day
    if ~isempty(time_idx)
        %search for a XX tab phrase
        try
            tab_idx = [1 find(contains(med_split, 'tab'))];
             %for each XX tab phrase, check that it is preceeded by a number
        total_tabs = 0;
        total_freq = 0;
        for i = 1:length(tab_idx)
            if any(num_inds(tab_idx(i):tab_idx(i+1)-1))
                num_tab_idx = find(num_inds(1:tab_idx(i)-1) ==1,1,'first');
                total_tabs = total_tabs + str2double(med_split(num_tab_idx));
                total_freq = total_freq + 1;
            end
        end
        catch
            total_freq = 0;
            fprintf('no tab statement')
        end
       
        if total_freq ~= 0
            spd = total_freq;
        else
            %search for XX unit
            %accepted units: mg
            unit_idx = find(contains(med_split, 'mg'));
            total_freq = length(unit_idx);
            %for each phrase with a unit, check if it contains numbers
            %inside. https://stackoverflow.com/questions/58979040/check-if-string-contains-any-numbers-in-matlab
            %If it doesn't then check if it is preceeded by a number
%             for i = 1:length(unit_idx)
%                 if cellfun(@any,regexp(med_split(unit_idx(i)), '[0-9]'))
%                     total_freq = total_freq + 1;
%                 elseif num_inds(unit_idx(i) - 1) == 1
%                     total_freq = total_freq + 1;
%                 end
%             end
            %attempt to find phrases of the format XX tab, or XX unit, or
            %XXunit
            if total_freq ~= 0
                spd = total_freq;
            end
            
        end
        if isnan(spd)
            fprintf(strcat("Unable to parse spd for SIG ", med_string))
        end
        %else, if a number preceeds the first time phrase, then calculate the frequency
    elseif num_inds(time_idx(1) - 1) == 1
        spd = str2double(med_split{time_idx(1)-1});
        if contains(lower(med_split{i}), 'hour')
            spd = 24 / spd;
        end
    else
        fprintf(strcat("Unable to parse spd for SIG ", med_string))
    end
end



%try to get the dose pps
num_idx = find(num_inds==1);
if ~isempty(time_idx)
    if num_idx(1)+1 ~= time_idx(1)
        pps = med_split{num_idx(1)};
    end
end
end

