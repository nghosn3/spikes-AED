% returns the normalized sum of all drugs, for an average AED load. Time
% series is given from start of EMU stay to end of EMU stay (based on drug data, NOT ieeg)

function [drug_sum]= get_drug_sum(med_names,max_dur,pt_drug_curves,pt_tHr,emu_dur)

drugs =zeros(length(med_names),max_dur*60); %450 hours of EMU stay in minutes
        for i =1:length(med_names)
            drug=pt_drug_curves{i};
            drug=drug./nanmax(drug); %normalize each drug curve
            if ~isempty(drug)
                dStart = round(pt_tHr{i}(1)*60)-1;
                drugs(i,dStart+1:dStart+length(drug))=drug;
            end
        end
        drug_sum=nansum(drugs,1);
        drug_sum=drug_sum./length(med_names);
        
        % cut off drug curve to only be length of emu stay(still includes pre-Ieeg recording)
        time =emu_dur*60;%number of minutes of emu stay
        drug_sum = drug_sum(1:time);
        drug_sum(drug_sum==0) =NaN;

end