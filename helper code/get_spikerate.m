function [rate] = get_spikerate(spikes) 


    nchs = length(spikes.file.block(1).chLabels);
    nblocks = length(spikes.file(1).block);
    
    % initialize spike rate matrix (default is nans so that I can
    % distinguish between a run where I skipped detection and a run when I
    % detected no spikes)
    rate = nan(nchs,nblocks); 
    block_inds = [];
    
    % Loop over blocks
    for h = 1:nblocks
        block = spikes.file(1).block(h);
        block_inds = [block_inds,h];
        if block.run_skip == 1
            continue;
        end
        
        gdf = spikes.file(1).block(h).gdf;
        
        %% Get spike rate
        if ~isempty(gdf)
            for ich = 1:nchs
                if ismember(ich,block.bad) ||...
                        ismember(ich,block.skip.all)
                    rate(ich,h) = nan; % leave it as a nan if I skipped detections on that channel
                    continue
                end
                rate(ich,h) = sum(gdf(:,1)==ich);
            end
        else
            rate(:,h) = 0;
        end
        
    end 

end 