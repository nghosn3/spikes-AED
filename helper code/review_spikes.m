function review_spikes(spikes)
%% Parameters
review_chs = 1;

%% Locations
locations = interictal_hub_locations;
addpath(genpath(locations.script_folder));

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get basic info about pt
name = spikes.name;
nfiles = length(spikes.file);

found_p = 0;
for p = 1:length(pt)
    if strcmp(pt(p).name,name)
        found_p = 1;
        break
    end
end
if found_p == 0
    error('Did not find patient in pt structure');
end

% Loop over files
for f = 1:nfiles

    chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
    dur = diff(pt(p).ieeg.file(f).block_times(1,:))/3600;
    ekg_chs = identify_ekg(chLabels);
    non_ekg_chs = find(~ekg_chs);
    non_ekg_labels = chLabels(~ekg_chs);
    
    nchs = length(chLabels);
    nblocks = length(spikes.file(f).block);
    
    % initialize spike rate matrix (default is nans so that I can
    % distinguish between a run where I skipped detection and a run when I
    % detected no spikes)
    rate = nan(nchs,nblocks); 
    bindices = [];
    
    % Loop over blocks
    for h = 1:nblocks
        block = spikes.file(f).block(h);
        bindices = [bindices,h];
        if block.run_skip == 1
            continue;
        end
        
        gdf = spikes.file(f).block(h).gdf;
        
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
    
    %% Raster plot of rate
    hf = figure;
    set(gcf,'position',[26 8 1001 797])
    h = turn_nans_gray(rate(~ekg_chs,:));
    
    set(h,'XData',[0:size(rate,2)*dur])
    xlim([0 size(rate,2)*dur])
    xlabel('Hour')
    

    yticks(1:length(non_ekg_chs))
    yticklabels(non_ekg_labels)
    title(sprintf('%s',name))
    
    if review_chs
        while 1
            try
                [x,y] = ginput;
            catch
                break
            end
            chLab = non_ekg_labels{round(y(end))};
            fidx = f;
            bidx = bindices(round(x(end)/dur));
            fprintf('\nShowing spikes for %s ch %s file %d block %d (hour %d)\n',...
                name,chLab,fidx,bidx,round(bidx*dur));

            plot_spikes_by_ch(p,chLab,fidx,bidx,spikes)
            
        end
    end
    
        
end



end


function h = turn_nans_gray(im)
    % white
    cmap = colormap;
    nanjet = [ 0.7,0.7,0.7; cmap  ];
    nanjetLen = length(nanjet); 
    pctDataSlotStart = 2/nanjetLen;
    pctDataSlotEnd   = 1;
    pctCmRange = pctDataSlotEnd - pctDataSlotStart;

    dmin = nanmin(im(:));
    dmax = nanmax(im(:));
    dRange = dmax - dmin;   % data range, excluding NaN

    cLimRange = dRange / pctCmRange;
    cmin = dmin - (pctDataSlotStart * cLimRange);
    cmax = dmax;
    h= imagesc(im);
    set(gcf,'colormap',nanjet);
    caxis([cmin cmax]);
end

function ekg_chs = identify_ekg(chLabels)

ekg_chs = zeros(length(chLabels),1);

for i = 1:length(chLabels)
    curr = chLabels{i};
    if contains(curr,'EKG') || contains(curr,'ECG')
        ekg_chs(i) = 1;
    end
end

end