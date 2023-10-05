function [] = reconstruct_bp_states(clusterCentroids,nClusters,nChannels)

nFrequencyBands = 6;
centroids = reshape(clusterCentroids, nChannels, nFrequencyBands, nClusters);
nRows = round(sqrt(nClusters));
nCols = floor(nClusters/nRows);
if nRows*nCols < nClusters
    nRows = nRows + 1;
end

for i = 1:nClusters
    subplot(nRows, nCols, i)
    imagesc((centroids(:, :, i)));
    if i == 1
        yticks(1:nChannels)
        %yticklabels(Regions)
    else
        yticks([])
    end
    xticks(1:nFrequencyBands)
    xticklabels(["Delta", "Theta", "Alpha", "Beta", "Gamma", "High Gamma"])
    xtickangle(45)
    colorbar
    
end
end