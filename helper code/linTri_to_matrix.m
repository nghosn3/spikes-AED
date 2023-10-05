function [out] = linTri_to_matrix(nChannels,data,plot,elecLabels)

% data is linear upper triangular matrix
% nchannels is the number of electrode channels used
% plot is logical of whether to plot the matrix or not


%upper triangular to connectivity matrix
a = triu(ones(nChannels), 1);
a(a > 0) = data; %allFeats(1,:)
out = (a + a')./(eye(nChannels)+1);
%plot
if plot
    imagesc(out/max(out, [], 'all')); colorbar;
    %labels = patient_localization(pt_ind).labels; elecLabels = labels(targetElectrodesRegionInds);
    xticks(1:nChannels);yticks(1:nChannels)
    xticklabels(elecLabels);yticklabels(elecLabels)
end

end