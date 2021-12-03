function [new_omma_centroids] = sizeThreshOmma(omma_centroids,omma_area,thresh)

%--------------------------------------------------------------------------
% screen out small centroids
%--------------------------------------------------------------------------

disp('Screening out small objects')

for i = 1:length(omma_area)
    count = 0;
    for j = 1:length(omma_area{i})
        if omma_area{i}(j) > thresh
            count = count + 1;
            new_omma_centroids{i}(count,1:2) = omma_centroids{i}(j,:);
        end
    end
end