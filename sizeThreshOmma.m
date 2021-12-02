function [new_omma_cent] = sizeThreshOmma(omma_cent,omma_area,thresh)

%--------------------------------------------------------------------------
% screen out small centroids
%--------------------------------------------------------------------------

disp('Screening out small objects')

for i = 1:length(omma_area)
    count = 0;
    for j = 1:length(omma_area{i})
        if omma_area{i}(j) > thresh
            count = count + 1;
            new_omma_cent{i}(count,1:2) = omma_cent{i}(j,:);
        end
    end
end