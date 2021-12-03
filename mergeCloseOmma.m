function [omma_centroids] = mergeCloseOmma(omma_centroids,kernal)

%--------------------------------------------------------------------------
% merge close centroids
%--------------------------------------------------------------------------

% create binary image with centroids painted in
omma_cent_disp = zeros(length(omma_centroids));
for j = 1:length(omma_centroids)
    for jj = 1:size(omma_centroids{j},1)
        omma_cent_disp(omma_centroids{j}(jj,2),omma_centroids{j}(jj,1),j) = 1;
    end
end

disp('Merging close objects')

% merge close centroids via dilation and recomputing centroid
se = strel('disk',kernal);
clear omma_centroids
for j = 1:size(omma_cent_disp,3)
    curr = omma_cent_disp(:,:,j);
    curr = imdilate(curr,se);
    curr2 = bwlabel(curr);
    temp_centroid = regionprops(curr2,'Centroid');
    omma_centroids{j} = round(cat(1,temp_centroid.Centroid));
    
end