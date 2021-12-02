function [omma_cent] = mergeCloseOmma(omma_cent,kernal)

%--------------------------------------------------------------------------
% merge close centroids
%--------------------------------------------------------------------------

% create binary image with centroids painted in
omma_cent_disp = zeros(size(omma_cent,3));
for j = 1:size(omma_cent,3)
    for jj = 1:size(omma_cent{j},1)
        omma_cent_disp(omma_cent{j}(jj,2),omma_cent{j}(jj,1),j) = 1;
    end
end

disp('Merging close objects')

% merge close centroids via dilation and recomputing centroid
se = strel('disk',kernal);
clear omma_cent
for j = 1:size(omma_cent_disp,3)
    
    curr = omma_cent_disp(:,:,j);
    curr = imdilate(curr,se);
    curr2 = bwlabel(curr);
    temp_centroid = regionprops(curr2,'Centroid');
    omma_cent{j} = round(cat(1,temp_centroid.Centroid));
    
end