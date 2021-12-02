function [omma_cent, omma_area] = initialSeg(ilastik_probabilities,threshold)


%--------------------------------------------------------------------------
% segment the ommatidia using a threshold
%--------------------------------------------------------------------------

disp('Thresholding ilastik probabilities to find ommatidia')

% store thresholded probabilites
segmented_ommatidia = zeros(size(ilastik_probabilities));

% loop through probability masks and threshold each one
for i = 1:size(ilastik_probabilities,3)
    
    % thresholding ilastik probability masks
    segmented_ommatidia(:,:,i) = im2bw(ilastik_probabilities(:,:,i),threshold);
    
end


%--------------------------------------------------------------------------
% find centroids of segmented ommatidia
%--------------------------------------------------------------------------

disp('Calculating ommatidia centroids')

labeled_ommatidia = zeros(size(segmented_ommatidia));

% bwlabel segmented ommatidia
for i = 1:size(labeled_ommatidia,3)
    labeled_ommatidia(:,:,i) = bwlabel(segmented_ommatidia(:,:,i));
end

% create cell array to hold centroid positions - each cell contains the
% centroids for one of the images
omma_cent = cell(size(labeled_ommatidia,3),1);

% compute centroid positions and screen out small ones
for i = 1:length(omma_cent)
    temp_centroid = regionprops(labeled_ommatidia(:,:,i),'Centroid');
    temp_area = regionprops(labeled_ommatidia(:,:,i),'Area');
    omma_cent{i} = round(cat(1,temp_centroid.Centroid));
    omma_area{i} = round(cat(1,temp_area.Area));
end



% for j = 1:size(omma_cent_disp,3)
%     omma_cent_disp(:,:,j) = zeros(size(omma_cent_disp(:,:,j)));
%     for jj = 1:size(omma_cent{j})
%         omma_cent_disp(omma_cent{j}(jj,2),omma_cent{j}(jj,1),j) = 1;
%     end
% end
% 
% dilate_omma_cent_disp = zeros(size(omma_cent_disp));
% se = strel('disk',2);
% for t = 1:size(omma_cent_disp,3)
%     dilate_omma_cent_disp(:,:,t) = imdilate(omma_cent_disp(:,:,t),se);
% end
