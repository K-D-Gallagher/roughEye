function [ROI_points] = customROI(raw_images,omma_cent)

for j = 1:size(raw_images,4)
    omma_cent_disp(:,:,j) = zeros(size(raw_images(:,:,:,j),1),size(raw_images(:,:,:,j),2));
    for jj = 1:size(omma_cent{j})
        omma_cent_disp(omma_cent{j}(jj,2),omma_cent{j}(jj,1),j) = 1;
    end
end

ROI = zeros(size(omma_cent_disp));

for t = 1:size(raw_images,4)
    
    % display image
    imshow(raw_images(:,:,:,t))
    hold on
    for ii = 1:length(omma_cent{t})
        plot(omma_cent{t}(ii,1),omma_cent{t}(ii,2),'om','Linewidth',2,'MarkerSize',3)
    end
    hold off
    
    % prompt user whether or not to draw custom ROI for this image
    answer = questdlg('Define ROI?', ...
	'Draw Custom ROI', ...
	'Yes','No','No');
    
    tf = strcmp( answer, 'No' );
    
    % if no, skip code below and go back to the top
    if tf
        continue 
    end
    
    % draw ROI and save coordinates
    p = drawpolygon('LineWidth',7,'Color','cyan');
    ROI_points{t} = p.Position;
    
    % create binary mask of ROI
    bw = double(poly2mask(ROI_points{t}(:,1),ROI_points{t}(:,2),size(raw_images,1),size(raw_images,2)));
    
    % apply mask to centroid map
    ROI(:,:,t) = omma_cent_disp(:,:,t) .* bw;
    
    % find new centroid positions
    % compute centroid positions
    [x,y] = find(ROI(:,:,t));
    temp_centroid = [y, x];
    ROI_centroids{t} = temp_centroid;
    
    % display new set of centroids within the ROI
    % display image
    imshow(raw_images(:,:,:,t))
    hold on
    for ii = 1:length(ROI_centroids{t})
        plot(ROI_centroids{t}(ii,1),ROI_centroids{t}(ii,2),...
            'om','Linewidth',2,'MarkerSize',3)
    end
    plot([ROI_points{t}(:,1);ROI_points{t}(1,1)], ...
        [ROI_points{t}(:,2);ROI_points{t}(1,2)], ...
        'LineWidth',7,'Color','cyan')
    hold off
    drawnow
    
    
    
end

