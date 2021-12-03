function [omma_triangles,clean_omma_centroids,delaunay_neighbors] ... 
    = triangulateAndFindNeighbors(omma_centroids)

for t = 1:length(omma_centroids)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delaunay triangulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % delaunay triangulation
    T = delaunay(omma_centroids{t}(:,1),omma_centroids{t}(:,2));
    TR = triangulation(T,omma_centroids{t});
    
    % find boundary
    boundary_cent = unique(boundary(omma_centroids{t}));
    boundary_cent = sort(boundary_cent,'descend');
    
    % find triangles that have vertices along the boundary and flip these 
    % to zeros in 'tempT'
    tempT = ones(size(T));
    for j = 1:size(T,1)
        for jj = 1:size(boundary_cent)
            if ismember(boundary_cent(jj),T(j,:))
                tempT(j,:) = [0 0 0];
            end
        end
    end
    % temp T is 0s and 1s, where 0s represent triangles that have
    % vertices on the edge - multiply by T to remove to flip these to 
    % zeros in T as well
    tempT = tempT .* T;
    
    % now we need to make a new matrix that doesn't have all these zero
    % entries because it won't work with the 'triangulation' function like this
    count = 0;
    new_length = find(not(tempT(:,1)));
    newT = zeros(length(new_length), 3);
    for j = 1:size(tempT,1)
        if tempT(j,1) ~= 0
            count = count + 1;
            newT(count,:) = tempT(j,:);
        end
    end
    
    % remove boundary centroids from our list
    newROIcentroids = zeros(size(omma_centroids{t},1) - length(unique(boundary_cent)),2);
    count = 0;
    count2 = 0;
    for j = 1:size(omma_centroids{t},1)
        
        % if its not a boundary centroid, then record it in our
        % new matrix
        if not(ismember(j,boundary_cent))
            count = count + 1;
            newROIcentroids(count,:) = omma_centroids{t}(j,:);
        
         end
        
    end
    
    % we need to update our triangular connectivity list according to the
    % centroids that have been removed (the boundary centroids). The
    % connectivity list refers to the rows of the centroid list. Therefore,
    % to account for the centroids we've removed, we need to subtract 1
    % from every connectivity reference GREATER than the one we just
    % removed
    for j = 1:length(boundary_cent)
        
        curr = boundary_cent(j);
        
        for i = 1:size(newT,1)
            for ii = 1:3
                if newT(i,ii) >= curr
                    newT(i,ii) = newT(i,ii) - 1;
                end
            end
        end
        
    end
    
    
    newTR = triangulation(newT,newROIcentroids);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define neighbors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x = newROIcentroids(:,1);
    y = newROIcentroids(:,2);

    % 1. Get all the triangles attached to a particular vertex in the
    % triangulation.  
    attachedTriangles = vertexAttachments(newTR);


    for i = 1: size(x,1)
        % 2. Use the connectivity list to get the vertex indices of all these
        % triangles
        verticesOfTI = newTR.ConnectivityList(attachedTriangles{i},:);
        % 3. Find all the unique vertices and remove the current vertex
        neighboursOfInternal{i} = setdiff(unique(verticesOfTI), i);
    end

    
    % store neighbor information for current time point in cell array
    % containing same info for all time points
    delaunay_neighbors{t} = neighboursOfInternal;
    
    original_delaunay{t} = T;
    edge_clean_delaunay{t} = newT;                  % definition of traingles - refers to indices of centroids in clean_omma_centroids
    original_triangulation{t} = TR;
    omma_triangles{t} = newTR;            % contains both centroid coordinates and definition of traingle connectivity
    clean_omma_centroids{t} = newROIcentroids;   % just a list of ommatidia centroids
    
    
    
%     figure(1)
%     
%     % plot triangulation
%     subplottight(1,2,1)
%     imshow(raw_images(:,:,:,t))
%     hold on
%     triplot(TR,'LineWidth',2,'Color','cyan')
%     hold off
%     
%     % after removing boundary triangles
%     subplottight(1,2,2)
%     imshow(raw_images(:,:,:,t))
%     hold on
%     triplot(newTR,'LineWidth',2,'Color','cyan')
%     hold off
%     
%     set(gcf,'PaperPosition',[0 0 (size(raw_images(:,:,1,1),2)/50)*2 (size(raw_images(:,:,1,1),1)/50)])
%     pbaspect([1 1 1])
%     filename = ['/Users/kevin/Documents/MATLAB/forSha/delaunay_edge_clean/T=' num2str(t,'%03i') '.png'];
%     print(gcf,'-dpng',filename)
    
    
end