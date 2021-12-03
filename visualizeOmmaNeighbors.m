function visualizeOmmaNeighbors(raw_images,clean_omma_centroids,...
    omma_triangles,delaunay_neighbors,time_point,save)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% go thru omma and visualize neighbors, making sure to avoid boundaries
%
% more or less a demonstration of how to index thru ommatidia while
% avoiding boundary ommatidia, which ensures that we will only be taking
% measurements of ommatidia with a full set of neighbors
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = time_point
    
    boundary_cent = unique(boundary(clean_omma_centroids{t},0.8));
    
    x = clean_omma_centroids{t}(:,1);
    y = clean_omma_centroids{t}(:,2);
    
        for j = 1:size(x,1)
            if not(ismember(j,boundary_cent))
                imshow(raw_images(:,:,:,t))
                hold on
                triplot(omma_triangles{t},'LineWidth',2,'Color','cyan')
                plot(x(j),y(j),'r*','MarkerSize',12, 'LineWidth', 4)
                plot(x(delaunay_neighbors{t}{j}),y(delaunay_neighbors{t}{j}),'r.', ...
                    'MarkerSize', 12)
                hold off
                drawnow
                if save
                    filename = ['/Users/kevin/Documents/MATLAB/forSha/visualizing_neighbors/T=' num2str(j,'%03i') '.png'];
                    print(gcf,'-dpng',filename)
                end
            end
        end
    
end