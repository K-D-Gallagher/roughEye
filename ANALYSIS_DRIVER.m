
% segment ommatidia from brightfield images of Drosophila adult eyes,
% define ommatidial lattice topology, and quantify roughness phenotype
% K-D-Gallagher https://github.com/K-D-Gallagher 2021

% Rough Eye Direct edge recovery system
% Rough-EyeDERS


%%
%--------------------------------------------------------------------------
% read in data
%--------------------------------------------------------------------------

% Here, we will read in our raw data and the corresponding pixel
% classification files exported from Ilastik.

% NOTE: the ilastik files (.h5 extension) should remain in the same folder
% as the raw images. This is the default behavior from Ilastik, so this
% should not be an issue.

filepath = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/cropped_images_768/';
[raw_images, ilastik_probabilities] = loadData(filepath);



%%
%--------------------------------------------------------------------------
% show ilastik classification probabilities on top of raw image
%--------------------------------------------------------------------------

% This is an optional block of code that will display a fusion of the raw
% images (green hue) with ilastik probabilities (purple huge)

% loop through number of images
for i = 1:size(raw_images,4)
    
    % show the ith image
    imshowpair(raw_images(:,:,:,i),ilastik_probabilities(:,:,i))
    pause(0.5)
    
end



%%
%--------------------------------------------------------------------------
% Segment ommatidia
%--------------------------------------------------------------------------

% Here, we will compute the centroids of the objects we detected via pixel
% classification in Ilastik. We achieve this by thresholding the pixel
% classification probabilities exported from Ilastik.
[omma_cent,omma_area] = initialSeg(ilastik_probabilities,0.3);

% remove out objects with less area than given threshold
[omma_cent] = sizeThreshOmma(omma_cent,omma_area,0.7);

% merge close ommatidia by dilating and then finding centroid again
[omma_cent] = mergeCloseOmma(omma_cent,4);



%%
%--------------------------------------------------------------------------
% hand correct using GUI
%--------------------------------------------------------------------------

% Here, we load the GUI for hand-correcting segmentation 

adultOmmatidiaSeg



%% 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% OPTIONAL: read in hand corrected centroids from file
% ONLY DO THIS IF YOU 
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

filepath = '/Users/kevin/Documents/MATLAB/forSha/omma_cent_disp_handCorrect_first75/*tif';
[omma_cent] = manualCorrectFromFile(filepath);


%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% OPTIONAL: draw ROI on each image and throw away centroids not in ROI
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

[omma_cent] = customROI(raw_images,omma_cent);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% triangulate points and find neighbors
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t = 1:length(ROI_centroids)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delaunay triangulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % delaunay triangulation
    T = delaunay(ROI_centroids{t}(:,1),ROI_centroids{t}(:,2));
    TR = triangulation(T,ROI_centroids{t});
    
    % find boundary
    boundary_cent = unique(boundary(ROI_centroids{t}));
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
    newROIcentroids = zeros(size(ROI_centroids{t},1) - length(unique(boundary_cent)),2);
    count = 0;
    count2 = 0;
    for j = 1:size(ROI_centroids{t},1)
        
        % if its not a boundary centroid, then record it in our
        % new matrix
        if not(ismember(j,boundary_cent))
            count = count + 1;
            newROIcentroids(count,:) = ROI_centroids{t}(j,:);
        
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
    edge_clean_delaunay{t} = newT;                  % definition of traingles - refers to indices of centroids in edge_clean_ROIcentroids
    original_triangulation{t} = TR;
    edge_clean_triangulation{t} = newTR;            % contains both centroid coordinates and definition of traingle connectivity
    edge_clean_ROIcentroids{t} = newROIcentroids;   % just a list of ommatidia centroids
    
    
    
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

%%

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

for t = 2
    
    boundary_cent = unique(boundary(edge_clean_ROIcentroids{t},0.8));
    
    x = edge_clean_ROIcentroids{t}(:,1);
    y = edge_clean_ROIcentroids{t}(:,2);
    
%     imshow(raw_images(:,:,:,t))
%     hold on
%     triplot(edge_clean_triangulation{t},'LineWidth',2,'Color','cyan')
%     for j = 1:length(boundary_cent)
%         plot(x(boundary_cent(j)),y(boundary_cent(j)),'ro','MarkerSize',12, 'LineWidth', 4)
%     end
%     filename = ['/Users/kevin/Documents/MATLAB/forSha/visualizing_neighbors/T=' num2str(0,'%03i') '.png'];
%     print(gcf,'-dpng',filename)
    
        for j = 1:size(x,1)
            if not(ismember(j,boundary_cent))
                imshow(raw_images(:,:,:,t))
                hold on
                triplot(edge_clean_triangulation{t},'LineWidth',2,'Color','cyan')
                plot(x(j),y(j),'r*','MarkerSize',12, 'LineWidth', 4)
                plot(x(delaunay_neighbors{t}{j}),y(delaunay_neighbors{t}{j}),'r.', ...
                    'MarkerSize', 12)
                hold off
    %             set(gcf,'PaperPosition',[0 0 (size(raw_images(:,:,1,1),2)/50)*2 (size(raw_images(:,:,1,1),1)/50)])
    %             pbaspect([1 1 1])
                filename = ['/Users/kevin/Documents/MATLAB/forSha/visualizing_neighbors/T=' num2str(j,'%03i') '.png'];
                print(gcf,'-dpng',filename)
            end
        end
    
end
    
%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% MAKE MEASUREMENTS
% first find fano facter per ommatidia, average per eye, then additionally
% average multiple eyes with the same genotype

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(edge_clean_ROIcentroids),1);

% initialize cell array for recording fano factor of omma
dist_fano = cell(length(edge_clean_ROIcentroids),1);
dist_fano_mean = zeros(length(edge_clean_ROIcentroids),1);

for t = 1:length(edge_clean_ROIcentroids)
    
    t
    
    % boundary centroids for current time
    boundary_cent{t} = unique(boundary(edge_clean_ROIcentroids{t},0.8));
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = edge_clean_ROIcentroids{t}(:,1);
    y = edge_clean_ROIcentroids{t}(:,2);
    
    % initialize list for storing fano for omma in current eye
    temp_dist_fano = nan(length(x),1);
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'edge_clean_ROIcentroids'
    for j = 1:size(x,1)
        
        % make sure its not a boundary point
        if not(ismember(j,boundary_cent{t}))
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %
            % variance in distance
            %
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            temp_distances = zeros(length(delaunay_neighbors{t}{j}),1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through current neighbors, calculate distance between
            % center point and each neighbor, and collect these
            % measurements in a list
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for jj = 1:length(delaunay_neighbors{t}{j})
                
                % store centroid components for current center position and
                % current neighbor
                curr_x = x(j);
                curr_y = y(j);
                neigh_x = x(delaunay_neighbors{t}{j}(jj));
                neigh_y = y(delaunay_neighbors{t}{j}(jj));
                
                % euclidean distance b/w two points
                temp_D = sqrt( ((neigh_x - curr_x)^2) + ((neigh_y - curr_y)^2) );
                
                % store in list of distances for current center point
                temp_distances(jj) = temp_D;
                
            end
            
            %--------------------------------------------
            % calculate index of dispersion (fano factor)
            %--------------------------------------------
            
            % calculate variance
            temp_dist_fano(j) = std(temp_distances) / mean(temp_distances);
            
%             %----------------------------------------------
%             % visualization of neighbors for each ommatidia
%             %----------------------------------------------
%             imshow(raw_images(:,:,:,t))
%             hold on
%             triplot(edge_clean_triangulation{t},'LineWidth',2,'Color','cyan')
%             plot(x(j),y(j),'r*','MarkerSize',12, 'LineWidth', 4)
%             plot(x(delaunay_neighbors{t}{j}),y(delaunay_neighbors{t}{j}),'r.', ...
%                 'MarkerSize', 12)
%             hold off
%             filename = ['/Users/kevin/Documents/MATLAB/forSha/visualizing_neighbors/T=' num2str(j,'%03i') '.png'];
%             print(gcf,'-dpng',filename)

        end
    end
    
    %------------------------------------
    % record fano of triangle edge length
    %------------------------------------
    dist_fano{t} = temp_dist_fano;
    dist_fano_mean(t) = nanmean(temp_dist_fano);
        
        
%         %---------------------------------------------------------
%         % visualization that we are properly finding all triangles
%         % associated with current R8
%         %---------------------------------------------------------
%         imshow(raw_images(:,:,:,t))
%         hold on
%         for q = 1:length(curr_triangles)
%             for qq = 1:3
%                 currR8 = edge_clean_delaunay{t}(curr_triangles(q),qq);
%                 currCent = edge_clean_ROIcentroids{t}(currR8,:);
%                 plot(currCent(1),currCent(2),'r*','MarkerSize',12, 'LineWidth', 4)
%             end
%         end
%         drawnow
%         hold off
        
    
    
%     %-----------------------------
%     % visualization triangle areas
%     %-----------------------------
%     
%     % colormap of triangle area
%     max_area = round(max(triangle_areas{t}));
%     min_area = round(min(triangle_areas{t}));
%     area_range = min_area:max_area;
%     colormap = parula(length(min_area:max_area));
%     colors = zeros(length(triangle_areas{t}),3);
%     for q = 1:length(triangle_areas{t})
%         [~,LOCB] = ismember(round(triangle_areas{t}(q)),area_range);
%         colors(q,:) = colormap(LOCB,:);
%     end
%     
%     imshow(raw_images(:,:,:,t))
%     hold on
%     triplot(edge_clean_triangulation{t},'LineWidth',2,'Color','cyan')
%     patch('vertices', edge_clean_ROIcentroids{t},'Faces', edge_clean_delaunay{t}, ...
%         'FaceColor','flat', 'FaceVertexCData', colors, 'CDataMapping', 'direct', ...
%         'FaceAlpha', 1);
%     hold off
%     drawnow
%     filename = ['/Users/kevin/Documents/MATLAB/forSha/media/triangle_areas/T=' num2str(t,'%03i') '.png'];
%     print(gcf,'-dpng',filename)
    
    
    
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% sort by mean fano factor
% record genotype of each image
% calculate average fano factor for each genotype
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

num_meas = length(dist_fano_mean);

%------------------------------
% sort according to fano factor
%------------------------------

%---------------------
% distance fano factor
%---------------------
sorted_distance_fano = zeros(num_meas,2);
% second column is the actual fano measurements
sorted_distance_fano(:,2) = dist_fano_mean;
% first column is just a list from 1:num_measurements so we can keep track
% of original index after we've sorted according to column 2
sorted_distance_fano(:,1) = linspace(1,num_meas,num_meas);
sorted_distance_fano = sortrows(sorted_distance_fano,2);


% define list of genotypes
genotypes = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];

% vector defining colors that correspond to each genotype
genotype_color = lines(length(genotypes));

%--------------------------------------------------------------------------
% create vectors and cell arrays that record filename, genotype, and
% plotting color for each image
%--------------------------------------------------------------------------
plt_color_vect = zeros(num_meas,3);
namestr = cell(num_meas,1);
genotype = cell(num_meas,1);
for t = 1:num_meas
    
    % parse out filename
    namestr{t} = Directory(sorted_distance_fano(t,1)).name;
    namestr{t} = namestr{t}(1:end-4);
    namestr{t} = strrep(namestr{t},'_',' ');
    
    % parse out the genotype specifically
    genotype{t} = strsplit(namestr{t});
    genotype{t} = genotype{t}{2};
    
    % record genotype color in our matrix for coloring
    [~,LOCB] = ismember(genotype{t},genotypes);
    plt_color_vect(t,:) = genotype_color(LOCB,:);
    
end

%--------------------------------------------------------------------------
% create average fano factor for each genotype - distance based fano
%--------------------------------------------------------------------------
sum_dist_fano_per_geno = cell(length(genotypes),1);
ave_dist_fano_per_geno = zeros(length(genotypes),1);
std_dist_fano_per_geno = zeros(length(genotypes),1);
% counts_per_geno_dist = zeros(length(genotypes),1);
all_dist_fano_per_geno = nan(100,length(genotypes));
sorted_dist_fano_per_geno = all_dist_fano_per_geno;
for t = 1:num_meas
    
    % figure out which genotype index
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    % add to sum for this genotype
    sum_dist_fano_per_geno{LOCB} = [sum_dist_fano_per_geno{LOCB}  sorted_distance_fano(t,2)];
    
%     % update the count for this genotype
%     counts_per_geno_dist(LOCB) = counts_per_geno_dist(LOCB) + 1;
    
end

for j = 1:length(genotypes)
    ave_dist_fano_per_geno(j) = mean(sum_dist_fano_per_geno{j});
    std_dist_fano_per_geno(j) = std(sum_dist_fano_per_geno{j});
    all_dist_fano_per_geno(1:length(sum_dist_fano_per_geno{j}),j) = sum_dist_fano_per_geno{j};
end

%--------------------------------------------------------------------------
% sort distributions of distance fano factor according to mean fano factor levels
%--------------------------------------------------------------------------

dist_sort_order = zeros(length(genotypes),3);
dist_sort_order(1:end,1) = 1:length(genotypes);         % entry 1: index from 1:num_genotypes
dist_sort_order(1:end,2) = ave_dist_fano_per_geno;      % entry 2: average fano factor for genotypes
dist_sort_order(1:end,3) = std_dist_fano_per_geno;      % entry 3: standard deviation for genotypes
dist_sort_order = sortrows(dist_sort_order,2);          % sort according to mean fano factor

for k = 1:length(genotypes)
    j = dist_sort_order(k,1);   % record genotype index
    sorted_dist_fano_per_geno(1:length(sum_dist_fano_per_geno{j}),k) = sum_dist_fano_per_geno{j};
    sorted_dist_genotype_name{k} = genotypes{j};
end


%------------------------------------------------------------------
%------------------------------------------------------------------
%------------------------------------------------------------------

% plot

%------------------------------------------------------------------
%------------------------------------------------------------------
%------------------------------------------------------------------



%------------------------
% plot mean and error bar
%------------------------

subplot(1,3,1)

er = errorbar(area_sort_order(:,2),area_sort_order(:,3)/2,'o','Linewidth',2);
xticks([1 2 3 4 5 6 7])
xticklabels(sorted_area_genotype_name)
xlim([0 8])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
ax = gca;
ax.FontSize = 16;

title(['mean Inter-R8-distance COV'],'FontSize',16)
xlabel(['Genotype'],'FontSize',16)
ylabel(['Coefficient of Variation'],'FontSize',16)

hold off


% boxplot

subplot(1,3,2)
boxplot(sorted_dist_fano_per_geno)
set(gca,'xticklabel',sorted_dist_genotype_name)
title('Inter-R8-distance COV boxplot')
xlabel(['Genotype'])
ylabel(['Coefficient of Variation'])
ax = gca;
ax.FontSize = 16;

% violin plot

subplot(1,3,3)
vs = violinplot(sorted_dist_fano_per_geno,sorted_dist_genotype_name,'ShowMean',true);
title('Inter-R8-distance COV violin plot')
xlabel(['Genotype'])
ylabel(['Coefficient of Variation'])
ax = gca;
ax.FontSize = 16;



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% MAKE MEASUREMENTS
% first find COV facter per ommatidia, aggregate by genotype (combine 
% ommatidia from all eyes of same genptype), then average

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(edge_clean_ROIcentroids),1);

% initialize cell array for recording COV factor of omma
dist_COV = cell(length(edge_clean_ROIcentroids),1);
dist_COV_mean = zeros(length(edge_clean_ROIcentroids),1);

for t = 1:length(edge_clean_ROIcentroids)
    
    t
    
    % boundary centroids for current time
    boundary_cent{t} = unique(boundary(edge_clean_ROIcentroids{t},0.8));
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = edge_clean_ROIcentroids{t}(:,1);
    y = edge_clean_ROIcentroids{t}(:,2);
    
    % initialize list for storing COV for omma in current eye
    temp_dist_COV = nan(length(x),1);
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'edge_clean_ROIcentroids'
    for j = 1:size(x,1)
        
        % make sure its not a boundary point
        if not(ismember(j,boundary_cent{t}))
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %
            % variance in distance
            %
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            temp_distances = zeros(length(delaunay_neighbors{t}{j}),1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through current neighbors, calculate distance between
            % center point and each neighbor, and collect these
            % measurements in a list
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for jj = 1:length(delaunay_neighbors{t}{j})
                
                % store centroid components for current center position and
                % current neighbor
                curr_x = x(j);
                curr_y = y(j);
                neigh_x = x(delaunay_neighbors{t}{j}(jj));
                neigh_y = y(delaunay_neighbors{t}{j}(jj));
                
                % euclidean distance b/w two points
                temp_D = sqrt( ((neigh_x - curr_x)^2) + ((neigh_y - curr_y)^2) );
                
                % store in list of distances for current center point
                temp_distances(jj) = temp_D;
                
            end
            
            %--------------------------------------------
            % calculate index of dispersion (COV)
            %--------------------------------------------
            
            % calculate variance
            temp_dist_COV(j) = std(temp_distances) / mean(temp_distances);
            

        end
    end
    
    %------------------------------------
    % record COV of triangle edge length
    %------------------------------------
    dist_COV{t} = temp_dist_COV;
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% determine genotype for each image

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% define list of genotypes
genotypes = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];

% vector defining colors that correspond to each genotype
genotype_color = lines(length(genotypes));

%--------------------------------------------------------------------------
% create vectors and cell arrays that record filename, genotype, and
% plotting color for each image
%--------------------------------------------------------------------------
plt_color_vect = zeros(num_meas,3);
namestr = cell(num_meas,1);
genotype = cell(num_meas,1);
for t = 1:num_meas
    
    % parse out filename
    namestr{t} = Directory(t).name;
    namestr{t} = namestr{t}(1:end-4);
    namestr{t} = strrep(namestr{t},'_',' ');
    
    % parse out the genotype specifically
    genotype{t} = strsplit(namestr{t});
    genotype{t} = genotype{t}{2};
    
    % record genotype color in our matrix for coloring
    [~,LOCB] = ismember(genotype{t},genotypes);
    plt_color_vect(t,:) = genotype_color(LOCB,:);
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% aggregate measurements based on genotype

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

aggregate_interR8_COV = cell(1,length(genotypes));

% loop through images
for t = 1:num_meas
    
    % find genotype of current image
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    aggregate_interR8_COV{LOCB} = [aggregate_interR8_COV{LOCB},dist_COV{t}'];
    
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% plot

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% find mean and std of each genotype
temp_mean_COV = [];
for t = 1:length(genotypes)
    temp_mean_COV(t) = nanmean(aggregate_interR8_COV{t});
end

% sort by mean fano factor
sort_order = zeros(length(genotypes),2);
sort_order(:,1) = [1:7];
sort_order(:,2) = temp_mean_COV;
sort_order = sortrows(sort_order,2);

sorted_aggregate_interR8_COV = cell(length(genotypes),1);
sorted_genotypes = cell(length(genotypes),1);
for t = 1:length(genotypes)
    sort_ind = sort_order(t,1);
    sorted_genotypes{t} = genotypes{sort_ind};
    sorted_aggregate_interR8_COV{t} = aggregate_interR8_COV{sort_ind};
end

% find mean and std of each genotype
mean_COV = [];
std_COV = [];
for t = 1:length(genotypes)
    mean_COV(t) = nanmean(sorted_aggregate_interR8_COV{t});
    std_COV(t) = nanstd(sorted_aggregate_interR8_COV{t});
end

% store in matrix
COV_matrix_form = nan(12000,7);
for t = 1:7
    COV_matrix_form(1:length(sorted_aggregate_interR8_COV{t}),t) = sorted_aggregate_interR8_COV{t};
end



%--------------------------------------------------------------------------
% plot inter-R8 distance
%--------------------------------------------------------------------------

subplot(1,3,1)

er = errorbar(mean_COV,std_COV/2,'o','Linewidth',2);
xticks([1 2 3 4 5 6 7])
xticklabels(sorted_genotypes)
xlim([0 8])
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ax = gca;
ax.FontSize = 16; 

title(['mean inter-R8-distance COV'],'FontSize',16)
xlabel(['genotype'],'FontSize',16)
ylabel(['coefficeint of variation'],'FontSize',16)



%--------------------------------------------------------------------------
% boxplot
%--------------------------------------------------------------------------

subplot(1,3,2)

boxplot(COV_matrix_form)
set(gca,'xticklabel',sorted_genotypes)
title('inter-R8-distance COV boxplot')
xlabel(['genotype'])
ylabel(['coefficient of variation'])
ax = gca;
ax.FontSize = 16;


%--------------------------------------------------------------------------
% violin plot
%--------------------------------------------------------------------------

subplot(1,3,3)

vs = violinplot(COV_matrix_form,sorted_genotypes,'ShowMean',true);
title('inter-R8-distance COV violin plot')
xlabel(['genotype'])
ylabel(['coefficient of variation'])
ax = gca;
ax.FontSize = 16;



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% MAKE MEASUREMENTS

% aggregate inter-R8_distances from multiple eyes

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measure inter-R8-distances for each individual image

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(edge_clean_ROIcentroids),1);

interR8_distances = {};
interR8_links = {};

% loop through images
for t = 1:length(edge_clean_ROIcentroids)
    
    t
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = edge_clean_ROIcentroids{t}(:,1);
    y = edge_clean_ROIcentroids{t}(:,2);
    
    % define list of unique centroid-centroid pairs for this current image
    unique_links = nan(1,4);
    link_count = 0;
    
    % define list of inter-R8 distances for current eye
    curr_distances = [];
    dist_count = 0;
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'edge_clean_ROIcentroids'
    for j = 1:size(x,1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % variance in distance
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp_distances = zeros(length(delaunay_neighbors{t}{j}),1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop through current neighbors, calculate distance between
        % center point and each neighbor, and collect these
        % measurements in a list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj = 1:length(delaunay_neighbors{t}{j})
            
            % store centroid components for current center position and
            % current neighbor
            curr_x = x(j);
            curr_y = y(j);
            neigh_x = x(delaunay_neighbors{t}{j}(jj));
            neigh_y = y(delaunay_neighbors{t}{j}(jj));
            
            % check to see if either combination of centroid pairs is
            % currently stored in 'unique_links'
            [LIA, ~] = ismember([curr_x curr_y neigh_x neigh_y], unique_links, 'rows');
            [LIA2,~] = ismember([neigh_x neigh_y curr_x curr_y], unique_links, 'rows');
            
            % if neither centroid pair is stored in 'unique_links', then we
            % can record the distance between these two points
            if not(LIA) && not(LIA2)
                
                % record current centroid-pair in 'unique_links'
                link_count = link_count + 1;
                unique_links(link_count,:) = [curr_x curr_y neigh_x neigh_y];
                
                % euclidean distance b/w two points
                temp_D = sqrt( ((neigh_x - curr_x)^2) + ((neigh_y - curr_y)^2) );
                
                % store in list of distances for current center point
                dist_count = dist_count + 1;
                curr_distances(dist_count) = temp_D;
                
            end
            
        end
        
    end
    
    interR8_distances{t} = curr_distances;
    interR8_links{t} = unique_links;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% determine genotype for each image

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% define list of genotypes
genotypes = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];

% vector defining colors that correspond to each genotype
genotype_color = lines(length(genotypes));

%--------------------------------------------------------------------------
% create vectors and cell arrays that record filename, genotype, and
% plotting color for each image
%--------------------------------------------------------------------------
plt_color_vect = zeros(num_meas,3);
namestr = cell(num_meas,1);
genotype = cell(num_meas,1);
for t = 1:num_meas
    
    % parse out filename
    namestr{t} = Directory(t).name;
    namestr{t} = namestr{t}(1:end-4);
    namestr{t} = strrep(namestr{t},'_',' ');
    
    % parse out the genotype specifically
    genotype{t} = strsplit(namestr{t});
    genotype{t} = genotype{t}{2};
    
    % record genotype color in our matrix for coloring
    [~,LOCB] = ismember(genotype{t},genotypes);
    plt_color_vect(t,:) = genotype_color(LOCB,:);
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% aggregate measurements based on genotype

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

aggregate_interR8_distances = cell(1,length(genotypes));

% loop through images
for t = 1:num_meas
    
    % find genotype of current image
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    aggregate_interR8_distances{LOCB} = [aggregate_interR8_distances{LOCB},interR8_distances{t}];
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% plot

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

sorted_aggregate_interR8_distances = cell(length(genotypes),1);
sorted_genotypes = cell(length(genotypes),1);
for t = 1:length(genotypes)
    sort_ind = sort_order(t,1);
    sorted_genotypes{t} = genotypes{sort_ind};
    sorted_aggregate_interR8_distances{t} = aggregate_interR8_distances{sort_ind};
end

% find mean and std of each genotype
mean_dist = [];
std_dist = [];
for t = 1:length(genotypes)
    mean_dist(t) = mean(sorted_aggregate_interR8_distances{t});
    std_dist(t) = std(sorted_aggregate_interR8_distances{t});
end

% store in matrix
dist_matrix_form = nan(12000,7);
for t = 1:7
    dist_matrix_form(1:length(sorted_aggregate_interR8_distances{t}),t) = sorted_aggregate_interR8_distances{t};
end

% find MEDIAN normalized distributions
normalized_interR8_distance_cellArray = {};
for t = 1:7
    
    temp_median = median(sorted_aggregate_interR8_distances{t});
    normalized_interR8_distance_cellArray{t} = sorted_aggregate_interR8_distances{t} - temp_median;
    
end

% store normalized distributions in matrix
normalized_interR8_dist_matrix = nan(12000,7);
for t = 1:7
    normalized_interR8_dist_matrix(1:length(normalized_interR8_distance_cellArray{t}),t) = normalized_interR8_distance_cellArray{t};
end

% find mean and std of each genotype or normalized data
norm_mean_dist = [];
norm_std_dist = [];
for t = 1:length(genotypes)
    norm_mean_dist(t) = mean(normalized_interR8_distance_cellArray{t});
    norm_std_dist(t) = std(normalized_interR8_distance_cellArray{t});
end

%--------------------------------------------------------------------------
% plot inter-R8 distance
%--------------------------------------------------------------------------

figure(1)

subplot(1,2,1)

er = errorbar(mean_dist,std_dist/2,'o','Linewidth',2);
xticks([1 2 3 4 5 6 7])
xticklabels(sorted_genotypes)
xlim([0 8])
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ax = gca;
ax.FontSize = 16; 

title(['raw inter-R8 distance'],'FontSize',16)
xlabel(['genotype'],'FontSize',16)
ylabel(['inter-R8-distance (pixels)'],'FontSize',16)

subplot(1,2,2)

er = errorbar(norm_mean_dist,norm_std_dist/2,'o','Linewidth',2);
xticks([1 2 3 4 5 6 7])
xticklabels(sorted_genotypes)
xlim([0 8])
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ax = gca;
ax.FontSize = 16; 

title(['MEDIAN NORMALIZED inter-R8 distance'],'FontSize',16)
xlabel(['genotype'],'FontSize',16)
ylabel(['normalized inter-R8-distance'],'FontSize',16)


%--------------------------------------------------------------------------
% boxplot
%--------------------------------------------------------------------------


figure(2)

subplot(1,2,1)

boxplot(dist_matrix_form)
set(gca,'xticklabel',sorted_genotypes)
title('raw inter-R8-distance')
xlabel(['genotype'])
ylabel(['inter-R8-distances (pixels)'])
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)

boxplot(normalized_interR8_dist_matrix)
set(gca,'xticklabel',sorted_genotypes)
title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['normalized inter-R8-distances'])
ax = gca;
ax.FontSize = 16;

%--------------------------------------------------------------------------
% violin plot
%--------------------------------------------------------------------------

figure(3)

subplot(1,2,1)

vs = violinplot(dist_matrix_form,sorted_genotypes,'ShowMean',true);
title('raw inter-R8-distance')
xlabel(['genotype'])
ylabel(['inter-R8-distances (pixels)'])
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)

vs = violinplot(normalized_interR8_dist_matrix,sorted_genotypes,'ShowMean',true);
title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['normalized inter-R8-distances'])
ax = gca;
ax.FontSize = 16;

%--------------------------------------------------------------------------
% ecdf
%--------------------------------------------------------------------------

figure(4)

subplot(1,2,1)

hold on

for j = 1:7
    
    [f,x] = ecdf(sorted_aggregate_interR8_distances{j});
    plot(x,f,'Linewidth',2)
    
end

title('raw inter-R8-distance')
legend(sorted_genotypes,'Location','Southeast')
ylabel('F(x)')
xlabel('x = inter-R8-distances (pixels)')
ax = gca;
ax.FontSize = 16;

hold off

subplot(1,2,2)

hold on

for j = 1:7
    
    [f,x] = ecdf(normalized_interR8_distance_cellArray{j});
    plot(x,f,'Linewidth',2)
    
end

title('MEDIAN NORMALIZED inter-R8-distance')
legend(sorted_genotypes,'Location','Southeast')
ylabel('F(x)')
xlabel('x = MEDIAN NORMALIZED inter-R8-distances')
ax = gca;
ax.FontSize = 16;

hold off

%% eCDF -20 to 20

% MANUALLY REMOVE MIR7
alt_genotypes = [sorted_genotypes(1:3); sorted_genotypes(5:7)];
alt_normalized_interR8_dist_matrix = [normalized_interR8_dist_matrix(:,1:3), normalized_interR8_dist_matrix(:,4:7)];
alt_normalized_interR8_distance_cellArray = {normalized_interR8_distance_cellArray{1}, ...
    normalized_interR8_distance_cellArray{2}, ...
    normalized_interR8_distance_cellArray{3}, ...
    normalized_interR8_distance_cellArray{5}, ...
    normalized_interR8_distance_cellArray{6}, ...
    normalized_interR8_distance_cellArray{7}};

subplot(1,2,1)

% vs = violinplot(normalized_interR8_dist_matrix,sorted_genotypes,'ShowMean',true);     % all genotypes
vs = violinplot(alt_normalized_interR8_dist_matrix,alt_genotypes,'ShowMean',true);      % MIR7 REMOVED
title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['normalized inter-R8-distances'])
ax = gca;
ax.FontSize = 16;
camroll(-90)

subplot(1,2,2)

coloz = lines(7);


coloz2 = zeros(size(coloz));
coloz2(2,:) = coloz(6,:);
coloz2(3,:) = coloz(4,:);
coloz2(4,:) = coloz(2,:);
coloz2(5,:) = coloz(7,:);
coloz2(6,:) = coloz(5,:);
coloz2(7,:) = coloz(3,:);

hold on

% for j = 1:length(normalized_interR8_distance_cellArray)     % all genotypes
for j = 1:length(alt_normalized_interR8_distance_cellArray)     % MIR7 REMOVED
    
    if j == 1
        
        [f,x] = ecdf(normalized_interR8_distance_cellArray{j});
        p  = patchline(x,f,'edgecolor',coloz(1,:),'linewidth',10,'edgealpha',1);
        
    else
        
%         [f,x] = ecdf(normalized_interR8_distance_cellArray{j});     % all genotypes
        [f,x] = ecdf(alt_normalized_interR8_distance_cellArray{j});     % MIR7 REMOVED
        plot(x,f,'Linewidth',3,'Color',coloz2(j,:))
        
    end
    
    
end

title('MEDIAN NORMALIZED inter-R8-distance')
% legend(sorted_genotypes,'Location','Southeast')    % MIR7 REMOVED
legend(alt_genotypes,'Location','Southeast')    % MIR7 REMOVED
ylabel('F(x)')
xlabel('x = MEDIAN NORMALIZED inter-R8-distances')
ax = gca;
ax.FontSize = 16;
xlim([-20 20])

hold off

%% eCDF -20 to 20 - SPECIFY WHICH GENOTYPES WE WANT

target_genotypes = [1 2 3 5 6 7];

alt_genotypes = {};
alt_normalized_interR8_dist_matrix = [];
alt_normalized_interR8_distance_cellArray = {};
for j = 1:length(target_genotypes)
    alt_genotypes{j} = sorted_genotypes{target_genotypes(j)};
    alt_normalized_interR8_dist_matrix(:,j) = normalized_interR8_dist_matrix(:,target_genotypes(j));
    alt_normalized_interR8_distance_cellArray{j} = normalized_interR8_distance_cellArray{target_genotypes(j)};
end

subplot(1,2,1)

vs = violinplot(alt_normalized_interR8_dist_matrix,alt_genotypes,'ShowMean',true);
% title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['MEDIAN NORMALIZED inter-R8-distances'])
ax = gca;
ax.FontSize = 16;
camroll(-90)

subplot(1,2,2)

coloz = lines(7);


coloz2 = zeros(size(coloz));
coloz2(2,:) = coloz(6,:);
coloz2(3,:) = coloz(4,:);
coloz2(4,:) = coloz(2,:);
coloz2(5,:) = coloz(7,:);
coloz2(6,:) = coloz(5,:);
coloz2(7,:) = coloz(3,:);

hold on

for j = 1:length(alt_normalized_interR8_distance_cellArray)
    
    if j == 1
        
        [f,x] = ecdf(normalized_interR8_distance_cellArray{j});
        p  = patchline(x,f,'edgecolor',coloz(1,:),'linewidth',10,'edgealpha',1);
        
    else
        
        [f,x] = ecdf(alt_normalized_interR8_distance_cellArray{j});
        plot(x,f,'Linewidth',3,'Color',coloz2(j,:))
        
    end
    
    
end

% title('MEDIAN NORMALIZED inter-R8-distance')
legend(alt_genotypes,'Location','Southeast')
ylabel('F(x)')
xlabel('x = MEDIAN NORMALIZED inter-R8-distances')
ax = gca;
ax.FontSize = 16;
xlim([-20 20])

hold off


% sgtitle('MEDIAN NORMALIZED inter-R8-distance')


%%

%--------------------------------------------------------------------------
%
% plot inter-R8-distance COV factor for each image separately
%
%--------------------------------------------------------------------------
for t = 1:75
    
    clf()
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot dispersion index vs. sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,1,'position',[0.1 0.1 0.4 0.82])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    hold on
    for j = 1:75
        plot(j,sorted_distance_COV(j,2),'.','Color',plt_color_vect(j,:),'MarkerSize',40,'Linewidth',1)
    end
    plot(t,sorted_distance_COV(t,2),'ok','MarkerSize',18,'Linewidth',6)
    ylabel('Mean index of dispersion (σ^2/μ)','FontSize',18)
    xlabel('Ranked sample order','FontSize',18)
    ax=gca;
    ax.FontSize = 14;
    
    %legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'.','Color',genotype_color(1,:),'MarkerSize',40);
    h(2) = plot(NaN,NaN,'.','Color',genotype_color(2,:),'MarkerSize',40);
    h(3) = plot(NaN,NaN,'.','Color',genotype_color(3,:),'MarkerSize',40);
    h(4) = plot(NaN,NaN,'.','Color',genotype_color(4,:),'MarkerSize',40);
    h(5) = plot(NaN,NaN,'.','Color',genotype_color(5,:),'MarkerSize',40);
    h(6) = plot(NaN,NaN,'.','Color',genotype_color(6,:),'MarkerSize',40);
    h(7) = plot(NaN,NaN,'.','Color',genotype_color(7,:),'MarkerSize',40);
    legend(h, genotypes(1),genotypes(2),genotypes(3),genotypes(4),...
        genotypes(5),genotypes(6),genotypes(7),'Location','northwest',...
        'FontSize',20);
    
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corresponding eye image with lattice drawn on top
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2,'position',[0.26 0.1 0.95 0.821])
    imshow(raw_images(:,:,:,sorted_distance_COV(t,1)))
    hold on
    triplot(edge_clean_triangulation{sorted_distance_COV(t,1)},'LineWidth',2,'Color','cyan')
    hold off
    annotation('textbox',[0.51 0.43 0.5 0.5],'string',namestr{t},...
   'linestyle','none','FontSize',30,'Color','r')
    
    sgtitle('Neighbor length - mean index of dispersion','FontSize',24)
    
    drawnow
    
    
    % print to file
    filename = ['/Users/kevin/Documents/MATLAB/forSha/media/individual_image_interR8/T=' num2str(t,'%03i') '.png'];
    print(gcf,'-dpng',filename)

end

%%

%---------------------------------------
%
%
% plot unsorted fano factor per genotype
%
%
%---------------------------------------

%-----------------------
% plot inter-R8 distance
%-----------------------

subplot(1,2,1)
bar(ave_dist_fano_per_geno)                
set(gca,'xticklabel',genotypes)
hold on

er = errorbar(ave_dist_fano_per_geno,std_dist_fano_per_geno/2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%---------------------------
% plot lattice triangle area
%---------------------------

subplot(1,2,2)
bar(ave_area_fano_per_geno)                
set(gca,'xticklabel',genotypes)
hold on

er = errorbar(ave_area_fano_per_geno,std_area_fano_per_geno/2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off


%%



