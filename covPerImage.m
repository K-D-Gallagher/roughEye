function [] = ...
    covPerImage(filepath,genotypes,edge_clean_ROIcentroids,delaunay_neighbors)


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


% vector defining colors that correspond to each genotype
genotype_color = lines(length(genotypes));

Directory = dir(strcat(filepath,'*.tif'));

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



% %------------------------
% % plot mean and error bar
% %------------------------
% 
% subplot(1,3,1)
% 
% er = errorbar(area_sort_order(:,2),area_sort_order(:,3)/2,'o','Linewidth',2);
% xticks([1 2 3 4 5 6 7])
% xticklabels(sorted_area_genotype_name)
% xlim([0 8])
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% ax = gca;
% ax.FontSize = 16;
% 
% title(['mean Inter-R8-distance COV'],'FontSize',16)
% xlabel(['Genotype'],'FontSize',16)
% ylabel(['Coefficient of Variation'],'FontSize',16)
% 
% hold off
% 
% 
% % boxplot
% 
% subplot(1,3,2)
% boxplot(sorted_dist_fano_per_geno)
% set(gca,'xticklabel',sorted_dist_genotype_name)
% title('Inter-R8-distance COV boxplot')
% xlabel(['Genotype'])
% ylabel(['Coefficient of Variation'])
% ax = gca;
% ax.FontSize = 16;

% violin plot

% subplot(1,3,3)
vs = violinplot(sorted_dist_fano_per_geno,sorted_dist_genotype_name,'ShowMean',true);
title('Inter-R8-distance COV violin plot')
xlabel(['Genotype'])
ylabel(['Coefficient of Variation'])
ax = gca;
ax.FontSize = 16;