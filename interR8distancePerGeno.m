function [] = ...
    interR8distancePerGeno(filepath,genotypes,edge_clean_ROIcentroids,delaunay_neighbors)

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

num_meas = length(edge_clean_ROIcentroids);

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

Directory = dir(strcat(filepath,'*.tif'));

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

% find mean and std of each genotype
temp_mean_distances = [];
for t = 1:length(genotypes)
    temp_mean_distances(t) = nanmean(aggregate_interR8_distances{t});
end

% sort by mean fano factor
sort_order = zeros(length(genotypes),2);
sort_order(:,1) = [1:7];
sort_order(:,2) = temp_mean_distances;
sort_order = sortrows(sort_order,2);

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