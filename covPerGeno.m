function [] = ...
    covPerGeno(filepath,genotypes,clean_omma_centroids,delaunay_neighbors)


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

num_meas = length(clean_omma_centroids);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(clean_omma_centroids),1);

% initialize cell array for recording COV factor of omma
dist_COV = cell(length(clean_omma_centroids),1);
dist_COV_mean = zeros(length(clean_omma_centroids),1);

for t = 1:length(clean_omma_centroids)
    
    t
    
    % boundary centroids for current time
    boundary_cent{t} = unique(boundary(clean_omma_centroids{t},0.8));
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = clean_omma_centroids{t}(:,1);
    y = clean_omma_centroids{t}(:,2);
    
    % initialize list for storing COV for omma in current eye
    temp_dist_COV = nan(length(x),1);
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'clean_omma_centroids'
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



% %--------------------------------------------------------------------------
% % plot inter-R8 distance
% %--------------------------------------------------------------------------
% 
% subplot(1,3,1)
% 
% er = errorbar(mean_COV,std_COV/2,'o','Linewidth',2);
% xticks([1 2 3 4 5 6 7])
% xticklabels(sorted_genotypes)
% xlim([0 8])
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% ax = gca;
% ax.FontSize = 16; 
% 
% title(['mean inter-R8-distance COV'],'FontSize',16)
% xlabel(['genotype'],'FontSize',16)
% ylabel(['coefficeint of variation'],'FontSize',16)
% 
% 
% 
% %--------------------------------------------------------------------------
% % boxplot
% %--------------------------------------------------------------------------
% 
% subplot(1,3,2)
% 
% boxplot(COV_matrix_form)
% set(gca,'xticklabel',sorted_genotypes)
% title('inter-R8-distance COV boxplot')
% xlabel(['genotype'])
% ylabel(['coefficient of variation'])
% ax = gca;
% ax.FontSize = 16;


%--------------------------------------------------------------------------
% violin plot
%--------------------------------------------------------------------------

% subplot(1,3,3)

vs = violinplot(COV_matrix_form,sorted_genotypes,'ShowMean',true);
title('inter-R8-distance COV violin plot')
xlabel(['genotype'])
ylabel(['coefficient of variation'])
ax = gca;
ax.FontSize = 16;







