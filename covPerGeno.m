function [] = ...
    covPerGeno(filepath,genotypes,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes,plot_style,ascending_mean,genotype_labels,...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size)


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
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            %
            %
            % variance in distance
            %
            %
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            temp_distances = zeros(length(delaunay_neighbors{t}{j}),1);
            
            %--------------------------------------------------------------
            % loop through current neighbors, calculate distance between
            % center point and each neighbor, and collect these
            % measurements in a list
            %--------------------------------------------------------------
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

Directory = dir(strcat(filepath,'*.tif'));

%--------------------------------------------------------------------------
% create vectors and cell arrays that record filename, genotype, and
% plotting color for each image
%--------------------------------------------------------------------------
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

% stats and prep for plotting (sort in order of ascending mean

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% find mean and std of each genotype in order to sort measurements in order
% to ascending mean COV
temp_mean_COV = [];
for t = 1:length(genotypes)
    temp_mean_COV(t) = nanmean(aggregate_interR8_COV{t});
end

%--------------------------------------------------------------------------
% sort by mean fano factor
%--------------------------------------------------------------------------
sort_order = zeros(length(genotypes),2);
sort_order(:,1) = [1:length(genotypes)];
sort_order(:,2) = temp_mean_COV;
sort_order = sortrows(sort_order,2);

sorted_aggregate_interR8_COV = cell(length(genotypes),1);
sorted_genotypes = cell(length(genotypes),1);
for t = 1:length(genotypes)
    sort_ind = sort_order(t,1);
    sorted_genotypes{t} = genotypes{sort_ind};
    sorted_genotype_labels(t) = genotype_labels(sort_ind);
    sorted_aggregate_interR8_COV{t} = aggregate_interR8_COV{sort_ind};
end
% convert list of sorted genotypes to string array
sorted_genotypes = string(sorted_genotypes);

%--------------------------------------------------------------------------
% find mean and std of each genotype
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

sorted_mean_COV = [];
sorted_std_COV = [];
for t = 1:length(genotypes)
    sorted_mean_COV(t) = nanmean(sorted_aggregate_interR8_COV{t});
    sorted_std_COV(t) = nanstd(sorted_aggregate_interR8_COV{t});
end

%---------
% unsorted
%---------

mean_COV = [];
std_COV = [];
for t = 1:length(genotypes)
    mean_COV(t) = nanmean(aggregate_interR8_COV{t});
    std_COV(t) = nanstd(aggregate_interR8_COV{t});
end

%--------------------------------------------------------------------------
% transfer cells into matrix
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

sorted_COV_matrix_form = nan(12000,length(genotypes));
for t = 1:length(genotypes)
    sorted_COV_matrix_form(1:length(sorted_aggregate_interR8_COV{t}),t) = sorted_aggregate_interR8_COV{t};
end

%---------
% unsorted
%---------

COV_matrix_form = nan(12000,length(genotypes));
for t = 1:length(genotypes)
    COV_matrix_form(1:length(aggregate_interR8_COV{t}),t) = aggregate_interR8_COV{t};
end


%--------------------------------------------------------------------------
% pull out target genotypes for plotting
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

alt_sorted_genotype_labels = [];
alt_sorted_COV_matrix_form = [];
alt_sorted_mean_COV = zeros(length(target_genotypes),1);
alt_sorted_std_COV = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(sorted_genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(sorted_genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_sorted_genotype_labels{count} = sorted_genotype_labels(j); % for plotting
        alt_sorted_COV_matrix_form(:,count) = sorted_COV_matrix_form(:,j);
        alt_sorted_mean_COV(count) = sorted_mean_COV(j);
        alt_sorted_std_COV(count) = sorted_std_COV(j);
        
    end
end


%---------
% unsorted
%---------

alt_genotype_labels = [];
alt_COV_matrix_form = [];
alt_mean_COV = zeros(length(target_genotypes),1);
alt_std_COV = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_genotypes{count} = genotypes(j);
        alt_genotype_labels{count} = genotype_labels(j); % for plotting
        alt_COV_matrix_form(:,count) = COV_matrix_form(:,j);
        alt_mean_COV(count) = mean_COV(j);
        alt_std_COV(count) = std_COV(j);
        
    end
end




figure_count = 0;



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% PLOT
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% mean + std
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'mean & std'))
    
    if ascending_mean
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)

        % plot
        er = errorbar(alt_sorted_mean_COV,alt_sorted_std_COV/2,'o','Linewidth',2);
        xticks(linspace(1,length(target_genotypes),length(target_genotypes)))
        xticklabels(alt_sorted_genotype_labels)
        xlim([0 length(target_genotypes)+1])
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        ax = gca;
        ax.FontSize = axes_label_size; 

        % title and axes labels
        title([plot_title],'FontSize',title_size)
        xlabel([x_label],'FontSize',axes_label_size)
        ylabel([y_label],'FontSize',axes_label_size)

        xtickangle(x_axis_text_angle)
        
    else
    
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)

        % plot
        er = errorbar(alt_mean_COV,alt_std_COV/2,'o','Linewidth',2);
        xticks(linspace(1,length(target_genotypes),length(target_genotypes)))
        xticklabels(alt_genotype_labels)
        xlim([0 length(target_genotypes)+1])
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        ax = gca;
        ax.FontSize = axes_label_size; 

        % title and axes labels
        title([plot_title],'FontSize',title_size)
        xlabel([x_label],'FontSize',axes_label_size)
        ylabel([y_label],'FontSize',axes_label_size)

        xtickangle(x_axis_text_angle)

    end
    
end


%--------------------------------------------------------------------------
% boxplot
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'box plot'))
    
    if ascending_mean
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)

        boxplot(alt_sorted_COV_matrix_form)
        set(gca,'xticklabel',alt_sorted_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel([x_label])
        ylabel([y_label])
        ax = gca;
        ax.FontSize = axes_label_size;

        xtickangle(x_axis_text_angle)
        
        
    else
    
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)

        boxplot(alt_COV_matrix_form)
        set(gca,'xticklabel',alt_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel([x_label])
        ylabel([y_label])
        ax = gca;
        ax.FontSize = axes_label_size;

        xtickangle(x_axis_text_angle)

    end
    
end


%--------------------------------------------------------------------------
% violin plot
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'violin plot'))
    
    if ascending_mean
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        vs = violinplot(alt_sorted_COV_matrix_form,alt_sorted_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel([x_label])
        ylabel([y_label])
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        vs = violinplot(alt_COV_matrix_form,alt_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel([x_label])
        ylabel([y_label])
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        
    end
    
end







