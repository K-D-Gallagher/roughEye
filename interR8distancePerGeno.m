function [] = ...
    interR8distancePerGeno(expInfo,genotypes,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes,plot_style,ascending_mean,median_normalized,genotype_labels,...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size, x_axis_lim,...
    conversion_factor,save_csv_to_file,y_axis_limit)

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

disp('\n')
disp("Measuring inter-ommatidial-distance and aggregating across images")
disp("based on genotype. Sample:   ")

num_meas = length(clean_omma_centroids);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measure inter-R8-distances for each individual image

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(clean_omma_centroids),1);

interR8_distances = {};
interR8_links = {};

% loop through images
for t = 1:length(clean_omma_centroids)
    
    % display counter
    if t > 1
        for j=0:log10(t-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',t)
    end
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = clean_omma_centroids{t}(:,1);
    y = clean_omma_centroids{t}(:,2);
    
    % define list of unique centroid-centroid pairs for this current image
    unique_links = nan(1,4);
    link_count = 0;
    
    % define list of inter-R8 distances for current eye
    curr_distances = [];
    dist_count = 0;
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'clean_omma_centroids'
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
                
                % convert to microns
                temp_D = conversion_factor * temp_D;
                
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

Directory = dir(strcat(expInfo.filepath_input,'*.tif'));

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
    genotype{t} = genotype{t}{1};
    
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

% make measurements and format for plotting

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

% find mean and std of each genotype
temp_mean_distances = [];
for t = 1:length(genotypes)
    temp_mean_distances(t) = nanmean(aggregate_interR8_distances{t});
end

% sort by mean fano factor
sort_order = zeros(length(genotypes),2);
sort_order(:,1) = 1:length(genotypes);
sort_order(:,2) = temp_mean_distances;
sort_order = sortrows(sort_order,2);

sorted_aggregate_interR8_distances = cell(length(genotypes),1);
sorted_genotypes = cell(length(genotypes),1);
for t = 1:length(genotypes)
    sort_ind = sort_order(t,1);
    sorted_genotypes{t} = genotypes{sort_ind};
    sorted_aggregate_interR8_distances{t} = aggregate_interR8_distances{sort_ind};
end

% find mean and std of each genotype, plus max length of measurements
sorted_mean_dist = [];
sorted_std_dist = [];
unsorted_mean_dist = [];
unsorted_std_dist = [];
max_length = 0;
for t = 1:length(genotypes)
    sorted_mean_dist(t) = mean(sorted_aggregate_interR8_distances{t});
    sorted_std_dist(t) = std(sorted_aggregate_interR8_distances{t});
    unsorted_mean_dist(t) = mean(aggregate_interR8_distances{t});
    unsorted_std_dist(t) = std(aggregate_interR8_distances{t});
    % find max length
    curr_length = length(aggregate_interR8_distances{t});
    if curr_length > max_length
        max_length = curr_length;
    end
end

% store in matrix
sorted_interR8_dist_matrix_form = nan(max_length,length(genotypes));
unsorted_interR8_dist_matrix_form = nan(max_length,length(genotypes));
for t = 1:length(genotypes)
    sorted_interR8_dist_matrix_form(1:length(sorted_aggregate_interR8_distances{t}),t) = sorted_aggregate_interR8_distances{t};
    unsorted_interR8_dist_matrix_form(1:length(aggregate_interR8_distances{t}),t) = aggregate_interR8_distances{t};
end


%--------------------------------------------------------------------------
% find MEDIAN normalized distributions
%--------------------------------------------------------------------------

sorted_median_normalized_interR8_distance_cellArray = {};
unsorted_median_normalized_interR8_distance_cellArray = {};
for t = 1:length(genotypes)
    
    sorted_temp_median = median(sorted_aggregate_interR8_distances{t});
    unsorted_temp_median = median(aggregate_interR8_distances{t});
    sorted_median_normalized_interR8_distance_cellArray{t} = sorted_aggregate_interR8_distances{t} - sorted_temp_median;
    unsorted_median_normalized_interR8_distance_cellArray{t} = aggregate_interR8_distances{t} - unsorted_temp_median;
    
end

% store normalized distributions in matrix
sorted_median_normalized_interR8_dist_matrix = nan(max_length,length(genotypes));
unsorted_median_normalized_interR8_dist_matrix = nan(max_length,length(genotypes));
for t = 1:length(genotypes)
    sorted_median_normalized_interR8_dist_matrix(1:length(sorted_median_normalized_interR8_distance_cellArray{t}),t) = sorted_median_normalized_interR8_distance_cellArray{t};
    unsorted_median_normalized_interR8_dist_matrix(1:length(unsorted_median_normalized_interR8_distance_cellArray{t}),t) = unsorted_median_normalized_interR8_distance_cellArray{t};

end

% find mean and std of each genotype of normalized data
sorted_norm_mean_dist = [];
sorted_norm_std_dist = [];
unsorted_norm_mean_dist = [];
unsorted_norm_std_dist = [];
for t = 1:length(genotypes)
    sorted_norm_mean_dist(t) = mean(sorted_median_normalized_interR8_distance_cellArray{t});
    sorted_norm_std_dist(t) = std(sorted_median_normalized_interR8_distance_cellArray{t});
    unsorted_norm_mean_dist(t) = mean(unsorted_median_normalized_interR8_distance_cellArray{t});
    unsorted_norm_std_dist(t) = std(unsorted_median_normalized_interR8_distance_cellArray{t});
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% save to file if desired
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if save_csv_to_file
    
    T1 = array2table(unsorted_interR8_dist_matrix_form,'VariableNames',...
        cellstr(genotypes));
    T2 = array2table(unsorted_median_normalized_interR8_dist_matrix,'VariableNames',...
        cellstr(genotypes));
    
    writetable(T1,fullfile(expInfo.filepath_output,strcat('/','interR8distance_unnormalized.xls')),...
        'WriteVariableNames',true);
    writetable(T2,fullfile(expInfo.filepath_output,strcat('/','interR8distance_MEDIANnormalized.xls')),...
        'WriteVariableNames',true);
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% pull out target genotypes for plotting
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% create new variables to house target genotype measurements
%--------------------------------------------------------------------------

target_sorted_mean_dist = [];
target_sorted_std_dist = [];
target_unsorted_mean_dist = [];
target_unsorted_std_dist = [];

target_sorted_norm_mean_dist = [];
target_sorted_norm_std_dist = [];
target_unsorted_norm_mean_dist = [];
target_unsorted_norm_std_dist = [];

target_sorted_interR8_dist_matrix_form = nan(size(sorted_interR8_dist_matrix_form,1),length(target_genotypes));
target_unsorted_interR8_dist_matrix_form = nan(size(sorted_interR8_dist_matrix_form,1),length(target_genotypes));
target_sorted_median_normalized_interR8_dist_matrix = nan(size(sorted_interR8_dist_matrix_form,1),length(target_genotypes));
target_unsorted_median_normalized_interR8_dist_matrix = nan(size(sorted_interR8_dist_matrix_form,1),length(target_genotypes));

target_sorted_aggregate_interR8_distances = {};
target_aggregate_interR8_distances = {};
target_sorted_median_normalized_interR8_distance_cellArray = {};
target_unsorted_median_normalized_interR8_distance_cellArray = {};

target_sorted_genotypes = {};

%--------------------------------------------------------------------------
% pull out measurements corresponding to target genotypes
%--------------------------------------------------------------------------

count1 = 0;
count2 = 0;
for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA1,~] = ismember(genotypes(j),target_genotypes);
    [LIA2,~] = ismember(sorted_genotypes(j),target_genotypes);
    
    if LIA1
        
        count1 = count1 + 1;

        target_unsorted_mean_dist(count1) = unsorted_mean_dist(j);
        target_unsorted_std_dist(count1) = unsorted_std_dist(j);

        target_unsorted_norm_mean_dist(count1) = unsorted_norm_mean_dist(j);
        target_unsorted_norm_std_dist(count1) = unsorted_norm_std_dist(j);
   
        target_unsorted_interR8_dist_matrix_form(:,count1) = unsorted_interR8_dist_matrix_form(:,j);
        target_unsorted_median_normalized_interR8_dist_matrix(:,count1) = unsorted_median_normalized_interR8_dist_matrix(:,j);

        target_aggregate_interR8_distances{count1} = aggregate_interR8_distances{j};
        target_unsorted_median_normalized_interR8_distance_cellArray{count1} = unsorted_median_normalized_interR8_distance_cellArray{j};
        
    end
        
    if LIA2
        
        count2 = count2 + 1;
        
        target_sorted_mean_dist(count2) = sorted_mean_dist(j);
        target_sorted_std_dist(count2) = sorted_std_dist(j);
        
        target_sorted_norm_mean_dist(count2) = sorted_norm_mean_dist(j);
        target_sorted_norm_std_dist(count2) = sorted_norm_std_dist(j);
        
        target_sorted_interR8_dist_matrix_form(:,count2) = sorted_interR8_dist_matrix_form(:,j);
        target_sorted_median_normalized_interR8_dist_matrix(:,count2) = sorted_median_normalized_interR8_dist_matrix(:,j);
        
        target_sorted_aggregate_interR8_distances{count2} = sorted_aggregate_interR8_distances{j};
        target_sorted_median_normalized_interR8_distance_cellArray{count2} = sorted_median_normalized_interR8_distance_cellArray{j};
        
        target_sorted_genotypes{count2} = sorted_genotypes{j};
        
    end
    
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% PLOT
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

close all

%--------------------------------------------------------------------------
% mean + std
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'mean & std'))
    
    % NON-MEDIAN NORMALIZED
    if not(median_normalized)
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            er = errorbar(target_sorted_mean_dist,target_sorted_std_dist/2,'o','Linewidth',2);
            xticks([1:length(target_genotypes)])
            xticklabels(target_sorted_genotypes)
            xlim([0 length(target_genotypes)+1])
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            ax = gca;
            ax.FontSize = axes_label_size; 

            title(plot_title,'FontSize',title_size)
            xlabel(x_label,'FontSize',axes_label_size)
            ylabel(y_label,'FontSize',axes_label_size)
            
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        % NON-MEDIAN NORMALIZED, PLOTTING IN UNALTERED ORDER
        else
            
            er = errorbar(target_unsorted_mean_dist,target_unsorted_std_dist/2,'o','Linewidth',2);
            xticks([1:length(target_genotypes)])
            xticklabels(target_genotypes)
            xlim([0 length(target_genotypes)+1])
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            ax = gca;
            ax.FontSize = axes_label_size; 

            title(plot_title,'FontSize',title_size)
            xlabel(x_label,'FontSize',axes_label_size)
            ylabel(y_label,'FontSize',axes_label_size)
            
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        end

    % MEDIAN NORMALIZED
    else
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            er = errorbar(target_sorted_norm_mean_dist,target_sorted_norm_std_dist/2,'o','Linewidth',2);
            xticks([1:length(target_genotypes)])
            xticklabels(target_sorted_genotypes)
            xlim([0 length(target_genotypes)+1])
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            ax = gca;
            ax.FontSize = axes_label_size; 

            title([strcat(plot_title, " (median normalized)")],'FontSize',title_size)
            xlabel(x_label,'FontSize',axes_label_size)
            ylabel(strcat(y_label," (median normalized)"),'FontSize',axes_label_size)

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        % MEDIAN NORMALIZED, PLOTTING IN UNALTERED ORDER
        else
            
            er = errorbar(target_unsorted_norm_mean_dist,target_unsorted_norm_std_dist/2,'o','Linewidth',2);
            xticks([1:length(target_genotypes)])
            xticklabels(target_genotypes)
            xlim([0 length(target_genotypes)+1])
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            ax = gca;
            ax.FontSize = axes_label_size; 

            title([strcat(plot_title, " (median normalized)")],'FontSize',title_size)
            xlabel(x_label,'FontSize',axes_label_size)
            ylabel(strcat(y_label," (median normalized)"),'FontSize',axes_label_size)

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        end
        
    end
end

%--------------------------------------------------------------------------
% boxplot
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'box plot'))
    
    % NON-MEDIAN NORMALIZED
    if not(median_normalized)
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            boxplot(target_sorted_interR8_dist_matrix_form)
            set(gca,'xticklabel',target_sorted_genotypes)
            title(plot_title,'FontSize',title_size)
            xlabel(x_label)
            ylabel(y_label)
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        % NON-MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            boxplot(target_unsorted_interR8_dist_matrix_form)
            set(gca,'xticklabel',target_genotypes)
            title(plot_title,'FontSize',title_size)
            xlabel(x_label)
            ylabel(y_label)
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        end

    % MEDIAN NORMALIZED
    else
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            boxplot(target_sorted_median_normalized_interR8_dist_matrix)
            set(gca,'xticklabel',target_sorted_genotypes)
            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            xlabel(x_label)
            ylabel(strcat(y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        % MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            boxplot(target_unsorted_median_normalized_interR8_dist_matrix)
            set(gca,'xticklabel',target_genotypes)
            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            xlabel(x_label)
            ylabel(strcat(y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        end
        
    end
end

%--------------------------------------------------------------------------
% violin plot
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'violin plot'))
    
    % NON-MEDIAN NORMALIZED
    if not(median_normalized)
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            vs = violinplot(target_sorted_interR8_dist_matrix_form,target_sorted_genotypes,'ShowMean',true);
            title(plot_title,'FontSize',title_size)
            xlabel(x_label)
            ylabel(y_label)
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        % NON-MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            vs = violinplot(target_unsorted_interR8_dist_matrix_form,target_genotypes,'ShowMean',true);
            title(plot_title,'FontSize',title_size)
            xlabel(x_label)
            ylabel(y_label)
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        end

    % MEDIAN NORMALIZED
    else
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            vs = violinplot(target_sorted_median_normalized_interR8_dist_matrix,target_sorted_genotypes,'ShowMean',true);
            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            xlabel(x_label)
            ylabel(strcat(y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
            
        % MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            vs = violinplot(target_unsorted_median_normalized_interR8_dist_matrix,target_genotypes,'ShowMean',true);
            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            xlabel(x_label)
            ylabel(strcat(y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)
        
        end
        
    end
end

%--------------------------------------------------------------------------
% ecdf
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'eCDF'))
    
    % NON-MEDIAN NORMALIZED
    if not(median_normalized)
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            hold on

            for j = 1:length(target_genotypes)

                [f,x] = ecdf(target_sorted_aggregate_interR8_distances{j});
                plot(x,f,'Linewidth',2)

            end

            title(plot_title,'FontSize',title_size)
            legend(target_sorted_genotypes,'Location','Southeast')
            ylabel('F(x)')
            xlabel(strcat("x = ",x_label))
            ax = gca;
            ax.FontSize = axes_label_size;
            
            xlim(x_axis_lim)
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)

            hold off
            
        % NON-MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            hold on

            for j = 1:length(target_genotypes)

                [f,x] = ecdf(target_aggregate_interR8_distances{j});
                plot(x,f,'Linewidth',2)

            end

            title(plot_title,'FontSize',title_size)
            legend(target_genotypes,'Location','Southeast')
            ylabel('F(x)')
            xlabel(strcat("x = ",x_label))
            ax = gca;
            ax.FontSize = axes_label_size;

            xlim(x_axis_lim)
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)

            hold off
            
        end

    % MEDIAN NORMALIZED
    else
        
        % PLOTTED IN ORDER OF ASCENDING MEAN
        if ascending_mean
            
            hold on

            for j = 1:length(target_genotypes)

                [f,x] = ecdf(target_sorted_median_normalized_interR8_distance_cellArray{j});
                plot(x,f,'Linewidth',2)

            end

            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            legend(target_sorted_genotypes,'Location','Southeast')
            ylabel('F(x)')
            xlabel(strcat("x = ",y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xlim(x_axis_lim)
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)

            hold off
            
        % MEDIAN NORMALIZED, PLOTTED IN UNALTERED ORDER
        else
            
            hold on

            for j = 1:length(target_genotypes)

                [f,x] = ecdf(target_unsorted_median_normalized_interR8_distance_cellArray{j});
                plot(x,f,'Linewidth',2)

            end

            title(strcat(plot_title, " (median normalized)"),'FontSize',title_size)
            legend(target_genotypes,'Location','Southeast')
            ylabel('F(x)')
            xlabel(strcat("x = ",y_label," (median normalized)"))
            ax = gca;
            ax.FontSize = axes_label_size;

            xlim(x_axis_lim)
            xtickangle(x_axis_text_angle)
            ylim(y_axis_limit)

            hold off
        
        end
        
    end
end