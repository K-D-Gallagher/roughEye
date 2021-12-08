function [] = ...
    covPerImage(filepath,genotypes,clean_omma_centroids,delaunay_neighbors, ...
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
% first find COV per ommatidia, average per eye, then additionally
% average multiple eyes with the same genotype

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

num_meas = length(clean_omma_centroids);

disp('\n')
disp("Measuring COV of inter-ommatidial-distance and averaging per image.")
disp("Then aggregating and averaging per genotype. Sample:   ")


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(clean_omma_centroids),1);

% initialize cell array for recording COV of omma
COV = cell(length(clean_omma_centroids),1);
COV_mean = zeros(length(clean_omma_centroids),1);

for t = 1:length(clean_omma_centroids)
    
    % display counter
    if t > 1
        for j=0:log10(t-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',t)
    end
    
    % boundary centroids for current time
    boundary_cent{t} = unique(boundary(clean_omma_centroids{t},0.8));
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = clean_omma_centroids{t}(:,1);
    y = clean_omma_centroids{t}(:,2);
    
    % initialize list for storing fano for omma in current eye
    temp_COV = nan(length(x),1);
    
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
            temp_COV(j) = std(temp_distances) / mean(temp_distances);

        end
    end
    
    %------------------------------------
    % record fano of triangle edge length
    %------------------------------------
    COV{t} = temp_COV;
    COV_mean(t) = nanmean(temp_COV);
        
        
%         %---------------------------------------------------------
%         % visualization that we are properly finding all triangles
%         % associated with current R8
%         %---------------------------------------------------------
%         imshow(raw_images(:,:,:,t))
%         hold on
%         for q = 1:length(curr_triangles)
%             for qq = 1:3
%                 currR8 = edge_clean_delaunay{t}(curr_triangles(q),qq);
%                 currCent = clean_omma_centroids{t}(currR8,:);
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
%     patch('vertices', clean_omma_centroids{t},'Faces', edge_clean_delaunay{t}, ...
%         'FaceColor','flat', 'FaceVertexCData', colors, 'CDataMapping', 'direct', ...
%         'FaceAlpha', 1);
%     hold off
%     drawnow
%     filename = ['/Users/kevin/Documents/MATLAB/forSha/media/triangle_areas/T=' num2str(t,'%03i') '.png'];
%     print(gcf,'-dpng',filename)
    
    
    
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% determine genotype for each image
%
%
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
% create average COV for each genotype - distance based fano
%--------------------------------------------------------------------------
COV_per_geno = cell(length(genotypes),1);
ave_COV_per_geno = zeros(length(genotypes),1);
std_COV_per_geno = zeros(length(genotypes),1);
all_COV_per_geno = nan(1000,length(genotypes));
sorted_COV_per_geno = all_COV_per_geno;

for t = 1:num_meas
    
    % figure out which genotype index
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    % add to sum for this genotype
    COV_per_geno{LOCB} = [COV_per_geno{LOCB}  COV_mean(t)];
    
end

% compute average and std for each geno
for j = 1:length(genotypes)
    ave_COV_per_geno(j) = mean(COV_per_geno{j});
    std_COV_per_geno(j) = std(COV_per_geno{j});
    all_COV_per_geno(1:length(COV_per_geno{j}),j) = COV_per_geno{j};
end



%--------------------------------------------------------------------------
% sort distributions of distance COV according to mean COV levels
%--------------------------------------------------------------------------

dist_sort_order = zeros(length(genotypes),3);
dist_sort_order(1:end,1) = 1:length(genotypes);         % entry 1: index from 1:num_genotypes
dist_sort_order(1:end,2) = ave_COV_per_geno;      % entry 2: average COV for genotypes
dist_sort_order(1:end,3) = std_COV_per_geno;      % entry 3: standard deviation for genotypes
dist_sort_order = sortrows(dist_sort_order,2);          % sort according to mean COV

for k = 1:length(genotypes)
    j = dist_sort_order(k,1);   % record genotype index
    sorted_COV_per_geno(1:length(COV_per_geno{j}),k) = COV_per_geno{j};
    sorted_genotypes{k} = genotypes{j};
    sorted_genotype_labels(k) = genotype_labels(j);     % for plotting
    sorted_ave_COV_per_geno(k) = ave_COV_per_geno(j);
    sorted_std_COV_per_geno(k) = std_COV_per_geno(j);
end




%--------------------------------------------------------------------------
% pull out target genotypes for plotting
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

alt_sorted_genotype_labels = [];
alt_sorted_COV_per_geno = [];
alt_sorted_ave_COV_per_geno = zeros(length(target_genotypes),1);
alt_sorted_std_COV_per_geno = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(sorted_genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(sorted_genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_sorted_genotype_labels{count} = sorted_genotype_labels(j); % for plotting
        alt_sorted_COV_per_geno(:,count) = sorted_COV_per_geno(:,j);
        alt_sorted_ave_COV_per_geno(count) = sorted_ave_COV_per_geno(j);
        alt_sorted_std_COV_per_geno(count) = sorted_std_COV_per_geno(j);
        
    end
    
end


%---------
% unsorted
%---------

alt_genotype_labels = {};
alt_COV_per_geno = [];
alt_ave_COV_per_geno = zeros(length(target_genotypes),1);
alt_std_COV_per_geno = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_genotype_labels{count} = genotype_labels(j); % for plotting
        alt_COV_per_geno(:,count) = all_COV_per_geno(:,j);
        alt_ave_COV_per_geno(count) = ave_COV_per_geno(j);
        alt_std_COV_per_geno(count) = std_COV_per_geno(j);
        
    end
    
end




%------------------------------------------------------------------
%------------------------------------------------------------------
%
%
% plot
%
%
%------------------------------------------------------------------
%------------------------------------------------------------------

figure_count = 0;

%--------------------------------------------------------------------------
% mean + std
%--------------------------------------------------------------------------

if any(strcmp(plot_style,'mean & std'))
    
    if ascending_mean
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)

        er = errorbar(alt_sorted_ave_COV_per_geno,alt_sorted_std_COV_per_geno/2,'o','Linewidth',2);
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
        
        er = errorbar(alt_ave_COV_per_geno,alt_std_COV_per_geno/2,'o','Linewidth',2);
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

        boxplot(alt_sorted_COV_per_geno)
        set(gca,'xticklabel',alt_sorted_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        boxplot(alt_COV_per_geno)
        set(gca,'xticklabel',alt_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
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

        vs = violinplot(alt_sorted_COV_per_geno,alt_sorted_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        vs = violinplot(alt_COV_per_geno,alt_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        
    end
    
end