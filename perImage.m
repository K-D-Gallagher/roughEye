function [] = ...
    perImage(expInfo,genotypes,clean_omma_centroids,delaunay_neighbors, ...
    target_genotypes,plot_style,ascending_mean,genotype_labels,...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size,...
    save_csv_to_file,y_axis_limit,...
    measurement_type,distance_cutoff)


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

% find ommatidia within distance cutoff

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% pass to new variable that we will modify based on distance_cutoff
omma_for_analysis = {};

for t = 1:num_meas
    
    omma_count = 0;
    
    % find the center of mass for the segmented ommatidia of the current
    % image
    center_x = sum(clean_omma_centroids{t}(:,1))/length(clean_omma_centroids{t}(:,1));
    center_y = sum(clean_omma_centroids{t}(:,2))/length(clean_omma_centroids{t}(:,2));
    
    for j = 1:length(clean_omma_centroids{t})
        
        term1 = (clean_omma_centroids{t}(j,1) - center_x)^2;
        term2 = (clean_omma_centroids{t}(j,2) - center_y)^2;
        dist = sqrt(term1 + term2);
        
        if dist < distance_cutoff
            
            omma_count = omma_count + 1;
            omma_for_analysis{t}(omma_count) = j;
            
        end
        
    end
    
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(clean_omma_centroids),1);

% initialize cell array for recording COV of omma
DATA = cell(length(clean_omma_centroids),1);
DATA_mean = zeros(length(clean_omma_centroids),1);

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
    temp_meas = nan(length(x),1);
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'clean_omma_centroids'
    for j = 1:size(x,1)
        
        % make sure its not a boundary point and within distance cutoff
        if not(ismember(j,boundary_cent{t})) && ismember(j,omma_for_analysis{t})
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            %
            %
            % variance in distance
            %
            %
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            temp_measurements = nan(length(delaunay_neighbors{t}{j}),1);
            
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
                temp_measurements(jj) = temp_D;
                
            end
            
            %--------------------------------------------
            % calculate index of dispersion (COV) or min/max
            %--------------------------------------------
            
            % calculate variance
            if strcmp('COV',measurement_type)
                
                temp_meas(j) = std(temp_measurements) / mean(temp_measurements);
                
            elseif strcmp('max_min',measurement_type)
                
                if not(isempty(temp_measurements))
                
                    temp_meas(j) = min(temp_measurements) / max(temp_measurements);
                
                end
                
            else
                
                fprintf('Measurement type not properly specified')
                
            end

        end
    end
    
    %------------------------------------
    % record fano of triangle edge length
    %------------------------------------
    DATA{t} = temp_meas;
    DATA_mean(t) = nanmean(temp_meas);
    
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
% create average COV for each genotype - distance based fano
%--------------------------------------------------------------------------
DATA_per_geno = cell(length(genotypes),1);
ave_DATA_per_geno = zeros(length(genotypes),1);
std_DATA_per_geno = zeros(length(genotypes),1);
all_DATA_per_geno = nan(num_meas,length(genotypes));
sorted_DATA_per_geno = all_DATA_per_geno;

for t = 1:num_meas
    
    % figure out which genotype index
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    % add to sum for this genotype
    DATA_per_geno{LOCB} = [DATA_per_geno{LOCB}  DATA_mean(t)];
    
end

% compute average and std for each geno
for j = 1:length(genotypes)
    ave_DATA_per_geno(j) = mean(DATA_per_geno{j});
    std_DATA_per_geno(j) = std(DATA_per_geno{j});
    all_DATA_per_geno(1:length(DATA_per_geno{j}),j) = DATA_per_geno{j};
end


%--------------------------------------------------------------------------
% sort distributions of distance COV according to mean COV levels
%--------------------------------------------------------------------------

dist_sort_order = zeros(length(genotypes),3);
dist_sort_order(1:end,1) = 1:length(genotypes);         % entry 1: index from 1:num_genotypes
dist_sort_order(1:end,2) = ave_DATA_per_geno;      % entry 2: average COV for genotypes
dist_sort_order(1:end,3) = std_DATA_per_geno;      % entry 3: standard deviation for genotypes
dist_sort_order = sortrows(dist_sort_order,2);          % sort according to mean COV

for k = 1:length(genotypes)
    j = dist_sort_order(k,1);   % record genotype index
    sorted_DATA_per_geno(1:length(DATA_per_geno{j}),k) = DATA_per_geno{j};
    sorted_genotypes{k} = genotypes{j};
    sorted_genotype_labels(k) = genotype_labels(j);     % for plotting
    sorted_ave_DATA_per_geno(k) = ave_DATA_per_geno(j);
    sorted_std_DATA_per_geno(k) = std_DATA_per_geno(j);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% save to file if desired
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if save_csv_to_file
    
    T1 = array2table(all_DATA_per_geno,'VariableNames',...
        cellstr(genotypes));
    
    writetable(T1,fullfile(expInfo.filepath_output,strcat('/','COVperImage.xls')),...
        'WriteVariableNames',true);
    
end


%--------------------------------------------------------------------------
% pull out target genotypes for plotting
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

alt_sorted_genotype_labels = [];
alt_sorted_DATA_per_geno = [];
alt_sorted_ave_DATA_per_geno = zeros(length(target_genotypes),1);
alt_sorted_std_DATA_per_geno = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(sorted_genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(sorted_genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_sorted_genotype_labels{count} = sorted_genotype_labels(j); % for plotting
        alt_sorted_DATA_per_geno(:,count) = sorted_DATA_per_geno(:,j);
        alt_sorted_ave_DATA_per_geno(count) = sorted_ave_DATA_per_geno(j);
        alt_sorted_std_DATA_per_geno(count) = sorted_std_DATA_per_geno(j);
        
    end
    
end


%---------
% unsorted
%---------

alt_genotype_labels = {};
alt_DATA_per_geno = [];
alt_ave_DATA_per_geno = zeros(length(target_genotypes),1);
alt_std_DATA_per_geno = zeros(length(target_genotypes),1);

count = 0;

for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        alt_genotype_labels{count} = genotype_labels(j); % for plotting
        alt_DATA_per_geno(:,count) = all_DATA_per_geno(:,j);
        alt_ave_DATA_per_geno(count) = ave_DATA_per_geno(j);
        alt_std_DATA_per_geno(count) = std_DATA_per_geno(j);
        
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

        er = errorbar(alt_sorted_ave_DATA_per_geno,alt_sorted_std_DATA_per_geno/2,'o','Linewidth',2);
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
        ylim(y_axis_limit)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        er = errorbar(alt_ave_DATA_per_geno,alt_std_DATA_per_geno/2,'o','Linewidth',2);
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
        ylim(y_axis_limit)
        
        
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

        boxplot(alt_sorted_DATA_per_geno)
        set(gca,'xticklabel',alt_sorted_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        ylim(y_axis_limit)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        boxplot(alt_DATA_per_geno)
        set(gca,'xticklabel',alt_genotype_labels)
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        ylim(y_axis_limit)
        
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

        vs = violinplot(alt_sorted_DATA_per_geno,alt_sorted_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        ylim(y_axis_limit)
        
    else
        
        % update figure count and generate new figure
        figure_count = figure_count + 1;
        figure(figure_count)
        
        vs = violinplot(alt_DATA_per_geno,alt_genotype_labels,'ShowMean',true);
        title(plot_title,'FontSize',title_size)
        xlabel(x_label)
        ylabel(y_label)
        ax = gca;
        ax.FontSize = axes_label_size;
        
        xtickangle(x_axis_text_angle)
        ylim(y_axis_limit)
        
    end
    
end