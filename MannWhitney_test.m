function [] = ...
    MannWhitney_test(expInfo,genotypes,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes, distance_cutoff, measurement_type, genotype_labels)


%--------------------------------------------------------------------------
% check to make sure we only have two target genotypes
%--------------------------------------------------------------------------

if length(target_genotypes) ~= 2
    disp('\n')
    disp("Please make sure you have selected two target genotypes for comparison")
    return
end


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

disp('\n')
disp("Measuring COV of inter-ommatidial-distance and aggregating across images")
disp("based on genotype. Sample:   ")

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

% initialize cell array for recording COV factor of omma
all_DATA = cell(length(clean_omma_centroids),1);
dist_COV_mean = zeros(length(clean_omma_centroids),1);

% loop thru eyes
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
    
    %----------------------------------------------------------------------
    % find average inter-ommatidial distance
    %----------------------------------------------------------------------
    all_distances = [];
    dist_count = 0;
    for j = 1:size(x,1)
        
        % make sure its not a boundary point and within distance cutoff
        if not(ismember(j,boundary_cent{t})) && ismember(j,omma_for_analysis{t})
            
            for jj = 1:length(delaunay_neighbors{t}{j})
                
                % store centroid components for current center position and
                % current neighbor
                curr_x = x(j);
                curr_y = y(j);
                neigh_x = x(delaunay_neighbors{t}{j}(jj));
                neigh_y = y(delaunay_neighbors{t}{j}(jj));
                
                % euclidean distance b/w two points
                temp_D = sqrt( ((neigh_x - curr_x)^2) + ((neigh_y - curr_y)^2) );
                
                dist_count = dist_count + 1;
                all_distances(dist_count) = temp_D;
            end
        end
    end
    
    % take mean distance
    mean_distance = mean(all_distances);
    
    
    %----------------------------------------------------------------------
    % apply metric for measuring disorder
    %----------------------------------------------------------------------
    
    % initialize list for storing COV for omma in current eye
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
            
            temp_distances = nan(length(delaunay_neighbors{t}{j}),1);
            
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
            % calculate index of dispersion (COV) or min/max
            %--------------------------------------------
            
            % calculate variance
            if strcmp('COV',measurement_type)
                
                temp_meas(j) = std(temp_distances) / mean(temp_distances);
                
            elseif strcmp('max_min',measurement_type)
                
                if not(isempty(temp_distances))
                
                    temp_meas(j) = (max(temp_distances) - min(temp_distances)) ...
                        / mean_distance;
                
                end
                
            else
                
                fprintf('Measurement type not properly specified')
                
            end
            

        end
    end
    
    %------------------------------------
    % record COV of triangle edge length
    %------------------------------------
    all_DATA{t} = temp_meas;
    
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
    namestr{t} = namestr{t}(1:end-4);           % trim off .tif 
    namestr{t} = strrep(namestr{t},'_',' ');    % remove any '_' characters
    
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

all_DATA_per_geno = cell(1,length(genotypes));

% loop through images
for t = 1:num_meas
    
    % find genotype of current image
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    all_DATA_per_geno{LOCB} = [all_DATA_per_geno{LOCB},all_DATA{t}'];
    
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
temp_mean_DATA = [];
for t = 1:length(genotypes)
    temp_mean_DATA(t) = mean(all_DATA_per_geno{t},'omitnan');
end

%--------------------------------------------------------------------------
% sort by mean fano factor
%--------------------------------------------------------------------------
sort_order = zeros(length(genotypes),2);
sort_order(:,1) = [1:length(genotypes)];
sort_order(:,2) = temp_mean_DATA;
sort_order = sortrows(sort_order,2);

sorted_aggregate_DATA = cell(length(genotypes),1);
sorted_genotypes = cell(length(genotypes),1);
for t = 1:length(genotypes)
    sort_ind = sort_order(t,1);
    sorted_genotypes{t} = genotypes{sort_ind};
    sorted_genotype_labels(t) = genotype_labels(sort_ind);
    sorted_aggregate_DATA{t} = all_DATA_per_geno{sort_ind};
end
% convert list of sorted genotypes to string array
sorted_genotypes = string(sorted_genotypes);

%--------------------------------------------------------------------------
% find mean and std of each genotype
%--------------------------------------------------------------------------

%---------------
% sorted by mean
%---------------

sorted_mean_DATA = [];
sorted_std_DATA = [];
for t = 1:length(genotypes)
    sorted_mean_DATA(t) = mean(sorted_aggregate_DATA{t},'omitnan');
    sorted_std_DATA(t) = std(sorted_aggregate_DATA{t},'omitnan');
end

%---------
% unsorted
%---------

mean_DATA = [];
std_DATA = [];
for t = 1:length(genotypes)
    mean_DATA(t) = mean(all_DATA_per_geno{t},'omitnan');
    std_DATA(t) = std(all_DATA_per_geno{t},'omitnan');
end

%--------------------------------------------------------------------------
% transfer cells into matrix
%--------------------------------------------------------------------------

%--------------------------------
% find max length of measurements
%--------------------------------
max_length = 0;
for t = 1:length(genotypes)
    % find max length
    curr_length = length(all_DATA_per_geno{t});
    if curr_length > max_length
        max_length = curr_length;
    end
end


%---------------
% sorted by mean
%---------------

sorted_DATA_matrix_form = nan(max_length,length(genotypes));
for t = 1:length(genotypes)
    sorted_DATA_matrix_form(1:length(sorted_aggregate_DATA{t}),t) = sorted_aggregate_DATA{t};
end

%---------
% unsorted
%---------

DATA_matrix_form = nan(max_length,length(genotypes));
for t = 1:length(genotypes)
    DATA_matrix_form(1:length(all_DATA_per_geno{t}),t) = all_DATA_per_geno{t};
end


%--------------------------------------------------------------------------
% pull out target genotypes for comparison
%--------------------------------------------------------------------------

count = 0;
target_meas = {};

for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        target_meas{count} = DATA_matrix_form(:,j);
        
    end
end

%--------------------------------------------------------------------------
% check to make sure we have successfully found two genotypes to compare
%--------------------------------------------------------------------------

if length(target_meas) ~= 2
    disp('\n')
    disp("Please double check you have selected two genotypes for comparison from the variable 'genotype_code'. They should be spelled them the same as they are in 'genotype_code'")
    return
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% perform statistical comparison

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

[p,h] = ranksum(target_meas{1},target_meas{2});

if h == 1
    decision = "significant";
else
    decision = "insignificant";
end

disp('\n')
disp(strcat("The difference between genotypes ",target_genotypes(1)," and ",target_genotypes(2)," is ",decision," with a p-value of ",num2str(p)))
disp(strcat("The -log10 p-value is: ", num2str(-log10(p))))

