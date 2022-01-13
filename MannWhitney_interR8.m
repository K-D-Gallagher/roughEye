function [] = ...
    MannWhitney_interR8(expInfo,genotypes,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes)

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
%
% pull out target genotypes for comparison
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

target_interR8_distances = {};

count1 = 0;
for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA1,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA1
        
        count1 = count1 + 1;

        target_interR8_distances{count1} = aggregate_interR8_distances{j};
        
    end
    
end

%--------------------------------------------------------------------------
% check to make sure we have successfully found two genotypes to compare
%--------------------------------------------------------------------------

if length(target_interR8_distances) ~= 2
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

[p,h] = ranksum(target_interR8_distances{1},target_interR8_distances{2});

if h == 1
    decision = "significant";
else
    decision = "insignificant";
end

disp('\n')
disp(strcat("The difference between genotypes ",target_genotypes(1)," and ",target_genotypes(2)," is ",decision," with a p-value of ",num2str(p)))





