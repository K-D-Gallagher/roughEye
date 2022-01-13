function [] = ...
    MannWhitney_COV(expInfo,genotypes,clean_omma_centroids,delaunay_neighbors,...
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

% measurements per ommatidia and then per eye

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

boundary_cent = cell(length(clean_omma_centroids),1);

% initialize cell array for recording COV factor of omma
dist_COV = cell(length(clean_omma_centroids),1);
dist_COV_mean = zeros(length(clean_omma_centroids),1);

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

aggregate_interR8_COV = cell(1,length(genotypes));

% loop through images
for t = 1:num_meas
    
    % find genotype of current image
    [~,LOCB] = ismember(genotype{t},genotypes);
    
    aggregate_interR8_COV{LOCB} = [aggregate_interR8_COV{LOCB},dist_COV{t}'];
    
end


%--------------------------------------------------------------------------
% pull out target genotypes for comparison
%--------------------------------------------------------------------------

count = 0;
target_COV = {};

for j = 1:length(genotypes)
    
    % check if current genotype in 'sorted_genotypes' is part of our list
    % of genotypes to plot
    [LIA,~] = ismember(genotypes(j),target_genotypes);
    
    if LIA
        
        count = count + 1;
        target_COV{count} = aggregate_interR8_COV{j};
        
    end
end

%--------------------------------------------------------------------------
% check to make sure we have successfully found two genotypes to compare
%--------------------------------------------------------------------------

if length(target_COV) ~= 2
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

[p,h] = ranksum(aggregate_interR8_COV{1},aggregate_interR8_COV{2});

if h == 1
    decision = "significant";
else
    decision = "insignificant";
end

disp('\n')
disp(strcat("The difference between genotypes ",target_genotypes(1)," and ",target_genotypes(2)," is ",decision," with a p-value of ",num2str(p)))

