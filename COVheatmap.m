
function COVheatmap(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors,...
    raw_images,save_individual_images,save_movies, use_local_scaling,use_global_scaling)


%--------------------------------------------------------------------------
% make directory
%--------------------------------------------------------------------------
close all
% make folder
base_path = expInfo.filepath_output;
masterpath_out = strcat(base_path,'/COV_heatmap/');
[ status, msg ] = mkdir(masterpath_out);
if status == 0
    msg
end


%--------------------------------------------------------------------------
% Set up the movie structure.
%--------------------------------------------------------------------------
% Preallocate movie, which will be an array of structures.
% First get a cell array with all the frames.
numberOfFrames = size(raw_images,4);
allTheFrames = cell(numberOfFrames,1);
vidHeight = size(raw_images,1);
vidWidth = size(raw_images,2);
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
% Next get a cell array with all the colormaps.
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
% Now combine these to make the array of structures.
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
set(gcf, 'renderer', 'zbuffer');


num_meas = size(raw_images,4);


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
%     genotype{t} = genotype{t}{1};
    genotype{t} = genotype{t}{2};
    
end




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% measurements per ommatidia

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



% find max and min COV for this dataset
local_max = [];
local_min = [];
for t = 1:size(raw_images,4)
    
    local_max = [local_max, max(dist_COV{t})];
    local_min = [local_min, min(dist_COV{t})];
    
    
end

% find max and min COV for this dataset
glob_max = 0;
glob_min = 10000;
for t = 1:size(raw_images,4)
    
    curr_max = max(dist_COV{t});
    curr_min = min(dist_COV{t});
    
    if curr_max > glob_max
        glob_max = curr_max;
    end
    
    if curr_min < glob_min
        glob_min = curr_min;
    end
    
end


% make heatmap of COV
global_conversion_heatmap = round(linspace(glob_min,glob_max,100),1);
global_heatmap = parula(100);


for t = 1:length(clean_omma_centroids)
    
    t
    % show raw
    imshow(imadjust(raw_images(:,:,1,t)))
    hold on
    
    % redefine centroid list as two separate lists, one for x and one for y
    x = clean_omma_centroids{t}(:,1);
    y = clean_omma_centroids{t}(:,2);
    
    
    % make heatmap of COV
    local_conversion_heatmap = round(linspace(local_min(t),local_max(t),100),1);
    local_heatmap = parula(100);
    
    % loop through ommatidia - the identity of ommatidia is defined by their
    % position within 'clean_omma_centroids'
    for j = 1:size(x,1)
        
        % make sure its not a boundary point
        if not(ismember(j,boundary_cent{t}))
            
            
            % find COV for current ommatidia
            curr_COV = dist_COV{t}(j);
            
            if not(isnan(curr_COV))
                
                if use_local_scaling
                    
                    [~,LOCB] = ismember(round(curr_COV,1),local_conversion_heatmap);
                    
                    plot(x(j),y(j),'Marker','.','Color',local_heatmap(LOCB,:),'MarkerSize',40)
                    
                    %
                elseif use_global_scaling
                    
                    [~,LOCB] = ismember(round(curr_COV,1,global_conversion_heatmap));
                    if LOCB > 0
                        plot(x(j),y(j),'Marker','.','Color',global_heatmap(LOCB,:),'MarkerSize',40)
                    end
                    
                    
                end
                
            end
            
            
        end
    end
    
    if use_local_scaling
    
        cmap = parula(11);
        lbl =  {num2str(local_conversion_heatmap(1)),...
        num2str(local_conversion_heatmap(10)),num2str(local_conversion_heatmap(20)),...
        num2str(local_conversion_heatmap(30)),num2str(local_conversion_heatmap(40)),...
        num2str(local_conversion_heatmap(50)),num2str(local_conversion_heatmap(60)),...
        num2str(local_conversion_heatmap(70)),num2str(local_conversion_heatmap(80)),...
        num2str(local_conversion_heatmap(90)),num2str(local_conversion_heatmap(100))};

        for ii = 1:size(cmap,1)
            p(ii) = patch(NaN, NaN, cmap(ii,:));
        end
        legend(p, lbl);
        
    elseif use_global_scaling
        
        cmap = parula(11);
        lbl =  {num2str(global_conversion_heatmap(1)),...
        num2str(global_conversion_heatmap(10)),num2str(global_conversion_heatmap(20)),...
        num2str(global_conversion_heatmap(30)),num2str(global_conversion_heatmap(40)),...
        num2str(global_conversion_heatmap(50)),num2str(global_conversion_heatmap(60)),...
        num2str(global_conversion_heatmap(70)),num2str(global_conversion_heatmap(80)),...
        num2str(global_conversion_heatmap(90)),num2str(global_conversion_heatmap(100))};

        for ii = 1:size(cmap,1)
            p(ii) = patch(NaN, NaN, cmap(ii,:));
        end
        legend(p, lbl);
    end
    
%     text(20,20,genotype{t},'Color','red','FontSize',20,'FontWeight','bold')
    
    hold off
    drawnow
    
    thisFrame = getframe(gca);
	tempFrame = uint8(zeros(vidHeight,vidWidth,3));
    tempFrame(:,:,1) = imresize(thisFrame.cdata(:,:,1),[vidHeight,vidWidth]);
    tempFrame(:,:,2) = imresize(thisFrame.cdata(:,:,2),[vidHeight,vidWidth]);
    tempFrame(:,:,3) = imresize(thisFrame.cdata(:,:,3),[vidHeight,vidWidth]);
    thisFrame.cdata = tempFrame;
	myMovie(t) = thisFrame;
    
    
    clf()
    
end

close all

%--------------------------------------------------------------------------
% save individual frames to file
%--------------------------------------------------------------------------
if save_individual_images
    
    disp('\n')
    disp("Saving frame:  ")
    
    % make folder
    filepath_out = strcat(masterpath_out,'/individualFrames/');
    [ status, msg ] = mkdir(filepath_out);
    if status == 0
        msg
    end
    
    % loop through and save each frame
    for j = 1:numberOfFrames
        
        %-------------------------------
        % assemble frame w/ text overlay
        %-------------------------------
        
        % assemble side-by-side image of raw + raw w/ ilastik probabilities
        new_frame = myMovie(j).cdata;
        
        % write genotype ontop of image
        position = [15 15];
        new_frame = insertText(new_frame,position,char(expInfo.genotypes_full{j}),'FontSize',30,...
            'BoxColor','white','TextColor','black');
        
        %-------------
        % save to file
        %-------------
        file_name = strcat(expInfo.filenames(j)," ",'latticeTriangulation.jpeg');
        imwrite(new_frame,fullfile(filepath_out,file_name));
        
        %----------------------------------
        % display counter in command window
        %----------------------------------
        if j > 1
            for q=0:log10(j-1)
              fprintf('\b'); % delete previous counter display
            end
            fprintf('%d',j)
        end
        
    end
    
end



%--------------------------------------------------------------------------
% save movies for individual genotypes
%--------------------------------------------------------------------------

% loop through list of genotypes to save
if save_movies
    
    for m = 1:length(genotype_code)
        
        if m == 1
            disp('\n')
        end
        
        disp(strcat("Saving movie for ", genotype_code(m), " genotype"))
        
        %-------------
        % render movie
        %-------------
        
        % pull out indices matching current target genotype
        ind = find(expInfo.genotypes_code == genotype_code(m));
        
        % make movie of this individual genotype
        temp_movie = uint8(zeros(vidHeight,vidWidth,3,length(ind)));
        
        % loop through indices that contain target genotype
        for z = 1:length(ind)
            
            j = ind(z);
            
            %-------------------------------
            % assemble frame w/ text overlay
            %-------------------------------
            
            % assemble side-by-side image of raw + raw w/ ilastik probabilities
            new_frame = myMovie(j).cdata;
            
            % write genotype ontop of image
            position = [15 15];
            new_frame = insertText(new_frame,position,char(expInfo.genotypes_full{j}),'FontSize',30,...
                'BoxColor','white','TextColor','black');
            
            %--------------------
            % store current frame
            %--------------------
            temp_movie(:,:,:,z) = new_frame;
        end
        
        %-----------
        % save movie
        %-----------
        
        % use genotype for movie name
        baseFileName = genotype_code(m);
        
        % assemble full file name
        fullFileName = fullfile(masterpath_out, baseFileName);
        
        % record video
        writerObj = VideoWriter(fullFileName,'MPEG-4');
        writerObj.FrameRate = 1;
        open(writerObj);
        writeVideo(writerObj,temp_movie)
        close(writerObj);
        
    end
end
