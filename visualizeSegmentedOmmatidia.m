function visualizeSegmentedOmmatidia(expInfo,raw_images,omma_centroids,...
    marker_type, marker_color, marker_size, line_width, ...
    save_individual_images, save_movies)

%--------------------------------------------------------------------------
% make directory
%--------------------------------------------------------------------------
close all
% make folder
base_path = expInfo.filepath_output;
masterpath_out = strcat(base_path,'/ommatidiaSegmentation/');
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



%--------------------------------------------------------------------------
% loop through number of images and record movie
%--------------------------------------------------------------------------

disp('\n')
disp("Rendering frame:  ")

for i = 1:numberOfFrames
    
    %----------------
    % display counter
    %----------------
    if i > 1
        for q=0:log10(i-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',i)
    end
    
    %-----------------------------------------
    % visualize triangles on top of raw images
    %-----------------------------------------
    imshow(raw_images(:,:,:,i))
    hold on
    for ii = 1:length(omma_centroids{i})
        plot(omma_centroids{i}(ii,1),omma_centroids{i}(ii,2),...
            marker_type,'Color',marker_color,'MarkerSize',marker_size, ...
            'LineWidth',line_width)
    end
    
	thisFrame = getframe(gca);
	myMovie(i) = imresize(thisFrame,[vidHeight,vidWidth]);
    hold off
    
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
        text_str = strcat(expInfo.genotypes_full(j), " ", expInfo.sex(j));
        position = [15 15];
        new_frame = insertText(new_frame,position,text_str,'FontSize',30,...
            'BoxColor','white','TextColor','black');
        
        %-------------
        % save to file
        %-------------
        file_name = strcat(expInfo.filenames(j)," ",'ommatidiaSegmentation.jpeg');
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
for m = 1:length(save_movies)
    
    disp(strcat("Saving movie for ", save_movies(m), " genotype"))
    
    %-------------
    % render movie
    %-------------
    
    % pull out indices matching current target genotype
    ind = find(expInfo.genotypes_code == save_movies(m));
    
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
        text_str = strcat(expInfo.genotypes_full(j), " ", expInfo.sex(j));
        position = [15 15];
        new_frame = insertText(new_frame,position,text_str,'FontSize',30,...
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
    baseFileName = save_movies(m);
    
    % assemble full file name
    fullFileName = fullfile(masterpath_out, baseFileName);
    
	% record video
	writerObj = VideoWriter(fullFileName,'MPEG-4');
    writerObj.FrameRate = 1;
	open(writerObj);
    writeVideo(writerObj,temp_movie)
    close(writerObj);
    
end