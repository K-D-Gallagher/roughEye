function visualizeTriangulation(expInfo,genotype_code,raw_images,omma_triangles,...
    save_individual_images,save_movies)

%--------------------------------------------------------------------------
% make directory
%--------------------------------------------------------------------------
close all
% make folder
base_path = expInfo.filepath_output;
masterpath_out = strcat(base_path,'/latticeTriangulation/');
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

for t = 1:numberOfFrames
    
    %----------------
    % display counter
    %----------------
    if t > 1
        for q=0:log10(t-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',t)
    end
    
    %-----------------------------------------
    % visualize triangles on top of raw images
    %-----------------------------------------
    imshow(raw_images(:,:,:,t))
    hold on
    triplot(omma_triangles{t},'LineWidth',2,'Color','cyan')
    hold off
    
	thisFrame = getframe(gca);
	tempFrame = uint8(zeros(vidHeight,vidWidth,3));
    tempFrame(:,:,1) = imresize(thisFrame.cdata(:,:,1),[vidHeight,vidWidth]);
    tempFrame(:,:,2) = imresize(thisFrame.cdata(:,:,2),[vidHeight,vidWidth]);
    tempFrame(:,:,3) = imresize(thisFrame.cdata(:,:,3),[vidHeight,vidWidth]);
    thisFrame.cdata = tempFrame;
	myMovie(t) = thisFrame;
    
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