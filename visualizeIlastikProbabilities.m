function visualizeIlastikProbabilities(expInfo,raw_images,ilastik_probabilities, ...
    save_individual_images, genotype_movies, genotype_display)

close all

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
    
    % display counter
    if i > 1
        for q=0:log10(i-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',i)
    end
    
    % show overlaid probabilities and raw image and record
    imshowpair(imcomplement(raw_images(:,:,:,i)),uint8(ilastik_probabilities(:,:,i)) * 255)
	thisFrame = getframe(gca);
	myMovie(i) = thisFrame;
    
end

close all



%--------------------------------------------------------------------------
% save individual frames to file
%--------------------------------------------------------------------------
if save_individual_images
    
    disp('\n')
    disp("Saving frame:  ")
    
    % make folder
    base_path = expInfo.filepath_output;
    filepath_out = strcat(base_path,'/IlastikClassification/');
    [ status, msg ] = mkdir(filepath_out);
    if status == 0
        msg
    end
    
    % loop through and save each frame
    for j = 1:numberOfFrames
        
        % assemble side-by-side image of raw + raw w/ ilastik probabilities
        new_frame = uint8(zeros(vidHeight,vidWidth*2,3));
        new_frame(:,1:vidWidth,:) = raw_images(:,:,:,j);
        new_frame(:,vidWidth+1:end,:) = myMovie(j).cdata;
        
        % write genotype ontop of image
        text_str = expInfo.genotypes_full(j);
        position = [15 15];
        new_frame = insertText(new_frame,position,text_str,'FontSize',30,...
            'BoxColor','white','TextColor','black');
        
        % save to file
        file_name = strcat(expInfo.filenames(j)," ",'ilastikClassification.jpeg');
        imwrite(new_frame,fullfile(strcat(filepath_out,'/individualFrames/'),file_name));
        
        % display counter in command window
        if j > 1
            for q=0:log10(j-1)
              fprintf('\b'); % delete previous counter display
            end
            fprintf('%d',j)
        end
        
    end
    
end






%--------------------------------------------------------------------------
% See if they want to save the movie to an avi file on disk.
%--------------------------------------------------------------------------
promptMessage = sprintf('Do you want to save this movie to disk?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'yes')
	% Get the name of the file that the user wants to save.
	% Note, if you're saving an image you can use imsave() instead of uiputfile().
	startingFolder = pwd;
	defaultFileName = {'*.avi';'*.mp4';'*.mj2'}; %fullfile(startingFolder, '*.avi');
	[baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
	if baseFileName == 0
		% User clicked the Cancel button.
		return;
	end
	fullFileName = fullfile(folder, baseFileName);
	% Create a video writer object with that file name.
	% The VideoWriter object must have a profile input argument, otherwise you get jpg.
	% Determine the format the user specified:
	[folder, baseFileName, ext] = fileparts(fullFileName);
	switch lower(ext)
		case '.jp2'
			profile = 'Archival';
		case '.mp4'
			profile = 'MPEG-4';
		otherwise
			% Either avi or some other invalid extension.
			profile = 'Uncompressed AVI';
	end
	writerObj = VideoWriter(fullFileName, profile);
	open(writerObj);
	% Write out all the frames.
	numberOfFrames = length(myMovie);
	for frameNumber = 1 : numberOfFrames 
	   writeVideo(writerObj, myMovie(frameNumber));
	end
	close(writerObj);
	% Display the current folder panel so they can see their newly created file.
	cd(folder);
	filebrowser;
	message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
	uiwait(helpdlg(message));
else
	uiwait(helpdlg('Done with demo!'));
end



