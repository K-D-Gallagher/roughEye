function [raw_images] = loadDataCrop(filepath)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% read in raw images
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% NOTE: this will read in all the files with the extension specificied
% after the * in the filed path. Therefore, only include images in this
% folder that you wish to load into matlab. They will be stored in a 3d
% matrix whose third dimension is 'number of images'

% NOTE: if you load in multiple images, they must all be the same size

% set directory for data
image_dir = dir(strcat(filepath,'*tif'));
folder = image_dir.folder;

%--------------------------------------------------------------------------
% read in first image so that we can check the dimensionality
%--------------------------------------------------------------------------
test_image = imread(fullfile(folder,image_dir(1).name));
x_resolution = size(test_image,1);  % size of x dimension
y_resolution = size(test_image,2);  % size of y dimension
num_files = length(image_dir);      % number of images

% initalize matrix to store raw images - needs to be unit 8
% raw_images = uint8(zeros(x_resolution,y_resolution,t_resolution));
raw_images = uint8(zeros(x_resolution,y_resolution,3,num_files));

disp('Reading in raw images  ')
       
% loop through files
for i = 1:num_files
    
    if i > 1
        for j=0:log10(i-1)
          fprintf('\b'); % delete previous counter display
        end
        fprintf('%d',i)
    end
    
    % read image
    curr_im = imread(fullfile(folder,image_dir(i).name));
    
    % store results of watershed segmentation
    raw_images(:,:,:,i) = curr_im;
    
end











