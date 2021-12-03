function [omma_centroids] = loadCentroidsFromFile(filepath)

% set directory for data
image_dir = dir(filepath);
folder = image_dir.folder;

% read in first image so that we can check the dimensionality
test_image = imread(fullfile(folder,image_dir(1).name));
x_resolution = size(test_image,1);  % size of x dimension
y_resolution = size(test_image,2);  % size of y dimension
num_files = length(image_dir);      % number of images

% initalize matrix to store raw images - needs to be unit 8
% raw_images = uint8(zeros(x_resolution,y_resolution,t_resolution));
omma_centroids_disp = uint8(zeros(x_resolution,y_resolution,num_files));

       
% loop through time
for i = 1:num_files
    
    % read image
    curr_im = imread(fullfile(folder,image_dir(i).name));
    
    % store results of watershed segmentation
    omma_centroids_disp(:,:,i) = curr_im;
    
end

dilate_omma_centroids_disp = zeros(size(omma_centroids_disp));
se = strel('disk',2);
for t = 1:size(omma_centroids_disp,3)
    dilate_omma_centroids_disp(:,:,t) = imdilate(omma_centroids_disp(:,:,t),se);
end

% find hand corrected centroids

omma_centroids = cell(size(dilate_omma_centroids_disp,3),1);

for t = 1:length(omma_cent)
    [x,y] = find(omma_centroids_disp(:,:,t));
    omma_centroids{t} = [y, x];
end
