function [hand_correct_omma_cent] = manualCorrectFromFile(filepath)

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
hand_correct_omma_cent_disp = uint8(zeros(x_resolution,y_resolution,num_files));

       
% loop through time
for i = 1:num_files
    
    % read image
    curr_im = imread(fullfile(folder,image_dir(i).name));
    
    % store results of watershed segmentation
    hand_correct_omma_cent_disp(:,:,i) = curr_im;
    
end

dilate_hand_correct_omma_cent_disp = zeros(size(hand_correct_omma_cent_disp));
se = strel('disk',2);
for t = 1:size(hand_correct_omma_cent_disp,3)
    dilate_hand_correct_omma_cent_disp(:,:,t) = imdilate(hand_correct_omma_cent_disp(:,:,t),se);
end

% find hand corrected centroids

hand_correct_omma_cent = cell(length(omma_cent),1);

for t = 1:length(omma_cent)
    [x,y] = find(hand_correct_omma_cent_disp(:,:,t));
    hand_correct_omma_cent{t} = [y, x];
end
