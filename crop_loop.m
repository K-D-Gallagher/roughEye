
% crop factor
crop_factor = 512;

% specify directories for input and output
input_dir = dir('/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/cropped_images_768 copy/*tif');
output_path = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/cropped_processed copy/';

% loop through directory
for j = 1:length(input_dir)
    
    % read current image
    curr_im = imread(fullfile(input_dir(j).folder, input_dir(j).name));
    
    % display current image
    f1 = figure(1);
    imshow(curr_im);
    
    % choose and record center pixel of cropped region
    [xTemp,yTemp] = getpts(f1);
    xTemp=round(xTemp); yTemp=round(yTemp);
    
    % create cropped image,make grayscale, histEq
    crop_im = curr_im(yTemp-crop_factor:yTemp+crop_factor,xTemp-crop_factor:xTemp+crop_factor,:);
%     bw_im = rgb2gray(crop_im);
%     eq_im = histeq(bw_im);
    
    imwrite(crop_im,strcat(output_path,input_dir(j).name));
    
end





