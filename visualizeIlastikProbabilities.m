function visualizeIlastikProbabilities(raw_images,ilastik_probabilities,save)

% loop through number of images
for i = 1:size(raw_images,4)
    
    % show the ith image
    imshowpair(raw_images(:,:,:,i),ilastik_probabilities(:,:,i),'ColorChannels','red-cyan' )
    pause(0.5)
    
    if save
        
    end
    
end