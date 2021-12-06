
% segment ommatidia from brightfield images of Drosophila adult eyes,
% define ommatidial lattice topology, and quantify roughness phenotype
% K-D-Gallagher https://github.com/K-D-Gallagher 2021

% Rough Eye Direct edge recovery system
% Rough-EyeDERS

% Eye structure from images containing lattice edges
% Eye Sicle

% REAP
% Rough Eye Analysis Program


%%
%--------------------------------------------------------------------------
% CROP LOOP?
%--------------------------------------------------------------------------





%%
%--------------------------------------------------------------------------
% record experiment metadata
%--------------------------------------------------------------------------

% NOTE: the ilastik files (.h5 extension) should remain in the same folder
% as the raw images. This is the default behavior from Ilastik, so this
% should not be an issue.
filepath = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/cropped_images_768/';
% [baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');

experimentID = '21.12.03 test'; % NO UNDERSCORES (_) please
genotype_code = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];
full_genotype = ...
    ["mir7" ...
    "UAS-Myc_control" ...
    "MycOverExpr" ...
    "UAS-Myc_control" ...
    "GMR>Gal4_control" ...
    "GMR>Gal4_control mir7delta" ...
    "MycOverExpr mir7delta"];

[expInfo] = initializeExperiment(experimentID,genotype_code,full_genotype,filepath);


%%
%--------------------------------------------------------------------------
% read in data
%--------------------------------------------------------------------------

% Here, we will read in our raw data and the corresponding pixel
% classification files exported from Ilastik.

[raw_images, ilastik_probabilities] = loadData(filepath);



%%
%--------------------------------------------------------------------------
% visualize ilastik classification overlaid on raw images
%--------------------------------------------------------------------------
% raw images = green hue
% ilastik classification = purple huge

% do you want to save each individual frame to file? (FOR ALL GENOTYPES)
save_individual_images = true;

% 
genotype_movies = [];

genotype_display = [];

visualizeIlastikProbabilities(expInfo,raw_images,ilastik_probabilities, ...
    save_individual_images, ...
    genotype_movies, ...
    genotype_display)



%%
%--------------------------------------------------------------------------
% Segment ommatidia
%--------------------------------------------------------------------------

% Compute the centroids of the ommatidia we've detected via pixel
% classification in Ilastik. We achieve this by thresholding the pixel
% classification probabilities exported from Ilastik.
[omma_centroids,omma_area] = initialSeg(ilastik_probabilities,0.3);

% remove objects with less area than given threshold (second argument)
[omma_centroids] = sizeThreshOmma(omma_centroids,omma_area,0.7);

% merge close ommatidia by dilating and then finding centroid again (second
% argument dilation kernal width)
[omma_centroids] = mergeCloseOmma(omma_centroids,4);



%% 
%--------------------------------------------------------------------------
% OPTIONAL: READ IN OMMATIDIA CENTROIDS FROM FILE
% This is an alternative to segmenting ommatidia from Ilastik pixel
% classification
%--------------------------------------------------------------------------

filepath = '/Users/kevin/Documents/MATLAB/forSha/omma_cent_disp_handCorrect_first75/*tif';
[omma_centroids] = loadCentroidsFromFile(filepath);



%%
%--------------------------------------------------------------------------
% hand correct using GUI
%--------------------------------------------------------------------------

% Load the GUI for hand-correcting segmentation. The GUI will allow you to
% export the 'omma_centroids' variable

omma_cent = omma_centroids;
adultOmmatidiaSeg



%%
%--------------------------------------------------------------------------
% OPTIONAL: define ROI for each image and throw away centroids not in ROI
%--------------------------------------------------------------------------

[omma_centroids] = customROI(raw_images,omma_centroids);



%%
%--------------------------------------------------------------------------
% triangulate points and find neighbors
%--------------------------------------------------------------------------


[omma_triangles,clean_omma_centroids,delaunay_neighbors] ... 
    = triangulateAndFindNeighbors(omma_centroids);



%%
%--------------------------------------------------------------------------
% visualize neighbors of each ommatidia at the specified time point
% option to save image to file
%--------------------------------------------------------------------------

which_file = 2;
save = false;

visualizeOmmaNeighbors(raw_images,clean_omma_centroids,...
    omma_triangles,delaunay_neighbors,which_file,save)


    
%%
%--------------------------------------------------------------------------
% calculate Coefficient of Variation across all ommatidia of each genotype
%--------------------------------------------------------------------------

covPerGeno(filepath,genotypes,clean_omma_centroids,delaunay_neighbors)


%%
%--------------------------------------------------------------------------
% calculate Coefficient of Variation across ommatidia for each image and
% then average multiple images according to genotype
%--------------------------------------------------------------------------

covPerImage(filepath,genotypes,clean_omma_centroids,delaunay_neighbors)



%%
%--------------------------------------------------------------------------
% analyze inter-R8-distance (not coefficient of variation) via raw
% distributions and eCDFs
%--------------------------------------------------------------------------

interR8distancePerGeno(filepath,genotypes,clean_omma_centroids,delaunay_neighbors)



%% eCDF -20 to 20

% MANUALLY REMOVE MIR7
alt_genotypes = [sorted_genotypes(1:3); sorted_genotypes(5:7)];
alt_normalized_interR8_dist_matrix = [normalized_interR8_dist_matrix(:,1:3), normalized_interR8_dist_matrix(:,4:7)];
alt_normalized_interR8_distance_cellArray = {normalized_interR8_distance_cellArray{1}, ...
    normalized_interR8_distance_cellArray{2}, ...
    normalized_interR8_distance_cellArray{3}, ...
    normalized_interR8_distance_cellArray{5}, ...
    normalized_interR8_distance_cellArray{6}, ...
    normalized_interR8_distance_cellArray{7}};

subplot(1,2,1)

% vs = violinplot(normalized_interR8_dist_matrix,sorted_genotypes,'ShowMean',true);     % all genotypes
vs = violinplot(alt_normalized_interR8_dist_matrix,alt_genotypes,'ShowMean',true);      % MIR7 REMOVED
title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['normalized inter-R8-distances'])
ax = gca;
ax.FontSize = 16;
camroll(-90)

subplot(1,2,2)

coloz = lines(7);


coloz2 = zeros(size(coloz));
coloz2(2,:) = coloz(6,:);
coloz2(3,:) = coloz(4,:);
coloz2(4,:) = coloz(2,:);
coloz2(5,:) = coloz(7,:);
coloz2(6,:) = coloz(5,:);
coloz2(7,:) = coloz(3,:);

hold on

% for j = 1:length(normalized_interR8_distance_cellArray)     % all genotypes
for j = 1:length(alt_normalized_interR8_distance_cellArray)     % MIR7 REMOVED
    
    if j == 1
        
        [f,x] = ecdf(normalized_interR8_distance_cellArray{j});
        p  = patchline(x,f,'edgecolor',coloz(1,:),'linewidth',10,'edgealpha',1);
        
    else
        
%         [f,x] = ecdf(normalized_interR8_distance_cellArray{j});     % all genotypes
        [f,x] = ecdf(alt_normalized_interR8_distance_cellArray{j});     % MIR7 REMOVED
        plot(x,f,'Linewidth',3,'Color',coloz2(j,:))
        
    end
    
    
end

title('MEDIAN NORMALIZED inter-R8-distance')
% legend(sorted_genotypes,'Location','Southeast')    % MIR7 REMOVED
legend(alt_genotypes,'Location','Southeast')    % MIR7 REMOVED
ylabel('F(x)')
xlabel('x = MEDIAN NORMALIZED inter-R8-distances')
ax = gca;
ax.FontSize = 16;
xlim([-20 20])

hold off

%% eCDF -20 to 20 - SPECIFY WHICH GENOTYPES WE WANT

target_genotypes = [1 2 3 5 6 7];

alt_genotypes = {};
alt_normalized_interR8_dist_matrix = [];
alt_normalized_interR8_distance_cellArray = {};
for j = 1:length(target_genotypes)
    alt_genotypes{j} = sorted_genotypes{target_genotypes(j)};
    alt_normalized_interR8_dist_matrix(:,j) = normalized_interR8_dist_matrix(:,target_genotypes(j));
    alt_normalized_interR8_distance_cellArray{j} = normalized_interR8_distance_cellArray{target_genotypes(j)};
end

subplot(1,2,1)

vs = violinplot(alt_normalized_interR8_dist_matrix,alt_genotypes,'ShowMean',true);
% title('MEDIAN NORMALIZED inter-R8-distance')
xlabel(['genotype'])
ylabel(['MEDIAN NORMALIZED inter-R8-distances'])
ax = gca;
ax.FontSize = 16;
camroll(-90)

subplot(1,2,2)

coloz = lines(7);


coloz2 = zeros(size(coloz));
coloz2(2,:) = coloz(6,:);
coloz2(3,:) = coloz(4,:);
coloz2(4,:) = coloz(2,:);
coloz2(5,:) = coloz(7,:);
coloz2(6,:) = coloz(5,:);
coloz2(7,:) = coloz(3,:);

hold on

for j = 1:length(alt_normalized_interR8_distance_cellArray)
    
    if j == 1
        
        [f,x] = ecdf(normalized_interR8_distance_cellArray{j});
        p  = patchline(x,f,'edgecolor',coloz(1,:),'linewidth',10,'edgealpha',1);
        
    else
        
        [f,x] = ecdf(alt_normalized_interR8_distance_cellArray{j});
        plot(x,f,'Linewidth',3,'Color',coloz2(j,:))
        
    end
    
    
end

% title('MEDIAN NORMALIZED inter-R8-distance')
legend(alt_genotypes,'Location','Southeast')
ylabel('F(x)')
xlabel('x = MEDIAN NORMALIZED inter-R8-distances')
ax = gca;
ax.FontSize = 16;
xlim([-20 20])

hold off


% sgtitle('MEDIAN NORMALIZED inter-R8-distance')


%%

%--------------------------------------------------------------------------
%
% plot inter-R8-distance COV factor for each image separately
%
%--------------------------------------------------------------------------
for t = 1:75
    
    clf()
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot dispersion index vs. sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,1,'position',[0.1 0.1 0.4 0.82])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    hold on
    for j = 1:75
        plot(j,sorted_distance_COV(j,2),'.','Color',plt_color_vect(j,:),'MarkerSize',40,'Linewidth',1)
    end
    plot(t,sorted_distance_COV(t,2),'ok','MarkerSize',18,'Linewidth',6)
    ylabel('Mean index of dispersion (σ^2/μ)','FontSize',18)
    xlabel('Ranked sample order','FontSize',18)
    ax=gca;
    ax.FontSize = 14;
    
    %legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'.','Color',genotype_color(1,:),'MarkerSize',40);
    h(2) = plot(NaN,NaN,'.','Color',genotype_color(2,:),'MarkerSize',40);
    h(3) = plot(NaN,NaN,'.','Color',genotype_color(3,:),'MarkerSize',40);
    h(4) = plot(NaN,NaN,'.','Color',genotype_color(4,:),'MarkerSize',40);
    h(5) = plot(NaN,NaN,'.','Color',genotype_color(5,:),'MarkerSize',40);
    h(6) = plot(NaN,NaN,'.','Color',genotype_color(6,:),'MarkerSize',40);
    h(7) = plot(NaN,NaN,'.','Color',genotype_color(7,:),'MarkerSize',40);
    legend(h, genotypes(1),genotypes(2),genotypes(3),genotypes(4),...
        genotypes(5),genotypes(6),genotypes(7),'Location','northwest',...
        'FontSize',20);
    
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corresponding eye image with lattice drawn on top
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2,'position',[0.26 0.1 0.95 0.821])
    imshow(raw_images(:,:,:,sorted_distance_COV(t,1)))
    hold on
    triplot(omma_triangles{sorted_distance_COV(t,1)},'LineWidth',2,'Color','cyan')
    hold off
    annotation('textbox',[0.51 0.43 0.5 0.5],'string',namestr{t},...
   'linestyle','none','FontSize',30,'Color','r')
    
    sgtitle('Neighbor length - mean index of dispersion','FontSize',24)
    
    drawnow
    
    
    % print to file
    filename = ['/Users/kevin/Documents/MATLAB/forSha/media/individual_image_interR8/T=' num2str(t,'%03i') '.png'];
    print(gcf,'-dpng',filename)

end

%%

%---------------------------------------
%
%
% plot unsorted fano factor per genotype
%
%
%---------------------------------------

%-----------------------
% plot inter-R8 distance
%-----------------------

subplot(1,2,1)
bar(ave_dist_fano_per_geno)                
set(gca,'xticklabel',genotypes)
hold on

er = errorbar(ave_dist_fano_per_geno,std_dist_fano_per_geno/2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%---------------------------
% plot lattice triangle area
%---------------------------

subplot(1,2,2)
bar(ave_area_fano_per_geno)                
set(gca,'xticklabel',genotypes)
hold on

er = errorbar(ave_area_fano_per_geno,std_area_fano_per_geno/2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off


%%



