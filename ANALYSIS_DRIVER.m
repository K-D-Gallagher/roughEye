
% segment ommatidia from brightfield images of Drosophila adult eyes,
% define ommatidial lattice topology, and quantify roughness phenotype
% K-D-Gallagher https://github.com/K-D-Gallagher 2021

%----------------------
% Software title ideas:
%----------------------
% Rough-EyeDERS
% Rough Eye Direct edge recovery system

% Eye Sicle
% Eye structure from images containing lattice edges

% REAP
% Rough Eye Analysis Program

% REALIZE
% Rough Eye Analysis via Lattice Orginization Estimate


%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% CROP DATA BEFORE USING ILASTIK
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% filepath = uigetdir('','Select directory containing data');
% filepath = strcat(filepath,'/');
filepath = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/males/';


% This will create a new folder called 'cropped/' within the same
% subdirectory as the folder containing your images
[cropInfo] = initializeCropping(filepath,'/');
[pre_crop_images] = loadDataCrop(filepath);

%--------------------
% launch cropTool GUI
%--------------------

% NOTE: Once you've finished cropping images, click the 'Save to File'
% button to save your cropped files

cropTool

%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% initialize experiment and read in data
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% NOTE: the ilastik files (.h5 extension) should remain in the same folder
% as the raw images. This is the default behavior from Ilastik, so this
% should not be an issue.

filepath = uigetdir('','Select directory containing data');
filepath = strcat(filepath,'/');
% filepath = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/cropped_images_768/';

experimentID = '21.12.03 test'; % NO UNDERSCORES (_) please
genotype_code = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];
full_genotype = genotype_code;
% full_genotype = ...
%     ["mir7delta positive control" ...
%     "UAS-Myc Control" ...
%     "MycOverExpr" ...
%     "UAS-Myc Control" ...
%     "GMR>Gal4 Control" ...
%     "GMR>Gal4 Control mir7delta" ...
%     "MycOverExpr mir7delta"];

[expInfo] = initializeExperiment(experimentID,genotype_code,full_genotype,filepath);
[raw_images, ilastik_probabilities] = loadData(filepath);



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% OPTIONAL: visualize ilastik classification overlaid on raw images
% raw images = green hue
% ilastik classification = purple huge
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%----------------------------------------------------------------------
% Do you want to save each individual frame to file? (FOR ALL GENOTYPES)
%----------------------------------------------------------------------
save_individual_images = true;
save_movies = true;

% call function
visualizeIlastikProbabilities(expInfo,genotype_code,raw_images,ilastik_probabilities, ...
    save_individual_images, save_movies)



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% Segment ommatidia from Ilastik classification
%
%
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%
%
% hand correct using GUI
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Load the GUI for hand-correcting segmentation. The GUI will allow you to
% export the 'omma_centroids' variable

ommatidiaSeg



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% SAVE WORKSPACE TO FILE
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

save(strcat(expInfo.filepath_output,experimentID,'_workspace','.mat'))



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% IF MATLAB CRASHES WHILE USING GUI:
% reload omma_centroids from backup file
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


load('backup_omma_centroids.mat')



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% visualize segmented ommatidia
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parameters for displaying segmented ommatidia centroids

% For additional options, see linespec:
% https://www.mathworks.com/help/matlab/creating_plots/specify-line-and-marker-appearance-in-plots.html

marker_type = 'o';         % alt options: o, +, *, x, ., etc
mark_color = 'cyan';     % alt options: red, green, blue, cyan, magenta, yellow, black, white
marker_size = 14;
line_width = 2;
%--------------------------------------------------------------------------

save_individual_images = true;
save_movies = true;

visualizeSegmentedOmmatidia(expInfo,genotype_code,raw_images,omma_centroids,...
    marker_type,mark_color,marker_size, line_width, ...
    save_individual_images,save_movies)



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% triangulate points and find neighbors
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

[omma_triangles,clean_omma_centroids,delaunay_neighbors] ... 
    = triangulateAndFindNeighbors(omma_centroids);



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% visualize triangulation
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


save_individual_images = true;
save_movies = true;

visualizeTriangulation(expInfo,genotype_code,raw_images,omma_triangles,...
    save_individual_images,save_movies)


%% visualize triangle areas

%     %-----------------------------
%     % visualization triangle areas
%     %-----------------------------
%     
%     % colormap of triangle area
%     max_area = round(max(triangle_areas{t}));
%     min_area = round(min(triangle_areas{t}));
%     area_range = min_area:max_area;
%     colormap = parula(length(min_area:max_area));
%     colors = zeros(length(triangle_areas{t}),3);
%     for q = 1:length(triangle_areas{t})
%         [~,LOCB] = ismember(round(triangle_areas{t}(q)),area_range);
%         colors(q,:) = colormap(LOCB,:);
%     end
%     
%     imshow(raw_images(:,:,:,t))
%     hold on
%     triplot(edge_clean_triangulation{t},'LineWidth',2,'Color','cyan')
%     patch('vertices', clean_omma_centroids{t},'Faces', edge_clean_delaunay{t}, ...
%         'FaceColor','flat', 'FaceVertexCData', colors, 'CDataMapping', 'direct', ...
%         'FaceAlpha', 1);
%     hold off
%     drawnow
%     filename = ['/Users/kevin/Documents/MATLAB/forSha/media/triangle_areas/T=' num2str(t,'%03i') '.png'];
%     print(gcf,'-dpng',filename)



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% visualize neighbors of each ommatidia at the specified time point
% option to save image to file
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% NOTE: this will make a movie of a single image (not multiple images from
% a given genotype), displaying neighbors of each ommatidia from your
% chosen image


which_file = 2;
save = false;

visualizeOmmaNeighbors(raw_images,clean_omma_centroids,...
    omma_triangles,delaunay_neighbors,which_file,save)


    
%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% COV per Genotype
%
% Coefficient of Variation (COV) in inter-ommatidial-distance aggregated
% according to genotype and then averaged
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% aggregate COV measurements of ommatidia from multiple images according to
% their genotype

%-------------------------
% WHICH GENOTYPES TO PLOT?
%-------------------------
target_genotypes = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];
% target_genotypes = ["q5" "q9" "q11" "q12" "q13" "q14"];

%-------------
% PLOT OPTIONS
%-------------
plot_style = "voilin plot";         % alt options: "mean & std", "box plot", "violin plot"
ascending_mean = true;              % do you want to sort genotypes by ascending mean?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?
x_axis_text_angle = 0;              % choose angle of x-axis text (rotate so they don't overlap)
plot_title = "Coefficient of Variation (COV) of inter-R8-distance";
title_size = 20;
x_label = "Genotypes";
y_label = "Coefficient of Variation (COV)";
axes_label_size = 16;

covPerGeno(filepath,genotype_code,clean_omma_centroids,delaunay_neighbors, ...
    target_genotypes, plot_style,ascending_mean,genotype_labels, ...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size)



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% COV per Image
%
% Average Coefficient of Variation (COV) in inter-ommatidial-distance
% calculated for each image and then averaged across images according
% to genotype
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%-------------------------
% WHICH GENOTYPES TO PLOT?
%-------------------------
% target_genotypes = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];
target_genotypes = ["q5" "q9" "q11" "q12" "q13" "q14"];

%-------------
% PLOT OPTIONS
%-------------
plot_style = "violin plot";
ascending_mean = true;              % do you want to sort genotypes by ascending mean?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?
x_axis_text_angle = 0;              % choose angle of x-axis text (rotate so they don't overlap)
plot_title = "Coefficient of Variation (COV) of inter-R8-distance";
title_size = 20;
x_label = "Genotypes";
y_label = "Coefficient of Variation (COV)";
axes_label_size = 16;

covPerImage(filepath,genotype_code,clean_omma_centroids,delaunay_neighbors, ...
    target_genotypes, plot_style,ascending_mean,genotype_labels, ...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size)



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% analyze inter-R8-distance (not coefficient of variation) via raw
% distributions and eCDFs
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



interR8distancePerGeno(filepath,genotype_code,clean_omma_centroids,delaunay_neighbors)



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



