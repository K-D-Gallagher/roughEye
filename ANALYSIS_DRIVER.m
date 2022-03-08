
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
% VERSION REQUIREMENTS
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% MATLAB Version 2018b
% Image Processing Toolbox
% Signal Processing Toolbox
% Statistics and Machine Learning Toolbox



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


%--------------------------------------------------------------------------
% define filepaths
%--------------------------------------------------------------------------

% NOTE: the ilastik files (.h5 extension) should remain in the same folder
% as the raw images. This is the default behavior from Ilastik, so this
% should not be an issue.

% filepath = uigetdir('','Select directory containing data');
% filepath = strcat(filepath,'/');
filepath = '/Users/kevin/Documents/MATLAB/forSha/raw_and_cropped_images/males/';

experimentID = '22.01.10 test'; % NO UNDERSCORES (_) please


%--------------------------------------------------------------------------
% define genotypes
%--------------------------------------------------------------------------

% NOTE: the way you write the genotypes must exactly mach the way its
% written in the filename

%---------
% old data
%---------
% genotype_code = ["mir7" "q5" "q9" "q11" "q12" "q13" "q14"];

%---------
% new data
%---------
genotype_code = ["q5" "q9" "q11" "q12" "q13" "q14"];
% genotype_code = ["Q5" "Q9" "Q11" "Q12" "Q13" "Q14"];
% genotype_code = ["Q23" "Q24" "Q29" "Q33" "Q34" "Q39"];
% genotype_code = ["Q46" "Q47" "Q50" "Q51" "Q52" "Q53"];

%----------
% OPTIONAL: write out full genotypes OR just make full_genotype == genotype_code
%----------

full_genotype = genotype_code;
% full_genotype = ...
%     ["mir7delta positive control" ...
%     "UAS-Myc Control" ...
%     "MycOverExpr" ...
%     "UAS-Myc Control" ...
%     "GMR>Gal4 Control" ...
%     "GMR>Gal4 Control mir7delta" ...
%     "MycOverExpr mir7delta"];


%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------

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
% IF MATLAB CRASHES WHILE USING GUI: reload workspace
%
% NOTE: the matlab workspace file was saved into the project directory
% located in the folder containing your input images
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

[file,path] = uigetfile;
load(fullfile(path,file),'-mat')



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
distance_cutoff = 200;
%--------------------------------------------------------------------------

save_individual_images = true;
save_movies = true;

visualizeSegmentedOmmatidia(expInfo,genotype_code,raw_images,omma_centroids,...
    marker_type,mark_color,marker_size, line_width, ...
    save_individual_images,save_movies,...
    distance_cutoff)



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
% COV heatmap
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

save_individual_images = true;
save_movies = true;

%-----------------------------------------------------------
% NOTE: choose only one (one must True and one must be false)
% this specifies whether to scale the COV colormap according
% to the current image or according to all images

use_local_scaling = true;
use_global_scaling = false;
%-----------------------------------------------------------

COVheatmap(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors,...
    raw_images,save_individual_images,save_movies, use_local_scaling, use_global_scaling)


    
%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% ANALYSIS PER GENOTYPES - data points are ommatidia
% i.e. measurements are made for individual ommatidia, which are then
% pooled according to genotype
%
% Two analysis options
%
% Option 1: 
% Coefficient of Variation (COV) in inter-ommatidial-distance
%
% Option 2:
% Min inter-ommatidial-distance / Max inter-ommatidial-distance
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Which genotypes to plot and the order in which we want to plot them:
%--------------------------------------------------------------------------
target_genotypes = ["q5" "q9" "q11" "q12" "q13" "q14"];
target_genotypes = ["q14" "q13" "q12" "q11" "q9" "q5"];

%--------------------------------------------------------------------------
% PLOT OPTIONS
%--------------------------------------------------------------------------
plot_style = "violin plot";         % alt options: "mean & std", "box plot", "violin plot"
ascending_mean = false;              % do you want to sort genotypes by ascending mean?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?

% plot_title = "Coefficient of Variation (COV) of inter-R8-distance";
% x_label = "Genotypes";
% y_label = "Coefficient of Variation (COV)";

plot_title = "Min/max inter-R8-distance";
x_label = "Genotypes";
y_label = "Min/max inter-R8-distance";

title_size = 20;
axes_label_size = 16;
x_axis_text_angle = 0;         % choose angle of x-axis text (rotate so they don't overlap)
y_axis_limit = [0 1];

distance_cutoff = 200;       % threshold, in pixels, within which we'll include ommatidia for measurement
measurement_type = 'max_min';     % OPTIONS: COV, max_min

save_csv_to_file = true;       % would you like to save a cvs document to file?
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


perGeno(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors, ...
    target_genotypes, plot_style,ascending_mean,genotype_labels, ...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size,save_csv_to_file,y_axis_limit,...
    measurement_type,distance_cutoff)



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% ANALYSIS PER EYE - data points are images
% i.e. measurements are first made for individual ommatidia, which are then
% averaged per sample (per image); averages for each image are then pooled
% according to genotype
%
% Two analysis options:
%
% Option 1: 
% Coefficient of Variation (COV) in inter-ommatidial-distance
%
% Option 2:
% Min inter-ommatidial-distance / Max inter-ommatidial-distance
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Which genotypes to plot and the order in which we want to plot them:
%--------------------------------------------------------------------------
target_genotypes = ["q5" "q9" "q11" "q12" "q13" "q14"];
target_genotypes = ["q14" "q13" "q12" "q11" "q9" "q5"];

%--------------------------------------------------------------------------
% PLOT OPTIONS
%--------------------------------------------------------------------------
plot_style = "violin plot";
ascending_mean = false;              % do you want to sort genotypes by ascending mean?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?

% plot_title = "Coefficient of Variation (COV) of inter-R8-distance";
% x_label = "Genotypes";
% y_label = "Coefficient of Variation (COV)";

plot_title = "Min/max inter-R8-distance";
x_label = "Genotypes";
y_label = "Min/max inter-R8-distance";

title_size = 20;
axes_label_size = 16;
x_axis_text_angle = 0;              % choose angle of x-axis text (rotate so they don't overlap)
y_axis_limit = [0 1];

distance_cutoff = 200;       % threshold, in pixels, within which we'll include ommatidia for measurement
measurement_type = 'max_min';     % OPTIONS: COV, max_min

save_csv_to_file = true;       % would you like to save a cvs document to file?
%--------------------------------------------------------------------------

perImage(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors, ...
    target_genotypes, plot_style,ascending_mean,genotype_labels, ...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size,save_csv_to_file,y_axis_limit,...
    measurement_type,distance_cutoff)



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

%--------------------------------------------------------------------------
% Which genotypes to plot and the order in which we want to plot them:
%--------------------------------------------------------------------------
target_genotypes = ["q5" "q9" "q11" "q12" "q13" "q14"];
target_genotypes = ["q14" "q13" "q12" "q11" "q9" "q5"];

%-------------
% PLOT OPTIONS
%-------------
plot_style = "violin plot";                % alt options: "mean & std", "box plot", "violin plot", and "eCDF"
ascending_mean = true;              % do you want to sort genotypes by ascending mean?
median_normalized = true;           % do you want to normalized the distributions based on the median?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?
x_axis_text_angle = 0;              % choose angle of x-axis text (rotate so they don't overlap)
plot_title = "Inter-R8-distances";
title_size = 20;
x_label = "Genotypes";
y_label = "Inter-R8-distances (microns)";
axes_label_size = 16;
x_axis_lim = [-30 30];
y_axis_limit = [-10 10];
conversion_factor = 0.454;      % pixels to um
save_csv_to_file = true;       % would you like to save a cvs document to file?

interR8distancePerGeno(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes,plot_style,ascending_mean,median_normalized,genotype_labels,...
    x_axis_text_angle, plot_title, title_size, ...
    x_label, y_label, axes_label_size, x_axis_lim,...
    conversion_factor,save_csv_to_file,y_axis_limit)


%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% use U-Mann-Whitley test to compare COV
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

target_genotypes = ["q5" "q11"];

MannWhitney_COV(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes)



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% use U-Mann-Whitley test to compare inter-R8-distances
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

target_genotypes = ["q9" "q11"];

MannWhitney_interR8(expInfo,genotype_code,clean_omma_centroids,delaunay_neighbors,...
    target_genotypes)




















%%

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

