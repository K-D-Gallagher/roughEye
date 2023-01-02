
% segment ommatidia from brightfield images of Drosophila adult eyes,
% define ommatidial lattice topology, and quantify roughness phenotype
% K-D-Gallagher https://github.com/K-D-Gallagher 2021



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
% Section1: Crop
% CROP DATA BEFORE USING ILASTIK
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% filepath = uigetdir('','Select directory containing data');
% filepath = strcat(filepath,'/');
filepath = '/your filepath/';


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
% For following sections,need to run
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
filepath = '/your filepath/';

experimentID = 'new folder name for output data'; % NO UNDERSCORES (_) please


%--------------------------------------------------------------------------
% define genotypes
%--------------------------------------------------------------------------

% NOTE: the way you write the genotypes must exactly mach the way its
% written in the filename

%---------
% genotype_code = ["Q5" "Q9" "Q11" "Q12" "Q13" "Q14"];
% genotype_code = ["Q23" "Q24" "Q29" "Q33" "Q34" "Q39"];
genotype_code = ["Q46" "Q47" "Q50" "Q51" "Q52" "Q53"];

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
% Section2-1 for hand correction
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
% Section2-2 for hand correction
% hand correct using GUI
%
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Load the GUI for hand-correcting segmentation. The GUI will allow you to
% export the 'omma_centroids' variable

ommatidiaSeg




%
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%
% OPTIONAL?visualize segmented ommatidia
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
% Section3-1 for Trangulation
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
%OPTIONAL:visualize triangulation
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
% OPTIONAL: visualize neighbors of each ommatidia at the specified time point
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
% Section3-2 for export .CSV DATA
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
%target_genotypes = ["Q5" "Q11" "Q13" "Q12" "Q9" "Q14"];
target_genotypes = ["Q47" "Q51" "Q46" "Q52" "Q53" "Q50"];
%target_genotypes = ["Q23" "Q29" "Q39" "Q33" "Q34" "Q24"];

%--------------------------------------------------------------------------
% PLOT OPTIONS
%--------------------------------------------------------------------------
plot_style = "violin plot";         % alt options: "mean & std", "box plot", "violin plot"
ascending_mean = false;              % do you want to sort genotypes by ascending mean?
genotype_labels = genotype_code;    % x-axis genotype labels - do you want to use the short hand code or full genotype?

% plot_title = "Coefficient of Variation (COV) of inter-R8-distance";
% x_label = "Genotypes";
% y_label = "Coefficient of Variation (COV)";

plot_title = "Value of D";
x_label = "Genotypes";
y_label = "Value of D";

title_size = 20;
axes_label_size = 16;
x_axis_text_angle = 0;         % choose angle of x-axis text (rotate so they don't overlap)
y_axis_limit = [0 1.5];

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

