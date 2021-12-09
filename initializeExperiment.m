function [expInfo] = initializeExperiment(experimentID,genotype_code,full_genotype,filepath)

% create new folder to save all plots, images, and database files within
filepath_out = strcat(filepath,'/',experimentID,'/');
[ status, msg ] = mkdir(filepath_out);
if status == 0
    msg
end

% create Directory based only on .tif files
Directory = dir(strcat(filepath,'*.tif'));
num_meas = length(Directory);

% initialize some data storage structures
filename = cell(num_meas,1);
genotype = cell(num_meas,1);
sex = cell(num_meas,1);

% string parsing based on filename to record the genotype and sex of each
% image
for t = 1:num_meas
    
    temp_namestr = Directory(t).name;
    temp_namestr = temp_namestr(1:end-4);               % trim of '.tif'
    temp_namestr = strrep(temp_namestr,'_',' ');        % replace '_' with ' '
    filename{t} = temp_namestr;                         % record filename
    split_str = strsplit(temp_namestr);                 % parse out the genotype specifically
    genotype{t} = split_str{2};                         % record genotype
    sex{t} = split_str{1};                              % record sex
    
end

% record full genotype
genotype_full = cell(num_meas,1);
for j = 1:length(genotype)
    [LIA,LOCB] = ismember(genotype{j},genotype_code);
    if LIA
        genotype_full{j} = full_genotype(LOCB);
    end
end

% convert to string arrays
filename = string(filename);
genotype = string(genotype);
% genotype_full = string(genotype_full);
sex = string(sex);

% make struct containing all the info we want
field1 = 'experimentID';    value1 = {experimentID};
field2 = 'filepath_input';  value2 = {filepath};
field3 = 'filepath_output'; value3 = {filepath_out};
field4 = 'filenames';       value4 = {filename};
field5 = 'genotypes_code';  value5 = {genotype};
field6 = 'genotypes_full';  value6 = {genotype_full};
field7 = 'sex';             value7 = {sex};

% return struct containing all the info we want
expInfo = struct(field1,value1,field2,value2,field3,value3,...
    field4,value4,field5,value5,field6,value6,field7,value7);