function [cropInfo] = initializeCropping(filepath,slash_direction)

%--------------------------------------------------------------------------
% determine path for output (make new 'crop/' folder in subdirectory of
% input directory
%--------------------------------------------------------------------------

filepath_parts = strsplit(filepath,slash_direction);

filepath_out = '';

for j = 1:length(filepath_parts)-2
    if j == 1
        filepath_out = strcat(filepath_parts{j},'/');
    else
        filepath_out = strcat(filepath_out,filepath_parts{j},'/');
    end
end

filepath_out = strcat(filepath_out,'cropped/');

%------------------------------------------------
% create new folder to save cropped images within
%------------------------------------------------

[ status, msg ] = mkdir(filepath_out);
if status == 0
    msg
end

% create Directory based only on .tif files
Directory = dir(strcat(filepath,'*.tif'));
num_meas = length(Directory);

%-----------------
% record filenames
%-----------------

filenames = {};
for t = 1:num_meas

    filenames{t} = Directory(t).name;
    
end

%--------------------------------------------------------------------------
% construct struct to store data
%--------------------------------------------------------------------------

% make struct containing all the info we want
field1 = 'filepath_in';     value1 = {filepath};
field2 = 'filepath_out';    value2 = {filepath_out};
field3 = 'filenames';       value3 = {filenames};

% return struct containing all the info we want
cropInfo = struct(field1,value1,field2,value2,field3,value3);