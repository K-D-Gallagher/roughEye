function varargout = ommatidiaSeg(varargin)
% OMMATIDIASEG MATLAB code for ommatidiaSeg.fig
%      OMMATIDIASEG, by itself, creates a new OMMATIDIASEG or raises the existing
%      singleton*.
%
%      H = OMMATIDIASEG returns the handle to a new OMMATIDIASEG or the handle to
%      the existing singleton*.
%
%      OMMATIDIASEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OMMATIDIASEG.M with the given input arguments.
%
%      OMMATIDIASEG('Property','Value',...) creates a new OMMATIDIASEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ommatidiaSeg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ommatidiaSeg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ommatidiaSeg

% Last Modified by GUIDE v2.5 14-Dec-2021 13:36:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ommatidiaSeg_OpeningFcn, ...
                   'gui_OutputFcn',  @ommatidiaSeg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before ommatidiaSeg is made visible.
function ommatidiaSeg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ommatidiaSeg (see VARARGIN)

% Choose default command line output for ommatidiaSeg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set intitial image to always be 1
handles.currentImage = 1;

% read in raw images
handles.raw_images = evalin('base','raw_images');

% read in ommatidia centroids
handles.omma_cent = evalin('base','omma_centroids');

% read in ommatidia centroids image
handles.omma_cent_disp = zeros(size(handles.raw_images,1),size(handles.raw_images,2),size(handles.raw_images,4));
for j = 1:length(handles.omma_cent)
    for jj = 1:size(handles.omma_cent{j},1)
        handles.omma_cent_disp(handles.omma_cent{j}(jj,2),handles.omma_cent{j}(jj,1),j) = 1;
    end
end

% read in ilastik pixel classification
handles.ilastik_probabilities = evalin('base','ilastik_probabilities');
handles.pix_classified = zeros(size(handles.ilastik_probabilities));

% loop through probability masks and threshold each one
for i = 1:size(handles.ilastik_probabilities,3)
    
    % thresholding ilastik probability masks
    handles.pix_classified(:,:,i) = im2bw(handles.ilastik_probabilities(:,:,i),0.3);
    
end

% find maximum number of time points and store in handles.max
handles.max = size(handles.raw_images,4);

% feed the max number of images to the 'num_images' function
handles.numImages = handles.max;

% feed the current image number to the 'current_image' callback
set(handles.current_image,'String', handles.currentImage);

% feed the max number of images to the 'num_images' callback
set(handles.num_images,'String',handles.max);

% set initialization values of toggle ilastik
set(handles.toggle_ilastik,'Value',0)

% set initialization values of toggle omma
set(handles.toggle_omma,'Value',1)

% set baseline x/y axes limits
handles.baseline_x_limit = [0.5,size(handles.raw_images,2)+0.5]; 
handles.baseline_y_limit = [0.5,size(handles.raw_images,1)+0.5]; 

% set current x/y axes limits to the baseline limits
handles.x_lim = handles.baseline_x_limit;
handles.y_lim = handles.baseline_y_limit;

% populate popupmenu marker_type options
strings = {'*', '+', '.', 'o'};
set(handles.marker_type, 'String', strings);
handles.marker_type_options = ['*', '+', '.', 'o'];
handles.marker_type_choice = handles.marker_type_options(1);

% populate popupmenu marker_color options
strings = {'cyan', 'blue', 'green', 'magenta', 'yellow', 'black', 'white', 'red'};
set(handles.marker_color, 'String', strings);
handles.marker_color_options = ['cyan', 'blue', 'green', 'magenta', 'yellow', 'black', 'white', 'red'];
handles.marker_color_choice = handles.marker_color_options(1);

% populate popupmenu marker_size options
strings = {'2', '4', '6', '8', '10', '12', '14'};
set(handles.marker_size, 'String', strings);
handles.marker_size_options = [2 4 6 8 10 12 14];
handles.marker_size_choice = handles.marker_size_options(3);

% populate popupmenu line_width options
strings = {'1', '2', '3', '4'};
set(handles.line_width, 'String', strings);
handles.line_width_options = [1 2 3 4];
handles.line_width_choice = handles.line_width_options(1);

% Update handles structure
update_arrows(handles)
guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)

% UIWAIT makes ommatidiaSeg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ommatidiaSeg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add ommatidia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.omma_cent = handles.backup_omma_cent;
handles.omma_cent_disp = handles.backup_omma_cent_disp;

guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add ommatidia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in add_omma.
function add_omma_Callback(hObject, eventdata, handles)
% hObject    handle to add_omma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% backup in case of undo
handles.backup_omma_cent = handles.omma_cent;
handles.backup_omma_cent_disp = handles.omma_cent_disp;

%  store time
time = str2num(handles.current_image.String);

%imageHandle = imshow(handles.L(:,:,handles.currentImage),'Parent',handles.axes1);
[temp_y, temp_x, ~] = impixel;

% make sure we're only adding a single ommatidia
for j = 1:length(temp_x)
    
    x = round(temp_x(j));
    y  = round(temp_y(j));
    
    % pull out current binary image of centroid positions
    curr_omma_disp = handles.omma_cent_disp(:,:,time);
    
    % flip this pixel white in 'omma_cent_display'
    curr_omma_disp(x,y) = 1;
    
    % compute centroid positions
    [x,y] = find(curr_omma_disp);
    
    temp_centroid = [y, x];
    
    handles.omma_cent_disp(:,:,time) = curr_omma_disp;
    handles.omma_cent{time} = temp_centroid;
    
end

guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete ommatidia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in delete_omma.
function delete_omma_Callback(hObject, eventdata, handles)
% hObject    handle to delete_omma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% backup in case of undo
handles.backup_omma_cent = handles.omma_cent;
handles.backup_omma_cent_disp = handles.omma_cent_disp;

%  store time
time = str2num(handles.current_image.String);

% draw region within which to delete ommatidia
h = imrect(handles.axes1);
mask = uint8(createMask(h));

% pull out current binary image of centroid positions
curr_omma_disp = handles.omma_cent_disp(:,:,time);

[x_mask,y_mask] = find(mask);

for j = 1:length(x_mask)
    curr_omma_disp(x_mask(j),y_mask(j)) = 0;
end

% compute centroid positions
[x,y] = find(curr_omma_disp);

temp_centroid = [y, x];

handles.omma_cent_disp(:,:,time) = curr_omma_disp;
handles.omma_cent{time} = temp_centroid;

guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define custom ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in define_ROI.
function define_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to define_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% backup in case of undo
handles.backup_omma_cent = handles.omma_cent;
handles.backup_omma_cent_disp = handles.omma_cent_disp;

%  store time
t = str2num(handles.current_image.String);

% draw ROI and save coordinates
p = drawpolygon('LineWidth',7,'Color','cyan');
ROI_points = p.Position;

% create binary mask of ROI
bw = double(poly2mask(ROI_points(:,1),ROI_points(:,2),size(handles.omma_cent_disp,1),size(handles.omma_cent_disp,2)));

% apply mask to centroid map
ROI = handles.omma_cent_disp(:,:,t) .* bw;

% find new centroid positions
% compute centroid positions
[x,y] = find(ROI);
temp_centroid = [y, x];
handles.omma_cent{t} = temp_centroid;

% update omma display
new_omma_disp = zeros(size(handles.omma_cent_disp(:,:,t)));
for j = 1:length(x)
    new_omma_disp(x(j),y(j)) = 0;
end
handles.omma_cent_disp(:,:,t) = new_omma_disp;


guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in export_labels.
function export_labels_Callback(hObject, eventdata, handles)
% hObject    handle to export_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','omma_centroids',handles.omma_cent);
    
guidata(hObject, handles)


% --- Executes on button press in save_to_file.
function save_to_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_to_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

omma_centroids = handles.omma_cent;

save('backup_omma_centroids.mat','omma_centroids');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Foward and back button display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Handles the display of the forward and back buttons.
function update_arrows(handles)
set(handles.current_image,'String', handles.currentImage);
if (handles.currentImage < handles.numImages) 
    set(handles.stack_forward_in_time, 'Enable', 'On');
    set(handles.stack_fast_forward_in_time,'Enable','On');
else
    set(handles.stack_forward_in_time, 'Enable', 'Off');
    set(handles.stack_fast_forward_in_time,'Enable','Off');
end
if (handles.currentImage > 1)
    set(handles.stack_backward_in_time, 'Enable', 'On');
    set(handles.stack_fast_backward_in_time,'Enable','On');
else
    set(handles.stack_backward_in_time, 'Enable', 'Off');
    set(handles.stack_fast_backward_in_time,'Enable','Off');
end

% --- Executes on button press in stack_backward_in_time.
function stack_backward_in_time_Callback(hObject, eventdata, handles)
% hObject    handle to stack_backward_in_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currentImage > 1)
    handles.currentImage = handles.currentImage - 1;
    update_arrows(handles); 
    stack_display_update(hObject, eventdata, handles)
end

% --- Executes on button press in stack_forward_in_time.
function stack_forward_in_time_Callback(hObject, eventdata, handles)
% hObject    handle to stack_forward_in_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currentImage < handles.numImages)
    handles.currentImage = handles.currentImage + 1;
    update_arrows(handles); 
    stack_display_update(hObject, eventdata, handles)
end

% --- Executes on button press in stack_fast_backward_in_time.
function stack_fast_backward_in_time_Callback(hObject, eventdata, handles)
% hObject    handle to stack_fast_backward_in_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currentImage > 10)
    handles.currentImage = handles.currentImage - 10;
    update_arrows(handles); 
    stack_display_update(hObject, eventdata, handles)
end

% --- Executes on button press in stack_fast_forward_in_time.
function stack_fast_forward_in_time_Callback(hObject, eventdata, handles)
% hObject    handle to stack_fast_forward_in_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currentImage < handles.numImages-10)
    handles.currentImage = handles.currentImage + 10;
    update_arrows(handles); 
    stack_display_update(hObject, eventdata, handles)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% draw rectangle and record mask of it

h = imrect(handles.axes1);
mask = uint8(createMask(h));

% get bounding box of the rectangular ROI drawn with the code above
zoom_dimensions = regionprops(mask,'BoundingBox');

% transfer from Struct to vector, Bounding Box is [x-start y-end length width]
% of the rectangle
zoom_dimensions = [round(zoom_dimensions.BoundingBox(1)) round(zoom_dimensions.BoundingBox(2)) ...
    round(zoom_dimensions.BoundingBox(3))+round(zoom_dimensions.BoundingBox(1)) ...
    round(zoom_dimensions.BoundingBox(4))+round(zoom_dimensions.BoundingBox(2))];

handles.x_lim = [zoom_dimensions(1) zoom_dimensions(3)];
handles.y_lim = [zoom_dimensions(2) zoom_dimensions(4)];


guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)


% --- Executes on button press in reset_zoom.
function reset_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to reset_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.x_lim = handles.baseline_x_limit;
handles.y_lim = handles.baseline_y_limit;

guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toggle layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in toggle_ilastik.
function toggle_ilastik_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_ilastik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_ilastik

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_ilastik.Value = current_state;

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in toggle_omma.
function toggle_omma_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_omma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_omma.Value = current_state;

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ommatidia marker choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in marker_size.
function marker_size_Callback(hObject, eventdata, handles)
% hObject    handle to marker_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns marker_size contents as cell array
%        contents{get(hObject,'Value')} returns selected item from marker_size

handles.marker_size_choice = handles.marker_size_options(get(handles.marker_size, 'Value'));

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes on selection change in marker_type.
function marker_type_Callback(hObject, eventdata, handles)
% hObject    handle to marker_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns marker_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from marker_type

handles.marker_type_choice = handles.marker_type_options(get(handles.marker_type, 'Value'));

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on selection change in marker_color.
function marker_color_Callback(hObject, eventdata, handles)
% hObject    handle to marker_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns marker_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from marker_color

handles.marker_color_choice = handles.marker_color_options(get(handles.marker_color, 'Value'));

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes on selection change in line_width.
function line_width_Callback(hObject, eventdata, handles)
% hObject    handle to line_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns line_width contents as cell array
%        contents{get(hObject,'Value')} returns selected item from line_width

handles.line_width_choice = handles.line_width_options(get(handles.line_width, 'Value'));

stack_display_update(hObject, eventdata, handles)
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stack_display_update(hObject, eventdata, handles)
% hObject    handle to watershed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cur = str2num(handles.current_image.String);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pixel classification and ommatidia seg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if handles.toggle_ilastik.Value ...
        && handles.toggle_omma.Value
    
    curr_raw = handles.raw_images(:,:,:,cur);
    curr_ilastik = handles.pix_classified(:,:,cur);
    imshowpair(curr_raw,curr_ilastik,'Parent',handles.axes1);
    hold on
    
    curr_omma = handles.omma_cent{cur};
    for ii = 1:length(curr_omma)
        plot(curr_omma(ii,1),curr_omma(ii,2),handles.marker_type_choice, ...
            'Color',handles.marker_color_choice, 'MarkerSize', handles.marker_size_choice, ...
            'Linewidth',handles.line_width_choice)
    end
    
    xlim(handles.x_lim)
    ylim(handles.y_lim);
    hold off
    
    update_arrows(handles);
    guidata(hObject, handles);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only ommatidia seg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_ilastik.Value) ...
        && handles.toggle_omma.Value
   
    curr_raw = handles.raw_images(:,:,:,cur);
    imshow(curr_raw,'Parent',handles.axes1)
    hold on
    
    curr_omma = handles.omma_cent{cur};
    for ii = 1:length(curr_omma)
        plot(curr_omma(ii,1),curr_omma(ii,2),handles.marker_type_choice, ...
            'Color',handles.marker_color_choice, 'MarkerSize', handles.marker_size_choice, ...
            'Linewidth',handles.line_width_choice)
    end
    
    xlim(handles.x_lim)
    ylim(handles.y_lim)
    hold off
    
    update_arrows(handles);    
    guidata(hObject, handles);

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only pixel classificatioin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_ilastik.Value ...
        && not(handles.toggle_omma.Value)
    
    curr_raw = handles.raw_images(:,:,:,cur);
    curr_ilastik = handles.pix_classified(:,:,cur);
    imshowpair(curr_raw,curr_ilastik,'Parent',handles.axes1);
    
    xlim(handles.x_lim)
    ylim(handles.y_lim)
    update_arrows(handles);
    guidata(hObject, handles);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neither pixel classification or ommatidia seg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_ilastik.Value) ...
        && not(handles.toggle_omma.Value)
    
    curr_raw = handles.raw_images(:,:,:,cur);
    imshow(curr_raw,'Parent',handles.axes1)
    
    xlim(handles.x_lim)
    ylim(handles.y_lim)
    
    update_arrows(handles);
    guidata(hObject, handles);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keyboard shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

 % determine the key that was pressed 
 keyPressed = eventdata.Key;
 
 % ROI shortcut
 if strcmpi(keyPressed,'r')
     % set focus to the button
     uicontrol(handles.define_ROI);
     % call the callback
     define_ROI(handles.define_ROI,[],handles);
     
 % add shortcut
 elseif strcmpi(keyPressed,'a')
     % set focus to the button
     uicontrol(handles.add_omma);
     % call the callback
     add_omma_Callback(handles.add_omma,[],handles);
     
 % delete shortcut
 elseif strcmpi(keyPressed,'d')
     % set focus to the button
     uicontrol(handles.delete_omma);
     % call the callback
     delete_omma_Callback(handles.delete_omma,[],handles);
     
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stuff you can't delete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Debug.
function Debug_Callback(hObject, eventdata, handles)
% hObject    handle to Debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard

% --- Executes during object creation, after setting all properties.
function num_images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function marker_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to marker_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function marker_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to marker_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function marker_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to marker_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function line_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to line_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
