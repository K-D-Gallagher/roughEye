function varargout = adultOmmatidiaSeg(varargin)
% ADULTOMMATIDIASEG MATLAB code for adultOmmatidiaSeg.fig
%      ADULTOMMATIDIASEG, by itself, creates a new ADULTOMMATIDIASEG or raises the existing
%      singleton*.
%
%      H = ADULTOMMATIDIASEG returns the handle to a new ADULTOMMATIDIASEG or the handle to
%      the existing singleton*.
%
%      ADULTOMMATIDIASEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADULTOMMATIDIASEG.M with the given input arguments.
%
%      ADULTOMMATIDIASEG('Property','Value',...) creates a new ADULTOMMATIDIASEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before adultOmmatidiaSeg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to adultOmmatidiaSeg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help adultOmmatidiaSeg

% Last Modified by GUIDE v2.5 15-Jan-2021 13:08:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @adultOmmatidiaSeg_OpeningFcn, ...
                   'gui_OutputFcn',  @adultOmmatidiaSeg_OutputFcn, ...
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

% --- Executes just before adultOmmatidiaSeg is made visible.
function adultOmmatidiaSeg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to adultOmmatidiaSeg (see VARARGIN)

% Choose default command line output for adultOmmatidiaSeg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set intitial image to always be 1
handles.currentImage = 1;

% read in raw images
handles.raw_images = evalin('base','raw_images');

% read in ommatidia centroids
handles.omma_cent = evalin('base','omma_cent');

% read in ommatidia centroids image
handles.omma_cent_disp = zeros(size(handles.omma_cent,3));
for j = 1:size(handles.omma_cent,3)
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


% Update handles structure
update_arrows(handles)
guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)

% UIWAIT makes adultOmmatidiaSeg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = adultOmmatidiaSeg_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in add_omma.
function add_omma_Callback(hObject, eventdata, handles)
% hObject    handle to add_omma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  store time
time = str2num(handles.current_image.String);

%imageHandle = imshow(handles.L(:,:,handles.currentImage),'Parent',handles.axes1);
[y, x, ~] = impixel;

% make sure we're only adding a single ommatidia
if length(x) == 1
    
    x = round(x);
    y  = round(y);
    
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

assignin('base','omma_cent',handles.omma_cent);
assignin('base','omma_cent_disp',handles.omma_cent_disp);
    
guidata(hObject, handles)

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
        plot(curr_omma(ii,1),curr_omma(ii,2),'+b','Linewidth',2)
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
        plot(curr_omma(ii,1),curr_omma(ii,2),'+b','Linewidth',2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shit you can't delete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% --- Executes on button press in Debug.
function Debug_Callback(hObject, eventdata, handles)
% hObject    handle to Debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard
