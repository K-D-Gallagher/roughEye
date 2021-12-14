function varargout = cropTool(varargin)
% CROPTOOL MATLAB code for cropTool.fig
%      CROPTOOL, by itself, creates a new CROPTOOL or raises the existing
%      singleton*.
%
%      H = CROPTOOL returns the handle to a new CROPTOOL or the handle to
%      the existing singleton*.
%
%      CROPTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROPTOOL.M with the given input arguments.
%
%      CROPTOOL('Property','Value',...) creates a new CROPTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cropTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cropTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cropTool

% Last Modified by GUIDE v2.5 09-Dec-2021 12:33:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cropTool_OpeningFcn, ...
                   'gui_OutputFcn',  @cropTool_OutputFcn, ...
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

% --- Executes just before cropTool is made visible.
function cropTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cropTool (see VARARGIN)

% Choose default command line output for cropTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set intitial image to always be 1
handles.currentImage = 1;

% read in uncropped images
handles.raw_images = evalin('base','pre_crop_images');

% populate popupmenu crop size options
strings = {'512 x 512','640 x 640', '768 x 768','896 x 896', '1024 x 1024', '1152 x 1152', '1280 x 1280'};
set(handles.crop_sizes, 'String', strings);
handles.crop_options = [512 640 768 896 1024 1152 1280];

% get current crop size
crop_value = handles.crop_options(get(handles.crop_sizes, 'Value'));
handles.crop_value = crop_value;

% generate cropped images matrix based on crop size
handles.cropped_images = uint8(zeros(crop_value,crop_value,3,size(handles.raw_images,4)));


% find maximum number of time points and store in handles.max
handles.max = size(handles.raw_images,4);

% feed the max number of images to the 'num_images' function
handles.numImages = handles.max;

% feed the current image number to the 'current_image' callback
set(handles.current_image,'String', handles.currentImage);

% feed the max number of images to the 'num_images' callback
set(handles.num_images,'String',handles.max);

% read in metainfo
handles.meta_info = evalin('base','cropInfo');

% create handle to record FOV
handles.crop_crop_coord = zeros(4,handles.max);

% Update handles structure
update_arrows(handles)
guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = cropTool_OutputFcn(hObject, eventdata, handles) 
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
% crop sizes dropdown menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in crop_sizes.
function crop_sizes_Callback(hObject, eventdata, handles)
% hObject    handle to crop_sizes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crop_sizes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crop_sizes

% get current crop size
crop_value = handles.crop_options(get(handles.crop_sizes, 'Value'));
handles.crop_value = crop_value;

% generate cropped images matrix based on crop size
handles.cropped_images = uint8(zeros(crop_value,crop_value,3,size(handles.raw_images,4)));

% Update handles structure
update_arrows(handles)
guidata(hObject, handles)
stack_display_update(hObject, eventdata, handles)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)
% hObject    handle to crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  store time
time = str2num(handles.current_image.String);
crop_value = handles.crop_value;
curr_im = handles.raw_images(:,:,:,time);

% %imageHandle = imshow(handles.L(:,:,handles.currentImage),'Parent',handles.axes1);
% [y, x, ~] = impixel;

% choose and record center pixel of cropped region
[xTemp,yTemp] = getpts(handles.axes1);
xTemp=round(xTemp); yTemp=round(yTemp);

plus = crop_value/2;
minus = crop_value/2 - 1;

% create cropped image,make grayscale, histEq
crop_im = curr_im(yTemp-minus:yTemp+plus,xTemp-minus:xTemp+plus,:);

% record crop coordinates
handles.crop_coord(:,time) = [yTemp-minus,yTemp+plus,xTemp-minus,xTemp+plus];

handles.cropped_images(:,:,:,time) = crop_im;

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

% --- Executes on button press in to_workspace.
function to_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to to_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','raw_images',handles.cropped_images);
assignin('base','crop_coord',handles.crop_coord);
    
guidata(hObject, handles)


% --- Executes on button press in save_to_file.
function save_to_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_to_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num_images = handles.max;

for j = 1:num_images
    fullname = fullfile(handles.meta_info.filepath_out,handles.meta_info.filenames{j});
    imwrite(handles.cropped_images(:,:,:,j),fullname)
end

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
% display update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function stack_display_update(hObject, eventdata, handles)
% hObject    handle to watershed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cur = str2num(handles.current_image.String);  
curr_raw = handles.raw_images(:,:,:,cur);
curr_crop = handles.cropped_images(:,:,:,cur);
imshow(curr_raw,'Parent',handles.axes1);
imshow(curr_crop,'Parent',handles.axes2);

update_arrows(handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other
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


% --- Executes during object creation, after setting all properties.
function crop_sizes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crop_sizes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

strings = {'512 x 512','640 x 640', '768 x 768','896 x 896', '1024 x 1024', '1152 x 1152', '1280 x 1280'};
handles.crop_sizes.Strings = strings;

guidata(hObject, handles);
