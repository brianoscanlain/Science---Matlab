function varargout = threshold_gui(varargin)
% THRESHOLD_GUI MATLAB code for threshold_gui.fig
%      THRESHOLD_GUI, by itself, creates a new THRESHOLD_GUI or raises the existing
%      singleton*.
%
%      H = THRESHOLD_GUI returns the handle to a new THRESHOLD_GUI or the handle to
%      the existing singleton*.
%
%      THRESHOLD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHOLD_GUI.M with the given input arguments.
%
%      THRESHOLD_GUI('Property','Value',...) creates a new THRESHOLD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before threshold_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to threshold_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help threshold_gui

% Last Modified by GUIDE v2.5 20-May-2014 16:36:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @threshold_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @threshold_gui_OutputFcn, ...
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


% --- Executes just before threshold_gui is made visible.
function threshold_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to threshold_gui (see VARARGIN)

% Choose default command line output for threshold_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes threshold_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = threshold_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assign a new callback function to the click event of the edit1 control
% It won’t work if you put it in your opening function since when 
% the opening function is called, the java objects have not been created yet!
set(findjobj(handles.edit1), ...
'MouseClickedCallback',{@MyButtonDownFcn , handles.edit1});
% Get default command line output from handles structure


keyboard;

function MyButtonDownFcn(hObject, eventdata,h)
% Manually assign the focus to the threshold_gui uicontrol
 
uicontrol(h)

function threshold_gui_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%       str2double(get(hObject,'String')) %returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_gui_WindowScrollWheelFcn(hObject, eventdata, handles)
 
h = gco;
if ~(h == handles.threshold_gui); return; end % exit if edit1 does not have the focus
 
s = str2double(get(h , 'String'));
if isnan(s); return; end; % exit if edit1 is empty
 
% VerticalScrollCount is positive or negative depending on 
% the direction of the rotation. It can also change in value 
% depending on the rotation speed.
if eventdata.VerticalScrollCount > 0
increment = -1;
else
increment = +1;
end
 
% Change the value of the edit1 string
set(h , 'String', num2str(s + increment));