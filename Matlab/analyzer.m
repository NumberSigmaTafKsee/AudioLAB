function varargout = analyzer(varargin)
% ANALYZER MATLAB code for analyzer.fig
%      ANALYZER, by itself, creates a new ANALYZER or raises the existing
%      singleton*.
%
%      H = ANALYZER returns the handle to a new ANALYZER or the handle to
%      the existing singleton*.
%
%      ANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZER.M_inserted with the given input arguments.
%
%      ANALYZER('Property','Value',...) creates a new ANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyzer

% Last Modified by GUIDE v2.5 16-Aug-2014 18:59:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzer_OutputFcn, ...
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


% --- Executes just before analyzer is made visible.
function analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyzer (see VARARGIN)

% Choose default command line output for analyzer
handles.output = hObject;
set(handles.A_slider,'Value',0.5);
set(handles.F_slider,'Value',5);
handles.A = get(handles.A_slider,'Value');
handles.f0 = get(handles.F_slider,'Value');
handles.N = str2double(get(handles.N_inserted, 'String'));
handles.M = str2double(get(handles.M_inserted, 'String'));
handles.W = rectwin(handles.N);
handles.signal = get(handles.Signal,'Value');
handles.show_W = 0;
handles.current_window = 'rectangular';
% Update handles structure
guidata(hObject, handles);
new_signal(handles);



set(handles.text7,'String',strcat(num2str(0.5),'V'));
set(handles.text9,'String',strcat(num2str(5),'Hz'));

% UIWAIT makes analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analyzer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Signal.
function Signal_Callback(hObject, eventdata, handles)
% hObject    handle to Signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.signal = get(handles.Signal,'Value');
% Update handles structure
guidata(hObject, handles);
new_signal(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns Signal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Signal


% --- Executes during object creation, after setting all properties.
function Signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function A_slider_Callback(hObject, eventdata, handles)
% hObject    handle to A_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text7,'String',strcat(num2str(single(get(hObject,'Value'))),'V'));
handles.A = get(handles.A_slider,'Value');
% Update handles structure
guidata(hObject, handles);
new_signal(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function A_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function F_slider_Callback(hObject, eventdata, handles)
% hObject    handle to F_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text9,'String',strcat(num2str(int8(get(hObject,'Value'))),'Hz'));
handles.f0 = get(handles.F_slider,'Value');
% Update handles structure
guidata(hObject, handles);
new_signal(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function F_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function N_inserted_Callback(hObject, eventdata, handles)
% hObject    handle to N_inserted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.N = str2double(get(handles.N_inserted, 'String'));
[handles.W,handles.current_window] = new_window(handles.current_window,handles);
% Update handles structure
guidata(hObject, handles);

new_signal(handles);
% Hints: get(hObject,'String') returns contents of N_inserted as text
%        str2double(get(hObject,'String')) returns contents of N_inserted as a double


% --- Executes during object creation, after setting all properties.
function N_inserted_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_inserted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_inserted_Callback(hObject, eventdata, handles)
% hObject    handle to M_inserted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.M = str2double(get(handles.M_inserted, 'String'));
% Update handles structure
guidata(hObject, handles);
new_signal(handles);
% Hints: get(hObject,'String') returns contents of M_inserted as text
%        str2double(get(hObject,'String')) returns contents of M_inserted as a double


% --- Executes during object creation, after setting all properties.
function M_inserted_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_inserted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_window.
function show_window_Callback(hObject, eventdata, handles)
% hObject    handle to show_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.show_W = mod(handles.show_W+1,2);
if handles.show_W == 1
   set(handles.show_window,'String','Show waveform');
else
    set(handles.show_window,'String','Show window');
end
guidata(hObject, handles);
new_signal(handles);
% Hint: get(hObject,'Value') returns toggle state of show_window

function new_signal(handles)

t0 = 0; %first time instant
tf = 1; %last time instant
%A = 0.5; % signal amplitude
%f0 = 5; % signal frequency
%N_inserted = 201; % # signal samples
%M_inserted = 256; %FFT points
%Signal generation
t = linspace(t0,tf,handles.N);
%Sampling rate and period
Ts = tf/handles.N;
fs = 1/Ts;
f = linspace(0,fs/2,handles.M/2);
if handles.show_W == 0
switch handles.signal
    case 1
        v = ((handles.A)*cos(2*pi*handles.f0*t)).*(handles.W)';
    case 2
        v = ((handles.A)*square(2*pi*handles.f0*t)).*(handles.W)';
    case 3
        v = ((handles.A)*sawtooth(2*pi*handles.f0*t)).*(handles.W)';
end

%Window over signal
v_windowed = v.*(handles.W)';
%FFT with dB representation
V = fft(v_windowed,handles.M); %We know here we have a M-N length 0 padding
V_dB = 20*log10(abs(V(1:length(V)/2)));
else
    v = handles.W;
    V = fft(handles.W,handles.M);
    V_dB = 20*log10(abs(V(1:length(V)/2)));

end
plot(handles.time, t, v','r');

plot(handles.freq, f, V_dB,'r');
axis([0 fs/2 -60 40]);
set(get(handles.freq,'XLabel'),'String','Frequency (Hz)');
set(get(handles.freq,'YLabel'),'String','Magnitude (dB)');
set(get(handles.time,'XLabel'),'String','Time (s)');
set(get(handles.time,'YLabel'),'String','Amplitude (V)');
% --------------------------------------------------------------------
function choose_window_Callback(hObject, eventdata, handles)
% hObject    handle to choose_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rectangular_Callback(hObject, eventdata, handles)
% hObject    handle to rectangular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Rectangular');
[handles.W,handles.current_window] = new_window('rectangular',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function triangular_Callback(hObject, eventdata, handles)
% hObject    handle to triangular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Triangular');
[handles.W,handles.current_window] = new_window('triangular',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);



% --------------------------------------------------------------------
function hamming_Callback(hObject, eventdata, handles)
% hObject    handle to hamming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Hamming');
[handles.W,handles.current_window] = new_window('hamming',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);

% --------------------------------------------------------------------
function hanning_Callback(hObject, eventdata, handles)
% hObject    handle to hanning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Hanning');
[handles.W,handles.current_window] = new_window('hanning',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);

% --------------------------------------------------------------------
function chebychev_Callback(hObject, eventdata, handles)
% hObject    handle to chebychev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function kaiser_Callback(hObject, eventdata, handles)
% hObject    handle to kaiser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function parameter_1_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Kaiser 1');
[handles.W,handles.current_window] = new_window('kaiser1',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function parameter_2_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Kaiser 4');
[handles.W,handles.current_window] = new_window('kaiser4',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function parameter_3_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Kaiser 9');
[handles.W,handles.current_window] = new_window('kaiser9',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Chebychev 20dB');
[handles.W,handles.current_window] = new_window('chebychev20',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Chebychev 30dB');
[handles.W,handles.current_window] = new_window('chebychev30',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.window_text,'String','Chebychev 40dB');
[handles.W,handles.current_window] = new_window('chebychev40',handles);
% Update handles structure
guidata(hObject, handles);
new_signal(handles);

function [W,current_window] = new_window(name,handles)
    current_window = name;
switch name

    case 'rectangular'
        W = rectwin(handles.N);
    case 'triangular'
        W = triang(handles.N);
    case 'hamming'
        W = hamming(handles.N);
    case 'hanning'
        W = hanning(handles.N);
    case 'chebychev20'
        W = chebwin(handles.N,20);
    case 'chebychev30'
        W = chebwin(handles.N,30);
    case 'chebychev40'
        W = chebwin(handles.N,40);
    case 'kaiser1'
        W = kaiser(handles.N,1);
    case 'kaiser4'
        W = kaiser(handles.N,4);
    case 'kaiser9'
        W = kaiser(handles.N,9);
end

