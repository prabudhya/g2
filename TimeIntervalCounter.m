% Editing notes at the end of the document

function varargout = TimeIntervalCounter(varargin)

%This code has been generated after designing the GUI interactively in
%GUIDE. Most comments have been automatically written by MATLAB.

% TIMEINTERVALCOUNTER MATLAB code for TimeIntervalCounter.fig
%      TIMEINTERVALCOUNTER, by itself, creates a new TIMEINTERVALCOUNTER or raises the existing
%      singleton*.
%
%      H = TIMEINTERVALCOUNTER returns the handle to a new TIMEINTERVALCOUNTER or the handle to
%      the existing singleton*.
%
%      TIMEINTERVALCOUNTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIMEINTERVALCOUNTER.M with the given input arguments.
%
%      TIMEINTERVALCOUNTER('Property','Value',...) creates a new TIMEINTERVALCOUNTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TimeIntervalCounter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TimeIntervalCounter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TimeIntervalCounter

% Last Modified by GUIDE v2.5 01-Aug-2017 18:23:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TimeIntervalCounter_OpeningFcn, ...
                   'gui_OutputFcn',  @TimeIntervalCounter_OutputFcn, ...
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
end

% --- Executes just before TimeIntervalCounter is made visible.
function TimeIntervalCounter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TimeIntervalCounter (see VARARGIN)


handles.output = hObject; % Choose default command line output for TimeIntervalCounter
TimeIntervalCounterFunctionPool('InitGUI',handles);
guidata(hObject, handles); % Update handles structure

% UIWAIT makes TimeIntervalCounter wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end

% --- Outputs from this function are returned to the command line.
function varargout = TimeIntervalCounter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%{
function varargout = TimeIntervalCounter_DeleteFcn(hObject, eventdata, handles)

TimeIntervalCounterFunctionPool('Cleanup', handles);
end
%}

% --- Executes on button press in measure_g2, which is for starting the measurement of g2.
function measure_g2_Callback(hObject, eventdata, handles)
% hObject    handle to measure_g2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TimeIntervalCounterFunctionPool('g2Measure',hObject,eventdata,handles);


end

% set_numSamples is where the total time taken for the measurement is entered.
function set_numSamples_Callback(hObject, eventdata, handles)
% hObject    handle to set_numSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_numSamples as text
%        str2double(get(hObject,'String')) returns contents of set_numSamples as a double

end

% --- Executes during object creation, after setting all properties.
function set_numSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_numSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in query, which is for sending
% queries.
function query_Callback(hObject, eventdata, handles)
% hObject    handle to query (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%global SR;
%fopen(SR.gpib);
TimeIntervalCounterFunctionPool('Query',handles);

end

% --- Executes on selection change in modeMenu.
function modeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to modeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TimeIntervalCounterFunctionPool('setMode', handles)
% Hints: contents = cellstr(get(hObject,'String')) returns modeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modeMenu
end

% --- Executes during object creation, after setting all properties.
function modeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function queryText_Callback(hObject, eventdata, handles)
% hObject    handle to queryText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of queryText as text
%        str2double(get(hObject,'String')) returns contents of queryText as a double
end

% --- Executes during object creation, after setting all properties.
function queryText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to queryText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in levelButton.
function levelButton_Callback(hObject, eventdata, handles)
% hObject    handle to levelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TimeIntervalCounterFunctionPool('setTriggerLevel',handles);
end

function voltageNum_Callback(hObject, eventdata, handles)
% hObject    handle to voltageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voltageNum as text
%        str2double(get(hObject,'String')) returns contents of voltageNum as a double
end

% --- Executes during object creation, after setting all properties.
function voltageNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voltageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in channelMenu.
function channelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to channelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channelMenu
end

% --- Executes during object creation, after setting all properties.
function channelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function binSize_Callback(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binSize as text
%        str2double(get(hObject,'String')) returns contents of binSize as a double
end

% --- Executes during object creation, after setting all properties.
function binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function minInterval_Callback(hObject, eventdata, handles)
% hObject    handle to minInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minInterval as text
%        str2double(get(hObject,'String')) returns contents of minInterval as a double
end

% --- Executes during object creation, after setting all properties.
function minInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function maxInterval_Callback(hObject, eventdata, handles)
% hObject    handle to maxInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxInterval as text
%        str2double(get(hObject,'String')) returns contents of maxInterval as a double
end

% --- Executes during object creation, after setting all properties.
function maxInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in plot_g2.
function plot_g2_Callback(hObject, eventdata, handles)
% hObject    handle to plot_g2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TimeIntervalCounterFunctionPool('g2Plot', handles);

end

% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, eventdata, handles)
% hObject    handle to stopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TimeIntervalCounterFunctionPool('g2Stop',handles);

end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global SR;
fclose(SR.gpib);
delete(SR.gpib);
clear SR.gpib;
delete(hObject);

end

%% EDIT NOTES
% 1. Check that the created object for SR620 is deleted in
% ImageFunctionPool
