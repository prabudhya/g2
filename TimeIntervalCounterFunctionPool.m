function TimeIntervalCounterFunctionPool(varargin)


switch varargin{1}
    case 'Init'                            
        Init(varargin{2}); 
    case 'InitGUI'
        %Expecting TimeIntervalCounter to call
        %TimeIntervalCounterFunctionPool('InitiGUI', handles) varargin{2}
        %is handles here
        InitGUI(varargin{2});
    case 'Query'                            
        Query(varargin{2}); 
    case 'setMode'
        setMode(varargin{2});
    case 'g2Measure'
        %Expecting TimeIntervalCounter to call
        %TimeIntervalCounterFunctionPool('g2Measure',hObject,eventdata,handles)
        g2Measure(varargin{2},varargin{3},varargin{4});
    case 'setTriggerLevel'
        setTriggerLevel(varargin{2});
    case'g2Plot'
        g2Calculate(varargin{2});
    case 'g2Stop'
        g2Stop(varargin{2});
end
        
end

function InitGUI(handles)

%Give some default parameters for histogram binning
set(handles.minInterval, 'String','-50');
set(handles.maxInterval, 'String','50');
set(handles.binSize, 'String','1');
%Call Init to handle initialization of GPIB object
Init(handles);


end

%Function for initializing the GPIB object.
function Init(handles)

global SR;

%Try to find a GPIB object

SR.gpib = instrfind('Type', 'gpib', 'BoardIndex', 7, 'PrimaryAddress', 16, 'Tag', '');

% Create the gpib object if it does not exist
% otherwise use the object that was found.
if isempty(SR.gpib)
    SR.gpib = gpib('AGILENT', 7, 16);
else
    fclose(SR.gpib);
    SR.gpib = SR.gpib(1);
end


 get(handles.modeMenu, 'String')
 disp('At end of Init')

IDN();
setClock();
end

%Function that sends query *IDN? as a test of communication
function IDN()
global SR;

if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end

fprintf(SR.gpib,'*IDN?');
identity=fscanf(SR.gpib);
if isempty(identity)
    warning('SR620 was not properly initialized.');
else
    disp(identity)
end

end

%Function for setting how many samples the SR620 takes in its measurement
function samples = set_numSamples(handles)

global SR;
if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end

entry = get(handles.set_numSamples, 'String');
%NaN is stored in a double array so isnumeric won't work
invalid = isnan(str2double(entry));

if (invalid) 
    %If a number is not provided put 10000 as a default value.
    entry = '1E4';
    disp('10000 samples chosen as a default value')
end

samples = str2double(entry);

end

%Function for running measurement of time intervals
function g2Measure(hObject,eventdata,handles)

%Setting a function that will execute whenever this function g2Measure ends
cleanupObj = onCleanup(@cleanupFcn);
 
global SR g2Go;

g2Go = true;

if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end

%Ensure SR620 is in Time mode
fprintf(SR.gpib, 'MODE 0');
pause(1)
%Also set arming mode. Use 1 for tests with cables. Use 0 for actual g2.
fprintf(SR.gpib, 'ARMM 0');
%fprintf(SR.gpib, 'ARMM 1');

%Pauses for 1s are to ensure commands execute. During testing, I found that
%they may not be executed quickly enough without the pause. Possible that this
%issue is only due to doing testing with a virtual machine that slows
%things down and pausing is unnecessary on RT1.
pause(1);
samples = set_numSamples(handles);

%Want to set input A as the source of START signals for real g2 measurement.
fprintf(SR.gpib, 'SRCE 0');

%fprintf(SR.gpib, 'SRCE 1');
%fprintf(SR.gpib,'SRCE 2');
pause(1)

global Enhanced_Factor;
Enhanced_Factor = 0;
warning = 0;

%Acquire ImageNVC's handles for use in tracking part of measurement. This
%is the same line as in TrackNV in the ExperimentYaoFunctionPool. 
ImageNVC_Handles=guidata(findobj(allchild(0), 'Tag', 'maxVx'));


tic
disp('Started Measurements')

%Binary dump mode can only be run with BDMP j, j ranging from 1 to 65535 so
%obtaining a larger number of samples requires breaking up the measurement.

collected = 0;
if samples > 20000
    
    filename = '';
    
    %Can only alter InputBufferSize when device is closed.
    fclose(SR.gpib);
    SR.gpib.Timeout = 60;
    SR.gpib.InputBufferSize = 8*20000+1;
    fopen(SR.gpib);
    itNum = 1;
    
    %Begin measurement by tracking
    ImageFunctionPool('NewTrackFast',hObject,eventdata,ImageNVC_Handles);
    while g2Go && samples > 0
        
       %Data is collected 20000 points at a time to reduce possibility of
       %timeouts if count rates are low.
        if samples > 20000
           
            fprintf(SR.gpib, 'BDMP%d',20000);
            %BDMP is run with 20000 but we collect twice that with fread
            %because half of it will be EOIs.
            raw = fread(SR.gpib,2*20000,'int');
            
        else 
            fprintf(SR.gpib, 'BDMP%d', samples);
            %The SR620 will get stuck in binary dump mode without the
            %factor of 2 in fread.
            raw = fread(SR.gpib,2*samples,'int');
        end
        
        
        data = convertData(raw);
        
        message = ['Iteration ',num2str(itNum), ':',' Collected 20000 points'];
        disp(message)
        
        if strcmp(filename, '')
            %Create a new file for storing first batch of data
            filename = g2Store(data,samples);
        else
            %Add the data to the existing file if not
            g2Append(data,filename);
        end
     
        samples = samples-20000;
        collected = collected + 20000;
        itNum = itNum +1;
        
        
        if ~g2Go
            disp('g2 measurement halted.');
            break; 
        end
        
        %Want to track after every 100000 intervals collected
        if mod(collected,100000) == 0
            Enhanced_Factor = 0;
            while Enhanced_Factor < .94
                ImageFunctionPool('NewTrackFast',hObject,eventdata,ImageNVC_Handles)
                if Enhanced_Factor<0.94
                    disp('Warning: Tracking optimized at 94% of the initial condition')
                    warning = warning+1;
                    if warning >=6
                        error('Lost NV')
                    end
                end
            end
            
        end
        %Need to put SR620 back in time mode after any tracking.
        fprintf(SR.gpib,'MODE 0');
        pause(1)
    end

else
    
    fclose(SR.gpib);

    SR.gpib.Timeout = 60;

    SR.gpib.InputBufferSize = 8*samples+1;
    fopen(SR.gpib);
    
    fprintf(SR.gpib, 'BDMP%d',samples);

    %Binary dump mode is supposed to send 8 byte two's complement integers.
    %fread should be able to preserve the sign of the time intervals
    raw = fread(SR.gpib,2*samples,'int');
    
    fprintf(SR.gpib,'MODE 0');
    data = convertData(raw);
    
    g2Store(data,samples);
end
disp('Finished Measurements')
toc
end

%Function for converting data from binary dump to meaningful times
function data = convertData(raw)
    
%A necessary conversion factor given on page 34 of SR620 manual.
timeFactor = 1.05963812934E-14;


%Removing the EOIs.
filtered = raw(1:2:end);

data = filtered*timeFactor;

end

%Function for storing time intervals in a new text file.
function file = g2Store(data,samples)

%Using code from ExperimentYaoFunctionPool to ensure g2 data is saved in
%the same location as other data on RT1's computer.

now = clock;
date = [num2str(now(1)),'-',num2str(now(2)),'-',num2str(round(now(3)))];
fullPath=fullfile('C:\Data\',date,'\');
if ~exist(fullPath,'dir')
    mkdir(fullPath);
end


%Create a filename for text file where g2 data will be stored.
time = datestr(datetime);
filename = strcat('g2_',time,'_',num2str(samples),'_samples','.txt');
%Remove whitespace and colons that cause problems in the file name
filename = strrep(filename,' ','_');
filename = strrep(filename,':','_');

file=strcat(fullPath,filename);
fid = fopen(file,'wt+');
fprintf(fid,'%e\n',data);
fclose(fid);
end

%Function for consolidating data into a single file.
function g2Append(data,file)

fid = fopen(file,'at+');
fprintf(fid,'%e\n',data);
fclose(fid);
end

%Function executed after stop button is pressed
function g2Stop(handles)
global g2Go;

%Set g2Go to false to break out of measurement loop.
g2Go = false;
disp('Stop button pressed');
end

%Function for creating a histogram from stored data.
function g2Calculate(handles)

%Want to select file with stored g2 data for use in histogram and other
%calculation.
now = clock;
date = [num2str(now(1)),'-',num2str(now(2)),'-',num2str(round(now(3)))];
fullPath=fullfile('C:\Data\',date,'\');
if ~exist(fullPath,'dir')
    mkdir(fullPath);
end
[baseName,folder] = uigetfile([fullPath,'*.*']);
filename = fullfile(folder,baseName);

%Read out values. Currently come out as a column vector.
fid = fopen(filename,'r');
intervals = fscanf(fid,'%f');
fclose(fid);

g2Hist = g2Bin(intervals,handles);

%This is code for producing a scatter plot of bin values along with the
%histogram if desired.
%{
values = g2Hist.Values;
edges = g2Hist.BinEdges;
width = g2Hist.BinWidth;
convertToScatter(values,edges,width)
%}

end

%Function for carrying out the binning of the data.
function hist = g2Bin(data, handles)

%Need to do binning in the time interval to get counts of number
%of times we have two photons on the A and B inputs separated by a certain
%amount of time.

edges = getBinEdges(handles);
edges = edges*(1E-9);

f = figure;
a = axes(f);

%Histogram will still be created even if user does not supply values for
%creating the bin edges.
if isnan(edges)
    hist = histogram(a,data);
    xlabel('Time (s)');
    ylabel('Counts');
else
    
    hist = histogram(a,data,edges);
    xlabel('Time (s)');
    ylabel('Counts');
end

end

%Function for finding bin edges based on user input to the GUI.
function edges = getBinEdges(handles)

minInterval = str2double(get(handles.minInterval, 'String'));

maxInterval = str2double(get(handles.maxInterval, 'String'));
binSize = str2double(get(handles.binSize, 'String'));

edges = minInterval:binSize:maxInterval-binSize;
end

%Function for representing bin values in a scatter plot
function convertToScatter(values,edges,width)

edges = edges(1:length(edges)-1);
centers = edges + width;
f2 = figure;
a2 = axes(f2);
scatter(centers,values,'filled');
end

%Function for setting mode of the SR620
function setMode(handles)
global SR;

%The Value property for the popup menu is a number indicating which string
%in the cell for the the String property was selected.
choiceIdx = get(handles.modeMenu, 'Value');
choices = cellstr(get(handles.modeMenu, 'String'));
disp(choiceIdx)
choice = choices(choiceIdx);
disp(choice)

if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end
SR.Mode = choice;

try
    switch choiceIdx
        %Index of the string in the cell array is one more than the number
        %used in the command to set the mode.
        case 1
            fprintf(SR.gpib, 'MODE 0');
        case 2
            fprintf(SR.gpib, 'MODE 1');
        case 3
            fprintf(SR.gpib, 'MODE 2');
        case 4
            fprintf(SR.gpib, 'MODE 3');
        case 5
            fprintf(SR.gpib, 'MODE 4');
        case 6
            fprintf(SR.gpib, 'MODE 5');
        case 7
            fprintf(SR.gpib, 'MODE 6');
    end   
catch ME
	fclose(SR.gpib);
	rethrow(ME);
end

fprintf(SR.gpib,'MODE?');
result = fscanf(SR.gpib);
disp('Mode is:')
disp(result);



end

%Function for creating and sending a query or command to the SR620
function Query(handles)

global SR;
if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
end

command = get(handles.queryText, 'String');

%If an invalid query is given, currently have to wait a bit for the time
%out message to appear on command line. Can then send correctly written
%queries.
fprintf(SR.gpib, command);
if(~isempty(strfind(command,'?')))
    response = fscanf(SR.gpib)

else
    disp('Command sent')
end

end

%Sets the clock source for the SR620.
function setClock()

global SR;
if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end


fprintf(SR.gpib, 'CLCK 0');
%fprintf(SR.gpib, 'CLCK 1');
fprintf(SR.gpib, 'CLCK?');clck = str2double(fscanf(SR.gpib));

if clck == 0
    disp('Clock set to internal source')
end

end

%Function for setting trigger levels of the channels.
function setTriggerLevel(handles)

global SR;
if(strcmp(SR.gpib.status,'closed'))
    fopen(SR.gpib);
    disp('Had to open device.')
end

voltage = get(handles.voltageNum, 'String');


choiceIdx = get(handles.channelMenu, 'Value');
choices = cellstr(get(handles.channelMenu, 'String'));
disp(choiceIdx);
choice = choices(choiceIdx);
disp(choice)

channelNum = num2str(choiceIdx-1);
command = ['LEVL',' ', channelNum,',',voltage]

fprintf(SR.gpib,command);
query = ['LEVL?',' ', channelNum];
fprintf(SR.gpib,query);
disp('Set Trigger Level');
result = fscanf(SR.gpib)

end

%Function that executes if measurement ends.
function cleanupFcn()
global g2Go SR;

g2Go = false;
disp('Measurement ended. Taking SR620 out of binary dump mode')
%Sending this command should prevent the SR620 from becoming stuck in
%binary dump mode if the measurement is ended earlier with Ctrl+C.
fprintf(SR.gpib, 'MODE 0');

end
