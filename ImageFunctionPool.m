function varargout = ImageFunctionPool(what,hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Written by Jeronimo Maze, July 2007 %%%%%%%%%%%%%%%%%%
%%%%%%%%%% Harvard University, Cambridge, USA  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch what 
    case 'NewTrackFast' 
        NewTrackFast(hObject, eventdata, handles);
    case {'ZScan','LargeScan','NormalScan'}
        PredefinedScans(what,hObject,eventdata,handles);
    case 'ManualChange' 
        ManChange(hObject, eventdata, handles);
    case 'UpdateCoords' 
        UpdateCoords(hObject, eventdata, handles);
    case 'GetCoords' 
        GetCoords(hObject, eventdata, handles);        
    case 'ShowMarker' 
        ShowMarker(hObject,eventdata,handles);
    case 'ApplyOffSet'
        ApplyOffSet(hObject, eventdata, handles);
    case 'Track'
        Track(hObject,eventdata,handles);
    case 'Start'
        ImageStart(hObject, eventdata, handles);
    case 'Scan'
        ImageScan('Scan',hObject, eventdata, handles);
    case 'Align'
        ImageScan('Align',hObject, eventdata, handles);
    case 'StopAlign'
        ImageScan('StopAlign',hObject, eventdata, handles);
    case 'StopScan'
        ImageScan('StopScan',hObject, eventdata, handles);
    case 'SetRange'
        ImageScan('SetRange',hObject, eventdata, handles);
    case 'ScanZ'
        ImageScan('ScanZ',hObject, eventdata, handles);
    case 'Save'
        ImageSaveImage('Save',hObject, eventdata, handles);
    case 'ChangeFileName'
        ImageSaveImage('ChangeFileName',hObject, eventdata, handles);
    case 'Edit'
        ImageSaveImage('Edit',hObject, eventdata, handles);
    case 'RePlot'
        ImageSaveImage('RePlot',hObject, eventdata, handles);
    case 'UploadPrevious'
        ImageUpload('UploadPrevious',hObject, eventdata, handles);
    case 'Upload'
        ImageUpload('Upload',hObject, eventdata, handles);
    case 'UploadNext'
        ImageUpload('UploadNext',hObject, eventdata, handles);
    case 'Delete'
        ImageUpload('Delete',hObject, eventdata, handles);
    case 'RunCPS'
        [varargout{1}]=ImageCPS('Run',hObject, eventdata, handles);
    case 'RunCPSOnce'
        [varargout{1}]=ImageCPS('RunOnce',hObject, eventdata, handles);
    case 'StopCPS'
        ImageCPS('Stop',hObject, eventdata, handles);
    case 'CursorType'
        ImageSetXYVoltage('CursorType',hObject, eventdata, handles);
    case 'Arrow'
        ImageSetXYVoltage('Arrow',hObject, eventdata, handles);
    case 'Fix'
        ImageSetXYVoltage('Fix',hObject, eventdata, handles);
    case 'SetZero'
        ImageSetXYVoltage('SetZero',hObject, eventdata, handles);
    case 'Cursor'
        ImageSetXYVoltage('Cursor',hObject, eventdata, handles);
    case 'Zoom On'
        ImageZoom('Zoom On',hObject, eventdata, handles);
    case 'ScreenCap'
        ImageScreenCap(hObject, eventdata, handles);
    case 'NVAdd'
        NVCatalog_Add(hObject, eventdata, handles);
    case 'NVRemove'
        NVCatalog_Update_Or_Remove(hObject, eventdata, handles, 1);
    case 'NVUpdate'
        NVCatalog_Update_Or_Remove(hObject, eventdata, handles, 0);
    otherwise
end


function NewTrackFast(hObject, eventdata, handles) 

global Enhanced_Factor gof countsExport dim
global scan1 scan2 scan3 count1 count2 count3 fit1 fit2 fit3
  
    dx = 0.005;
    dy = 0.005;
    dz = 4.0;
    dt = 0.01;
    NN = 21;
    
counter = '';
if handles.sr620Button.Value
    counter = 'SR620';
end
    global gScan gbManChange gConfocal; 

    % measure initial count
    WriteVoltage('Dev1/ao0',gScan.FixVx + gConfocal.XOffSet);
    WriteVoltage('Dev1/ao1',gScan.FixVy + gConfocal.YOffSet);
    WriteVoltage('Obj_Piezo',gScan.FixVz); 
    
    %Edit 7/19/2017--Leon Otis
    if ~strcmp(counter,'SR620')
        
        CPS_initial = measureCountRate(dt);
    else
       
        CPS_initial = measureSR620CountRate(dt);
    end
       
    % tracking start
    disp('Tracking NV...');
    
    x0 = gScan.FixVx; % this is in voltage
    y0 = gScan.FixVy;
    z0 = gScan.FixVz; 
    
    xscan = [x0-dx/2:dx/(NN-1):x0+dx/2];
    yscan = [y0-dy/2:dy/(NN-1):y0+dy/2];
%     if handles.Random_Walk.Value==1
%         Nn=round(NN/2);
%         zscan = [z0-dz/2:dz/(Nn-1):z0+dz/2];
%     else
    zscan = [z0-dz/2:dz/(NN-1):z0+dz/2];
%     end
    
    offset=[gConfocal.XOffSet gConfocal.YOffSet 0];
    
    %Edit Thomas Mittiga 5/1/2017
    gof=cell(1,3);
    scanS=zeros(6,NN);
    for i=1:dim % iterate x, y, and z direction
        %End Edit
        % scanning device selection
        if i == 1
            %disp('Tracking an ideal X position...');
            device = 'Dev1/ao0';
            scan = xscan; % change x
            WriteVoltage('Dev1/ao1', y0 + offset(2)); % fixed y
            WriteVoltage('Obj_Piezo', z0); % fixed z
            NNit=NN;
        elseif i == 2
            %disp('Tracking an ideal Y position...');
            device = 'Dev1/ao1';
            scan = yscan; % change y
            WriteVoltage('Dev1/ao0', x0 + offset(1)); % fixed x
            WriteVoltage('Obj_Piezo', z0); % fixed z
            NNit=NN;
        else
            %disp('Tracking an ideal Z position...');
            device = 'Obj_Piezo';
            scan = zscan; % change y
            WriteVoltage('Dev1/ao0', x0 + offset(1)); % fixed x
            WriteVoltage('Dev1/ao1', y0 + offset(2)); % fixed y
            pause(1)
%             if handles.Random_Walk.Value==1
%                 NNit=Nn;
%             else
%                 NNit=NN;
%             end
        end
            
        count =[];
        for j=1:NNit

            if get(handles.StopTracking, 'Value'), break; end
            
            %get counts
            WriteVoltage(device, scan(j) + offset(i));
            if i==3
                pause(1);
                if j==1
                    pause(7); %First time setting z, wait longer
                end
            end
            
             %Edit 7/19/2017
             
            if ~strcmp(counter,'SR620')
                CPS0 = measureCountRate(dt);
            else
                
                CPS0 = measureSR620CountRate(dt);
            end
            count = [count CPS0];
            %End Edit
        end
     

        g = fittype('a*exp(-4*log(2)*((x-b)/c).^2)+d');
        
        if i == 3
            bguess = 3;
            [f,gof{i}] = fit(scan',count',g,'StartPoint',[max(count) scan(round(NN/2)) bguess min(count)]);
        else
            rSq=-1*10^-30; %sometime rsquared can be negative
            for bguess = 0.001:0.001:0.005;
                [f,gof{i}] = fit(scan',count',g,'StartPoint',[max(count) scan(round(NN/2)) bguess min(count)]);
                if gof{i}.rsquare>rSq
                    rSq=gof{i}.rsquare;
                    fFinal=f;
                end
            end
            f=fFinal;
        end
        disp(gof{i})
        scanS(2*i-1,:)=scan'; %save temporarily
        scanS(2*i,:)=count;
        
        maxloc = f.b;
        best=f.a;
        RandWalk=true;
        if maxloc > max(scan) || maxloc < min(scan)
            maxind = find(f(scan) == max(f(scan)));
            maxloc = scan(maxind);
            best=count(maxind);
            RandWalk=false;
        end
        
        if i == 1
            x0 = maxloc;
            gScan.FixVx = x0;            
        elseif i == 2
            y0 = maxloc;
            gScan.FixVy = y0; 
        else
            z0 = maxloc;
            gScan.FixVz = z0;
        end
        
        if get(handles.bTrackPlot, 'Value')
            figure(1000);
            subplot(2,4,i)
%             cla;
            if i~=3
                scan=scan/gConfocal.Vx_per_um;
                plot(scan,count,'.r',scan,f(scan*gConfocal.Vx_per_um),'-b')
                xlim([min(scan) max(scan)])
            else
                plot(scan,count,'.r',scan,f(scan),'-b')
                xlim([min(scan) max(scan)])
            end          
            subtitle = 'xyz';
           % title(subtitle(i))
            %xlim([min(scan) max(scan)])
            xlabel([subtitle(i) ' (um)'])
            if i==1
            ylabel('Counts (C/s)')
            end
        end
        
        %Edit 3/20/2017
        if i==3 && RandWalk && handles.Random_Walk.Value==1 %Don't rand walk if at edge
            mcounts=[];
            %%Rectify z piezo hysteresis with random walk track
            WriteVoltage('Obj_Piezo',z0); % go to best z position
            pause(3)
            zCurr=z0;
            zStep = dz/(NN-1);
            bestN=0; %number of times this point was best
            while bestN<4 && handles.Random_Walk.Value==1
                if ~strcmp(counter,'SR620')  
                    Cmax=measureCountRate(dt);
                else
                    Cmax = measureSR620CountRate(dt);
                end
                if rand>0.5
                    pm=1;
                else
                    pm=-1;
                end
                
                if zCurr+pm*zStep>max(scan) || zCurr+pm*zStep<min(scan)
                    pm=-1*pm; %don't go out of range
                end
                
                WriteVoltage('Obj_Piezo',zCurr+pm*zStep);
                pause(0.3)
                if ~strcmp(counter,'SR620')
                    newVal=measureCountRate(dt);
                else
                    newVal=measureSR620CountRate(dt);
                end
                if newVal > Cmax
                    zCurr = zCurr + pm*zStep;
                    Cmax = newVal;
                end
                
                WriteVoltage('Obj_Piezo',zCurr);
                mcounts=[mcounts Cmax];
                
                if Cmax>0.95*best
                    bestN=bestN+1;
                end
                pause(0.3)
            end
            figure(1000);
            subplot(2,4,4);
            plot(mcounts)
            z0=zCurr;
            gScan.FixVz = z0;
        end

        %End Edit
        
        if i == 3 && handles.Random_Walk.Value==1
            CPS_final = Cmax;
            Enhanced_Factor = (CPS_final/CPS_initial);
            countsExport = CPS_final;
        elseif i==3
            maxind = find(f(scan) == max(f(scan)));
            CPS_final = count(maxind);
            Enhanced_Factor = (CPS_final/CPS_initial);
            countsExport = CPS_final;
        end
        
        
        if i==1
            scan1=scan;
            count1=count;
            fit1=f;
        elseif i==2
            scan2=scan;
            count2=count;
            fit2=f;
        else
            scan3=scan;
            count3=count;
            fit3=f;
        end
     
    end   
    
    
%Check goodness of fit criterion
saveQ=false;
for i=1:3
    if isstruct(gof{i})
       if gof{i}.rsquare <0.5
          saveQ=true ;
       end
    end
end

if saveQ
    global NVCat_path
    fid = fopen([NVCat_path 'trackinglog.txt'],'a'); %Open as appending 
    for i=1:6
        for j=1:length(scanS(i,:))
            fprintf(fid,'%d\t', scanS(i,j));
        end
    end
    fclose(fid);
end
set(handles.FixVx,'String',num2str(x0)); 
set(handles.FixVy,'String',num2str(y0)); 
set(handles.FixVz,'String',num2str(z0)); 

%go to new position
WriteVoltage('Dev1/ao0',x0 + gConfocal.XOffSet);
WriteVoltage('Dev1/ao1',y0 + gConfocal.YOffSet);
WriteVoltage('Obj_Piezo',z0); 

set(handles.StopTracking,'Value',0); 
%disp('Tracking finished!')
disp(['Done! Signal Enhancement Factor : ' num2str(Enhanced_Factor)])
    
    
function countRate = measureCountRate(dt)
%%%%%% measure count rate %%%%%%
global hCPS;

NRead = 10;
DT = dt;
TimeOut = DT * NRead * 1.1;
Freq = 1/DT;

[~, hPulse] = DigPulseTrainCont(Freq,0.5,10000);

hCounter = SetCounter(NRead+1);

hCPS.hCounter = hCounter;
hCPS.hPulse = hPulse;

status = DAQmxStartTask(hCounter);  DAQmxErr(status);
status = DAQmxStartTask(hPulse);    DAQmxErr(status);

DAQmxWaitUntilTaskDone(hCounter,TimeOut);
DAQmxStopTask(hPulse);

A = ReadCounter(hCounter,NRead+1);
DAQmxStopTask(hCounter);
A = diff(A);

countRate = sum(A)/(NRead * DT);

DAQmxClearTask(hPulse);
DAQmxClearTask(hCounter); 
    
 %Function for measuring count rate when tracking and scanning are done with SR620
function countRate = measureSR620CountRate(dt)

%TimeIntervalCounter.m must be run first so that global variable SR has
%been initialized.
global SR;

%Change SR620 to Count mode
fprintf(SR.gpib,'MODE 6');

%Set SR620 to take counts for time dt and do a measurement.
gate = ['GATE',' ',num2str(dt)];
fprintf(SR.gpib, gate);
fprintf(SR.gpib, 'SIZE 1');

%Do measurements first on channel A and then on channel B, combining
%results to give one count rate.
fprintf(SR.gpib, 'SRCE 0');
fprintf(SR.gpib, 'MEAS? 0');
cr1 = fscanf(SR.gpib);

fprintf(SR.gpib, 'SRCE 1');
fprintf(SR.gpib, 'MEAS? 0');
cr2 = fscanf(SR.gpib);

cr1 = str2double(cr1);
cr2 = str2double(cr2);
cr = cr1+cr2;
%Return an averaged rate that is in counts per second
countRate = cr/dt;
   
function DrawCrossHairs(hObject,eventdata,handles)
global gScan

% get the current limits of the axes1
XLimits = get(handles.axes1,'XLim');
YLimits = get(handles.axes1,'YLim');

CP = get(handles.axes1,'CurrentPoint');
p.x = CP(1,1);
p.y = CP(1,2);

% draw two lines
if isfield(gScan,'xLineHandle'),
    if ishandle(gScan.xLineHandle),
        set(gScan.xLineHandle,'XData',[XLimits(1),XLimits(2)]);
        set(gScan.xLineHandle,'YData',[p.y p.y]);
    else
        gScan.xLineHandle = line([XLimits(1),XLimits(2)],[p.y p.y],[1,1]);
    end
else
    gScan.xLineHandle = line([XLimits(1),XLimits(2)],[p.y p.y],[1,1]);
end


if isfield(gScan,'yLineHandle'),
    if ishandle(gScan.yLineHandle),
        set(gScan.yLineHandle,'XData',[p.x p.x]);
        set(gScan.yLineHandle,'YData',[YLimits(1),YLimits(2)]);
    else
        gScan.yLineHandle = line([p.x p.x],[YLimits(1),YLimits(2)],[1,1]);
    end
else
    gScan.yLineHandle = line([p.x p.x],[YLimits(1),YLimits(2)],[1,1]);
end

function PredefinedScans(what,hObject,eventdata,handles)

switch what
    case 'LargeScan'
        set(handles.minVx,'String',-0.4);
        set(handles.maxVx,'String',+0.4);
        set(handles.minVy,'String',-0.4);
        set(handles.maxVy,'String',+0.4);
        set(handles.bFixVx,'Value',0);
        set(handles.bFixVy,'Value',0);
        set(handles.bFixVz,'Value',1);
        set(handles.bFixDT,'Value',1);        
        set(handles.NVz,'String',10);
        set(handles.FixDT,'String',0.001);
    case 'ZScan'
        set(handles.bFixVx,'Value',1);
        set(handles.bFixVy,'Value',1);
        set(handles.bFixVz,'Value',0);
        set(handles.bFixDT,'Value',1); 
        set(handles.FixDT,'String',0.01);
        set(handles.NVz,'String',100);
    case 'NormalScan'
        set(handles.bFixVx,'Value',0);
        set(handles.bFixVy,'Value',0);
        set(handles.bFixVz,'Value',1);
        set(handles.bFixDT,'Value',1);        
        set(handles.NVz,'String',10);
        set(handles.FixDT,'String',0.001);
end

function Track(hObject,eventdata,handles)
global gScan gConfocal gTracking;

%% added new functionality for Track Continuously checkbox
TrackCenter(1);
bTrackContinuously = get(handles.cbTrackCont,'Value');
while bTrackContinuously,
    pause(10);
    TrackCenter(1);
    bTrackContinuously = get(handles.cbTrackCont,'Value');
    L=size(gTracking.Trajectory,1)
    if L>1,
        gTracking.Trajectory(L-1,:);
        gTracking.Trajectory(L,:);
        sum(abs(gTracking.Trajectory(end-1,:) - gTracking.Trajectory(end,:)));
    end
end;

set(handles.FixVx,'String',gScan.FixVx);
set(handles.FixVy,'String',gScan.FixVy);
set(handles.FixVz,'String',gScan.FixVz);

function ApplyOffSet(hObject, eventdata, handles)
global gConfocal;

gConfocal.XOffSet = eval(get(handles.XOffSet,'String'));
gConfocal.YOffSet = eval(get(handles.YOffSet,'String'));

ImageSetXYVoltage('Fix',hObject, eventdata, handles);

SaveConfocalFile;

function SaveConfocalFile
global gConfocal;
fid = fopen([gConfocal.path gConfocal.filename],'wt');
fprintf(fid,'XOffSet= %f\n',gConfocal.XOffSet);
fprintf(fid,'YOffSet= %f',gConfocal.YOffSet);
fclose(fid);

function ImageStart(hObject, eventdata, handles)
global gScan gSaveImg;
global gConfocal gmSEQ dim NVCat_path;

dim=3;

ImageFillUpForm('Start', hObject, eventdata, handles);

LoadNIDAQmx
% if ~libisloaded('mynidaqmx')
%     disp('Matlab: Load nicaiu.dll')
%     warning('off','MATLAB:loadlibrary:TypeNotFound')
%     
%     % This is to load ni library
%     loadlibrary('nicaiu.dll','nidaqmx.h','alias','mynidaqmx');
%     
%     warning('on','MATLAB:loadlibrary:TypeNotFound')
%     %funclist = libfunctions('myni','-full')
%     %libfunctionsview('myni')
% end

DAQmx_Val_Volts= 10348; % measure volts

%Read Starting File
gConfocal.filename = 'StartConfocal.txt';
gConfocal.path = 'C:\MATLAB Scripts\MATLAB_Code\Sets\';
gConfocal.positionfile='CurrentNV_Position.txt';
gConfocal.V_per_um='V_per_um.txt';
ReadStartingFile(handles);

NVCat_path = 'C:\MATLAB Scripts\MATLAB_Code\Sets\';

% #EO
% Load EO-Drive Piezo libraries 
warning('off', 'MATLAB:loadlibrary:ClassIsLoaded')
if libisloaded('C:\Program Files\Edmund Optics\EO-Drive\EO-Drive.dll') == 0
    loadlibrary('C:\Program Files\Edmund Optics\EO-Drive\EO-Drive.dll','C:\Program Files\Edmund Optics\EO-Drive\eo-drive.h','alias','EO-Drive');    
end
calllib('EO-Drive','EO_ReleaseAllHandles');
warning('on', 'MATLAB:loadlibrary:ClassIsLoaded')

global EO_handle;
EO_handle = int32(calllib('EO-Drive','EO_InitHandle')); %handle for EO Piezo
if EO_handle==0
    warning('Piezo was not properly initialized.')
end

set(handles.XOffSet,'String',gConfocal.XOffSet);
set(handles.YOffSet,'String',gConfocal.YOffSet);
set(handles.bTrackPlot,'Value',1);
set(handles.TrackTimes,'Value',2);
set(handles.bSaveImg,'Value',1);
gmSEQ.meas='SPCM';

function ReadStartingFile(handles)
global gConfocal;

file = [gConfocal.path gConfocal.filename];
fid = fopen(file,'rt');

gConfocal.XOffSet = fscanf(fid,'XOffSet= %f'); fscanf(fid,'\n');
gConfocal.YOffSet = fscanf(fid,'YOffSet= %f');

fclose(fid);

file = [gConfocal.path gConfocal.positionfile];
fid = fopen(file,'rt');
set(handles.FixVx,'String',num2str(fscanf(fid,'X= %f')-gConfocal.XOffSet,4)); fscanf(fid,'\n');
set(handles.FixVy,'String',num2str(fscanf(fid,'Y= %f')-gConfocal.YOffSet,4)); fscanf(fid,'\n');
set(handles.FixVz,'String',num2str(fscanf(fid,'Z= %f'),4));
fclose(fid);

file = [gConfocal.path gConfocal.V_per_um];
fid = fopen(file,'rt');
gConfocal.Vx_per_um = fscanf(fid,'Vx= %f'); fscanf(fid,'\n');
gConfocal.Vy_per_um = fscanf(fid,'Vy= %f');
fclose(fid);

function ImageFillUpForm(what, hObject, eventdata, handles)
global gScan;

switch what
    case 'Start'
        gScan = DefaultScan;
        UpdateScan(gScan,hObject, eventdata, handles);
    case 'UpdateScan'
        UpdateScan(gScan,hObject, eventdata, handles);
    otherwise
end

function Scan = DefaultScan

Scan.minVx = -0.025; Scan.maxVx = +0.025; Scan.NVx = 51;  
Scan.minVy = -0.025; Scan.maxVy = +0.025; Scan.NVy = 51;
Scan.minVz = +0; Scan.maxVz = +100; Scan.NVz = 51;
Scan.minDT = 1e-3; Scan.maxDT = +3e-3; Scan.NDT = 10;

Scan.bFixVx = 0;  Scan.FixVx = 0;
Scan.bFixVy = 0;  Scan.FixVy = 0;
Scan.bFixVz = 1;  Scan.FixVz = 50;
Scan.bFixDT = 1;  Scan.FixDT = .005;


function UpdateScan(Scan,hObject, eventdata, handles)
global gSaveImg;

set(handles.minVx,'String',num2str(Scan.minVx,4));
set(handles.maxVx,'String',num2str(Scan.maxVx,4));
set(handles.NVx,'String',num2str(Scan.NVx,4));
set(handles.bFixVx,'Value',Scan.bFixVx);
set(handles.FixVx,'String',num2str(Scan.FixVx,4));

set(handles.minVy,'String',num2str(Scan.minVy,4));
set(handles.maxVy,'String',num2str(Scan.maxVy,4));
set(handles.NVy,'String',num2str(Scan.NVy,4));
set(handles.bFixVy,'Value',Scan.bFixVy);
set(handles.FixVy,'String',num2str(Scan.FixVy,4));

set(handles.minVz,'String',Scan.minVz);
set(handles.maxVz,'String',Scan.maxVz);
set(handles.NVz,'String',Scan.NVz);
set(handles.bFixVz,'Value',Scan.bFixVz);
set(handles.FixVz,'String',Scan.FixVz);

set(handles.minDT,'String',Scan.minDT);
set(handles.maxDT,'String',Scan.maxDT);
set(handles.NDT,'String',Scan.NDT);
set(handles.bFixDT,'Value',Scan.bFixDT);
set(handles.FixDT,'String',Scan.FixDT);

xlabel(handles.axes1,'distance [um]');
ylabel(handles.axes1,'distance [um]');

set(handles.bSaveImgEPS,'Value',0);
set(handles.bSaveImgData,'Value',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%           ImageScan             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageScan(what,hObject, eventdata, handles)
global gScan;

gScan = GetScanValues(gScan, hObject, eventdata, handles);
%PrepareAxes(gScan, hObject, eventdata, handles);
switch what
    case 'Scan'
        CallMakeScan(gScan,hObject, eventdata, handles);
    case 'StopScan'
        StopScan;
    case 'Align'
        Align(gScan,hObject, eventdata, handles);
    case 'StopAlign'
        StopAlign;
    case 'ScanZ'
        ScanZ(gScan,hObject, eventdata, handles);
    case 'SetRange'
        SetRange(hObject, eventdata, handles);
    case 'GetScanValues'
        return;
    otherwise
end

function SetRange(hObject, eventdata, handles)
global gScan;

axes(handles.axes1);
sAxes = get(gca);
gScan.minVx = sAxes.XLim(1);
gScan.maxVx = sAxes.XLim(2);
gScan.minVy = sAxes.YLim(1);
gScan.maxVy = sAxes.YLim(2);

RX = gScan.maxVx - gScan.minVx;
RY = gScan.maxVy - gScan.minVy;

if RX>RY
    gScan.NVy = round( gScan.NVx * RY / RX ); 
else
    gScan.NVx = round( gScan.NVy * RX / RY ); 
end

ImageFillUpForm('UpdateScan', hObject, eventdata, handles);



function ScanZ(Scan,hObject, eventdata, handles)
% global bGo hTasks;
% 
% bGo = true;
% 
% VZ = Scan.minVz:(Scan.maxVz-Scan.minVz)/(Scan.NVz-1):Scan.maxVz;
% VZinv = Scan.maxVz:(Scan.minVz-Scan.maxVz)/(Scan.NVz-1):Scan.minVz;
% V = [VZ VZinv];
% N = length(V);
% 
% plot(handles.axes1,V);
% 
% Freq = 1/Scan.FixDT;
% TimeOut = Scan.FixDT * Scan.NVz * 2;
% 
% status = DAQmxResetDevice('Dev1');
% disp(['NI: Reset Devie          :' num2str(status)]);
% 
% [status, hPulse] = DigPulseTrainCont(Freq,0.5,10000);
% DAQmxErr(status);
% 
% hAlign = SetScanZ(V,N,Freq);
% 
% hTasks.hAlign = hAlign;
% hTasks.hPulse = hPulse;
% N=0;
% while N<1000
%     N=N+1;
%     axes(handles.axes1);
%     status = DAQmxStartTask(hAlign);    disp(['NI: Start AO task        :' num2str(status)]);
%     status = DAQmxStartTask(hPulse);    disp(['NI: Start Pulse          :' num2str(status)]);
% 
%     DAQmxWaitUntilTaskDone(hAlign,TimeOut);
% 
%     DAQmxStopTask(hPulse);
%     DAQmxStopTask(hAlign);
% 
%     if ~bGo, break; end
% end
%     
% DAQmxClearTask(hPulse);
% DAQmxClearTask(hAlign);


function CallMakeScan(gScan,hObject, eventdata, handles);
global bGo
bScanCont = get(handles.bScanCont,'Value');
bGo = true;
if bScanCont
    R = 100;
else
    R = 0;
end
r = 0;
while bGo & r<=R
    r = r+1;
    MakeScan(gScan,hObject, eventdata, handles);
    bScanCont = get(handles.bScanCont,'Value');
    if ~bScanCont,break; end
end



function StopAlign
global hTasks bGo;

bGo = false;

hAlign = hTasks.hAlign;
hPulse = hTasks.hPulse ;

DAQmxStopTask(hPulse);
DAQmxStopTask(hAlign);

DAQmxClearTask(hPulse);
DAQmxClearTask(hAlign);


function Align(Scan,hObject, eventdata, handles)
global hTasks bGo gConfocal;

bGo = true;

ds = (Scan.maxVx-Scan.minVx)/(Scan.NVx-1);
R = sqrt(Scan.maxVx^2 + Scan.maxVy^2);
%2piR/N =ds
N = 10 * 2*pi*R/ds;
THETA = 0:10 * 2*pi/(N-1):10 * 2*pi;
VX =  gConfocal.XOffSet + R*cos(THETA);
VY =  gConfocal.YOffSet + R*sin(THETA);

%plot(THETA,VX,'b',THETA,VY,'r')
%return;

V = [VX VY];
N = length(VX);

Freq = 1/Scan.FixDT;
TimeOut = Scan.FixDT * N * 1.1;

% status = DAQmxResetDevice('Dev1');
% disp(['NI: Reset Devie          :' num2str(status)]);

[status, hPulse] = DigPulseTrainCont(Freq,0.5,10000);
if status
disp(['NI: Set Pulse            :' num2str(status)]);
end

hAlign = SetXYAlign(V,N,Freq);

hTasks.hAlign = hAlign;
hTasks.hPulse = hPulse;
N=0;
while N<1000
    N=N+1;
    axes(handles.axes1);
    status = DAQmxStartTask(hAlign);     disp(['NI: Start AO task        :' num2str(status)])
    status = DAQmxStartTask(hPulse);    disp(['NI: Start Pulse          :' num2str(status)])

    DAQmxWaitUntilTaskDone(hAlign,TimeOut);

    DAQmxStopTask(hPulse);
    DAQmxStopTask(hAlign);

    if ~bGo, break; end
end
    
DAQmxClearTask(hPulse);
DAQmxClearTask(hAlign);

%Come back to the originalposition
WriteVoltage(1,Scan.FixVx + gConfocal.XOffSet);
WriteVoltage(2,Scan.FixVy  + gConfocal.YOffSet);



function StopScan
global hTasks bGo;

bGo = false;

counter = '';
if handles.sr620Button.Value
    counter = 'SR620';
end

hScan = hTasks.hScan;
hCounter = hTasks.hCounter;
hPulse = hTasks.hPulse ;

if ~strcmp(counter, 'SR620')
    DAQmxStopTask(hPulse);
    DAQmxStopTask(hScan);
    DAQmxStopTask(hCounter);
    DAQmxClearTask(hPulse);
    DAQmxClearTask(hScan);
    DAQmxClearTask(hCounter);
end


function PrepareAxes(Scan, hObject, eventdata, handles)
axes(handles.axes1);
xlim([Scan.minVx Scan.maxVx]);
ylim([Scan.minVy Scan.maxVy]);

function MakeScan(Scan,hObject, eventdata, handles)
global hTasks bGo gScan gConfocal gmSEQ;

bGo = true;

counter = '';
if handles.sr620Button.Value
    counter = 'SR620';
end

planScan = PlanScan(Scan,hObject, eventdata, handles);
%planScan.ContRange
%planScan.FLRange
%status = DAQmxResetDevice('Dev1');
%disp(['NI: Reset Devie          :' num2str(status)]);

%Turn On Laser : JM 2008-07-27
% PBFunctionPool('PBON',0);

I = zeros(planScan.SizeImg);

itl = 0;

%Allocate array for the data of a row of the scan if using SR620
rowdata = zeros(1,planScan.NRead*2+1);
dt = str2double(handles.FixDT.String);

if ~strcmp(counter, 'SR620')
    for tl = planScan.TL
        itl = itl + 1;
        
        [planScan] = SetPulseAndContScan(3,tl,planScan,hObject, eventdata, handles);
        
        WriteVoltage(planScan.WhatTL,tl);
        isl = 0;
        for sl = planScan.SL
            isl = isl + 1;
            
            [planScan] = SetPulseAndContScan(2,sl,planScan,hObject, eventdata, handles);
            
            WriteVoltage(planScan.WhatSL,sl);
            ifl = 0;
            temp = 0;
            for fl = planScan.FL
                ifl = ifl + 1;
                [planScan] = SetPulseAndContScan(1,fl,planScan,hObject, eventdata, handles);
                WriteVoltage(planScan.WhatFL,fl);
                
                    
                    status = DAQmxStartTask(planScan.hCounter);
                    DAQmxErr(status);
                    status = DAQmxStartTask(planScan.hScan);
                    DAQmxErr(status);
                    status = DAQmxStartTask(planScan.hPulse);
                    DAQmxErr(status);

                    if ~get(handles.bFastScan,'Value')
                        DAQmxWaitUntilTaskDone(planScan.hScan,planScan.TimeOut);
                    else
                        DAQmxWaitUntilTaskDone(planScan.hScan,tl*planScan.NRead);
                    end

                    DAQmxStopTask(planScan.hPulse);
                    DAQmxStopTask(planScan.hScan);
                    A = ReadCounter( planScan.hCounter, planScan.NRead*2+1);
                    A = ProcessDataVector(A, planScan.DT);
                    I(ifl,:) = (A(1:planScan.NRead)+fliplr(A(planScan.NRead+1:planScan.NRead*2)))/2;


                    DAQmxStopTask(planScan.hCounter);
               
              
                % progress description
                Scan_Prog = round(100*ifl/length(planScan.FL));
                if Scan_Prog == 1 || Scan_Prog == 25 || Scan_Prog == 50 || Scan_Prog == 75 || Scan_Prog == 100 
                    if temp ~= Scan_Prog
                        disp(['Scanning Progress   :  ' num2str(Scan_Prog) ' %']);
                        temp = Scan_Prog;
                    end
                end
                
                %plotting takes ~200 ms
                if ~get(handles.bFastScan,'Value') || fl==planScan.FL(end)
                PlotScan(I,planScan,hObject, eventdata, handles,'Quick');

                CPS = ProcessDataCPS(I(ifl,:),planScan.NRead,1);
                set(handles.CPS,'String',CPS);
                end
                drawnow;
                if ~bGo, break; end
            end
            
            
            if get(handles.bSaveImg,'Value')
                if planScan.WhatSL == 3
                    gScan.FixVz = sl;
                    set(handles.FixVz,'String',num2str(gScan.FixVz));
                    disp(['z = ' num2str(sl) ' um' '  &  ' num2str(isl) ' / ' num2str(length(planScan.SL)) ' images done']);
                end
                ImageSaveImage('Save',hObject, eventdata, handles);
            end
            drawnow;
            if ~bGo, break; end
        end
        drawnow;
        if ~bGo, break; end
    end
    
 else
     %Use SR620 instead of DAQ for this block of loops. 
     %Need to put count rate data in the variable I in the right way so it can 
     %be used to plot the scan.
     for tl = planScan.TL
        
        WriteVoltage(planScan.WhatTL,tl);
        isl = 0;
        for sl = planScan.SL
            isl = isl + 1;
            WriteVoltage(planScan.WhatSL,sl);
            ifl = 0;
            temp = 0;
            for fl = planScan.FL
                ifl = ifl + 1;
                WriteVoltage(planScan.WhatFL,fl);
                for i = 1:length(planScan.Cont)
                    
                    %This line should allow scan to end more quickly if stop
                    %button is pressed.
                    if ~bGo, break; end
                    
                    %Count rate data for constructing the scan is collected
                    %point by point with the SR620 while the DAQ is able to
                    %quickly collect a row of pixels.
                    WriteVoltage(planScan.WhatCont,planScan.Cont(i));
                    rowdata(i) = measureSR620CountRate(dt);
                end
                I(ifl,:) = (rowdata(1:planScan.NRead)+fliplr(rowdata(planScan.NRead+1:planScan.NRead*2)))/2;
               
                
                 %Repeat code needed to plot the scan as it progresses so
                 %that scanning with SR620 behaves the same way as with DAQ.
                 
                % progress description
                    Scan_Prog = round(100*ifl/length(planScan.FL));
                    if Scan_Prog == 1 || Scan_Prog == 25 || Scan_Prog == 50 || Scan_Prog == 75 || Scan_Prog == 100 
                        if temp ~= Scan_Prog
                            disp(['Scanning Progress   :  ' num2str(Scan_Prog) ' %']);
                            temp = Scan_Prog;
                        end
                    end

                    %plotting takes ~200 ms
                    if ~get(handles.bFastScan,'Value') || fl==planScan.FL(end)
                    PlotScan(I,planScan,hObject, eventdata, handles,'Quick');

                    CPS = ProcessDataCPS(I(ifl,:),planScan.NRead,1);
                    set(handles.CPS,'String',CPS);
                    end
                    drawnow;
                    if ~bGo, break; end
            end
            if get(handles.bSaveImg,'Value')
                    if planScan.WhatSL == 3
                        gScan.FixVz = sl;
                        set(handles.FixVz,'String',num2str(gScan.FixVz));
                        disp(['z = ' num2str(sl) ' um' '  &  ' num2str(isl) ' / ' num2str(length(planScan.SL)) ' images done']);
                    end
                    ImageSaveImage('Save',hObject, eventdata, handles);
            end
                drawnow;
                if ~bGo, break; end
        end
        drawnow;
        if ~bGo, break; end
     end     
end


if ~get(handles.bScanCont,'Value')
    PlotScan(I,planScan,hObject, eventdata, handles,'Final'); 
end


% if ~get(handles.bFastScan,'Value')
if ~strcmp(counter, 'SR620')
    DAQmxClearTask(planScan.hPulse);
    DAQmxClearTask(planScan.hScan);
    DAQmxClearTask(planScan.hCounter);
end
% else
%     DAQmxClearTask(planScan.hPulse(1));
%     DAQmxClearTask(planScan.hPulse(2));
%     DAQmxClearTask(planScan.hScan(1));
%     DAQmxClearTask(planScan.hScan(2));
%     DAQmxClearTask(planScan.hCounter);
% end

% added 13 Aug 2008
% if Marker is checked, make cross hairs
if get(handles.cbMarker,'Value') && ~get(handles.bScanCont,'Value')
    DrawCrossHairs(hObject,eventdata,handles);
end
% go to initial point
if ~get(handles.bScanCont,'Value')
    WriteVoltage('Dev1/ao0',gConfocal.XOffSet);
    WriteVoltage('Dev1/ao1',gConfocal.YOffSet);   
    WriteVoltage('Obj_Piezo',50);
end



function PlotScan(I,planScan,hObject, eventdata, handles, mode)
global gRawImg gConfocal;

if nargin == 5
    mode = 'Quick';
end

switch planScan.WhatCont
    case 1
        multiplier(1)=1/gConfocal.Vx_per_um;
    case 2
        multiplier(1)=1/gConfocal.Vy_per_um;
    otherwise
        multiplier(1)=1;
end

if planScan.ND==1
    switch mode
        case 'Quick'
%             plot(handles.axes1,planScan.Cont(1:planScan.NRead),I,'b');
%             xlabel(handles.axes1,planScan.LabelCont);
            plot(handles.axes1,planScan.Cont(1:planScan.NRead)*multiplier(1),I,'b');
            xlabel(handles.axes1,'distance [um]');
            drawnow;
        case 'Final'
%             plot(handles.axes1,planScan.Cont(1:planScan.NRead),I,'b');
%             xlabel(handles.axes1,planScan.LabelCont);
            plot(handles.axes1,planScan.Cont(1:planScan.NRead)*multiplier(1),I,'b');
            xlabel(handles.axes1,'distance [um]');
            set(handles.axes1,'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles});
    end
elseif planScan.ND==2
    switch planScan.WhatFL
        case 1
            multiplier(2)=1/gConfocal.Vx_per_um;
        case 2
            multiplier(2)=1/gConfocal.Vy_per_um;
        otherwise
            multiplier(2)=1;
    end
    switch mode
        case 'Quick'
            %axes(handles.axes1);

            imagesc(I,'XData',planScan.ContRange*multiplier(1),'YData',planScan.FLRange*multiplier(2),'Parent',handles.axes1);
%            daspect(handles.axes1,[1 1 1]);
            colormap(handles.axes1,pink);
            RePlot(hObject, eventdata, handles);
            hc = colorbar('peer',handles.axes1,'location','EastOutside'); 
            
            xlabel(handles.axes1,'distance [um]');
            ylabel(handles.axes1,'distance [um]');
            drawnow;
        case 'Final'
            %axes(handles.axes1);
            
            cla(handles.axes1);
            imagesc(I,'XData',planScan.ContRange*multiplier(1),'YData',planScan.FLRange*multiplier(2),...
                'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles},...
                'Parent',handles.axes1);
            RePlot(hObject, eventdata, handles);
            
            xlabel(handles.axes1,'distance [um]');
            ylabel(handles.axes1,'distance [um]');
            
            colormap(handles.axes1,pink);
            hc = colorbar('peer',handles.axes1,'location','EastOutside');
            drawnow;
            
%             xlabel(handles.axes1,planScan.LabelCont);
%             ylabel(handles.axes1,planScan.LabelFL);
            
%             daspect(handles.axes1,[1 1 1]);

            set(handles.axes1,'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles});
    end
end



function planScan = SetPulseAndContScan(cloop,cvalue,planScan,hObject, eventdata, handles)
global hTasks;

if cloop ~= planScan.DTloop
    return;
end
planScan.DT = cvalue;

%if ~get(handles.bFastScan,'Value') %%% Normal Scan
    planScan.TimeOut = cvalue * planScan.NRead * 4;

    %planScan.hCounter = SetCounter( planScan.NRead + 1 );
    planScan.hCounter = SetCounter( planScan.NRead*2 +1, 1/cvalue);
    [status, planScan.hPulse] = DigPulseTrainCont(1/cvalue,0.5,planScan.NRead*2+1);
    DAQmxErr(status);
    planScan.hScan = SetXYOutPut(planScan,1/cvalue,hObject, eventdata, handles);

    hTasks.hScan = planScan.hScan;
    hTasks.hCounter = planScan.hCounter;
    hTasks.hPulse = planScan.hPulse;

%    DupCount = DAQmxGet(planScan.hCounter, 'CI.DupCountPrevent', 'Dev1/ctr0');
%    DAQmxSet(planScan.hCounter, 'CI.DupCountPrevent', 'Dev1/ctr0',1);

% else %%% Fast Scan
%     planScan.TimeOut = cvalue * planScan.NRead * 4;
% 
%     Ntot = (planScan.NRead*2)*length(planScan.FL);
%     planScan.hCounter = SetCounter(Ntot);
%     
%     [status, planScan.hPulse(1)] = DigPulseTrainCont(1/cvalue,0.5,10000);
%     [statusFL, planScan.hPulse(2)] = DigPulseTrainCont(1/cvalue/(planScan.NRead*2),0.5,10000);
%     if status
%     disp(['NI: Set Pulse            :' num2str(status)]);
%     end
%     if statusFL
%     disp(['NI: Set PulseFL            :' num2str(statusFL)]);
%     end
%     planScan.hScan = SetXYOutPut(planScan,1/cvalue,hObject, eventdata, handles);
%     
%     hTasks.hCounter = planScan.hCounter;
%     hTasks.hScan = planScan.hScan;
%     % planScan.hScan(1),(2)
%     hTasks.hPulse = planScan.hPulse;
%     % planScan.hScan(1),(2)
% 
%     DupCount = DAQmxGet(planScan.hCounter, 'CI.DupCountPrevent', 'Dev1/ctr0');
%     DAQmxSet(planScan.hCounter, 'CI.DupCountPrevent', 'Dev1/ctr0',1);
% end
    

function planScan = PlanScan(Scan,hObject, eventdata, handles)
global gConfocal;

bScanX = ~get(handles.bFixVx,'Value');
bScanY = ~get(handles.bFixVy,'Value');
bScanZ = ~get(handles.bFixVz,'Value');
bScanDT = ~get(handles.bFixDT,'Value');

if bScanX
    VV{1} = Scan.minVx:(Scan.maxVx-Scan.minVx)/(Scan.NVx-1):Scan.maxVx;
    VVRange{1}= [min(VV{1}) max(VV{1})];
    VV{1} = VV{1} + gConfocal.XOffSet;

else
    VV{1} =  gConfocal.XOffSet + Scan.FixVx; 
    VVRange{1}= [min(VV{1}) max(VV{1})] - gConfocal.XOffSet;
end
if bScanY
    VV{2} =  Scan.minVy:(Scan.maxVy-Scan.minVy)/(Scan.NVy-1):Scan.maxVy;
    VVRange{2}= [min(VV{2}) max(VV{2})];
    VV{2} = VV{2} + gConfocal.YOffSet;

else
    VV{2} =  gConfocal.YOffSet + Scan.FixVy; 
    VVRange{2}= [min(VV{2}) max(VV{2})] - gConfocal.YOffSet;
end
if bScanZ
    VV{3} = Scan.minVz:(Scan.maxVz-Scan.minVz)/(Scan.NVz-1):Scan.maxVz;
    VVRange{3}= [min(VV{3}) max(VV{3})];
else
    VV{3} = Scan.FixVz; 
    VVRange{3}= [min(VV{3}) max(VV{3})];
end
if bScanDT
    VV{4} = Scan.minDT:(Scan.maxDT-Scan.minDT)/(Scan.NDT-1):Scan.maxDT;
    VVRange{4}= [min(VV{4}) max(VV{4})];
elseif bScanZ && ~bScanX && ~bScanY
    error('Please scan in 2D (XZ or YZ) rather than only in Z.');
else
    VV{4} = Scan.FixDT/2; 
    VVRange{4}= [min(VV{4}) max(VV{4})];
end

PS = [1 bScanX 1; 2 bScanY 2; 3 bScanZ 3; 4 bScanDT 4];
PS = sortrows(PS,[-2 3]);

Label{1} = 'V_x';
Label{2} = 'V_y';
Label{3} = 'V_z';
Label{4} = 'DT';

%Set Continuous
planScan.Cont = [VV{PS(1,1)}(1) VV{PS(1,1)} fliplr(VV{PS(1,1)})];
planScan.ContRange = VVRange{PS(1,1)};
planScan.WhatCont = PS(1,1);
planScan.NRead = (length(planScan.Cont)-1)/2;
planScan.LabelCont = Label{PS(1,1)};

planScan.FL = VV{PS(2,1)};
planScan.WhatFL = PS(2,1);
planScan.FLRange = VVRange{PS(2,1)};
planScan.LabelFL = Label{PS(2,1)};
planScan.Ntot = planScan.NRead*length(planScan.FL);

planScan.SL = VV{PS(3,1)};
planScan.WhatSL = PS(3,1);
planScan.SLRange = VVRange{PS(3,1)};
planScan.LabelSL = Label{PS(3,1)};

planScan.TL = VV{PS(4,1)};
planScan.WhatTL = PS(4,1);
planScan.TLRange = VVRange{PS(4,1)};
planScan.LabelTL = Label{PS(4,1)};

planScan.ND = min([2 sum(PS(:,2))]); % Number of Dimensions to plot

planScan.SizeImg = [length(planScan.FL) planScan.NRead ];

planScan.DTloop = find(PS(:,1)==4)-1;

planScan.hCounter = [];
planScan.hPulse = [];
planScan.hScan = [];

function CallImageSetXYVoltage(scr,eventdata,arg1,arg2)
ImageSetXYVoltage('Cursor',arg1,eventdata,arg2);

%WriteVoltage('Dev1/ao0',Scan.FixVx);
%WriteVoltage('Dev1/ao1',Scan.FixVy);
function readArray = ReadCounter(task,N)
numSampsPerChan = N;
timeout = 0;
readArray = zeros(1,N);
global gmSEQ
if strcmp(gmSEQ.meas,'SPCM')
    [status, readArray]= DAQmxReadCounterU32(task, numSampsPerChan, timeout, zeros(1,numSampsPerChan), numSampsPerChan, libpointer('int32Ptr',0) ); 
elseif strcmp(gmSEQ.meas,'APD')
    [status, readArray] = DAQmxFunctionPool('ReadAnalogVoltage',task, numSampsPerChan, timeout);
end
DAQmxErr(status);



% function task = SetCounter(N)
% DAQmx_Val_Volts= 10348; % measure volts
% DAQmx_Val_Rising = 10280; % Rising
% DAQmx_Val_FiniteSamps = 10178; % Finite Samples
% DAQmx_Val_CountUp = 10128; % Count Up
% DAQmx_Val_CountDown = 10124; % Count Down
% DAQmx_Val_GroupByChannel = 0; % Group per channel
% DAQmx_Val_ContSamps =10123; % Continuous Samples
% 
% [ status, TaskName, task ] = DAQmxCreateTask([]);
% 
% if status,
% disp(['NI: Create Counter Task  :' num2str(status)]);
% end
% status = DAQmxCreateCICountEdgesChan(task,'Dev1/ctr0','',...
%     DAQmx_Val_Rising , 0, DAQmx_Val_CountUp);
% 
% if status,
% disp(['NI: Create Counter       :' num2str(status)]);
% end
% status = DAQmxCfgSampClkTiming(task,'/Dev1/PFI13',1.0,...
%     DAQmx_Val_Rising,DAQmx_Val_FiniteSamps ,N);
% if status,
% disp(['NI: Cofigure the Clk     :' num2str(status)]);
% end

function [task] = SetCounter(varargin)
%varargin(1) is the number of total samples
%varargin(2) is the frequency of the gating to expect
% Initialize DAQ
global gmSEQ
if strcmp(gmSEQ.meas,'SPCM')
    [status, task ] = DAQmxFunctionPool('SetCounter',varargin{1});
elseif strcmp(gmSEQ.meas,'APD')
    [status, task ] = DAQmxFunctionPool('CreateAIChannel','/Dev1/PFI13',varargin{1},varargin{2});
end
DAQmxErr(status);

function task = SetXYOutPut(planScan,rate,hObject, eventdata, handles)
DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
DAQmx_Val_FiniteSamps = 10178; % Finite Samples
DAQmx_Val_CountUp = 10128; % Count Up
DAQmx_Val_CountDown = 10124; % Count Down
DAQmx_Val_GroupByChannel = 0; % Group per channel

switch planScan.WhatCont
    case 1
        Device = 'Dev1/ao0';
    case 2
        Device = 'Dev1/ao1';
    case 3
        Device = 'Obj_Piezo';
    case 4
        disp('This functionality is not supported! 2007-08-23');
        return;
end

V = planScan.Cont;
N = length(planScan.Cont);

%Setting the clock
[ status, TaskName, task ] = DAQmxCreateTask([]);
DAQmxErr(status);
    
status = DAQmxCreateAOVoltageChan(task,Device,-10,10,DAQmx_Val_Volts);
DAQmxErr(status);
    
status = DAQmxCfgSampClkTiming(task,'/Dev1/PFI13',rate,...
        DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,N);
DAQmxErr(status);

    
empty_ptr = libpointer('doublePtr',[]);
zero_ptr = libpointer('int32Ptr',zeros(1,N));
%    zero_ptr = libpointer('uint32Ptr',0);   %The number of events counted by the counter
status = DAQmxWriteAnalogF64(task, N, 0, 10,...
        DAQmx_Val_GroupByChannel, V, zero_ptr);
DAQmxErr(status);

% %%%% Fast Scan 
% if get(handles.bFastScan,'Value')
%     switch planScan.WhatFL
%     case 1
%         DeviceFL = 'Dev1/ao0';
%     case 2
%         DeviceFL = 'Dev1/ao1';
%     case 3
%         DeviceFL = 'Obj_Piezo';
%     case 4
%         disp('This functionality is not supported! 2014-05-24');
%         return;
%     end
%     
%     VFL = [planScan.FL fliplr(planScan.FL)];
%     NFL = 2*length(planScan.FL);
%     rateFL = rate/N;
%     
%     [ statusFL, TaskNameFL, taskFL ] = DAQmxCreateTask([]);
%     if statusFL,
%         disp(['NI: Create AO TaskFL       :' num2str(statusFL)]);
%     end
%     
%     statusFL = DAQmxCreateAOVoltageChan(taskFL,DeviceFL,-10,10,DAQmx_Val_Volts);
%     if statusFL,
%         disp(['NI: Create AO ChannelFL    :' num2str(statusFL)]);
%     end
%     
%     statusFL = DAQmxCfgSampClkTiming(taskFL,'/Dev1/PFI13',rateFL,...
%             DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,NFL);
%     if statusFL,
%         disp(['NI: Config Sample ClockFL  :' num2str(statusFL)]);
%     end
%     
%     empty_ptr = libpointer('doublePtr',[]);
%     zero_ptr = libpointer('int32Ptr',zeros(1,NFL));
%     %    zero_ptr = libpointer('uint32Ptr',0);   %The number of events counted by the counter
%     statusFL = DAQmxWriteAnalogF64(taskFL, NFL, 0, 10,...
%             DAQmx_Val_GroupByChannel, VFL, zero_ptr);
%     if statusFL,
%     disp(['NI: Write Analog F64     :' num2str(statusFL)])
%     end
%     task = [task taskFL];
% end


function task = SetXYAlign(V,N,rate)
DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
DAQmx_Val_FiniteSamps = 10178; % Finite Samples
DAQmx_Val_CountUp = 10128; % Count Up
DAQmx_Val_CountDown = 10124; % Count Down
DAQmx_Val_GroupByChannel = 0; % Group per channel

%Setting the clock
[ status, TaskName, task ] = DAQmxCreateTask([]);
DAQmxErr(status);
status = DAQmxCreateAOVoltageChan(task,'Dev1/ao0:1',-10,10,DAQmx_Val_Volts);
DAQmxErr(status);

status = DAQmxCfgSampClkTiming(task,'/Dev1/PFI13',rate,...
    DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,N);
DAQmxErr(status);

empty_ptr = libpointer('doublePtr',[]);
zero_ptr=libpointer('int32Ptr',zeros(1,2*N));
%    zero_ptr = libpointer('uint32Ptr',0);   %The number of events counted by the counter
status = DAQmxWriteAnalogF64(task, N, 0, 10,...
    DAQmx_Val_GroupByChannel, V, zero_ptr);
DAQmxErr(status);

function task = SetScanZ(V,N,rate)
DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
DAQmx_Val_FiniteSamps = 10178; % Finite Samples
DAQmx_Val_CountUp = 10128; % Count Up
DAQmx_Val_CountDown = 10124; % Count Down
DAQmx_Val_GroupByChannel = 0; % Group per channel

%Setting the clock
[ status, TaskName, task ] = DAQmxCreateTask([]);
disp(['NI: Create AO Task       :' num2str(status)]);

status = DAQmxCreateAOVoltageChan(task,'Dev1/ao2',-10,10,DAQmx_Val_Volts);
disp(['NI: Create AO Channel    :' num2str(status)]);

status = DAQmxCfgSampClkTiming(task,'/Dev1/PFI13',rate,...
    DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,N);
disp(['NI: Config Sample Clock  :' num2str(status)]);

empty_ptr = libpointer('doublePtr',[]);
zero_ptr=libpointer('int32Ptr',zeros(1,2*N));
%    zero_ptr = libpointer('uint32Ptr',0);   %The number of events counted by the counter
status = DAQmxWriteAnalogF64(task, N, 0, 10,...
    DAQmx_Val_GroupByChannel, V, zero_ptr);
disp(['NI: Write Analog F64     :' num2str(status)])

%%%% Joonhee (05/24/2014)
function task = ScanVoltage(what,minVoltage,maxVoltage,freq)
DAQmx_Val_Volts = 10348; % measure volts

switch what
    case {1,'Dev1/ao0'}
        Device = 'Dev1/ao0';
    case {2,'Dev1/ao1'}
        Device = 'Dev1/ao1';
    case {3,'Obj_Piezo'}
        % Device = 'Dev1/ao2'; #PI
        
        % #EO
        global EO_handle;

        if Voltage > 100
            Voltage = 100;
        end
        if Voltage < 0
            Voltage = 0;
        end

        calllib('EO-Drive','EO_Move', EO_handle, Voltage); 
        % #EO
       
        
        return;
        
    case 4
        return;
%     case 'ao_2' 
%         Device = 'Dev1/ao2';
%     case 'ao_3'
%         Device = 'Dev1/ao3';
    otherwise
        disp('Error in Write Voltage: I dont get it!');
end

[ status, TaskName, task ] = DAQmxCreateTask([]);
status = status + DAQmxCreateAOVoltageChan(task,Device,-10,10,DAQmx_Val_Volts);
status = status + DAQmxWriteAnalogScalarF64(task,1,0,Voltage);
if status ~= 0
    disp(['Error in writing voltage in Device ' Device]);
end
DAQmxClearTask(task);


function task = WriteVoltage(what,Voltage)
DAQmx_Val_Volts= 10348; % measure volts

switch what
    case {1,'Dev1/ao0'}
        Device = 'Dev1/ao0';
    case {2,'Dev1/ao1'}
        Device = 'Dev1/ao1';
    case {3,'Obj_Piezo'}
        % Device = 'Dev1/ao2'; #PI
        
        % #EO
        global EO_handle;

        if Voltage > 100
            Voltage = 100;
        end
        if Voltage < 0
            Voltage = 0;
        end

        calllib('EO-Drive','EO_Move', EO_handle, Voltage); 
        % #EO
       
        
        return;
        
    case 4
        return;
%     case 'ao_2' 
%         Device = 'Dev1/ao2';
%     case 'ao_3'
%         Device = 'Dev1/ao3';
    otherwise
        disp('Error in Write Voltage: I dont get it!');
end

[ status, TaskName, task ] = DAQmxCreateTask([]);
status = status + DAQmxCreateAOVoltageChan(task,Device,-10,10,DAQmx_Val_Volts);
status = status + DAQmxWriteAnalogScalarF64(task,1,0,Voltage);
if status ~= 0
    disp(['Error in writing voltage in Device ' Device]);
end
DAQmxClearTask(task);

function Scan = GetScanValues(Scan,hObject, eventdata, handles)

Scan.minVx = eval(get(handles.minVx,'String'));
Scan.maxVx = eval(get(handles.maxVx,'String'));
Scan.NVx = eval(get(handles.NVx,'String'));
Scan.bFixVx = get(handles.bFixVx,'Value');
Scan.FixVx = eval(get(handles.FixVx,'String'));

Scan.minVy = eval(get(handles.minVy,'String'));
Scan.maxVy = eval(get(handles.maxVy,'String'));
Scan.NVy = eval(get(handles.NVy,'String'));
Scan.bFixVy = get(handles.bFixVy,'Value');
Scan.FixVy = eval(get(handles.FixVy,'String'));

Scan.minVz = eval(get(handles.minVz,'String'));
Scan.maxVz = eval(get(handles.maxVz,'String'));
Scan.NVz = eval(get(handles.NVz,'String'));
Scan.bFixVz = get(handles.bFixVz,'Value');
Scan.FixVz = eval(get(handles.FixVz,'String'));

Scan.minDT = eval(get(handles.minDT,'String'));
Scan.maxDT = eval(get(handles.maxDT,'String'));
Scan.NDT = eval(get(handles.NDT,'String'));
Scan.bFixDT = get(handles.bFixDT,'Value');
Scan.FixDT = eval(get(handles.FixDT,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   IageSaveImage                    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImageSaveImage(what,hObject, eventdata, handles)
global gSaveImg
%Determine what to do
switch what
    case 'Save'
        SaveImage('Save',hObject, eventdata, handles);
    case 'ChangeFileName'
        ChangeFileName(hObject, eventdata, handles);
    case 'Edit'
        SaveImage('Edit',hObject, eventdata, handles);
    case 'RePlot'
        RePlot(hObject, eventdata, handles);
    otherwise
end


function RePlot(hObject, eventdata, handles)

if ~get(handles.bMinThresh,'Value')&&~get(handles.bMaxThresh,'Value')
    caxis auto
    return
end

v=caxis;
if get(handles.bMinThresh,'Value')
    MinThresh = eval(get(handles.MinThresh,'String'));
else MinThresh=v(1);
end

if get(handles.bMaxThresh,'Value')
    MaxThresh = eval(get(handles.MaxThresh,'String'));
else MaxThresh=v(2);
end

caxis([MinThresh MaxThresh]);
xlabel(handles.axes1,'distance [um]');
ylabel(handles.axes1,'distance [um]');
    
% 
% %h = imagesc(ISC,'XData',sImage.XData,'YData',sImage.YData);
% h = imagesc(ISC,'XData',sImage.XData,'YData',sImage.YData,...
%     'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles},'Parent',handles.axes1);
% daspect(handles.axes1,[1 1 1]);
% colormap(handles.axes1,pink);
% hc = colorbar('location','eastoutside');%,'peer',handles.axes1);
% xlabel(handles.axes1,'distance [um]');
% ylabel(handles.axes1,'distance [um]');

function ChangeFileName(hObject, eventdata, handles)
global gSaveImg;

file = [gSaveImg.path gSaveImg.file];
[file, path, filterindex] = uiputfile('Image*', 'Save File As',file);

% Added if statement in case user presses "Cancel"
% in dialog box
if (file),
    gSaveImg.path = path;
    gSaveImg.file = file;

    set(handles.strFileName,'String',gSaveImg.file);
end

function SaveImage(what,hObject, eventdata, handles)
global gSaveImg gScan gConfocal;

%get(handles.bFixVx,'Value') 

switch what
    case 'Save'
        bSave = true;
    case 'Edit'
        bSave = false;
end

if bSave 
    now = clock;
    date = [num2str(now(1)),'-',num2str(now(2)),'-',num2str(round(now(3)))];
    fullPath=fullfile('C:\Data\',date,'\');
    if ~exist(fullPath,'dir')
        mkdir(fullPath);
    end
gSaveImg.path = fullPath;

gSaveImg.file = ['Image_' date '.txt'];
set(handles.strFileName,'String',gSaveImg.file);
    ND = 4 - ( gScan.bFixVx + gScan.bFixVy + gScan.bFixVz + gScan.bFixDT );
else 
    ND = 2; 
end 

bEPS = get(handles.bSaveImgEPS,'Value');
bData = get(handles.bSaveImgData,'Value');

%Get image from plot
sAxes = get(handles.axes1);
sImage = get(sAxes.Children(end));
if ND == 1 
    Data = [ sImage.XData' sImage.YData' ];
else
    Data = sImage.CData;
end

%File name and prompt
file = [gSaveImg.path gSaveImg.file];
raw_filename = [gSaveImg.path gSaveImg.file];
bOverWrite = get(handles.bOverwrite,'Value');

if ~bOverWrite
    mfile = strrep(file,'.txt','*');
    mfilename = strrep(gSaveImg.file,'.txt','');
    A = ls(mfile);
    ImgN = 0;
    for f = 1:size(A,1)
        sImgN = sscanf(A(f,:),[mfilename '_Img%d.txt']);
        if ~isempty(sImgN)
            if sImgN > ImgN
                ImgN = sImgN;
            end
        end
    end
    ImgN = ImgN + 1;
    file = strrep(file,'.txt',sprintf('_Img%03d.txt',ImgN));
end

%Save File as Data
if bData & bSave
    %CurrentVz = DAQmxFunctionPool('ReadVoltage',);
    fid = fopen(file,'wt');
%     fprintf(fid,'VxRange: [%f %f] NVx:%d bFixVx:%d FixVx:%f\n',...
%         gScan.minVx,gScan.maxVx,gScan.NVx,gScan.bFixVx,gScan.FixVx);
%     fprintf(fid,'VyRange: [%f %f] NVy:%d bFixVy:%d FixVy:%f\n',...
%         gScan.minVy,gScan.maxVy,gScan.NVy,gScan.bFixVy,gScan.FixVy);
%     fprintf(fid,'VzRange: [%f %f] NVz:%d bFixVz:%d FixVz:%f\n',...
%         gScan.minVz,gScan.maxVz,gScan.NVz,gScan.bFixVz,gScan.FixVz);
%     fprintf(fid,'DTRange: [%f %f] NDT:%d bFixVDT:%d FixVDT:%f\n',...
%         gScan.minDT,gScan.maxDT,gScan.NDT,gScan.bFixDT,gScan.FixDT);
    fprintf(fid,'VxRange in distance [um]: [%f %f] NVx:%d bFixVx:%d FixVx:%f\n',...
        gScan.minVx/gConfocal.Vx_per_um,gScan.maxVx/gConfocal.Vx_per_um,gScan.NVx,gScan.bFixVx,gScan.FixVx);
    fprintf(fid,'VyRange in distance [um]: [%f %f] NVy:%d bFixVy:%d FixVy:%f\n',...
        gScan.minVy/gConfocal.Vy_per_um,gScan.maxVy/gConfocal.Vy_per_um,gScan.NVy,gScan.bFixVy,gScan.FixVy);
    fprintf(fid,'VzRange: [%f %f] NVz:%d bFixVz:%d FixVz:%f\n',...
        gScan.minVz,gScan.maxVz,gScan.NVz,gScan.bFixVz,gScan.FixVz);
    fprintf(fid,'DTRange: [%f %f] NDT:%d bFixVDT:%d FixVDT:%f\n',...
        gScan.minDT,gScan.maxDT,gScan.NDT,gScan.bFixDT,gScan.FixDT);

    fprintf(fid,'Size: [%.0f %.0f]\n',size(Data));
    fprintf(fid,'%.3f\t',Data');
    fprintf(fid,'\nNotes: %s',get(handles.Notes,'String'));
    fprintf(fid,'\nTime: %s',datestr(datenum(now), 'yyyy-mm-dd_HH-MM-SS'));
    fclose(fid);
end
%UpdateName in form
auxstr = strrep(file,gSaveImg.path,'');
set(handles.strUploadedFile,'String',auxstr);

global gDataPlot 

if get(handles.bExtPlot,'Value')
    bEPS = 1;
else
    bEPS = 0;
end

if bEPS
    file = strrep(file,'.txt','.eps');
    fig = figure(120);
    hca = gca;
    figure(120)
    imagesc(Data,'XData',sImage.XData,'YData',sImage.YData,'Parent',hca);
    daspect(hca,[1 1 1]);
    colormap(hca,pink);
    hc = colorbar('peer',hca,'location','EastOutside');
%     xlabel(hca,'V_x [volts]');
%     ylabel(hca,'V_y [volts]');
    if ~get(handles.bMinThresh,'Value')&&~get(handles.bMaxThresh,'Value')
        caxis auto
        return
    end
    
    v=caxis;
    if get(handles.bMinThresh,'Value')
        MinThresh = eval(get(handles.MinThresh,'String'));
    else MinThresh=v(1);
    end
    
    if get(handles.bMaxThresh,'Value')
        MaxThresh = eval(get(handles.MaxThresh,'String'));
    else MaxThresh=v(2);
    end
    caxis([MinThresh MaxThresh]);
    xlabel(hca,'distance [um]');
    ylabel(hca,'distance [um]');

%     title(hca,strcat(strrep(strrep(file,gSaveImg.path,''),'_',' '),' - ',...
%         get(handles.UploadedNotes,'String')));
    title(hca,raw_filename);
    
    %%%%%% Convert to Axes %%%%%%
    AxisPoints = size(Data);
    Nx = AxisPoints(1);
    Ny = AxisPoints(2);
    
    x = (sImage.XData(1) : (sImage.XData(2)-sImage.XData(1)) / (Nx-1) : sImage.XData(2));
    y = (sImage.YData(1) : (sImage.YData(2)-sImage.YData(1)) / (Ny-1) : sImage.YData(2));
    
    %figure(1)
    if length(x) > 2 && length(y) > 2
        CutData = Data;
%         title('Data saved as gDataPlot') 
    else
        if length(x) > 2
            CutData = (Data(:,2))';
            plot(x, CutData, '-b'), hold on, plot(x, CutData, '.r'), hold off
%             title('Data saved as gDataPlot') 
        else
            CutData = Data(2,:);
            plot(y, CutData, '-b'), hold on, plot(y, CutData, '.r'), hold off
%             title('Data saved as gDataPlot') 
        end
    end
    
    gDataPlot.X1 = x;
    gDataPlot.X2 = y;
    gDataPlot.Y = CutData;
    if bSave, saveas(gcf,file,'eps');
        close(fig);
    end
end

%Save File as EPS
% if bEPS
%     file = strrep(file,'.txt','.eps');
%     fig = figure;
%     hca = gca;
%     if ND==1
%         plot(hca,sImage.XData,sImage.YData,'r');
%         
%         
%         gDataPlot.X = sImage.XData; 
%         gDataPlot.Y = sImage.YData; 
%         
%         
%         %xlabel(sAxes.XLabel.);
%     else
%         h = imagesc(Data,'XData',sImage.XData,'YData',sImage.YData,'Parent',hca);
%         daspect(hca,[1 1 1]);
%         colormap(hca,pink);
%         hc = colorbar('peer',hca,'location','EastOutside');
%         xlabel(hca,'V_x [volts]');
%         ylabel(hca,'V_y [volts]');
%         title(hca,strcat(strrep(strrep(file,gSaveImg.path,''),'_',' '),' - ',...
%             get(handles.UploadedNotes,'String')));
%         
%         gDataPlot.X = sImage.XData; 
%         gDataPlot.Y = sImage.YData; 
%         gDataPlot.Z = Data; 
%         
%     end
%     if bSave, saveas(gcf,file,'eps'); 
%         close(fig);
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     ImageUpload                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageUpload(what, hObject, eventdata, handles)

switch what
    case 'Upload'
        Upload(hObject, eventdata, handles);
    case 'UploadNext'
        UploadNext(1,hObject, eventdata, handles);
    case 'UploadPrevious'
        UploadNext(-1,hObject, eventdata, handles);
    case 'Delete'
        DeleteImg(hObject, eventdata, handles);
    otherwise
end

function DeleteImg(hObject, eventdata, handles)
global gSaveImg;

path = [gSaveImg.path];
file = get(handles.strUploadedFile,'String');

filename = [path file];
filename2 = strrep(filename,'.txt','.eps');

Answer=input(['Do you want to delete ' filename '?(Yes:Y No:N)'],'s');
switch Answer
    case 'Y'
        delete(filename);
        delete(filename2);
    otherwise
end


function UploadNext(Inc,hObject, eventdata, handles)
global gSaveImg gUploadedImgPath;

if isempty(gUploadedImgPath)
    path = [gSaveImg.path];
else
    path = gUploadedImgPath;
end
file = get(handles.strUploadedFile,'String');

if Inc == 1
    file = NextImg(path,file);
else
    file = PreviousImg(path,file);
end

fid = fopen([path file]);
if fid >0
    fclose(fid);
    set(handles.strUploadedFile,'String',file);
    [IMG, Info]=ReadImageFile([path file]);

    PlotIMG(IMG,Info,hObject, eventdata, handles);
end

%To plot externally
if get(handles.bExtPlot,'Value')
    SaveImage('Edit',hObject, eventdata, handles);
end


function Upload(hObject, eventdata, handles)
global gSaveImg gUploadedImgPath;

file = [gSaveImg.path];
[file, path, filterindex] = uigetfile('Image*.txt', 'Upload Image',file);
gUploadedImgPath = path;

set(handles.strUploadedFile,'String',file);

[IMG, Info]=ReadImageFile([path file]);

PlotIMG(IMG,Info,hObject, eventdata, handles);

function PlotIMG(IMG,Info,hObject, eventdata, handles)
global gRawImg;

imagesc(IMG.Data,'XData',[IMG.ContRange(1) IMG.ContRange(2)],'YData',[IMG.FLRange(1) IMG.FLRange(2)],...
                'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles},'Parent',handles.axes1);

daspect(handles.axes1,[1 1 1]);
colormap(handles.axes1,pink);
RePlot(hObject,eventdata,handles);
colorbar('peer',handles.axes1,'location','EastOutside');
% 
% Info2 =char(Info{:});
% set(handles.strInfo,'String',Info2);
set(handles.UploadedNotes,'String',IMG.Notes);
gRawImg = IMG.Data;
xlabel(handles.axes1,strcat(IMG.LabelCont,' distance [um]'));
ylabel(handles.axes1,strcat(IMG.LabelFL,' distance [um]'));
set(handles.axes1,'ButtonDownFcn',{@CallImageSetXYVoltage,hObject,handles});

function [IMG, Info]=ReadImageFile(file)
%file
fid = fopen(file,'r');
%sline = fgetl(fid)
A = fscanf(fid,'VxRange in distance [um]: [%f %f] NVx:%d bFixVx:%d FixVx:%f');
Scan.minVx = A(1); Scan.maxVx = A(2); Scan.NVx = A(3); Scan.bFixVx = A(4); Scan.FixVx = A(5); 
fscanf(fid,'\n');
ScanRange{1}=[Scan.minVx Scan.maxVx];
A = fscanf(fid,'VyRange in distance [um]: [%f %f] NVy:%d bFixVy:%d FixVy:%f');
Scan.minVy = A(1); Scan.maxVy = A(2); Scan.NVy = A(3); Scan.bFixVy = A(4); Scan.FixVy = A(5); 
fscanf(fid,'\n');
ScanRange{2}=[Scan.minVy Scan.maxVy];
A = fscanf(fid,'VzRange: [%f %f] NVz:%d bFixVz:%d FixVz:%f');
Scan.minVz = A(1); Scan.maxVz = A(2); Scan.NVz = A(3); Scan.bFixVz = A(4); Scan.FixVz = A(5); 
fscanf(fid,'\n');
ScanRange{3}=[Scan.minVz Scan.maxVz];
A = fscanf(fid,'DTRange: [%f %f] NDT:%d bFixVDT:%d FixVDT:%f');
Scan.minDT = A(1); Scan.maxDT = A(2); Scan.NDT = A(3); Scan.bFixDT = A(4); Scan.FixDT = A(5); 
fscanf(fid,'\n');
ScanRange{4}=[Scan.minDT Scan.maxDT];
PS = [1 Scan.bFixVx; 2 Scan.bFixVy; 3 Scan.bFixVz; 4 Scan.bFixDT];
PS = sortrows(PS,[2 1]);
Label(1)='x';
Label(2)='y';
Label(3)='z';
Label(4)='t';
IMG.ContRange=ScanRange{PS(1,1)};
IMG.LabelCont=Label(PS(1,1));
IMG.FLRange=ScanRange{PS(2,1)};
IMG.LabelFL=Label(PS(2,1));


% IMG.minVx = A(1); IMG.maxVx = A(2); IMG.NVx = A(3); IMG.bFixVx = A(4); IMG.FixVx = A(5); 
% fscanf(fid,'\n');
% [A] = fscanf(fid,'VyRange in distance [um]: [%f %f] NVy:%d bFixVy:%d FixVy:%f');
% IMG.minVy = A(1); IMG.maxVy = A(2); IMG.NVy = A(3); IMG.bFixVy = A(4); IMG.FixVy = A(5); 
% fscanf(fid,'\n');
% [A] = fscanf(fid,'VzRange: [%f %f] NVz:%d bFixVz:%d FixVz:%f');
% IMG.minVz = A(1); IMG.maxVz = A(2); IMG.NVz = A(3); IMG.bFixVz = A(4); IMG.FixVz = A(5); 
% fscanf(fid,'\n');
% [A] = fscanf(fid,'DTRange: [%f %f] NDT:%d bFixVDT:%d FixVDT:%f');
% IMG.minDT = A(1); IMG.maxDT = A(2); IMG.NDT = A(3); IMG.bFixDT = A(4); IMG.FixDT = A(5); 
% fscanf(fid,'\n');
[A] = fscanf(fid,'Size: [%d %d]');
IMG.Nhor = A(1); IMG.Nvert = A(2);
fscanf(fid,'\n');
IMG.Data = fscanf(fid,'%f\t',[IMG.Nhor IMG.Nvert]);
IMG.Data = IMG.Data';
%IMG
fscanf(fid,'\nNotes: ');
IMG.Notes = fread(fid,'*char');
IMG.Notes = IMG.Notes';
fclose(fid);

% fid = fopen(file,'r');
% Info{1} = sprintf('Vx = [%+.3f %+.3f] N = %.0f IsFix = %.0f Fix = %+.3f',...
%     IMG.minVx,IMG.maxVx,IMG.NVx,IMG.bFixVx,IMG.FixVx);
% Info{2} = sprintf('Vy = [%+.3f %+.3f] N = %.0f IsFix = %.0f Fix = %+.3f',...
%     IMG.minVy,IMG.maxVy,IMG.NVy,IMG.bFixVy,IMG.FixVy);
% Info{3} = sprintf('Vz = [%+.3f %+.3f] N = %.0f IsFix = %.0f Fix = %+.3f',...
%     IMG.minVz,IMG.maxVz,IMG.NVz,IMG.bFixVz,IMG.FixVz);
% Info{4} = sprintf('DT = [%+.3f %+.3f] N = %.0f IsFix = %.0f Fix = %+.3f',...
%     IMG.minDT,IMG.maxDT,IMG.NDT,IMG.bFixDT,IMG.FixDT);
% fclose(fid);
Info=0;

function [NextFileName, NImgN ]= NextImg(path,file)
str2 = strrep(file,'.txt','');
mfilename = str2(1:length(str2)-7);
cImgN = eval(str2(length(str2)-2:length(str2)));
mfile = strcat(path,mfilename,'*');
A = ls(mfile);
ImgN = 1000;
for f = 1:size(A,1)
    sImgN = sscanf(A(f,:),[mfilename '_Img%d.txt']); 
    if ~isempty(sImgN)
        if sImgN > cImgN & sImgN < ImgN
            ImgN = sImgN;
        end
    end
end
if ImgN == 1000
    NextFileName = file;
    disp(['Next Image: same file : ' file]);
else
    NextFileName = strcat(mfilename,sprintf('_Img%03d.txt',ImgN));
    disp(['Next Image: ' NextFileName]);
end

function [PreviousFileName, NImgN ]= PreviousImg(path,file)
str2 = strrep(file,'.txt','');
mfilename = str2(1:length(str2)-7);
cImgN = eval(str2(length(str2)-2:length(str2)));
mfile = strcat(path,mfilename,'*');
A = ls(mfile);
ImgN = 1;
for f = 1:size(A,1)
    sImgN = sscanf(A(f,:),[mfilename '_Img%d.txt']);
    if ~isempty(sImgN)
        if sImgN < cImgN & sImgN > ImgN
            ImgN = sImgN;
        end
    end
end
if ImgN == 0
    PreviousFileName = file;
    disp(['Previous Image: same file : ' file]);
else
    PreviousFileName = strcat(mfilename,sprintf('_Img%03d.txt',ImgN));
    disp(['Previous Image: ' PreviousFileName]);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        ImageCPS                     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout=ImageCPS(what,hObject, eventdata, handles)
global gScan;

switch what
    case 'Run'
        %Turn On Laser : JM 2008-07-27
%        PBFunctionPool('PBON',2^5);
        varargout{1} = RunCPS(gScan,hObject, eventdata, handles);
    case 'RunOnce'
        varargout{1} = RunCPSOnce(gScan,hObject, eventdata, handles);
    case 'Stop'
        %Turn On Laser : JM 2008-07-27
%        PBFunctionPool('PBON',0);
        StopCPS;
    otherwise
end


function StopCPS
global bGo hCPS;

bGo = false;

hCounter = hCPS.hCounter;
hPulse = hCPS.hPulse ;
counter = '';
if handles.sr620Button.Value
    counter = 'SR620';
end

if ~strcmp(counter, 'SR620')
    DAQmxStopTask(hPulse);
    DAQmxStopTask(hCounter);

    DAQmxClearTask(hPulse);
    DAQmxClearTask(hCounter);
end


function CPS = RunCPS(Scan,hObject, eventdata, handles)
global bGo hCPS;

bGo = true;

NRead = 20;

DT = 0.01;
TimeOut = DT * NRead * 1.1;
Freq = 1/DT;

counter = '';
if handles.sr620Button.Value
    counter = 'SR620';
end

% status = DAQmxResetDevice('Dev1');
% disp(['NI: Reset Devie          :' num2str(status)]);
if ~strcmp(counter, 'SR620')
    try
    hCounter = SetCounter(NRead, Freq);
    [~, hPulse] = DigPulseTrainCont(Freq,0.5,NRead);

    hCPS.hCounter = hCounter;
    hCPS.hPulse = hPulse;

    strt=tic;
    time=toc(strt);
    while bGo && time < 3*60


        status = DAQmxStartTask(hCounter);  DAQmxErr(status);
        status = DAQmxStartTask(hPulse);    DAQmxErr(status);

        DAQmxWaitUntilTaskDone(hCounter,TimeOut);

        DAQmxStopTask(hPulse);

        A = ReadCounter(hCounter,NRead);
        DAQmxStopTask(hCounter);
        %axes(handles.axes1)
        A=ProcessDataVector(A,1);
        CPS = ProcessDataCPS(A, NRead, DT);

        set(handles.CPS,'String',CPS);
        drawnow;
        time=toc(strt);
    end

    DAQmxClearTask(hPulse);
    DAQmxClearTask(hCounter);
    catch ME
        KillAllTasks;
        rethrow(ME);
    end

else
    strt=tic;
    time=toc(strt);
    while bGo && time < 3*60
        
        CPS = measureSR620CountRate(DT);
        set(handles.CPS,'String',CPS);
        drawnow;
        time=toc(strt);
    end
end

function CPS = RunCPSOnce(Scan,hObject, eventdata, handles)
global bGo hCPS;

bGo = true;

NRead = 20;

DT = 0.01;
TimeOut = DT * NRead * 1.1;
Freq = 1/DT;

% status = DAQmxResetDevice('Dev1');
% disp(['NI: Reset Devie          :' num2str(status)]);
try
hCounter = SetCounter(NRead, Freq);
[~, hPulse] = DigPulseTrainCont(Freq,0.5,NRead);

hCPS.hCounter = hCounter;
hCPS.hPulse = hPulse;

    status = DAQmxStartTask(hCounter);  DAQmxErr(status);
    status = DAQmxStartTask(hPulse);    DAQmxErr(status);

    DAQmxWaitUntilTaskDone(hCounter,TimeOut);

    DAQmxStopTask(hPulse);

    A = ReadCounter(hCounter,NRead);
    DAQmxStopTask(hCounter);
    A=ProcessDataVector(A,1);
    CPS = ProcessDataCPS(A, NRead, DT);
    set(handles.CPS,'String',CPS);
    drawnow;


DAQmxClearTask(hPulse);
DAQmxClearTask(hCounter);
catch ME
    KillAllTasks;
    rethrow(ME);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    ImageSetXYVoltage      %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImageSetXYVoltage(what,hObject, eventdata, handles)
%global gScan;
switch what
    case 'Cursor'
        SetXYFromCursor(hObject, eventdata, handles);
    case 'CursorType'
        CursorType(hObject, eventdata, handles);
    case 'Arrow'
        ArrowFunction(hObject, eventdata, handles);
    case 'Fix'
        FixOutputs(hObject, eventdata, handles);
    case 'SetZero'
        SetZeroOutputs(hObject, eventdata, handles);
    otherwise
end

function ManChange(hObject, eventdata, handles) 
global gScan gbManChange gConfocal;

% p.x = eval(get(handles.FixVx,'String')); 
% p.y = eval(get(handles.FixVy,'String')); 
% p.z = eval(get(handles.FixVz,'String')); 
dx = (gScan.maxVx - gScan.minVx)/(gScan.NVx-1); 
dy = (gScan.maxVy - gScan.minVy)/(gScan.NVy-1); 
dz = (gScan.maxVz - gScan.minVz)/(gScan.NVz-1); 

gScan.FixVx = gScan.FixVx + gbManChange.X * dx; 
gScan.FixVy = gScan.FixVy + gbManChange.Y * dy; 
gScan.FixVz = gScan.FixVz + gbManChange.Z * dz; 

disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(gScan.FixVz) ')'])


WriteVoltage('Dev1/ao0',gScan.FixVx + gConfocal.XOffSet);
WriteVoltage('Dev1/ao1',gScan.FixVy + gConfocal.YOffSet);
WriteVoltage('Obj_Piezo',gScan.FixVz);

UpdateCoords(hObject, eventdata, handles) 

function GetCoords(hObject, eventdata, handles) 
global EO_handle gScan;
[err_code obj_piz_pos]=calllib('EO-Drive','EO_GetCommandPosition', EO_handle,0);
disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(obj_piz_pos) ')'])


function UpdateCoords(hObject, eventdata, handles) 
global gScan;

set(handles.FixVx,'String',num2str(gScan.FixVx));
set(handles.FixVy,'String',num2str(gScan.FixVy));
set(handles.FixVz,'String',num2str(gScan.FixVz));
disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(gScan.FixVz) ')'])


function FixOutputs(hObject, eventdata, handles)
global gScan gConfocal;
ImageScan('GetScanValues',hObject, eventdata, handles);


WriteVoltage('Dev1/ao0',gScan.FixVx + gConfocal.XOffSet);
WriteVoltage('Dev1/ao1',gScan.FixVy + gConfocal.YOffSet);    
WriteVoltage('Obj_Piezo',gScan.FixVz);
disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(gScan.FixVz) ')'])

temp_filename = 'CurrentNV_Position.txt';
temp_path = 'C:\MATLAB Scripts\MATLAB_Code\Sets\';
fid = fopen([temp_path temp_filename],'wt');
fprintf(fid,'X= %f\n', gScan.FixVx + gConfocal.XOffSet);
fprintf(fid,'Y= %f\n', gScan.FixVy + gConfocal.YOffSet);
fprintf(fid,'Z= %f\n', gScan.FixVz);
fclose(fid);

function SetZeroOutputs(hObject, eventdata, handles)
global gScan gConfocal;
gScan.FixVx = 0;
gScan.FixVy = 0;


WriteVoltage('Dev1/ao0',gScan.FixVx + gConfocal.XOffSet);
WriteVoltage('Dev1/ao1',gScan.FixVy + gConfocal.YOffSet);    
disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(gScan.FixVz) ')'])

%set(handles.FixVx,'String',num2str(gScan.FixVx));
%set(handles.FixVy,'String',num2str(gScan.FixVy));

temp_filename = 'CurrentNV_Position.txt';
temp_path = 'C:\MATLAB Scripts\MATLAB_Code\Sets\';
fid = fopen([temp_path temp_filename],'wt');
fprintf(fid,'X= %f\n', gScan.FixVx + gConfocal.XOffSet);
fprintf(fid,'Y= %f\n', gScan.FixVy + gConfocal.YOffSet);
fprintf(fid,'Z= %f\n', gScan.FixVz);
fclose(fid);

% % #EO
% global EO_handle;
% 
% if gScan.FixVz > 100
%     gScan.FixVz = 100;
% end
% if gScan.FixVz < 0
%     gScan.FixVz = 0;
% end
% 
% calllib('EO-Drive','EO_Move', EO_handle, double(gScan.FixVz)); 
% % #EO


function ArrowFunction(hObject, eventdata, handles)
global gScan;
%if ~get(handles.bCursorXY,'Value'), return; end

Object = get(hObject);
c= Object.CurrentCharacter;
cCode = unicode2native(c);
if isempty(cCode), return; end

switch cCode
    case 29
        disp('Move Right')
        MoveCursorAndSetVoltage('Right',hObject, eventdata, handles);
    case 28
        disp('Move Left')
        MoveCursorAndSetVoltage('Left',hObject, eventdata, handles);
    case 30
        disp('Move Up')
        MoveCursorAndSetVoltage('Up',hObject, eventdata, handles);
    case 31
        disp('Move Down')
        MoveCursorAndSetVoltage('Down',hObject, eventdata, handles);
    otherwise
end

function MoveCursorAndSetVoltage(where,hObject, eventdata, handles)
global gScan;

%CP = get(handles.axes1,'CurrentPoint');
p.x = eval(get(handles.FixVx,'String'));
p.y = eval(get(handles.FixVy,'String'));
p
dx = (gScan.maxVx - gScan.minVx)/(gScan.NVx-1);
dy = (gScan.maxVy - gScan.minVy)/(gScan.NVy-1);

switch where
    case 'Right'
        p.x = p.x + dx; %Right
    case 'Left'
        p.x = p.x - dx; %Left
    case 'Up'
        p.y = p.y - dy; %Up
    case 'Down'
        p.y = p.y + dy; %Down
    otherwise
end
p
SetVoltage(p,hObject, eventdata, handles)
return;

%set(0,'PointerLocation',[CP(1,1) CP(1,2)]) 
Axes = get(handles.axes1);
CP1 = get(handles.axes1,'CurrentPoint')
P1 = get(0,'PointerLocation');
P2(1)=P1(1) + 1;
P2(2)=P1(2) + 1;
set(0,'PointerLocation',P2);
pause(0.1)
CP2 = get(gca,'CurrentPoint')
Pos = get(0,'PointerLocation');
DCP = CP2-CP1
ScreenSize = get(0,'ScreenSize');
XRange = (Axes.XLim(2)-Axes.XLim(1))
 T = LocalAxis2Screen(handles.axes1);
OT = 1./T
% 
% CurrentPoint=Axes.CurrentPoint
 OuterPosition=Axes.OuterPosition;
 XRangeOP = OuterPosition(3)-OuterPosition(1)
 Position=Axes.Position;
 XRangeP = Position(3) - Position(1)
% Parent = get(Axes.Parent);
% ParentPosition=Parent.Position
% Parent2 = get(Parent.Parent);
% Parent2Position=Parent2.Position
%moveptr(handles.axes1,'init');
%p.x,p.y
%moveptr(handles.axes1,'move',p.x,p.y)

function CursorType(hObject, eventdata, handles)

axes(handles.axes1);
if get(handles.bCursorXY,'Value')
    set(gcf,'Pointer','crosshair');
else
    set(gcf,'Pointer','arrow');
end


function CallImageSetXYVoltageKey(scr,hObject,eventdata,handles)
disp('Display Key Press');

function Scan = SetXYFromCursor(hObject, eventdata, handles)
global gScan gConfocal;
CP = get(handles.axes1,'CurrentPoint');
p.x = CP(1,1)* gConfocal.Vx_per_um;
p.y = CP(1,2)* gConfocal.Vy_per_um;
if get(handles.bCursorXY,'Value')

%     DrawArrow(p,hObject, eventdata, handles);
    gScan.FixVx = p.x;
    if ~get(handles.bFixVy,'Value')
        gScan.FixVy = p.y;
        gScan.bFixVy = 0;
        gScan.bFixVz = 1;
    else
        gScan.FixVz = p.y;
        gScan.bFixVz = 0;
        gScan.bFixVy = 1;
    end
    
    ImageFillUpForm('UpdateScan',hObject, eventdata, handles);
    hX = WriteVoltage('Dev1/ao0',p.x + gConfocal.XOffSet);
    hY = WriteVoltage('Dev1/ao1',p.y + gConfocal.YOffSet);
    disp(['Galvo 1 Position (X, Y, Z) = (' num2str(gScan.FixVx) ', ' num2str(gScan.FixVy) ', ' num2str(gScan.FixVz) ')'])
end
% added 10 Aug 2008, jhodges
% draw crosshairs of the XY location
DrawCrossHairs(hObject,eventdata,handles)


function SetVoltage(p,hObject, eventdata, handles)
global gScan;


if get(handles.bCursorXY,'Value')
%     DrawArrow(p,hObject, eventdata, handles);
    gScan.FixVx = p.x;
    gScan.FixVy = p.y;
    ImageFillUpForm('UpdateScan',hObject, eventdata, handles);
    hX = WriteVoltage('Dev1/ao0',p.x);
    hY = WriteVoltage('Dev1/ao1',p.y);
end

% 
% function DrawArrow(p,hObject, eventdata, handles)
% sAxes = get(handles.axes1);
% sImage = get(sAxes.Children);
% s=0:1/49:1;
% one = ones(size(s));
% Ly = abs(sImage.YData(2) - sImage.YData(1));
% Lx = abs(sImage.XData(2) - sImage.XData(1));
% Lx,Ly
% HLine.X = (p.x-Lx/2) + (p.x+Lx/2)*s;
% HLine.Y = (p.y)*one;
% VLine.Y = (p.y-Ly/2) + (p.y+Ly/2)*s;
% VLine.X = (p.x)*one;
% hold on;
% hH = plot(HLine.X,HLine.Y,'r');
% hV = plot(VLine.X,VLine.Y,'r');
% hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = LocalAxis2Screen(ax)
% Axis to screen coordinate transformation.
%
%   T = AXIS2SCREEN(AX) computes a coordinate transformation 
%       T = [xo yo rx ry] 
%   that relates the normalized axes coordinates [xa;ya] of a 
%   given point to its screen coordinate [xs;ys] (in the root 
%   units) by
%       xs = xo + rx * xa
%       ys = yo + ry * ya
%
%   See also SISOTOOL.

% Get axes normalized position in figure
AxisPos = LocalGetNormPos(ax);

% Get figure's normalized position in screen
FigPos = LocalGetNormPos(get(ax,'Parent'));

% Transformation norm. axis coord -> norm. fig. coord.
T = AxisPos;

% Transformation norm. axis coord -> norm. screen coord.
T(1:2) = FigPos(1:2) + FigPos(3:4) .* T(1:2);
T(3:4) = FigPos(3:4) .* T(3:4);

% Transform to screen units
ScreenSize = get(0,'ScreenSize');
T(1:2) = ScreenSize(1:2) + ScreenSize(3:4) .* T(1:2);
T(3:4) = ScreenSize(3:4) .* T(3:4);


%%%%%%%%%%%%%%%%%%%
% LocalGetNormPos %
%%%%%%%%%%%%%%%%%%%
function Pos = LocalGetNormPos(Handle)
%   Determine if mag. plot is slipping out of focus.

Units = get(Handle,'Units');
if strcmp(Units,'normalized')
   Pos = get(Handle,'Position');
else
   set(Handle,'Units','normalized')
   Pos = get(Handle,'Position');
   set(Handle,'Units',Units)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%          ImageZoom             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageZoom(what,hObject, eventdata, handles)

switch what
    case 'Zoom On'
        axes(handles.axes1);
        if get(handles.bZoom,'Value')
            disp('Zoom')
            zoom on;
        else
            zoom off;
        end
    otherwise
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ShowMarker 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ShowMarker(hObject,eventdata,handles)

X = str2num(get(handles.FixVx,'String'));
Y = str2num(get(handles.FixVy,'String'));

dX = .05;
dY = .05;
bbox = dsxy2figxy(handles.axes1,[X-dX,Y-dY,X,Y]);
%annotation('ellipse',bbox);

function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from axis to figure coords
% Transforms [axx axy] or [xypos] from axes hAx (data) coords into coords
% wrt GCF for placing annotation objects that use figure coords into data
% space. The annotation objects this can be used for are
%    arrow, doublearrow, textarrow
%    ellipses (coordinates must be transformed to [x, y, width, height])
% Note that line, text, and rectangle anno objects already are placed
% on a plot using axes coordinates and must be located within an axes.
% Usage: Compute a position and apply to an annotation, e.g.,
%   [axx axy] = ginput(2);
%   [figx figy] = getaxannopos(gca, axx, axy);
%   har = annotation('textarrow',figx,figy);
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')'])

%% Obtain arguments (only limited argument checking is performed).
% Determine if axes handle is specified
if length(varargin{1})== 1 && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end);
else
	hAx = gca;
end;
% Parse either a position vector or two 2-D point tuples
if length(varargin)==1	% Must be a 4-element POS vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Two tuples (start & end points)
end
%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));
%% Transform data from figure space to data space
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
	varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
	pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
	pos(3) = pos(3)*axpos(3)/axwidth;
	pos(4) = pos(4)*axpos(4)/axheight;
	varargout{1} = pos;
end
%% Restore axes units
set(hAx,'Units',axun)


function ImageScreenCap(hObject, eventdata, handles)
screencapture(gcf,[390 160 405 380],'clipboard');  

function data= ProcessDataCPS(RawData, NRead, intTime)
global gmSEQ
if strcmp(gmSEQ.meas,'SPCM')
    data=sum(RawData)/(NRead * intTime);
elseif strcmp(gmSEQ.meas,'APD')
    data=mean(RawData);
end

function data= ProcessDataVector(RawData, dt)
global gmSEQ
if strcmp(gmSEQ.meas,'SPCM')
    data=diff(RawData);
    data=data/dt;
elseif strcmp(gmSEQ.meas,'APD')
    data=RawData(2:length(RawData));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%          NV Catalog             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The catalog is intended to store the relative positions of NVs. When
%adding a new NV to a preexisting catalog, update the position of any NV in
%the catalog, then add the new NV. If one wishes to change the position of
%an NV without affecting the others, remove that NV, update the position of
%another NV in the catalog, and then add the NV with its new position.


function NVCatalog_Add(hObject, eventdata, handles)
%Adds to or creates a catalog

global NVCat_path

%Get current NV position
currX=str2double(get(handles.FixVx,'String'));
currY=str2double(get(handles.FixVy,'String'));

%Get catalog number
catN=get(handles.Catalog_Number, 'String');

%Get NV number
nvN=str2double(get(handles.nv_nr, 'String'));

%Don't overwrite unless given permission
NVCatalog_Check_Names(hObject, eventdata, handles)

temp_filename = ['NV_Catalog_' catN '.txt'];
fid = fopen([NVCat_path temp_filename],'a'); %Open as appending 
fprintf(fid,'%d\t', nvN);
fprintf(fid,'%d\t', currX);
fprintf(fid,'%d\t', currY);
fclose(fid);

function NVCatalog_Update_Or_Remove(hObject, eventdata, handles, UpdateQ)
%When updating an NV position, translates all other NV positions as well.

global NVCat_path

%Get current NV position
currX=str2double(get(handles.FixVx,'String'));
currY=str2double(get(handles.FixVy,'String'));


%Get catalog number
catN=get(handles.Catalog_Number, 'String');

%Get NV number
nvN=get(handles.nv_nr, 'String');

temp_filename = ['NV_Catalog_' catN '.txt'];
fid = fopen([NVCat_path temp_filename],'r');
s=textscan(fid, '%s', 'delimiter', '\t'); %Read text into variable delimited by tabs
fclose(fid);

idx=find(strcmp(s{1}, nvN), 1,'first');

if isempty(idx)
    disp('No NV of that number cataloged. Abort.')
    return
else
    lenS=length(s{1});
    storeS=zeros(1,lenS);
        
    %Begin overwrite
    fid = fopen([NVCat_path temp_filename],'a');
    
    if UpdateQ
        %Find displacement vector
        disX=currX-storeS(idx+1);
        disY=currY-storeS(idx+2);
        
        for i=1:lenS %Store and update current values
            storeS(i)=str2double(s{1}(i));
            
            %first index is NV #, second is X position, third is Y
            if mod(i,3)==1
                fprintf(fid,'%d\t',storeS(i));
            elseif mod(i,3)==2
                fprintf(fid,'%d\t',storeS(i)+disX);
            else
                fprintf(fid,'%d\t',storeS(i)+disY);
            end
        end
    else
        nvN=str2double(nvN);
        i=1;
        while i<=lenS %Store and update current values
            storeS(i)=str2double(s{1}(i));
            
            %first index is NV #, second is X position, third is Y
            if mod(i,3)==1
                if storeS(i)==nvN
                    i=i+2;
                else
                    fprintf(fid,'%d\t',storeS(i));
                end
            elseif mod(i,3)==2
                fprintf(fid,'%d\t',storeS(i));
            else
                fprintf(fid,'%d\t',storeS(i));
            end
            i=i+1;
        end
    end
    fclose(fid);
    
end

function NVCatalog_Check_Names(hObject, eventdata, handles)
%Check whether NV is already in catalog

global NVCat_path

%Get catalog number
catN=get(handles.Catalog_Number, 'String');

%Get NV number
nvN=get(handles.nv_nr, 'String');

temp_filename = ['NV_Catalog_' catN '.txt'];
try
    fid = fopen([NVCat_path temp_filename],'r');
    s=textscan(fid, '%s', 'delimiter', '\t'); %Read text into variable delimited by line breaks
    fclose(fid);
    
    idx=find(strcmp(s{1}, nvN), 1,'first');
    
    if ~isempty(idx)
        NVCatalog_OverwriteQ()
        waitfor(1010)
        if ~yn
            disp('Abort')
            return
        end
    end
catch MExc
    getReport(MExc)
    disp('No Catalog Exists... Probably')
end

function NVCatalog_OverwriteQ()
% ask for permission to overwrite NV

% Create a figure and axes
f = figure(1010);

txt='NV already exists in catalog. Do you want to overwrite it?';

txtFig = uicontrol('Style', 'text', 'String', txt , 'Position', [140 180 300 50]) ;

% Create push buttons
btnY = uicontrol('Style', 'pushbutton', 'String', 'Yes',...
    'Position', [220 160 50 20],...
    'Callback', 'global yn; yn=1; close(1010)');

btnN = uicontrol('Style', 'pushbutton', 'String', 'No',...
    'Position', [320 160 50 20],...
    'Callback', 'global yn; yn=0; close(1010)');
