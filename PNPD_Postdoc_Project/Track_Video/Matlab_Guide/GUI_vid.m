%% Guide to synchronize behavior with the electrophysiological record.

% The code relies on the following functions:

% --> showFrameOnAxis_2.m --> Display input video frame on axis
%     Computer Vision Toolbox - Copyright 2004-2010 The MathWorks, Inc.

% --> load_open_ephys_data.m --> Open *.continuos from Openephys
%     https://github.com/open-ephys/analysis-tools


% by Vinicius Carvalho. Nucleo de Neurociencias NNC.
% email: vrcarva@gmail.com 

% by Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com

% Universidade Federal de Minas Gerais
% Started in:  02/2019
% Last update: 04/2020

%%
function varargout = GUI_vid(varargin)
% GUI_VID MATLAB code for GUI_vid.fig
%      GUI_VID, by itself, creates a new GUI_VID or raises the existing
%      singleton*.
%
%      H = GUI_VID returns the handle to a new GUI_VID or the handle to
%      the existing singleton*.
%
%      GUI_VID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_VID.M with the given input arguments.
%
%      GUI_VID('Property','Value',...) creates a new GUI_VID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_vid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_vid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_vid

% Last Modified by GUIDE v2.5 10-Apr-2020 00:00:03
global StopF %stop button flag
global uTScount%user timestamp(user selects onset and offset according to the video - i.e freezing) counter
global LFPvar %if .mat, the LFP data variable name
global LFPfiles
global uTSframe %user timestamps (frame)
global uTSsec   %user timestamps (video time)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_vid_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_vid_OutputFcn, ...
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


% --- Executes just before GUI_vid is made visible.
function GUI_vid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_vid (see VARARGIN)

% Choose default command line output for GUI_vid
handles.output = hObject;

handles.vidStr = '';
handles.LFPStr = '';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_vid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_vid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%defines sampling period of videoTimestamps - one video timestamp is
%captured every "periodo_vTS" frames
function h = periodo_vTS
h = 10; 

% --- Executes on button press in tbPlayPause.
function tbPlayPause_Callback(hObject, eventdata, handles)
% hObject    handle to tbPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbPlayPause
global StopF

%loops frame display
if get(hObject,'Value') %play
    StopF = 0;
    while (handles.currFrameNum<handles.totalFrames) && ~StopF %check if stop was pressed, and if the current frame is greater than the total number of frames
        tic
        handles.currFrameNum = 1+handles.currFrameNum;%updates number of frames
        currFrame = read(handles.videoH,handles.currFrameNum);
        showFrameOnAxis_2(handles.axes1, currFrame);
        axes(handles.axes1);
        imagesc(currFrame);%shows current frame on figure
        set(handles.sliderplay,'Value',handles.currFrameNum );%update slider
        set(handles.FrameStr,'String',sprintf('Frames: %.0f/%.0f',handles.currFrameNum,handles.totalFrames));%update counters
        set(handles.FrameStrSec,'String',sprintf('Time: %.1f/%.1f',handles.videoH.CurrentTime,handles.videoH.Duration));
        set (gca,'visible','off')
        
        %handles.currFrameNum
        if mod(handles.currFrameNum-1,periodo_vTS)== 0 && isfield(handles,'Ephys')
            updateLFPaxis_Callback(handles) %one video timestamp for every 10 video frames - updates the LFP plot
        end
        ttemp = toc;
        
        if ttemp < 1/handles.FR
            pause(((1/handles.FR)-ttemp)) 
            %java.lang.Thread.sleep(((1/handles.FR)-ttemp)*1000);
        end
        guidata(hObject, handles);
    end
else % pause video display
    StopF = 1;
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function load_vid_Callback(hObject, eventdata, handles)
% hObject    handle to load_vid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format long

[fileStr,pathStr] = uigetfile({'*.asf;*.asx;*.avi;*.m4v;*.mj2;*.mdl;*.mov;*.mp4;*.mpg;*.wmv',...
    'Video Files (*.mp4,*.asf,*.asx,*.avi,*.m4v,*.mj2,*.mdl,*.mov,*.mp4,*.mpg,*.wmv)';...
    '*.*',  'All Files (*.*)'}, ...
   'Select a Video File');
handles.vidStr = [pathStr fileStr];
set(handles.fileStr,'String',sprintf('Video:%s\nLFP:%s',handles.vidStr,handles.LFPStr));


%open video
handles.videoH = VideoReader([pathStr fileStr]);
handles.totalFrames = handles.videoH.NumberOfFrames; %get number of frames
handles.videoH = VideoReader([pathStr fileStr]);
handles.FR = handles.videoH.FrameRate;%get frame rate

%get first frame
currFrame = read(handles.videoH,1);
handles.currFrameNum = 1;

%show frame on axes
axes(handles.axes1)
%set(handles.axes1,'XLim',[1 handles.videoH.Width]);
%set(handles.axes1,'YLim',[1 handles.videoH.Height]);
showFrameOnAxis_2(handles.axes1, currFrame);

%update slider
ssteps = [1, 1]/(handles.totalFrames-1);
set(handles.sliderplay,'Max',handles.totalFrames);
set(handles.sliderplay,'Min',1);
set(handles.sliderplay,'Value',1);
set(handles.sliderplay,'SliderStep',ssteps);

%updates counter string
set(handles.FrameStr,'String',sprintf('Frames: %.0f/%.0f',handles.currFrameNum,handles.totalFrames));
set(handles.FrameStrSec,'String',sprintf('Time: %.1f/%.1f',0,handles.videoH.Duration));

guidata(hObject, handles);


%abre registro openEphys
% --------------------------------------------------------------------
function OpenEphys_Callback(hObject, eventdata, handles)
% hObject    handle to OpenEphys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LFPfiles
%select 1 open ephys channel
[fileStr,pathStr] = uigetfile({'*.continuous',...
    'Recording (*.continuous)';...
    '*.*',  'All Files (*.*)'}, ...
   'Select OpenEphys channel recording','MultiSelect','On');
disp('Loading...please wait');

if ~iscell(fileStr)
    fileStr = {fileStr};
end

[handles.Ephys.data, handles.Ephys.dataTime, Ephysinfo] = load_open_ephys_data([pathStr fileStr{1}]);%if >1 channels selected, shows the first
%handles.Ephys has electrophysiology data - recording (data), sample
%time points (dataTime), info, as well as event timestamps, given next...
%load events
[tsLabels,tsTimes, ~] = load_open_ephys_data([pathStr 'all_channels.events']);
assignin('base','Ephysinfo',Ephysinfo);% TODO - conferir se n�o precisa mesmo do campo info do open ephys 
handles.FsOrig = Ephysinfo.header.sampleRate;
%updates static text with filenames
handles.LFPStr = [pathStr fileStr{1}];
set(handles.fileStr,'String',sprintf('Video:%s\nLFP:%s',handles.vidStr,handles.LFPStr));
set(handles.popCanais,'String',cellfun(@(x) x(1:find(x=='.',1,'first')-1), fileStr,'UniformOutput',false));%atualiza menu popup com os canais selecionados

handles.Ephys.vTSlabels = tsLabels(tsLabels==0);
handles.Ephys.vTStime = tsTimes(tsLabels==0);
%
handles.Ephys.trialTSlabels = tsLabels(tsLabels~=0);
handles.Ephys.trialTStime = tsTimes(tsLabels~=0);
%trialTS = trial time stamps (TODO: conferir como ser�o os protocolos - ex: funciona ok para trigger de estimula��o, mas nao para os que usam ASSR)
%vTS = video Time stamps used for sync (label = 1) - every 10 frames, one vTS is saved

%filtra e decima o sinal
%interface
%abre caixa de dialogo para selecionar alguns parametros
strDialog = {'Fator de Decima��o',...
     'Filtro Passa-altas (cutoff 1 Hz) S/N',...
     'Filtro Notch (60 e 120 Hz) S/N',...
     'Filtro Passa-baixas (0-450 Hz) S/N'};
handles.sPreProc = inputdlg(strDialog,'Sele��o de Par�metros',[1 1 1 1],{'1','S','S','S'});%selected preprocessing 
decrate = str2double(handles.sPreProc{1});

[handles.Ephys.data handles.Ephys.dataTime] = filtraSinal(handles);
handles.Ephys.Fs = handles.FsOrig/decrate;
%LFP type and file
LFPfiles = fileStr;
handles.LFPtype = 'OpenEphys';

disp('Ok, done!');
%plot LFP
axes(handles.LFPaxes)
cla(handles.LFPaxes);
handles.LFP_plot = plot(handles.Ephys.dataTime,handles.Ephys.data,'Color','[0.019607843137255, 0.384313725490196, 0.631372549019608]');

box off

%plot trial timestamps
hold on
handles.trialTSplot = plot(0,0,'Color','[0.6350, 0.0780, 0.1840]');
hold off
%hides trial TS plots, if box is not checked
if ~get(handles.trialTSbox,'Value')
    set(handles.trialTSplot,'visible','off')
end

handles.LFP_plot.XDataSource = 'handles.Ephys.dataTime(idxs)';
handles.LFP_plot.YDataSource = 'handles.Ephys.data(idxs)';

%selecionar soh labels 1 - frames  - label 2 aparece como linha no plot
guidata(hObject, handles);

%processamento do sinal - decimacao e filtragem
function [dataOut, timeOut] = filtraSinal(handles)

decrate = str2double(handles.sPreProc{1});

if decrate > 1
    handles.Ephys.data = decimate(handles.Ephys.data,decrate);
    handles.Ephys.dataTime = handles.Ephys.dataTime(1:decrate:end);
    handles.Ephys.Fs = handles.FsOrig/decrate;
else
    handles.Ephys.Fs = handles.FsOrig;
end

%filters
if strcmpi(handles.sPreProc{2},'s')%highpass
    handles.Ephys.data = fun_myfilters(handles.Ephys.data,handles.Ephys.Fs ,[1 0],'iir',0);
end
if strcmpi(handles.sPreProc{3},'s') %notch
    Wo = 60/( handles.Ephys.Fs/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    
    Wo = 120/( handles.Ephys.Fs/2);  BW = Wo/35;
    [b2,a2] = iirnotch(Wo,BW);
    handles.Ephys.data = filtfilt(b,a,handles.Ephys.data);
    handles.Ephys.data = filtfilt(b2,a2,handles.Ephys.data);
end
if strcmpi(handles.sPreProc{4},'s')%lowpass
    handles.Ephys.data = fun_myfilters(handles.Ephys.data,handles.Ephys.Fs,[0 450],'iir',0);
end

dataOut = handles.Ephys.data;
timeOut = handles.Ephys.dataTime;


%abre registro .mat
% --------------------------------------------------------------------
function openMat_Callback(hObject, eventdata, handles)
% hObject    handle to openMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LFPfiles
global LFPvar

%Seleciona registro .mat
[fileStr,pathStr] = uigetfile({'*.mat'},'Select .mat recording');
%seleciona qual variavel eh o registro
varsMat = who('-file', [pathStr fileStr]);
sel3 = listdlg('PromptString','Selecione o nome da vari�vel:',...
    'SelectionMode','single',...
    'ListString',varsMat);
LFPvar = varsMat{sel3};
if find(strcmp(varsMat,'Fs'))%se o .mat tiver uma variavel Fs, essa eh a frequencia de amostragem
    %abre registro
    load([pathStr fileStr],LFPvar,'Fs');
else%caso contrario, usuario digita no console
    Fs = input('\nDigite a Frequ�ncia de Amostragem (em Hz):');
end

handles.Ephys.Fs = Fs;

LFPmat = eval(LFPvar);

%mostra 1o canal
if size(LFPmat,1)>1
    handles.Ephys.data = LFPmat(1,:);
end
%cria vetor de tempo
handles.Ephys.dataTime = [0:length(LFPmat)-1]/handles.Ephys.Fs;
handles.FsOrig = handles.Ephys.Fs;

%atualiza campos
handles.LFPStr = [pathStr fileStr];
set(handles.fileStr,'String',sprintf('Video:%s\nLFP:%s',handles.vidStr,handles.LFPStr));
%atualiza popup canais
set(handles.popCanais,'String',cellfun(@num2str,num2cell(1:size(LFPmat,1)),'UniformOutput',false));%atualiza menu popup com os canais selecionados

%filtra e decima o sinal
%interface
%abre caixa de dialogo para selecionar alguns parametros
strDialog = {'Fator de Decima��o',...
     'Filtro Passa-altas (cutoff 1 Hz) S/N',...
     'Filtro Notch (60 e 120 Hz) S/N',...
     'Filtro Passa-baixas (0-450 Hz) S/N',...
     'LFP offset (quantos segundos o LFP est� adiantado (>0) ou atrasado(<0) em rela��o ao v�deo)'};
handles.sPreProc = inputdlg(strDialog,'Sele��o de Par�metros',[1 1 1 1 1],{'1','S','S','S','0'});%selected preprocessing 
handles.sPreProc{5} = str2double(handles.sPreProc{5});
decrate = str2double(handles.sPreProc{1});
%LFP preprocessing - filter, decimate
[handles.Ephys.data handles.Ephys.dataTime] = filtraSinal(handles);
handles.Ephys.Fs = handles.FsOrig/decrate;%new sampling frequency
handles.LFPtype = 'mat';

%plot LFP
axes(handles.LFPaxes)
cla(handles.LFPaxes);
handles.LFP_plot = plot(handles.Ephys.dataTime,handles.Ephys.data);
%plot trial timestamps
hold on
handles.trialTSplot = plot(0,0,'r.');
hold off

handles.LFP_plot.XDataSource = 'handles.Ephys.dataTime(idxs)';
handles.LFP_plot.YDataSource = 'handles.Ephys.data(idxs)';

guidata(hObject, handles);

% --- Executes on slider movement.
function sliderplay_Callback(hObject, eventdata, handles)
% hObject    handle to sliderplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the frame selected with the slider
handles.currFrameNum = round(get(hObject,'Value'));
%update axes with selected frame
currFrame = read(handles.videoH,handles.currFrameNum);
showFrameOnAxis_2(handles.axes1, currFrame);
%update counter string
set(handles.FrameStr,'String',sprintf('Frames: %.0f/%.0f',handles.currFrameNum,handles.totalFrames));
set(handles.FrameStrSec,'String',sprintf('Time: %.1f/%.1f',handles.videoH.CurrentTime,handles.videoH.Duration));

%updates LFP axes
if isfield(handles,'Ephys')
    updateLFPaxis_Callback(handles)
end

guidata(hObject, handles);

% --- Executes on button press in pbTS.
function pbTS_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to pbTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global uTScount
global uTSframe
global uTSsec
%Add selected timestamp
uTScount = uTScount+1;

if length(get(handles.list1,'string')) == length(get(handles.list2,'string'))%onset TS
    %get current frame and adds to TS variable (frame and seconds)
    uTSframe(ceil(uTScount/2),:) = [handles.currFrameNum 0];
    uTSsec(ceil(uTScount/2),:) = [handles.videoH.CurrentTime 0];
    
    %get current list and adds new TS to the bottom
    a = get(handles.list1,'string');
    temp = char(a);
    c = strvcat(temp,sprintf('%.0f',handles.currFrameNum));
    set(handles.list1,'string',cellstr(c));

else
    %offset TS
    %get current frame and adds to TS variable (frame and seconds)
    uTSframe(ceil(uTScount/2),2) = handles.currFrameNum;
    uTSsec(ceil(uTScount/2),2) = handles.videoH.CurrentTime;
    
    %get current list and adds new TS to the bottom
    a = get(handles.list2,'string');
    temp = char(a);
    c = strvcat(temp,sprintf('%.0f',handles.currFrameNum));
    set(handles.list2,'string',cellstr(c));     
end

guidata(hObject, handles);

% --- Executes on button press in pbTSdel.
function pbTSdel_Callback(hObject, eventdata, handles)
% hObject    handle to pbTSdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global uTScount
global uTSframe
global uTSsec
%Deletes previous timestamp

if uTScount >= 1
    if length(get(handles.list1,'string')) ~= length(get(handles.list2,'string'))%onset TS
        %deletes last TS recorded on variables TS and TSsec
        uTSframe(ceil(uTScount/2),:) = [];
        uTSsec(ceil(uTScount/2),:) = [];
        %get current list and deletes last TS     
        a = get(handles.list1,'string');
        a(end) = [];
        set(handles.list1,'string',a);
    else%offset TS
        %deletes last TS recorded on variables TS and TSsec
        uTSframe(ceil(uTScount/2),2) = NaN;
        uTSsec(ceil(uTScount/2),2) = NaN;
        %get current list and deletes last TS     
        a = get(handles.list2,'string');
        a(end) = [];
        set(handles.list2,'string',a);
    end
end

%decrease number of selected TS
uTScount = uTScount-1;
guidata(hObject, handles);

% --------------------------------------------------------------------
function TSexport_Callback(hObject, eventdata, handles)
% hObject    handle to TSexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global uTSframe
global uTSsec

uTS_LFPsec = zeros(size(uTSframe));
uTS_LFPidx = zeros(size(uTSframe));

if isfield(handles,'Ephys')%user selected timestamps relative to LFP - if OpenEphys video Timestamps, every 10 frames = 1 timestamp...
    %...(video and LFP share the same time, but not the same indexes)
    if strcmpi(handles.LFPtype, 'OpenEphys') 
        for ii = 1:size(uTSframe,1)
            videoTStemp = handles.Ephys.vTStime(ceil(uTSframe(ii,1)/periodo_vTS));%timestamp of selected frame
            uTS_LFPsec(ii,1) = videoTStemp;
            uTS_LFPidx(ii,1) = find(handles.Ephys.dataTime>=videoTStemp,1,'first');
            videoTStemp2 = handles.Ephys.vTStime(ceil(uTSframe(ii,2)/periodo_vTS));%timestamp of selected frame
            uTS_LFPsec(ii,2) = videoTStemp2;
            uTS_LFPidx(ii,2) = find(handles.Ephys.dataTime>=videoTStemp2,1,'first');
        end
    else
        %registro .mat
        for ii = 1:size(uTSframe,1)
            uTS_LFPsec(ii,1) = uTSsec(ii,1) + handles.sPreProc{5};%timestamp of selected frame
            uTS_LFPidx(ii,1) = find(handles.Ephys.dataTime>=uTS_LFPsec(ii,1),1,'first');
            uTS_LFPsec(ii,2) = uTSsec(ii,2) + handles.sPreProc{5};
            uTS_LFPidx(ii,2) = find(handles.Ephys.dataTime>=uTS_LFPsec(ii,2),1,'first');    
        end       
            
    end
    assignin('base','TS_LFPsec', uTS_LFPsec);
    assignin('base','TS_LFPindex', uTS_LFPidx);
       
end

assignin('base','TSframes', uTSframe);
assignin('base','TSseconds', uTSsec);


function updateLFPaxis_Callback(handles)

%gets axes scales from textboxes
xScale = get(handles.xaxisScale,'String');
yScale = get(handles.yaxisScale,'String');


if strcmpi(handles.LFPtype,'OpenEphys') || strcmpi(handles.LFPtype,'TDT') %open Ephys recording with videoTimestamps
    %current frame corresponds to which LFP time point?
    frameTime = handles.Ephys.vTStime(ceil(handles.currFrameNum/periodo_vTS));
else %.mat file
    frameTime = handles.videoH.CurrentTime + handles.sPreProc{5}; 
end

xScale = str2double(xScale);
xAxisLim = [frameTime-xScale/2 frameTime+xScale/2];%x axis limits
yScale = str2double(yScale);
yAxisLim = [-yScale +yScale];

idxs = handles.Ephys.dataTime>xAxisLim(1) & handles.Ephys.dataTime<xAxisLim(2);

%axes(handles.LFPaxes);
%set(handles.LFP_plot, 'XData', handles.Ephys.dataTime(idxs), 'YData', handles.Ephys.data(idxs));
refreshdata(handles.LFP_plot,'caller')
%drawnow

%set LFP plot
%set(handles.LFP_plot, 'XData', handles.Ephys.dataTime(idxs), 'YData', handles.Ephys.data(idxs));

if get(handles.trialTSbox,'Value') && strcmpi(handles.LFPtype,'OpenEphys')
%set timestamps plot
    tempTSplot = handles.Ephys.trialTStime(handles.Ephys.trialTStime>xAxisLim(1) & handles.Ephys.trialTStime<xAxisLim(2));
    set(handles.trialTSplot, 'XData',tempTSplot , 'YData',0*tempTSplot,'Marker','.');
end

xlim(handles.LFPaxes,xAxisLim);
ylim(handles.LFPaxes,yAxisLim)

%axis(handles.LFPaxes,[xAxisLim yAxisLim]);

guidata(handles.LFPaxes, handles);

function xaxisScale_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxisScale as text
%        str2double(get(hObject,'String')) returns contents of xaxisScale as a double
if isfield(handles,'Ephys')
updateLFPaxis_Callback(handles)
end
guidata(hObject, handles);


function yaxisScale_Callback(hObject, eventdata, handles)
% hObject    handle to yaxisScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaxisScale as text
%        str2double(get(hObject,'String')) returns contents of yaxisScale as a double
if isfield(handles,'Ephys')
    updateLFPaxis_Callback(handles)
end
guidata(hObject, handles);


% --- Executes on button press in trialTSbox.
function trialTSbox_Callback(hObject, eventdata, handles)
% hObject    handle to trialTSbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%hides trial TS plots, if box is not checked
if ~get(hObject,'Value')
    set(handles.trialTSplot,'visible','off')
else
    set(handles.trialTSplot,'visible','on')
end
guidata(hObject, handles);


% --- Executes on selection change in popCanais.
function popCanais_Callback(hObject, eventdata, handles)
% hObject    handle to popCanais (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popCanais contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popCanais
global LFPfiles % selected LFP file(s)
global LFPvar %if .mat, the LFP data variable name
selCh = get(hObject,'Value');

switch lower(handles.LFPtype)
    case 'openephys'
    pathStr = handles.LFPStr(1:find(handles.LFPStr=='\',1,'last'));
    
    handles.LFPStr = [pathStr LFPfiles{selCh}];
    [handles.Ephys.data, handles.Ephys.dataTime, ~] = load_open_ephys_data([pathStr LFPfiles{selCh}]);%if >1 channels selected, shows the first
    set(handles.fileStr,'String',sprintf('Video:%s\nLFP:%s',handles.vidStr,handles.LFPStr));
    %filtra o sinal
    decrate = str2double(handles.sPreProc{1});
    handles.Ephys.Fs = handles.FsOrig/decrate;
    
    case 'tdt'
        handles.Ephys.data = evalin('base', sprintf('tdtLFPdata(%d,:)',selCh)); 
    case 'mat' %TODO use base workspace for other channels?
    load(handles.LFPStr);
    tempMat = eval(LFPvar);
    handles.Ephys.data = tempMat(selCh,:);        
end

[handles.Ephys.data, handles.Ephys.dataTime] = filtraSinal(handles);

updateLFPaxis_Callback(handles)

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popCanais_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCanais (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function sliderplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function FrameStr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function FrameTotStr_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to FrameTotStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function FrameStrSec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameStrSec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in list1.
function list1_Callback(hObject, eventdata, handles)  %list with onset timestamps 
% hObject    handle to list1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list1

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to TSexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in list2.
function list2_Callback(hObject, eventdata, handles) %#ok<*DEFNU> %list with offset timestamps 
% hObject    handle to list2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list2

% --- Executes during object creation, after setting all properties.
function list2_CreateFcn(hObject, eventdata, handles) %list with offset timestamps 
% hObject    handle to list2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function xaxisScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxisScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function yaxisScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxisScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function list1_CreateFcn(hObject, eventdata, handles)  %list with onset timestamps 
% hObject    handle to list1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global uTScount
global uTSframe
global uTSsec
%create timestamp variables
uTSframe = zeros(1,2);%user defined timestamps, in frames
uTSsec = zeros(1,2);%user defined timestamps, in seconds (relative do VIDEO)
uTScount = 0;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tbPlayPause_CreateFcn(hObject, ~, ~)
% hObject    handle to tbPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'cdata',imread('PP.png'))


% --- Executes on button press in pbDebug.
function pbDebug_Callback(hObject, eventdata, handles)
% hObject    handle to pbDebug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','handles',handles)

guidata(hObject, handles);

% --------------------------------------------------------------------
function openTDT_Callback(hObject, eventdata, handles)
% hObject    handle to openTDT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%opens TDT recording + video trigger
teste = figure();
TDT_header = HeaderBloco();%select tank-block and read TDT header

%read video trigger
teste = figure();
[ TDT_videoOut ] = LerTrigger2( TDT_header,'Select video trigger (default VIDF)' );
if ~isempty(TDT_videoOut)
    handles.Ephys.vTStime = TDT_videoOut.time(TDT_videoOut.sort==0);
    %vTS = video Time stamps used for sync (label = 1) - every 10 frames, one vTS is saved
end

%read stream - LFP
teste = figure();
[ TDT_streamOut ] = LerStream( TDT_header,0,0 );

handles.FsOrig = TDT_streamOut{1}.FreqAmos;%sampling frequency
 
Ncanais = TDT_streamOut{1}.nChannels;
tdtLFPdata = double(TDT_streamOut{1}.Dados)';
assignin('base','tdtLFPdata',tdtLFPdata); %all channel data to base workspace
handles.Ephys.data = tdtLFPdata(1,:);
handles.Ephys.dataTime = [0:size(tdtLFPdata,2)-1]/handles.FsOrig;

handles.LFPStr = sprintf('Tank %s , Block %s',TDT_header.TANK_in,TDT_header.BLOCK_in);
set(handles.fileStr,'String',sprintf('Video:%s\nLFP:%s',handles.vidStr,handles.LFPStr));
set(handles.popCanais,'String',cellfun(@num2str,num2cell(1:Ncanais),'UniformOutput',false));%atualiza menu popup com os canais selecionados

%TODO load (trial/stim) trigger labels

%filtra e decima o sinal
%interface
%abre caixa de dialogo para selecionar alguns parametros
strDialog = {'Fator de Decima��o',...
     'Filtro Passa-altas (cutoff 1 Hz) S/N',...
     'Filtro Notch (60 e 120 Hz) S/N',...
     'Filtro Passa-baixas (0-450 Hz) S/N'};
handles.sPreProc = inputdlg(strDialog,'Sele��o de Par�metros',[1 1 1 1],{'1','S','S','S'});%selected preprocessing 
decrate = str2double(handles.sPreProc{1});

[handles.Ephys.data handles.Ephys.dataTime] = filtraSinal(handles);
handles.Ephys.Fs = handles.FsOrig/decrate;
%LFP type and file
handles.LFPtype = 'TDT';

%plot LFP
axes(handles.LFPaxes)
cla(handles.LFPaxes);
handles.LFP_plot = plot(handles.Ephys.dataTime,handles.Ephys.data);
%plot trial timestamps
hold on
handles.trialTSplot = plot(0,0,'r.');
hold off

handles.LFP_plot.XDataSource = 'handles.Ephys.dataTime(idxs)';
handles.LFP_plot.YDataSource = 'handles.Ephys.data(idxs)';

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fileStr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in mice.
function mice_Callback(hObject, eventdata, handles)
% hObject    handle to mice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'cdata',imread('mouse.png'))

% Hint: get(hObject,'Value') returns toggle state of mice


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------

%% last update 10/04/2020 - 01:43am
%  listening: Early Day Miners - Sterling Provisions