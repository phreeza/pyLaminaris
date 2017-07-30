function varargout = playRipple(varargin)
% PLAYRIPPLE MATLAB code for playRipple.fig
%      PLAYRIPPLE, by itself, creates a new PLAYRIPPLE or raises the existing
%      singleton*.
%
%      H = PLAYRIPPLE returns the handle to a new PLAYRIPPLE or the handle to
%      the existing singleton*.
%
%      PLAYRIPPLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYRIPPLE.M with the given input arguments.
%
%      PLAYRIPPLE('Property','Value',...) creates a new PLAYRIPPLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before playRipple_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to playRipple_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help playRipple

% Last Modified by GUIDE v2.5 05-Feb-2016 16:18:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @playRipple_OpeningFcn, ...
    'gui_OutputFcn',  @playRipple_OutputFcn, ...
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

function settings = do_init_param(handles,whichrig)

switch whichrig
    case 'Ji_desktop'
        settings.stim_path = 'C:\Users\Ji Liu\Dropbox\RippleStim';
    case 'Ji_laptop'
        settings.stim_path = 'C:\Users\Ji Liu\Dropbox\RippleStim';
    case 'Cheetah'
        settings.stim_path = 'C:\Documents and Settings\i suck\My Documents\Google Drive\Stim_cheetah';
    case 'Prairie'
        settings.stim_path = 'C:\GoogleDrive\Stim';
    case 'BScope'
        settings.stim_path = 'D:\Ji\Google Drive\Stim';
        
end

% settings.nCycle = 1; % cycle of continuous repeats
settings.nRepeat = 1; % repreat each stimulus for how many times

%% adapted from Paul Watkins

switch whichrig
    case 'Ji_desktop'
        settings.TDT.circuit_fpath = fullfile('C:\Users\Ji Liu\Dropbox\CheetahRig',...
            'ripple_test.rcx');
    case 'Ji_laptop'
        settings.TDT.circuit_fpath = fullfile('C:\Users\Ji Liu\Dropbox\CheetahRig',...
            'ripple_test.rcx');
    case 'Cheetah'
        %         settings.TDT.circuit_fpath_otherStim = fullfile('C:\Documents and Settings\i suck\My Documents\Google Drive\Code\CheetahRig',...
        %             'ripple_test12.rcx');
        % the above circuit is working!!!
        %         settings.TDT.circuit_fpath_longStim = fullfile('C:\Documents and Settings\i suck\My Documents\Google Drive\Code\CheetahRig',...
        %             'ripple_Escabi_test1.3.rcx');
        
        root = 'C:\Documents and Settings\i suck\My Documents\Google Drive\Code\CheetahRig';
        
        settings.TDT.circuit_fpath_otherStim_UserP = fullfile(root,...
            'otherStim_UserP(0.0)_test1.1.rcx');
        settings.TDT.circuit_fpath_otherStim_frameClk = fullfile(root,...
            'otherStim_frameClk_working_softTrg_1.4.rcx');
        
        settings.TDT.circuit_fpath_longStim_UserP = fullfile(root,...
            'ripple_Escabi_UserP(0.0)_test1.5.rcx');
        settings.TDT.circuit_fpath_longStim_frameClk = fullfile(root,...
            'ripple_Escabi_frameClk_test1.3.rcx');
        
        settings.TDT.circuit_fpath_otherStimWithITD_UserP = fullfile(root,...
            'otherStim_UserP(0.0)_ITD_v2.0.rcx');


    case 'Prairie'
        root = 'C:\GoogleDrive\Code\CheetahRig';
        settings.TDT.circuit_fpath_otherStim_UserP = fullfile(root,...
            'otherStim_UserP(0.0)_test1.1.rcx');
        settings.TDT.circuit_fpath_otherStim_frameClk = fullfile(root,...
            'otherStim_frameClk_working_softTrg_1.4.rcx');
        settings.TDT.circuit_fpath_longStim_UserP = fullfile(root,...
            'ripple_Escabi_UserP(0.0)_test1.5.rcx');
        settings.TDT.circuit_fpath_longStim_frameClk = fullfile(root,...
            'ripple_Escabi_frameClk_test1.3.rcx');
        
    case 'BScope'
        root = 'D:\Ji\Google Drive\Code\CheetahRig';
        settings.TDT.circuit_fpath_otherStim_UserP = fullfile(root,...
            'otherStim_UserP(0.0)_test1.1.rcx');
        settings.TDT.circuit_fpath_otherStim_frameClk = fullfile(root,...
            'otherStim_frameClk_working_softTrg_1.4.rcx');
        settings.TDT.circuit_fpath_longStim_UserP = fullfile(root,...
            'ripple_Escabi_UserP(0.0)_test1.5.rcx');
        settings.TDT.circuit_fpath_longStim_frameClk = fullfile(root,...
            'ripple_Escabi_frameClk_test1.3.rcx');
        
end

% specify mapping from channels to TDT circuit buffers
% settings.TDT.CHANNELS = {'play_1' 'play_2'};
switch whichrig
    case 'Cheetah'
        settings.TDT.nPA5 = 2; % we have
    case 'Prairie'
        settings.TDT.nPA5 = 1; % we have
    case 'BScope'
        settings.TDT.nPA5 = 1; % we have
end
settings.TDT.range = 9.9999; % +- 10V, but prevent clipping
settings.TDT.callback_period = 0.050; % in seconds, period for the buffer update timer
% specify mapping from trigger type to TDT mux select
%settings.stim_triggers = strvcat('User (P0.0)','Frame Clk','Line Clk');
% settings.TDT.TRIGGER_MUX_SELECT = [0 1 1];
% buffer size in circuit MUST be divisible by number of chunks.
% must have at least two chunks.
% settings.TDT.NCHUNKS = 4;
settings.TDT.sampling_rate = 195312.50;
settings.TDT.fs_n = 5; % second argument for RX6.LoadCOFsf, the first being the circuit file path, 5->200kHz, w/ realizable rate 195312.50Hz
settings.TDT.maxSPL = 80; % maximum sound level speaker is calibrated to
settings.TDT.atten = 0;  % default attenuation
settings.TDT.timer_name = 'playRippleTimer';


return

function do_init_gui(handles)

set(handles.stim_path,'string',handles.settings.stim_path);
set(handles.atten,'string',num2str(handles.settings.TDT.atten));
set(handles.nRepeat,'string',num2str(handles.settings.nRepeat));
findStim(handles);
setappdata(handles.stim_uselist,'n_stim_use',0); % when initialize, no TORC selected
setappdata(handles.playRipple,'armed',0);
return

function findStim(handles)
folder = get(handles.stim_path,'string');
files = dir(folder);
stim = cellstr('N/A'); idx = 0;
for i = 1:length(files)
    identifier = 'stim';
    if files(i).isdir && strncmp(identifier,files(i).name,length(identifier))
        idx = idx + 1;
        stim{idx} = files(i).name;
    end
end

if idx == 0
    % no stim folder is found
    set(handles.stim_idlelist,'string',stim);
    set(handles.stim_idlelist,'enable','off');
else
    set(handles.stim_idlelist,'value',1);
    set(handles.stim_idlelist,'enable','on');
    set(handles.stim_idlelist,'string',stim);
    setappdata(handles.stim_idlelist,'n_stim_idle',idx);
    set(handles.stim_uselist,'value',1);
    set(handles.stim_uselist,'string',cellstr('stim to use'));
    setappdata(handles.stim_idlelist,'n_stim_use',0);
end

return

function findWaveform(handles,subfolder,hObject)
% this function finds all the .waveform files under subfolder
allStim = getappdata(handles.playRipple,'allStim');
nStim = getappdata(handles.playRipple,'nStim');
folder = get(handles.stim_path,'string');
path = fullfile(folder,subfolder);
files = dir(path);
for i = 1:length(files)
    [~,~,ext] = fileparts(files(i).name);
    if strcmp('.waveform',ext)
        nStim = nStim + 1;
        allStim{nStim} = fullfile(path,...
            files(i).name);
    end
end
setappdata(handles.playRipple,'allStim',allStim);
setappdata(handles.playRipple,'nStim',nStim);
return

function findWaveform_mat(handles,subfolder)
% this function finds all the .mat files under subfolder
allStim = getappdata(handles.playRipple,'allStim');
nStim = getappdata(handles.playRipple,'nStim');
folder = get(handles.stim_path,'string');
path = fullfile(folder,subfolder);

if get(handles.sel_long_stim,'value') == 1
    setappdata(handles.playRipple,'root',folder);
    setappdata(handles.playRipple,'folder_name',subfolder)
    param_fpath = fullfile(path,sprintf('%s_param.mat',subfolder));
    hMat_param = matfile(param_fpath,'writable',true);
    setappdata(handles.playRipple,'hMat_param',hMat_param);    
    nStim = 1;
    allStim = cellstr(param_fpath);
    hMat = [];
    
end
if get(handles.sel_other,'value') == 1
    files = what(path);
%     hMat = cell(1);
    for i = 1:length(files.mat)
        fullname = fullfile(path,files.mat{i});
        h = matfile(fullname);
        if ~isempty(whos(h,'wave'))
            nStim = nStim + 1;
            allStim{nStim} = fullname;
%             hMat{nStim} = h;
        end
    end
end
setappdata(handles.playRipple,'allStim',allStim);
setappdata(handles.playRipple,'nStim',nStim);
setappdata(handles.playRipple,'hMat_stim',[]);
return

% --- Executes just before playRipple is made visible.
function playRipple_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to playRipple (see VARARGIN)

% Choose default command line output for playRipple
handles.output = hObject;
movegui(handles.playRipple,'center');
if isempty(varargin)
    fprintf('not enough input arguments')
else
    whichrig = varargin{1};
    handles.whichrig = whichrig;
    % Update handles structure
    switch whichrig
        case 'Ji_desktop'
            
        case 'Ji_laptop'
            
        case 'Cheetah'
            set(handles.channel_sel,'value',1);
            set(handles.cheetahMode,'value',1);
            set(handles.trgSrc,'value',1);
        case 'Prairie'
            set(handles.channel_sel,'value',2);
            set(handles.imagingMode,'value',1);
        case 'BScope'
            set(handles.channel_sel,'value',2);
            set(handles.imagingMode,'value',1);
    end
    
    handles.settings = do_init_param(handles,whichrig);
    
    do_init_gui(handles); % putting setting values into gui interface
    
    do_init_TDT(handles);
end
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes playRipple wait for user response (see UIRESUME)
% uiwait(handles.playRipple);


% --- Outputs from this function are returned to the command line.
function varargout = playRipple_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function stim_path_Callback(hObject, eventdata, handles)
% hObject    handle to stim_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_path as text
%        str2double(get(hObject,'String')) returns contents of stim_path as a double


% --- Executes during object creation, after setting all properties.
function stim_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in loadStimPath.
function loadStimPath_Callback(hObject, eventdata, handles)
% hObject    handle to loadStimPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.whichrig
    case 'Ji_desktop'
        start_path = 'C:\Users\Ji Liu\Dropbox\RippleStim';
    case 'Ji_laptop'
        start_path = 'C:\Users\Ji Liu\Dropbox\RippleStim';
    case 'Cheetah'
        start_path = handles.settings.stim_path;
    case 'Prairie'
        start_path = handles.settings.stim_path;
    case 'BScope'
        start_path = handles.settings.stim_path;
end
stim_path = uigetdir(start_path,...
    'Please Select a Stim Path');
if stim_path
    set(handles.stim_path,'string',stim_path);
    findStim(handles);
    guidata(hObject,handles);
end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in stim_uselist.
function stim_uselist_Callback(hObject, eventdata, handles)
% hObject    handle to stim_uselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stim_uselist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stim_uselist


% --- Executes during object creation, after setting all properties.
function stim_uselist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_uselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on stim_idlelist and none of its controls.
function stim_idlelist_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to stim_idlelist (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'm')
    if getappdata(handles.stim_idlelist,'n_stim_idle')>=1
        idx = get(handles.stim_idlelist,'value');
        subfolder = get(handles.stim_idlelist,'string');
        if ~iscell(subfolder)
            subfolder = cellstr(subfolder);
        end
        putStimtoUse(handles,subfolder{idx});
        %         removeStimfromIdle(handles,subfolder,idx);
    else
        % do nothing?
    end
end


function putStimtoUse(handles,subfolder)
uselist = get(handles.stim_uselist,'string');
n = getappdata(handles.stim_uselist,'n_stim_use');
if ~iscell(uselist)
    uselist = cellstr(uselist);
end
uselist{n+1} = subfolder;
set(handles.stim_uselist,'string',uselist);
setappdata(handles.stim_uselist,'n_stim_use',n+1);
return

function removeStimfromIdle(handles,subfolder,idx)
n = getappdata(handles.stim_idlelist,'n_stim_idle');
setappdata(handles.stim_idlelist,'n_stim_idle',n-1);
subfolder(idx) = [];
tmp = min(idx,size(subfolder,1));
set(handles.stim_idlelist,'value',tmp);
set(handles.stim_idlelist,'string',subfolder);
return

% --- Executes on button press in arm_disarm.
function arm_disarm_Callback(hObject, eventdata, handles)
% hObject    handle to arm_disarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~getappdata(handles.playRipple,'armed')
    succeeded = do_arm(handles,hObject);
    if succeeded
        setappdata(handles.playRipple,'armed',1);
        update_arm_disarm_toggle(handles);
    else
        
    end
else
    succeeded = do_disarm(handles,hObject);
    if succeeded
        setappdata(handles.playRipple,'armed',0);
        update_arm_disarm_toggle(handles);
        report_status(handles,'Disarmed!');
    end
end

function succeeded = loadWaveform_other(handles)
allStim = getappdata(handles.playRipple,'allStim');
nStim = getappdata(handles.playRipple,'nStim');
allStimWav = cell(nStim,1);nSamps = zeros(nStim,1);
nErr=0;
softGain = 10^(str2double(get(handles.softGain,'string'))/20);
fs = handles.settings.TDT.sampling_rate;
nCmp = round(0.2*fs);  % test if stim is shorter than 200ms, which is chosen
%to aid circuit status detection, then pad with zeros to increase length

skipStim = zeros(nStim,1);
for i = 1:nStim
    f = load(allStim{i});
    try
        %wave = f.wave;
        %wave = wave*handles.settings.TDT.range; % normalize waveform to span RX6 output range, which is 10volts
        
        if max(abs(f.wave))*softGain > 1 % handles.settings.TDT.range
            report_status(handles,'software gain causing signal amplitude over RX6 input range!')
            skipStim(i) = 1;
        end

        allStimWav{i,1} = 0; % changed 1/25/2016 by JL, we load each stim during stim presentation not now
%         nSamps(i) = length(wave);
        nSamps(i) = 0;
        if i == 1
            if ~isempty(f.param)
                param{1} = f.param;
            else
                param{1} = 0;
            end
        else
            if ~isempty(f.param)
                param{i} = f.param;
            else
                param{i} = 0;
            end
        end
    catch
        nErr = nErr+1;
        fprintf('\nfile format error:\n%s',allStim{i});
        allStim(i) = [];
        continue
    end
end
if nErr>=nStim
    succeeded = 0;
else
    setappdata(handles.playRipple,'allStimWav',allStimWav);
    setappdata(handles.playRipple,'skipStim',skipStim);
    setappdata(handles.playRipple,'nSamps',nSamps);
    setappdata(handles.playRipple,'nStim',nStim-nErr);
    setappdata(handles.playRipple,'allStim',allStim);
    setappdata(handles.playRipple,'param',param);
    succeeded = 1;
end

return

function succeeded = loadWaveform_DRC(handles)
allStim = getappdata(handles.playRipple,'allStim');
nStim = getappdata(handles.playRipple,'nStim');
allStimWav = cell(nStim,1);nSamps = zeros(nStim,1);
nErr=0;
softGain = 10^(str2double(get(handles.softGain,'string'))/20);
fs = handles.settings.TDT.sampling_rate;

for i = 1:nStim
    f = load(allStim{i});
    try
        wave = f.wave;
        wave = wave*handles.settings.TDT.range; % normalize waveform to span RX6 output range, which is 10volts
        
        if max(abs(wave))*softGain > handles.settings.TDT.range
            report_status(handles,'software gain causing signal amplitude over RX6 input range!')
        else
            wave = wave*softGain;
        end

        allStimWav{i,1} = wave;
        nSamps(i) = length(wave);
        if i == 1
            param = f.param;
        else
            param(i) = f.param;
        end
    catch
        nErr = nErr+1;
        fprintf('\nfile format error:\n%s',allStim{i});
        allStim(i) = [];
        continue
    end
end
if nErr>=nStim
    succeeded = 0;
else
    setappdata(handles.playRipple,'allStimWav',allStimWav);
    setappdata(handles.playRipple,'nSamps',nSamps);
    setappdata(handles.playRipple,'nStim',nStim-nErr);
    setappdata(handles.playRipple,'allStim',allStim);
    setappdata(handles.playRipple,'param',param);
    succeeded = 1;
end

return

function succeeded = loadWaveform_long(handles)
allStim = getappdata(handles.playRipple,'allStim');
nStim = getappdata(handles.playRipple,'nStim');
allStimWav = cell(nStim,1);nSamps = zeros(nStim,1);
nSeg = zeros(nStim,1); nSamps_seg = zeros(nStim,1);
nErr=0;
softGain = 10^(str2double(get(handles.softGain,'string'))/20);
for i = 1:nStim
    f = matfile(allStim{i});
    try

        if abs(f.wav_min) * softGain  > 1 || abs(f.wav_max) * softGain  > 1
            report_status(handles,'software gain causing signal amplitude over RX6 input range!')
            succeeded = 0;
            return
        else
            allStimWav = 0;
            nSeg(i) = f.nSeg;
            nSamps_seg(i) = f.nSamps_seg;
        end

    catch
        nErr = nErr+1;
        fprintf('\nfile format error:\n%s',allStim{i});
        allStim(i) = [];
        continue
    end
end
if nErr>=nStim
    succeeded = 0;
else
    setappdata(handles.playRipple,'allStimWav',allStimWav);
    if get(handles.sel_other,'value')
        setappdata(handles.playRipple,'nSamps',nSamps);
    else
        setappdata(handles.playRipple,'nSamps_seg',nSamps_seg);
        setappdata(handles.playRipple,'nSeg',nSeg);
    end
    setappdata(handles.playRipple,'nStim',nStim-nErr);
    setappdata(handles.playRipple,'allStim',allStim);
    succeeded = 1;
end

return


function succeeded = do_arm(handles,hObject)
if get(handles.sel_other,'value')  % if want to play other stim than Escabi's DMR, which is hugely long
    %% find waveform files
    % each entry in TORC uselist is a folder, so containing possibly multiple
    % stimuli
    n = getappdata(handles.stim_uselist,'n_stim_use');
    if n>=1
        subfolder = get(handles.stim_uselist,'string');
        if ~iscell(subfolder)
            subfolder = cellstr(subfolder);
        end
        setappdata(handles.playRipple,'allStim',cellstr('n/a'));
        setappdata(handles.playRipple,'nStim',0);
        for i = 1:n % find all .waveform files in the folder
            findWaveform_mat(handles,subfolder{i}); % findWaveform will put all waveform files
        end
        
    else
        report_status(handles,'No stim selected!');
        succeeded = 0;
        return
    end
    if getappdata(handles.playRipple,'nStim')<=0
        report_status(handles,'No waveform files found!');
        succeeded = 0;
        return
    else
        report_status(handles,'waveform files found!')
        
    end
    succeeded = loadWaveform_other(handles);
    if ~succeeded
        report_status(handles,'loading waveform failure!');
        return
    else
        report_status(handles,'succeeded loading waveform!');
        
    end
    %% randomize sequence of repetitins of all stimuli
    atten = eval(get(handles.atten,'string'));
    if ~isempty(find(atten<0))
        report_status(handles,'Negative attenuation not allowed!');
        succeeded = 0;
        return
    end
    nStim = getappdata(handles.playRipple,'nStim');
    nAtten = length(atten);
    nRepeat = str2double(get(handles.nRepeat,'string'));
    stimSeq = zeros(nStim*nAtten*nRepeat,2);
    k=1;
    for i = 1:nStim
        for j = 1:nAtten
            stimSeq(k:k+nRepeat-1,:) = repmat([i,atten(j)],nRepeat,1);
            k=k+nRepeat;
        end
    end
    rndIdx = randperm(nStim*nAtten*nRepeat);
    stimSeq = stimSeq(rndIdx,:); % randomize stim presentation
    setappdata(handles.playRipple,'stimSeq',stimSeq);
    setappdata(handles.playRipple,'nRepeat',nRepeat);
    setappdata(handles.playRipple,'nPresent',size(stimSeq,1)); % save number of all presentations of stimuli
    setappdata(handles.playRipple,'testMode',get(handles.testMode,'value'));
    setappdata(handles.playRipple,'cheetahMode',get(handles.cheetahMode,'value'));
    %% initialize TDT circuit
    report_status(handles,'Arming...');
    succeeded = do_start_TDT_hardware_otherStims(handles);
    if ~succeeded
        return
    end
    
    %% start timer
    % use a matlab timer instead of polling on the buffer pointer.
    % this allows use of the matlab console while the playback is occurring,
    % and also frees up some cycles on the CPU. the period must be low enough,
    % however, so that the buffer does not overrun. the TDT circuit contains a
    % hardware flag for buffer overrun, which results in a fatal error.
    % remove any previous timers created here and create a new one.
    
    t = timerfindall('Tag',handles.settings.TDT.timer_name);
    if ~isempty(t), stop(t); delete(t); end
    %     if get(handles.imagingMode,'value') % different timer_callback
    %         TDT_TIMER = timer('Period',handles.settings.TDT.callback_period,...
    %             'TasksToExecute',inf,...
    %             'TimerFcn',{@playRipple_timer_callback_imaging,handles},...
    %             'StartFcn',{@playRipple_timer_startfcn,handles},...
    %             'StopFcn',{@playRipple_timer_stopfcn,handles},...
    %             'ExecutionMode','fixedRate',...
    %             'BusyMode','drop',...
    %             'Tag',handles.settings.TDT.timer_name);
    %     else
    %         TDT_TIMER = timer('Period',handles.settings.TDT.callback_period,...
    %             'TasksToExecute',inf,...
    %             'TimerFcn',{@playRipple_timer_callback_ephys,handles},...
    %             'StartFcn',{@playRipple_timer_startfcn,handles},...
    %             'StopFcn',{@playRipple_timer_stopfcn,handles},...
    %             'ExecutionMode','fixedRate',...
    %             'BusyMode','drop',...
    %             'Tag',handles.settings.TDT.timer_name);
    %     end
    TDT_TIMER = timer('Period',handles.settings.TDT.callback_period,...
        'TasksToExecute',inf,...
        'TimerFcn',{@playRipple_timer_callback_workForAll,handles},...
        'StartFcn',{@playRipple_timer_startfcn,handles},...
        'StopFcn',{@playRipple_timer_stopfcn,handles},...
        'ExecutionMode','fixedRate',...
        'BusyMode','drop',...
        'Tag',handles.settings.TDT.timer_name);
    setappdata(handles.playRipple,'TDT_TIMER',TDT_TIMER);
    start(TDT_TIMER);
    %% now wait for ephus trigger
    
    %%
    
end
if get(handles.sel_long_stim,'value') %  if want to play long continuous stim
    n = getappdata(handles.stim_uselist,'n_stim_use');
    if n==1
        subfolder = get(handles.stim_uselist,'string');
        if ~iscell(subfolder)
            subfolder = cellstr(subfolder);
        end
        setappdata(handles.playRipple,'allStim',cellstr('n/a'));
        setappdata(handles.playRipple,'nStim',0);
        findWaveform_mat(handles,subfolder{1}); % findWaveform will put all waveform files
        
    else
        report_status(handles,'No stim selected! Or one stim at a time!');
        succeeded = 0;
        return
    end
    if getappdata(handles.playRipple,'nStim')<=0
        report_status(handles,'No waveform files found!');
        succeeded = 0;
        return
    else
        report_status(handles,'Waveforms found!')
    end
    succeeded = loadWaveform_long(handles);
    if ~succeeded
        report_status(handles,'loading waveform failure!');
        return
    else
        report_status(handles,'succeeded loading waveform!');
        
    end
    % no need to randomize presentation of stimuli, usually only one DMR
    % should be selected, so what's left is to start TDT system, w/o the
    % timer this time. Stim play out can be monitored through while loop
    % this time!!!
    atten = eval(get(handles.atten,'string'));
    if ~isempty(find(atten<0))
        report_status(handles,'Negative attenuation not allowed!');
        succeeded = 0;
        return
    end
    nStim = getappdata(handles.playRipple,'nStim');
    nAtten = length(atten);
    nRepeat = str2double(get(handles.nRepeat,'string'));
    stimSeq = zeros(nStim*nAtten*nRepeat,2);
    k=1;
    for i = 1:nStim
        for j = 1:nAtten
            stimSeq(k:k+nRepeat-1,:) = repmat([i,atten(j)],nRepeat,1);
            k=k+nRepeat;
        end
    end
    rndIdx = randperm(nStim*nAtten*nRepeat);
    stimSeq = stimSeq(rndIdx,:); % randomize stim presentation
    setappdata(handles.playRipple,'stimSeq',stimSeq);
    setappdata(handles.playRipple,'nRepeat',nRepeat);
    setappdata(handles.playRipple,'nPresent',size(stimSeq,1)); % save number of all presentations of stimuli
    setappdata(handles.playRipple,'testMode',get(handles.testMode,'value'));
    setappdata(handles.playRipple,'cheetahMode',get(handles.cheetahMode,'value'));
    
    %% initialize TDT circuit
    report_status(handles,'Arming...');
    if get(handles.sel_long_stim,'value')
        succeeded = do_start_TDT_hardware_rippleEscabi(handles);
    end

    if ~succeeded
        return
    end
    
    
    t = timerfindall('Tag',handles.settings.TDT.timer_name);
    if ~isempty(t), stop(t); delete(t); end
    
    
    if get(handles.sel_long_stim,'value')
        TDT_TIMER = timer('Period',handles.settings.TDT.callback_period,...
            'TasksToExecute',inf,...
            'TimerFcn',{@playRipple_timer_callback_workForAll,handles},...
            'StartFcn',{@playRipple_timer_startfcn,handles},...
            'StopFcn',{@playRipple_timer_stopfcn,handles},...
            'ExecutionMode','fixedRate',...
            'BusyMode','drop',...
            'Tag',handles.settings.TDT.timer_name);
    end
    setappdata(handles.playRipple,'TDT_TIMER',TDT_TIMER);
    start(TDT_TIMER);
end



return

function succeeded = do_disarm(handles,hObject)
succeeded = 1;
TDT_TIMER = getappdata(handles.playRipple,'TDT_TIMER');
if isvalid(TDT_TIMER)
    stop(TDT_TIMER);delete(TDT_TIMER);
end
setappdata(handles.playRipple,'armed',0);
return

function update_arm_disarm_toggle(handles)

if getappdata(handles.playRipple,'armed')
    set(handles.arm_disarm,'String','Disarm');
    set(handles.arm_disarm,'backgroundColor',[0 70/256 112/256]);
    %   set(handles.arm_disarm,'Value',1);
    set(handles.arm_disarm,'FontAngle','Normal');
    set(handles.arm_disarm,'FontName','Times New Roman');
    set(handles.arm_disarm,'foregroundColor',[1 1 1]);
    
else
    set(handles.arm_disarm,'String','Arm');
    set(handles.arm_disarm,'FontAngle','italic');
    set(handles.arm_disarm,'backgroundColor',[0 114/256 189/256]);
    set(handles.arm_disarm,'foregroundColor',[1 1 1]);
    set(handles.arm_disarm,'FontName','MS Sans Serif');
    %   set(handles.arm_disarm,'Value',0);
end

function succeeded = do_start_TDT_hardware_rippleEscabi(handles)
RX6 = getappdata(handles.playRipple,'RX6');
RX6.Halt;
RX6.ClearCOF;
trgSrcStr = get(handles.trgSrc,'String');
sel = get(handles.trgSrc,'value');
if ~iscell(trgSrcStr)
    trgSrcStr = cellstr(trgSrcStr);
end
trgSrc = trgSrcStr{sel};

if strncmp(trgSrc, 'User P(0.0)', length('User P(0.0)'))
    COF_path = handles.settings.TDT.circuit_fpath_longStim_UserP;
    useFrameClk = 0;
end
if strncmp(trgSrc, 'Frame Clk', length('Frame Clk'));
    COF_path = handles.settings.TDT.circuit_fpath_longStim_frameClk;
    useFrameClk = 1;
    nFrame = str2double(get(handles.cntFrames,'string'));
    setappdata(handles.playRipple,'nFrame',nFrame);
    setappdata(handles.playRipple,'nBaseFrame',nFrame);
end
setappdata(handles.playRipple,'useFrameClk',useFrameClk);
if ~exist(COF_path,'file')
    handles.errorMSG = 'TDT circuit file does not exist!';
    report_status(handles,handles.errorMSG);
    succeeded = 0;
    return
end
% RX6.LoadCOFsf(COF_path,handles.settings.TDT.fs_n);
switch handles.tdt_fs.Value
    case 1  % 100kHz fs
        n_fs = 4;
    case 2  % 200kHz fs
        n_fs = 5;
end
RX6.LoadCOFsf(COF_path,n_fs);
fprintf('circuit path: %s\n',COF_path);
% second argument in the above line(5, in this case) specifies sampling
% rate of RX6(see RX6 manual). Actual rate realizable by RX6 corrosponding
% to 200kHz is 195312.50Hz

%% load first stim waveform into MATLAB memory

succeeded = playRipple_set_TDT_tags_rippleEscabi(handles,useFrameClk); % 0: source0; 1 source1; set tags specific to each source
if ~succeeded
    return
end

% succeeded = playRipple_set_TDT_other_tags(handles);% set other tags that are same for both source
% if ~succeeded
%     return
% end
%% set attenuation for the first stimulus
succeeded = playRipple_set_atten(handles,1);
if ~succeeded
    return
end
RX6.Run;

return
%%
function succeeded = do_start_TDT_hardware_otherStims(handles)
RX6 = getappdata(handles.playRipple,'RX6');
RX6.Halt;
RX6.ClearCOF;
trgSrcStr = get(handles.trgSrc,'String');
sel = get(handles.trgSrc,'value');
if ~iscell(trgSrcStr)
    trgSrcStr = cellstr(trgSrcStr);
end
trgSrc = trgSrcStr{sel};

% select different circuit according to trigger source
if strncmp(trgSrc, 'User P(0.0)', length('User P(0.0)'))
    if ~get(handles.withITD,'value')
        COF_path = handles.settings.TDT.circuit_fpath_otherStim_UserP;
    else
        COF_path = handles.settings.TDT.circuit_fpath_otherStimWithITD_UserP;
    end
    useFrameClk = 0;
end

if strncmp(trgSrc, 'Frame Clk', length('Frame Clk'));
    COF_path = handles.settings.TDT.circuit_fpath_otherStim_frameClk;
    useFrameClk = 1;
    nFrame = str2double(get(handles.cntFrames,'string'));
    setappdata(handles.playRipple,'nFrame',nFrame);
end


setappdata(handles.playRipple,'useFrameClk',useFrameClk);

% load circuit 
if ~exist(COF_path,'file')
    handles.errorMSG = 'TDT circuit file does not exist!';
    report_status(handles,handles.errorMSG);
    succeeded = 0;
    return
end

% print circuit and load circuit into RX6
fprintf('circuit file: %s\n',COF_path);
switch handles.tdt_fs.Value
    case 1  % 100kHz fs
        n_fs = 4;
    case 2  % 200kHz fs
        n_fs = 5;
end
fs = RX6.GetSFreq;
fprintf('current fs: %fHz\n',fs);
setappdata(handles.playRipple,'curr_fs',fs);
RX6.LoadCOFsf(COF_path,n_fs);
if RX6.LoadCOFsf(COF_path,n_fs);
    fprintf('\nPASSED loading circuit\n');
else
    fprintf('\nFAILED loading circuit\n');
end

% second argument in the above line(5, in this case) specifies sampling
% rate of RX6(see RX6 manual). Actual rate realizable by RX6 corrosponding
% to 200kHz is 195312.50Hz

%% load first stim waveform into MATLAB memory
[wav,nSamps,ITD,gain2] = playRipple_load_next_stim(handles,1);
%% set TDT circuit tags for next stimulus
succeeded = playRipple_set_TDT_tags(handles,wav,nSamps,useFrameClk,ITD,gain2); % 0: source0; 1 source1; set tags specific to each source

if ~succeeded
    return
end
% succeeded = playRipple_set_TDT_other_tags(handles);% set other tags that are same for both source
% if ~succeeded
%     return
% end
%% set attenuation for the first stimulus
succeeded = playRipple_set_atten(handles,1);
if ~succeeded
    return
end
RX6.Run;
return

function succeeded = do_start_TDT_hardware_DRC(handles,hObject)
RX6 = getappdata(handles.playRipple,'RX6');
RX6.Halt;
RX6.ClearCOF;
trgSrcStr = get(handles.trgSrc,'String');
sel = get(handles.trgSrc,'value');
if ~iscell(trgSrcStr)
    trgSrcStr = cellstr(trgSrcStr);
end
trgSrc = trgSrcStr{sel};
if strncmp(trgSrc, 'UserP(0.0)', length('UserP(0.0)'))
    
    COF_path = handles.settings.TDT.circuit_fpath_DRC_UserP;
    useFrameClk = 0;
end
if strncmp(trgSrc, 'Frame Clk', length('Frame Clk'));
    
    COF_path = handles.settings.TDT.circuit_fpath_DRC_frameClk;
    useFrameClk = 1;

end
setappdata(handles.playRipple,'useFrameClk',useFrameClk);
if ~exist(COF_path,'file')
    handles.errorMSG = 'TDT circuit file does not exist!';
    report_status(handles,handles.errorMSG);
    succeeded = 0;
    return
end
fprintf('circuit file: %s\n',COF_path);
if RX6.LoadCOFsf(COF_path,handles.settings.TDT.fs_n)
    fprintf('\nPASSED loading circuit\n');
else
    fprintf('\nFAILED loading circuit\n');
end
% second argument in the above line(5, in this case) specifies sampling
% rate of RX6(see RX6 manual). Actual rate realizable by RX6 corrosponding
% to 200kHz is 195312.50Hz

%% load first stim waveform into MATLAB memory
%% set TDT circuit tags for next stimulus
succeeded = playRipple_set_TDT_tags_DRC(handles); % 0: source0; 1 source1; set tags specific to each source
if ~succeeded
    return
end
% succeeded = playRipple_set_TDT_other_tags(handles);% set other tags that are same for both source
% if ~succeeded
%     return
% end
%% set attenuation for the first stimulus
succeeded = playRipple_set_atten(handles,1);
if ~succeeded
    return
end
RX6.Run;
return

%% ADAPTED FROM PAUL WATKIN, auditoryStim GUI
function succeeded = do_init_TDT(handles)
succeeded = 1;
% iniate communication with RX6.
report_status(handles,'Connecting to RX6');

% Create an Active X control that can talk to a TDT RX6
TDTRX6 = actxcontrol('RPco.x', [5 5 0 0], handles.playRipple);

% Tell the Active X control to establish comm. link with the RX6
%if ~TDTRX6.ConnectRX6('USB', 1) % xxx - make parameter in handles.TDT

if ~TDTRX6.ConnectRX6('GB', 1)
    handles.errorMSG = 'FAILED connecting to RX6 through GB!';
    report_status(handles,handles.errorMSG);
    if ~TDTRX6.ConnectRX6('USB', 1)
        handles.errorMSG = 'FAILED connecting to RX6 through USB!';
        report_status(handles,handles.errorMSG);
        succeeded = 0;return
    else
        report_status(handles,'PASSED connecting to RX6 through USB!');
    end
else
    report_status(handles,'PASSED connecting to RX6 through GB!');
end


% handles.TDT.RX6 = TDTRX6;
setappdata(handles.playRipple,'RX6',TDTRX6);
% initialize the attenuators.
report_status(handles,'Connecting to PA5');

% Create Active X controls to talk to PA5s
%TDTPA5 = zeros(1,handles.TDT.NUMPA5);

clear TDTPA5
for i=1:handles.settings.TDT.nPA5
    TDTPA5(i) = actxcontrol('PA5.x',[5 5 0 0], handles.playRipple);
    
    % Tell the Active X control to establish comm. link with the PA5s
    %if ~TDTPA5(i).ConnectPA5('USB',i) % xxx - make parameter in handles.TDT
    if ~TDTPA5(i).ConnectPA5('GB',i)
        
        report_status(handles,sprintf('FAILED connecting to PA5 %d/%d through GB!',...
            i,handles.settings.TDT.nPA5));
        if ~TDTPA5(i).ConnectPA5('USB',i)
            handles.errorMSG = sprintf('FAILED connecting to PA5 %d/%d through USB!',...
                i,handles.settings.TDT.nPA5);
            report_status(handles,handles.errorMSG);
            
            succeeded = 0; return
        else
            report_status(handles,sprintf('PASSED connecting to PA5 %d/%d through USB!',...
                i,handles.settings.TDT.nPA5));
        end
        %     succeeded = 0; return
    else
        report_status(handles,sprintf('PASSED connecting to PA5 %d/%d through GB!',...
            i,handles.settings.TDT.nPA5));
    end
end

% handles.TDT.PA5 = TDTPA5;
setappdata(handles.playRipple,'PA5',TDTPA5);

% update hw information application data for callbacks
% setappdata(handles.playRipple,'TDT',handles.TDT);
report_status(handles,'TDT initiation complete!');


function report_status(handles,message)
set(handles.status_str,'string',message,'visible','on');
display(sprintf('\n%s\n',message));
drawnow;




% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.TORC1,'value')
    set(handles.stim_idlelist,'enable','on');
    set(handles.stim_uselist,'enable','on');
    set(handles.TORC2_list,'enable','off');
    set(handles.TORC2_uselist,'enable','off');
else
    set(handles.TORC2_list,'enable','on');
    set(handles.TORC2_uselist,'enable','on');
    set(handles.stim_idlelist,'enable','off');
    set(handles.stim_uselist,'enable','off');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over stim_uselist.
function stim_uselist_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to stim_uselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function putStimtoIdle(handles,subfolder)
idleList = get(handles.stim_idlelist,'string');
n = getappdata(handles.stim_idlelist,'n_stim_idle');
if ~iscell(idleList)
    idleList = cellstr(idleList);
end
idleList{n+1} = subfolder;
set(handles.stim_idlelist,'string',idleList);
setappdata(handles.stim_idlelist,'n_stim_idle',n+1);


return

function removeStimfromUse(handles,subfolder,idx)
n = getappdata(handles.stim_uselist,'n_stim_use');
setappdata(handles.stim_uselist,'n_stim_use',n-1);
subfolder(idx) = [];
tmp = min(idx,size(subfolder,1));
set(handles.stim_uselist,'value',tmp);
set(handles.stim_uselist,'string',subfolder);
return


% --- Executes on key press with focus on stim_uselist and none of its controls.
function stim_uselist_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to stim_uselist (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% if strcmp(eventdata.Key, 'n')
%     if getappdata(handles.stim_uselist,'n_stim_use')>=1
%         idx = get(handles.stim_uselist,'value');
%         subfolder = get(handles.stim_uselist,'string');
%         if ~iscell(subfolder)
%             subfolder = cellstr(subfolder);
%         end
%         putStimtoIdle(handles,subfolder{idx});
%         removeStimfromUse(handles,subfolder,idx);
%     else
%         % do nothing?
%     end
% end









function nCycle_Callback(hObject, eventdata, handles)
% hObject    handle to nCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nCycle as text
%        str2double(get(hObject,'String')) returns contents of nCycle as a double


% --- Executes during object creation, after setting all properties.
function nCycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nRepeat_Callback(hObject, eventdata, handles)
% hObject    handle to nRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nRepeat as text
%        str2double(get(hObject,'String')) returns contents of nRepeat as a double


% --- Executes during object creation, after setting all properties.
function nRepeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function atten_Callback(hObject, eventdata, handles)
% hObject    handle to atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atten as text
%        str2double(get(hObject,'String')) returns contents of atten as a double


% --- Executes during object creation, after setting all properties.
function atten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stim_idlelist.
function stim_idlelist_Callback(hObject, eventdata, handles)
% hObject    handle to stim_idlelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stim_idlelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stim_idlelist


% --- Executes during object creation, after setting all properties.
function stim_idlelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_idlelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function stim_idlelist_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to stim_idlelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in modeChange.
function modeChange_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in modeChange
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.testMode,'value')
    looptest(handles.playRipple);
    
end
if get(handles.cheetahMode,'value')
end


% --- Executes on key press with focus on arm_disarm and none of its controls.
function arm_disarm_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to arm_disarm (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toggle_arm_disarm.
function toggle_arm_disarm_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_arm_disarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~getappdata(handles.playRipple,'armed')
    succeeded = do_arm(handles,hObject);
    if succeeded
        setappdata(handles.playRipple,'armed',1);
        %         update_arm_disarm_toggle(handles);
        set(handles.toggle_arm_disarm,'String','Disarm');
    else
        set(handles.toggle_arm_disarm,'value',0);
    end
else
    succeeded = do_disarm(handles,hObject);
    if succeeded
        setappdata(handles.playRipple,'armed',0);
        %         update_arm_disarm_toggle(handles);
        report_status(handles,'Disarmed!');
        set(handles.toggle_arm_disarm,'String','Arm');
    else
        set(handles.toggle_arm_disarm,'value',1);
        
    end
end
% Hint: get(hObject,'Value') returns toggle state of toggle_arm_disarm


% --- Executes on button press in SoftTrg_1.
function SoftTrg_1_Callback(hObject, eventdata, handles)
% hObject    handle to SoftTrg_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RX6 = getappdata(handles.playRipple,'RX6');
RX6.SoftTrg(1);


% --- Executes when user attempts to close playRipple.
function playRipple_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to playRipple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isappdata(handles.playRipple,'RX6')
    RX6 = getappdata(handles.playRipple,'RX6');
    RX6.Halt;
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over stim_idlelist.
function stim_idlelist_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to stim_idlelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function softGain_Callback(hObject, eventdata, handles)
% hObject    handle to softGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of softGain as text
%        str2double(get(hObject,'String')) returns contents of softGain as a double


% --- Executes during object creation, after setting all properties.
function softGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to softGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trgSrc.
function trgSrc_Callback(hObject, eventdata, handles)
% hObject    handle to trgSrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trgSrc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trgSrc


% --- Executes during object creation, after setting all properties.
function trgSrc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trgSrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cntFrames_Callback(hObject, eventdata, handles)
% hObject    handle to cntFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cntFrames as text
%        str2double(get(hObject,'String')) returns contents of cntFrames as a double


% --- Executes during object creation, after setting all properties.
function cntFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cntFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SoftTrg_3.
function SoftTrg_3_Callback(hObject, eventdata, handles)
% hObject    handle to SoftTrg_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RX6 = getappdata(handles.playRipple,'RX6');
RX6.SoftTrg(3);


% --- Executes on selection change in channel_sel.
function channel_sel_Callback(hObject, eventdata, handles)
% hObject    handle to channel_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel_sel


% --- Executes during object creation, after setting all properties.
function channel_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in imagingMode.
function imagingMode_Callback(hObject, eventdata, handles)
% hObject    handle to imagingMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of imagingMode


function nBaseFrame_Callback(hObject, eventdata, handles)
% hObject    handle to nBaseFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nBaseFrame as text
%        str2double(get(hObject,'String')) returns contents of nBaseFrame as a double


% --- Executes during object creation, after setting all properties.
function nBaseFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nBaseFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loopTrg3.
function loopTrg3_Callback(hObject, eventdata, handles)
% hObject    handle to loopTrg3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nPresent = getappdata(handles.playRipple,'nPresent');
RX6 = getappdata(handles.playRipple,'RX6');
% nFrame = get(handles.cntFrames,'string'); nFrame = str2double(nFrame);
% nBaseFrame = get(handles.nBaseFrame,'string');  nBaseFrame = str2double(nBaseFrame);
prompt = {'loop how many times?','interval in between loops'};
dlg_title = ' ';
num_lines = 1;
def = {'500','0.5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    for i = 1:str2double(answer{1})
        RX6.SoftTrg(3);
        pause(str2double(answer{2}));
    end
end

% --- Executes on button press in clear_uselist.
function clear_uselist_Callback(hObject, eventdata, handles)
% hObject    handle to clear_uselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stim_uselist,'value',1);
set(handles.stim_uselist,'string',{[]});
setappdata(handles.stim_uselist,'n_stim_use',0);


% --- Executes on selection change in sel_trg_src.
function sel_trg_src_Callback(hObject, eventdata, handles)
% hObject    handle to sel_trg_src (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sel_trg_src contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_trg_src


% --- Executes during object creation, after setting all properties.
function sel_trg_src_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_trg_src (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sound_only.
function sound_only_Callback(hObject, eventdata, handles)
% hObject    handle to sound_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sound_only


% --- Executes on button press in withITD.
function withITD_Callback(hObject, eventdata, handles)
% hObject    handle to withITD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of withITD


% --- Executes on selection change in tdt_fs.
function tdt_fs_Callback(hObject, eventdata, handles)
% hObject    handle to tdt_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tdt_fs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tdt_fs


% --- Executes during object creation, after setting all properties.
function tdt_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tdt_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
