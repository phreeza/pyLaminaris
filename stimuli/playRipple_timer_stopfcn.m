function playRipple_timer_stopfcn(obj,event,handles)

RX6 = getappdata(handles.playRipple,'RX6');
RX6.Halt;
playRipple_reset_atten(handles);
hMat = getappdata(handles.playRipple,'hMat');
nStimPlayed = getappdata(handles.playRipple,'nStimPlayed');
if nStimPlayed == 0
    delete(fullfile(getappdata(handles.playRipple,'stimfile_dir'),...
        getappdata(handles.playRipple,'stimfile_name')));
    return
end

stimPlayedID = zeros(nStimPlayed,1);
stimPlayedAtten = zeros(nStimPlayed,1);

i = 1;
name = sprintf('stimPlayedID%d',i);
if isappdata(handles.playRipple,name)
    stimPlayedID(i) = getappdata(handles.playRipple,name);
    rmappdata(handles.playRipple,name);  % remove appdata
end
name = sprintf('stimPlayedAtten%d',i);
if isappdata(handles.playRipple,name)
    stimPlayedAtten(i) = getappdata(handles.playRipple,name);
    rmappdata(handles.playRipple,name);  % remove appdata
end
cheetahMode = getappdata(handles.playRipple,'cheetahMode');
if cheetahMode
    name = sprintf('stimPlayedEvent%d',i);
    if isappdata(handles.playRipple,name)
        stimPlayedEvent = getappdata(handles.playRipple,name);
        rmappdata(handles.playRipple,name);  % remove appdata
        % preallocate memory for stimPlayedEvent
        stimPlayedEvent = repmat(stimPlayedEvent,nStimPlayed,1);
    end
end

for i = 2:nStimPlayed
    name = sprintf('stimPlayedID%d',i);
    if isappdata(handles.playRipple,name)
        stimPlayedID(i) = getappdata(handles.playRipple,name);
        rmappdata(handles.playRipple,name);  % remove appdata
    end
    name = sprintf('stimPlayedAtten%d',i);
    if isappdata(handles.playRipple,name)
        stimPlayedAtten(i) = getappdata(handles.playRipple,name);
        rmappdata(handles.playRipple,name);  % remove appdata
    end
    
    if cheetahMode
        name = sprintf('stimPlayedEvent%d',i);
        if isappdata(handles.playRipple,name)
            stimPlayedEvent(i) = getappdata(handles.playRipple,name);
            rmappdata(handles.playRipple,name);  % remove appdata
        end
    end
end
% hMat.nCycle = getappdata(handles.playRipple,'nCycle');
hMat.skipStim = getappdata(handles.playRipple,'skipStim');
hMat.softGain = str2double(get(handles.softGain,'string'));
hMat.nStimPlayed = nStimPlayed;
hMat.stimPlayedID = stimPlayedID;
hMat.stimPlayedAtten = stimPlayedAtten;
nRepeat = getappdata(handles.playRipple,'nRepeat');
hMat.reps_total = nRepeat;
nStim = length(hMat.allStim);
hMat.reps_completed = floor(nStimPlayed/nRepeat/nStim);
hMat.fs = handles.settings.TDT.sampling_rate;
hMat.randomized = 1; % playRipple GUI always randonmize stim presentation
if get(handles.sel_other,'value')
    hMat.param = getappdata(handles.playRipple,'param');
end
% try
%     if get(handles.sel_DRC,'value')
%         hMat.param = getappdata(handles.playRipple,'param');
%     end
% catch
% end
if cheetahMode
    hMat.stimPlayedEvent = stimPlayedEvent;
end
sel = get(handles.trgSrc,'value');
str = get(handles.trgSrc,'string');
trgSrc = str{sel};
if strncmp(trgSrc, 'User P(0.0)', length('User P(0.0)'))
    useFrameClk = 0;
end
if strncmp(trgSrc, 'Frame Clk', length('Frame Clk'))
    useFrameClk = 1;
end
if get(handles.imagingMode,'value')
    if useFrameClk
        hMat.nBaseFrame = str2double(get(handles.nBaseFrame,'string'));
        hMat.framesPerStim = str2double(get(handles.cntFrames,'string'));
    else
    end
end
fprintf('done saving sound file data to: \n%s\n',hMat.Properties.Source);
return