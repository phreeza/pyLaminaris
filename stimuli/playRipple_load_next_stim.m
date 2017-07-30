function [wave,nSamps,ITD,gain2] = playRipple_load_next_stim(handles,idx)
% idx: indicates position of stimSeq
% wav: actual waveform in .waveform file
% nSamps: number of samples in waveform

%------ old version before 1/25/2016 ------------
% stimSeq = getappdata(handles.playRipple,'stimSeq');
% allStimWav = getappdata(handles.playRipple,'allStimWav');
% nSAMPS = getappdata(handles.playRipple,'nSamps');
% wave = allStimWav{stimSeq(idx,1)};  
% nSamps = nSAMPS(stimSeq(idx,1));

%------- new verstion after 1/25/2016 ------------
skipStim = getappdata(handles.playRipple,'skipStim');
stimSeq = getappdata(handles.playRipple,'stimSeq');
fs = handles.settings.TDT.sampling_rate;
nCmp = round(0.2*fs);  % test if stim is shorter than 200ms, which is chosen
%to aid circuit status detection, then pad with zeros to increase length
tic
if ~skipStim(stimSeq(idx)) % if this stim is not skipped then it means after softGain it remains within TDT output range
    
    allStim = getappdata(handles.playRipple,'allStim');
    f = load(allStim{stimSeq(idx)});
    softGain = 10^(str2double(get(handles.softGain,'string'))/20);
    wave = f.wave*handles.settings.TDT.range*softGain;
    if length(wave) < nCmp
        wave = [wave zeros(1,nCmp - length(wave))];
    end
    nSamps = length(wave);
    if get(handles.withITD,'value')
        ITD = f.ITD; % ITD in microseconds;
        try
            gain2 = f.gain2;
            fprintf('succeeded load gain2\n');
        catch
            gain2 = 0;
        end
    else
        ITD = [];
        gain2 = [];
    end
else
    wave = zeros(1,nCmp);
    nSamps = nCmp;
    ITD = [];
    gain2 = [];
end
t = toc;
fprintf('new version! time to load data: %g sec\n',t);

return
