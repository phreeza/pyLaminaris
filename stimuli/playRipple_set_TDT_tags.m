function succeeded = playRipple_set_TDT_tags(varargin)
succeeded = 1;
switch nargin
    case  3
        handles = varargin{1};
        wav = varargin{2};
        nSamps = varargin{3};
        
    case 4
        handles = varargin{1};
        wav = varargin{2};
        nSamps = varargin{3};
        useFrameClk = varargin{4};
    case 5
        handles = varargin{1};
        wav = varargin{2};
        nSamps = varargin{3};
        useFrameClk = varargin{4};
        ITD = varargin{5};
    case 6
        handles = varargin{1};
        wav = varargin{2};
        nSamps = varargin{3};
        useFrameClk = varargin{4};
        ITD = varargin{5};
        gain2 = varargin{6};
        
end

% switch get(handles.trgSrc,'value')
%     case 1 % UserP(0.0)
%         useFrameClk = 0;
%     case 2 % Frame Clk
%         useFrameClk = 1;
% end

RX6 = getappdata(handles.playRipple,'RX6');


% ------------ commented JL 1/25/2016 --------------------------
% if get(handles.cheetahMode,'value')
%     if ~RX6.SetTagVal('bufSz',nSamps);
%         report_status(handles,'error setting tag bufSz!');
%         succeeded = 0;
%         return
%     end
%     if ~RX6.SetTagVal('nCycle',getappdata(handles.playRipple,'nCycle'))
%         report_status(handles,'error setting tag nCycle!');
%         succeeded = 0;
%         return
%     end
% end

if get(handles.imagingMode,'value')
    if useFrameClk
        if ~RX6.SetTagVal('nBaseFrame',str2double(get(handles.nBaseFrame,'string')) +1 )
            report_status(handles,'error setting tag nBaseFrame!');
            succeeded = 0;
            return
        end
        if ~RX6.SetTagVal('nFrame',getappdata(handles.playRipple,'nFrame'))
            report_status(handles,'error setting tag nFrame!');
            succeeded = 0;
            return
        end
    else
        
        
    end
    
    if get(handles.sound_only,'value') % if checked, set trg source to SoftTrg such that do not require frame clk to play sound
        RX6.SetTagVal('trg_sel',1);
    else
        RX6.SetTagVal('trg_sel',0);
    end
end

if size(wav,1)~=1 % if wav not a row vector
    wav = wav';
end

if get(handles.withITD,'value')
    fs = getappdata(handles.playRipple,'curr_fs');
    if ITD >=0
        nBlank = ceil(fs*ITD*1e-3); % ITD should be in mseconds
        wave2write1 = [wav zeros(1,nBlank)];
        wave2write2 = [zeros(1,nBlank),wav];
    else
        nBlank = ceil(fs*abs(ITD)*1e-3);
        wave2write1 = [zeros(1,nBlank),wav];
        wave2write2 = [wav zeros(1,nBlank)];
    end
    % apply additional gain on channel2
    wave2write2 = 10^(gain2/20) * wave2write2;
    fprintf('succeeded setting gain2\n');
    
    N = length(wave2write1);
    if ~RX6.SetTagVal('bufSzCmp',N-1);
        report_status(handles, 'error setting tag bufSzCmp!');
        succeeded = 0;
        return
    end
    
    if ~RX6.WriteTagV('data',0,wave2write1); % WriteTagV takes row vectors
        report_status(handles,'error setting tag data!');
        succeeded = 0;
        return
    end
    
    if ~RX6.WriteTagV('data2',0,wave2write2); % WriteTagV takes row vectors
        report_status(handles,'error setting tag data!');
        succeeded = 0;
        return
    end
    
else
    if ~RX6.SetTagVal('bufSzCmp',nSamps-1);
        report_status(handles, 'error setting tag bufSzCmp!');
        succeeded = 0;
        return
    end
    
    if ~RX6.WriteTagV('data',0,wav); % WriteTagV takes row vectors
        report_status(handles,'error setting tag data!');
        succeeded = 0;
        return
    end
    
    
end

if ~RX6.SetTagVal('channel_sel',get(handles.channel_sel,'value'))
    report_status(handles,'error setting tag channel_sel!');
    succeeded = 0;
    return
end




return

function report_status(handles,message)
set(handles.status_str,'string',message,'visible','on');
display(sprintf('\n%s\n',message));
drawnow;