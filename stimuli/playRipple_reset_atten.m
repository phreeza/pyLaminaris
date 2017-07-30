function succeeded = playRipple_reset_atten(handles)
% idx: which stim in stimSeq to play
succeeded = 1;
atten = 120;
PA5 = getappdata(handles.playRipple,'PA5');
for i = 1:length(PA5)
    PA5(i).SetAtten(atten);
    error = PA5(i).GetError();
    if ~isempty(error)
        PA5(i).Display(error,0);
        report_status(handles,'error setting atten!');
        succeeded = 0;
    end
end

return

function report_status(handles,message)
set(handles.status_str,'string',message,'visible','on');
display(sprintf('\n%s\n',message));
drawnow;
return