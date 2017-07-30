function playRipple_timer_startfcn(obj,event,handles)

%     
%         setappdata(handles.playRipple,'waitForTrg',1);
%         setappdata(handles.playRipple,'nStimPlayed',0);
%         nPresent = getappdata(handles.playRipple,'nPresent');
%         report_status(handles,sprintf('Waiting for trigger, stim 1 / %d',nPresent));
%         
% 
%         
%         setappdata(handles.playRipple,'stimfile_dir',...
%                 stimfile_dir);
%         allStim = getappdata(handles.playRipple,'allStim');
%         nStimPlayed = 0;
%         stimfile_name = sprintf('stim_file_%s.mat',datestr(now,'yyyy-mm-dd_HH-MM-SS'));
%         setappdata(handles.playRipple,'stimfile_name',stimfile_name);
%         save(fullfile(stimfile_dir,stimfile_name),'allStim','nStimPlayed','-v7.3');
%         hMat = matfile(fullfile(stimfile_dir,stimfile_name),'Writable',true); % create a handle for .mat file just created
%         setappdata(handles.playRipple,'hMat',hMat); % save handle to .mat file as appdata
%         % let's assume by now cheetah has been successfully connected
%         if cheetahMode
%         end
        

        
        %%
        setappdata(handles.playRipple,'waitForTrg',1);
        setappdata(handles.playRipple,'nStimPlayed',0);
        nPresent = getappdata(handles.playRipple,'nPresent');
        report_status(handles,sprintf('Waiting for trigger, stim 1 / %d',nPresent)); 
        switch handles.whichrig
            case 'Cheetah'
                cheetahMode = getappdata(handles.playRipple,'cheetahMode');
                if cheetahMode
                    %hCheetah = findall(0,'type','figure','tag','cheetah_helper'); % find handle of cheetah_helper
                    hCheetah = findobj('tag','cheetah_helper');
                    if isempty(hCheetah)
                        error('cheetah_helper is not running!');
                        return
                    end
                    stimfile_dir = getappdata(hCheetah,'cheetah_dir');
                    GetCheetahEvent; % flush all previous events
                else
                    stimfile_dir = 'D:\sound_files_playRipple';
                    if ~isdir(stimfile_dir)
                        mkdir(stimfile_dir);
                    end
                end
                
            case 'Prairie'
                stimfile_dir = 'D:\sound_files_playRipple';
                if ~isdir(stimfile_dir)
                    mkdir(stimfile_dir);
                end
            case 'BScope'
                stimfile_dir = 'D:\Ji\sound_files_playRipple';
                if ~isdir(stimfile_dir)
                    mkdir(stimfile_dir);
                end
        end

        setappdata(handles.playRipple,'stimfile_dir',stimfile_dir);
        allStim = getappdata(handles.playRipple,'allStim');
        nStimPlayed = 0;
        stimfile_name = sprintf('sound_file_%s.mat',datestr(now,'yyyy-mm-dd_HH-MM-SS'));
        setappdata(handles.playRipple,'stimfile_name',stimfile_name);
        save(fullfile(stimfile_dir,stimfile_name),'allStim','nStimPlayed','-v7.3');   
        hMat = matfile(fullfile(stimfile_dir,stimfile_name),'Writable',true); % create a handle for .mat file just created
        setappdata(handles.playRipple,'hMat',hMat); % save handle to .mat file as appdata
        if get(handles.sel_long_stim,'value')
            setappdata(handles.playRipple,'nSegPlayed',0);
            setappdata(handles.playRipple,'stimInfoSaved',0);
        end
        RX6 = getappdata(handles.playRipple,'RX6');
        RX6.SoftTrg(1);
        RX6.SoftTrg(2); % reset all counters in TDT circuit to 0

return


function report_status(handles,message)
set(handles.status_str,'string',message,'visible','on');
display(sprintf('\n%s\n',message));
drawnow;