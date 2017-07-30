function playRipple_reset_TDT_circuit(handles)
RX6 = getappdata(handles.playRipple,'RX6');
RX6.SoftTrg(2); % software trigger resets circuit
RX6.Halt;
return