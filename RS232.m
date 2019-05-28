%% Instrument Connection

% Find a serial port object.
rs232 = instrfind('Type', 'serial', 'Port', 'COM9', 'Tag', '');

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(rs232)
    rs232 = serial('COM9');
else
    fclose(rs232);
    rs232 = rs232(1);
end

% Connect to instrument object
fopen(rs232);


%% Configure instrument object
set(rs232, 'BaudRate', 115200);
set(rs232, 'DataTerminalReady', 'off');
set(rs232, 'RequestToSend', 'off');
set(rs232, 'Timeout', 3.0);
set(rs232, 'DataBits', 8);
rs232.StopBits = 1;
% rs232.Terminator = 'CR'
% 
% byte1 = uint8(63)
% value = [uint8(63), uint8(31), uint8(15), uint8(7), uint8(1)]
% ed_ID = [uint8(1), uint8(1)]
% bindata = typecast([byte1,value,ed_ID], 'double')
bindata = 0;
strdata = 7;
% 
% typecast(fread(handle, samples, 'double'), 'int64')
% Communicating with instrument object, obj1.
flushinput(rs232)
% fwrite(rs232, uint32(5000000), 'uint32')
for i = 1:10
    fwrite(rs232, uint8(0));
    fwrite(rs232, uint8(9));
    fwrite(rs232, uint8(55));
end

fread(rs232, 30)
% fread(rs232,2)
% fread(rs232,4)
% fwrite(rs232, 13)
% fwrite(rs232, 10)
% flushinput(rs232)
% fread(rs232,4)
% 
% % fwrite(rs232, 13); %CR
% for i = 0:15
% %     fwrite(rs232, i);
%     fwrite(rs232, i);
%     pause(0.01)
%     flushinput(rs232);
%     fread(rs232, 1)
% %     pause(1) 
% end
% flushinput(rs232);
% fread(rs232,16)

% Close and delete
fclose(rs232);
delete(rs232);