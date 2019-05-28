
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
set(rs232, 'BaudRate', 9600);
set(rs232, 'DataTerminalReady', 'off');
set(rs232, 'RequestToSend', 'off');
set(rs232, 'Timeout', 2.0);
set(rs232, 'DataBits', 8);
rs232.StopBits = 1;
% rs232.Terminator = 'CR'

CR = 13;
NL = 10;

% fwrite(rs232, uint32(25000000), 'uint32');
% data = uint16(2);
% de2bi(data)
% fwrite(rs232, CR);
% fwrite(rs232, NL);
% fwrite(rs232, uint16(data), 'uint16');
% fread(rs232, 1)

%Info type ID bytes
state0 = uint8(0);
period = uint8(1);
ed = uint8(2);
outer_period = uint8(3);
clr = uint8(5);

%Default values
period_a = int32(125000000);
state0_a = uint8(1); %0000_0001
ed_a0 = uint64(bitshift(25000000, 4)); %Start with the edge value, then shift over 4
ed_a0 = bitset(ed_a0, 1, 1); %Turn on this edge
% ed_a0 = bitrevorder(ed_a0)
de2bi(ed_a0);

flushinput(rs232)

% Clear edges
fwrite(rs232, clr);
fwrite(rs232, CR);
fwrite(rs232, NL);
fread(rs232, 1)

%% Set initial state
fwrite(rs232, state0);
fwrite(rs232, state0_a);
fwrite(rs232, CR);
fwrite(rs232, NL);
fread(rs232, 1)

%% Set period
fwrite(rs232, period);
fwrite(rs232, period_a, 'uint32');
fwrite(rs232, CR);
fwrite(rs232, NL);
fread(rs232, 1)

%% Set outer period
fwrite(rs232, outer_period);
fwrite(rs232, int32(10), 'uint32');
fwrite(rs232, CR);
fwrite(rs232, NL);
fread(rs232, 1)

%% Set first edge
fwrite(rs232, ed);
fwrite(rs232, uint8(0)); % Edge ID 
fwrite(rs232, int32(62.5000000), 'int32'); % x = 125,000,000
fwrite(rs232, int32(0), 'int32'); % dx = 0
fwrite(rs232, uint8(8), 'uint8'); % Channel 0, enabled on
fwrite(rs232, CR);
fwrite(rs232, NL);
fread(rs232, 1)

    %% Set edge
    fwrite(rs232, ed);
    fwrite(rs232, uint8(1)); % Edge ID 
    fwrite(rs232, int32(200000000), 'int32'); % x = 200,000,000
    fwrite(rs232, int32(10000000), 'int32'); % dx = 20,000,000
    fwrite(rs232, uint8(9)); % Channel 1, enabled on
    fwrite(rs232, CR);
    fwrite(rs232, NL);
    fread(rs232, 1)
%     
%         %% Set edge
%     fwrite(rs232, ed);
%     fwrite(rs232, uint8(2)); % Edge ID 
%     fwrite(rs232, int32(30000), 'int32'); % x = 20,000,000
%     fwrite(rs232, int32(0), 'int32'); % dx = 0
%     fwrite(rs232, uint8(10)); % Channel 2, enabled on
%     fwrite(rs232, CR);
%     fwrite(rs232, NL);
%     fread(rs232, 1)
%     
%             %% Set edge
%     fwrite(rs232, ed);
%     fwrite(rs232, uint8(3)); % Edge ID 
%     fwrite(rs232, int32(35000), 'int32'); % x = 20,000,000
%     fwrite(rs232, int32(0), 'int32'); % dx = 0
%     fwrite(rs232, uint8(10)); % Channel 2, enabled on
%     fwrite(rs232, CR);
%     fwrite(rs232, NL);
%     fread(rs232, 1)
%     
%                 %% Set edge
%     fwrite(rs232, ed);
%     fwrite(rs232, uint8(4)); % Edge ID 
%     fwrite(rs232, int32(40000 ), 'int32'); % x = 20,000,000
%     fwrite(rs232, int32(0), 'int32'); % dx = 0
%     fwrite(rs232, uint8(11)); % Channel 1, enabled on
%     fwrite(rs232, CR);
%     fwrite(rs232, NL);
%     fread(rs232, 1)
%    

% fwrite(rs232, period);
% fwrite(rs232, period_a);
% fwrite(rs232, CR);
% 
% fwrite(rs232, ed);
% fwrite(rs232, ed_a0);
% fwrite(rs232, edID_a0);
% fwrite(rs232, CR);

% flushinput(rs232);
% fread(rs232,16)

% Close and delete
fclose(rs232);
delete(rs232);