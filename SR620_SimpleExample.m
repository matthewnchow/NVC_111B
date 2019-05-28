% Example communication with SR620
% Find a GPIB object. Call it obj1
obj1 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = gpib('NI', 0, 16);
else
    fclose(obj1);
    obj1 = obj1(1);
end
fopen(obj1);
% Communicating with instrument object, obj1.
fprintf(obj1, 'MODE6'); %Set to count mode
fprintf(obj1, 'SRCE0'); %Set source to channel A
% fprintf(obj1, 'GATE 0.001'); %Set gate to 1 ms
fprintf(obj1, 'SIZE1');
fprintf(obj1, 'TSLP1,0'); %Set ChA trigger slope to positive
fprintf(obj1, 'LEVL1,.48'); %Set gate trigger level to 1V
data1 = query(obj1, 'XALL?'); %Measure and put data for mean, rel, jitter, max, min of measurement in output buffer

% Simple, but prob too slow:
% data2 = []
% for i = 1:1000
%     data2 = [data2, query(obj1, 'MEAS?0')];
% end 
% data2

fprintf(obj1, 'BDMP2');
A = fread(obj1, 4, 'double');
b = de2bi(A);
fprintf(obj1, 'MODE6'); %Reset to escape binary dump
size(A)
% Disconnect from instrument object, obj1.
fclose(obj1);
