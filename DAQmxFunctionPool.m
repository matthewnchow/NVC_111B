function varargout = DAQmxFunctionPool(varargin)

switch varargin{1}
    case 'GetCounts'
        %[meanCounts, stdCOunts] = GetCounts(SamplingFreq,Samples);
        [varargout{1},varargout{2}] = GetCounts(varargin{2},varargin{3});
    case 'WriteVoltage'
        %[varargout{1}] = WriteVoltage(Device,Voltage);
        [varargout{1}] = WriteVoltage(varargin{2},varargin{3});
    case 'WriteVoltages'
        [varargout{1}] = WriteVoltages(varargin{2},varargin{3});
    case 'CreateAIChannel'
        [varargout{1},varargout{2}] = CreateAIChannel(varargin{2},varargin{3},varargin{4});
    case 'ReadVoltageScalar'
        [varargout{1}] = ReadVoltageScalar(varargin{2});
    case 'WriteDigitalChannel'
        WriteDigitalChannel(varargin{2},varargin{3});
    case 'ReadDigitalChannel'
        ReadDigitalChannel(varargin{2});
    case 'SetGatedCounter'
        SetGatedCounter(varargin{2},varargin{3},varargin{4},varargin{5});
    case 'SetGatedNCounter'
        [varargout{1},varargout{2}] = SetGatedNCounter(varargin{2});
    case 'ReadCounterScalar'
        [varargout{1},varargout{2}] = ReadCounterScalar(varargin{2});
    case 'WriteAnalogVoltage'
        [varargout{1},varargout{2}] = WriteAnalogVoltage(varargin{2},varargin{3},varargin{4},varargin{5});
    case 'ReadAnalogVoltage'
        [varargout{1},varargout{2}] = ReadAnalogVoltage(varargin{2},varargin{3},varargin{4});
    case 'SetCounter'
        [varargout{1},varargout{2}] = SetCounter(varargin{2});
    case 'ReadCounter'
        [varargout{1}] = ReadCounter(varargin{2}, varargin{3}, varargin{4});
    case 'SetStartStopCounter'
        switch length(varargin)
            case 1
                N = 2;
                startstop = PortMap('StartStop');
                armstart = PortMap('ArmStart');
            case 2
                N = varargin{2};
                startstop = PortMap('StartStop');
                armstart = PortMap('ArmStart');
            case 3
                N = varargin{2};
                startstop = varargin{3};
                armstart = PortMap('ArmStart');
            case 4
                N = varargin{2};
                startstop = varargin{3};
                armstart = varargin{4};
        end
        %[status, task] = SetStartStopcounter(N, sampclk, armstart)
        [varargout{1}, varargout{2}] = SetStartStopCounter(N, startstop, armstart);
    case 'resetDAQ_StartStopSeq'
        [varargout{1}, varargout{2}, varargout{3}] = resetDAQ_StartStopSeq(varargin{2});
    otherwise
        disp('DAQmxFunctionPool - I dont get it!');
end

function [status, task] = SetStartStopCounter(N, startstop, armstart)
% N = starts + stops, replacement for N gated counter
% Hardware triggered
% Added by M Chow on 6-1-19
    DAQmx_Val_Rising = 10280; % Rising
    DAQmx_Val_Falling =10171; % Falling
    DAQmx_Val_FiniteSamps = 10178; % Finite Samples
    DAQmx_Val_CountUp = 10128; % Count Up
    DAQmx_Val_CountDown = 10124; % Count Down
    DAQmx_Val_DigEdge = 10150;

    [ stat, ~, task ] = DAQmxCreateTask([]);
    DAQmxErr(stat);

    % If you switch to a non-inverting peripheral for this channel out of
    % the FPGA, change from falling to rising edge! ~MChow 6/6/19
    stat = DAQmxCreateCICountEdgesChan(task,PortMap('Ctr in'),'',... 
        DAQmx_Val_Falling, 0, DAQmx_Val_CountUp);
    DAQmxErr(stat);
    
    stat = DAQmxCfgSampClkTiming(task,startstop,1.0,... 
        DAQmx_Val_Falling,DAQmx_Val_FiniteSamps ,N);
    DAQmxErr(stat)

    % This is the hardware trigger, called "Arm Start". Ignores sampleclk
    % edges until ArmStart goes high
    DAQmxErr(calllib('mynidaqmx', 'DAQmxSetArmStartTrigType', task, DAQmx_Val_DigEdge))
    DAQmxErr(calllib('mynidaqmx', 'DAQmxSetDigEdgeArmStartTrigSrc', task, armstart));
    status = calllib('mynidaqmx', 'DAQmxSetDigEdgeArmStartTrigEdge', task, DAQmx_Val_Rising);
    DAQmxErr(status)

    
function [status, hCounter] = resetDAQ_StartStopSeq(N)
    % Added by MChow 6-1-19
    global gScan gConfocal
    DAQmxResetDevice('Dev1')
    DAQmxResetDevice('myDAQ1')
    %Fix Galvos in base position
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo x'), gScan.FixVx + gConfocal.XOffSet);
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo y'), gScan.FixVy + gConfocal.YOffSet);
    DAQmxFunctionPool('WriteVoltage',PortMap('Obj_Piezo'), gScan.FixVz);    
    
    [status, hCounter] = SetStartStopCounter(N, PortMap('StartStop')); %Normally '/Dev1/PFI9'
    DAQmxErr(status)


function ReadDigitalChannel(Device, Data)
DAQmx_Val_ChanPerLine =0; % One Channel For Each Line
DAQmx_Val_ChanForAllLines =1; % One Channel For All Lines
DAQmx_Val_GroupByChannel = 0; % Group per channel

data = [0 0 0 0 0 0 0 0];
read = [0 0 0 0 0 0 0 0];
bytesPerSamp = 0;
% DAQmx Configure Code
[ status, ~, task ] = DAQmxCreateTask([]);
DAQmxErr(status);
DAQmxErr(DAQmxCreateDIChan(task,'Dev1/port0/line0:7','',DAQmx_Val_ChanForAllLines));
% DAQmx Start Code
DAQmxErr(DAQmxStartTask(task));
% DAQmx Read Code
DAQmxErr(DAQmxReadDigitalLines(task,1,10.0,DAQmx_Val_GroupByChannel,data,100,read,bytesPerSamp));
% DAQmx Stop Code
DAQmxStopTask(task);
DAQmxClearTask(task);


function WriteDigitalChannel(Device,Data)
DAQmx_Val_ChanPerLine =0; % One Channel For Each Line
DAQmx_Val_ChanForAllLines =1; % One Channel For All Lines
DAQmx_Val_GroupByChannel = 0; % Group per channel

sampsPerChanWritten = 0;

% DAQmx Configure Code
[ status, TaskName, task ] = DAQmxCreateTask([]);
DAQmxErr(status);
DAQmxErr(DAQmxCreateDOChan(task,Device,'',DAQmx_Val_ChanForAllLines));
% DAQmx Start Code
DAQmxErr(DAQmxStartTask(task));
% DAQmx Write Code
DAQmxErr(DAQmxWriteDigitalLines(task,1,1,10.0,DAQmx_Val_GroupByChannel,Data,sampsPerChanWritten));
% DAQmx Stop Code
DAQmxStopTask(task);
DAQmxClearTask(task);


function [meanCounts, stdCounts] = GetCounts(SamplingFreq,Samples)

    TimeOut = 1.1*Samples/SamplingFreq;
    hCounter = SetCounter( Samples+1 );
    [status, hPulse] = DigInfPulse(SamplingFreq,0.5,Samples); %Changed to Inf by MChow
    DAQmxErr(status)

    status = DAQmxStartTask(hCounter);
    if status ~= 0,
        disp(['NI: Start Counter        :' num2str(status)])
    end
    status = DAQmxStartTask(hPulse);    
    if status ~= 0,
        disp(['NI: Start Pulse          :' num2str(status)])
    end

    A = ReadCounter( hCounter, Samples , TimeOut);

    DAQmxStopTask(hPulse);
    DAQmxStopTask(hCounter);
    DAQmxClearTask(hPulse);
    DAQmxClearTask(hCounter);

    A = diff(A);
    meanCounts = mean(A);
    stdCounts = std(A);


function [status, task] = SetCounter(N)
    DAQmx_Val_Volts= 10348; % measure volts
    DAQmx_Val_Rising = 10280; % Rising
    DAQmx_Val_FiniteSamps = 10178; % Finite Samples
    DAQmx_Val_CountUp = 10128; % Count Up
    DAQmx_Val_CountDown = 10124; % Count Down
    DAQmx_Val_GroupByChannel = 0; % Group per channel
    DAQmx_Val_ContSamps =10123; % Continuous Samples

    [ status, TaskName, task ] = DAQmxCreateTask([]);
    DAQmxErr(status);

    status = DAQmxCreateCICountEdgesChan(task,PortMap('Ctr in'),'',... %
        DAQmx_Val_Rising , 0, DAQmx_Val_CountUp);
    DAQmxErr(status);

    status = DAQmxCfgSampClkTiming(task,PortMap('Ctr sampclk'),1.0,... %'myDAQ1/PFI2';PortMap('Ctr Trig')
        DAQmx_Val_Rising,DAQmx_Val_FiniteSamps ,N);
    DAQmxErr(status);

function readArray = ReadCounter(task,N,TimeOut)
    numSampsPerChan = N;
    %readArray = libpointer('int64Ptr',zeros(1,N));
    readArray = zeros(1,N);
    arraySizeInSamps = N;
    sampsPerChanRead = libpointer('int32Ptr',0);

    [status, readArray]= DAQmxReadCounterF64(task, numSampsPerChan,...
        TimeOut, readArray, arraySizeInSamps, sampsPerChanRead );
    DAQmxErr(status)


function status = WriteVoltages(Devices,Voltages)
    DAQmx_Val_Volts= 10348; % measure volts
    status = -1;
    for k=1:length(Devices)
        switch Devices{k}
            case {PortMap('Galvo x'),PortMap('Galvo y')}
                if abs(Voltages(k)) > 3
                    disp('Error in WriteVoltage (Voltage exceeds 3 volts)');
                    return;
                end
        end
        [ status, TaskName, task ] = DAQmxCreateTask([]);
        status = status + DAQmxCreateAOVoltageChan(task,Devices{k},-10,10,DAQmx_Val_Volts);
        status = status + DAQmxWriteAnalogScalarF64(task,1,0,Voltages(k));
        if status ~= 0
            disp(['Error in writing voltage in Device ' Devices{k}]);
        end
        DAQmxClearTask(task);
    end

function status = WriteVoltage(Device,Voltage)
    DAQmx_Val_Volts= 10348; % measure volts
    status = -1;
    switch Device
        case {PortMap('Galvo x'),PortMap('Galvo y')}
            if abs(Voltage) > 3
                disp('Error in WriteVoltage (Voltage exceeds 3 volts)');
                return;
            end
    end

    [ status, TaskName, task ] = DAQmxCreateTask([]);
    status = status + DAQmxCreateAOVoltageChan(task,Device,-10,10,DAQmx_Val_Volts);
    status = status + DAQmxWriteAnalogScalarF64(task,1,0,Voltage);
    if status ~= 0
        disp(['Error in writing voltage in Device ' Device]);
    end
    DAQmxClearTask(task);

function [data, status] = ReadCounterScalar(task)
% added by Satcher 10/19/2016
global gTimeOut
data=0;
timeout=gTimeOut;
[status, data] = calllib('mynidaqmx','DAQmxReadCounterScalarU32',...
task, timeout,data,[]); 
DAQmxErr(status);
DAQmxStopTask(task);
%DAQmxClearTask(task);

function [status, hScan] = WriteAnalogVoltage(chan,vec, samps, freq)
DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
% DAQmx_Val_Falling = 10171; % Falling
DAQmx_Val_FiniteSamps = 10178; % Finite Samples
DAQmx_Val_GroupByChannel = 0; % Group per channel

[ status, ~, hScan ] = DAQmxCreateTask([]);
DAQmxErr(status);
status = DAQmxCreateAOVoltageChan(hScan,chan,-5,5,DAQmx_Val_Volts);
DAQmxErr(status);
status = DAQmxCfgSampClkTiming(hScan,PortMap('Ctr Trig'),freq,...
    DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,samps);
DAQmxErr(status);
zero_ptr = libpointer('int32Ptr',zeros(1,samps));
status = DAQmxWriteAnalogF64(hScan, samps, 0, 10,...
    DAQmx_Val_GroupByChannel, vec, zero_ptr);
DAQmxErr(status);

function [status, hRead] = CreateAIChannel(trig, samps, freq)
DAQmx_Val_RSE=10083;
DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
%DAQmx_Val_Falling = 10171; % Falling
DAQmx_Val_FiniteSamps = 10178; % Finite Samples

[ status, ~, hRead ] = DAQmxCreateTask([]);

DAQmxErr(status);

status = DAQmxCreateAIVoltageChan(hRead,PortMap('APD in'),'',...
    DAQmx_Val_RSE, -5,5, DAQmx_Val_Volts,[]);

DAQmxErr(status);
status = DAQmxCfgSampClkTiming(hRead,trig,freq,DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,samps); %ctr1 out
DAQmxErr(status);

function [status, RawData] = ReadAnalogVoltage(task, samps, timeout)
DAQmx_Val_GroupByChannel = 0; % Group per channel
sampsPerChanRead = libpointer('int32Ptr',0);
RawData = zeros(1,samps);
[status, RawData]=DAQmxReadAnalogF64(task,samps,timeout,DAQmx_Val_GroupByChannel,RawData,samps,sampsPerChanRead,[]);
DAQmxErr(status);

function [Answer] = ReadVoltageScalar(Device)
%This does not work
DAQmx_Val_Volts= 10348; % measure volts
[ status, TaskName, task ] = DAQmxCreateTask([]);
status = DAQmxCreateAOVoltageChan(task,Device,-10,10,DAQmx_Val_Volts);
[status, output] = DAQmxReadAnalogScalarF64(task);
if status ~= 0
    disp(['Error while reading voltage in Device ' Device]);
    DAQmxErr(status);
end
DAQmxClearTask(task);

Answer.status = status;
Answer.Voltage = output;