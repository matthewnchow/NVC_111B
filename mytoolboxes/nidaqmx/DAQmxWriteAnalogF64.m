function [status, a, b]=DAQmxWriteAnalogF64(taskHandle, numSampsPerChan, autoStart,...
    timeout, dataLayout, writeArray, sampsPerChanWritten)

[status, a, b]=calllib('mynidaqmx','DAQmxWriteAnalogF64',taskHandle, numSampsPerChan,...
    autoStart,timeout, dataLayout, writeArray, sampsPerChanWritten,[]);
