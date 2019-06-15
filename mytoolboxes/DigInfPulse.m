%% Continuous version of DigPulseTrainCont
%  Creates task to put out pulses at freq and duty cycle continuously
%  Written by Matthew Chow 5/24/19 to patch compatibility 
%    with 2 counter, M series devices (NI 6251)
function [status, task] = DigInfPulse(Freq,DutyCycle)

DAQmx_Val_Volts= 10348; % measure volts
DAQmx_Val_Rising = 10280; % Rising
DAQmx_Val_FiniteSamps = 10178; % Finite Samples
DAQmx_Val_CountUp = 10128; % Count Up
DAQmx_Val_CountDown = 10124; % Count Down
DAQmx_Val_GroupByChannel = 0; % Group per channel
DAQmx_Val_ContSamps =10123; % Continuous Samples
DAQmx_Val_Hz = 10373; % Hz
DAQmx_Val_Low = 10214; % Low

[status,~,task] = DAQmxCreateTask('');
DAQmxErr(status);
try
    status = DAQmxCreateCOPulseChanFreq(task,PortMap('Ctr out'),'',DAQmx_Val_Hz,DAQmx_Val_Low,0.0,Freq, DutyCycle);
    DAQmxErr(status);
    status = DAQmxCfgImplicitTiming(task,DAQmx_Val_ContSamps,0);
    DAQmxErr(status);
catch ME
    DAQmxStopTask(task);
    DAQmxClearTask(task);
    rethrow(ME);
end

end