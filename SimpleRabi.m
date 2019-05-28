%% Run Parameters
MWpower = -5; %MW power(dBm);
MWfreq = 2.853; %(GHz) -1 branch, low field, 4 degenerate directions;

steps = 100; %Number of steps to take in MW on time
dt = 20 * 10^-9; %20 ns step size
counttime = 300 * 10^-9; %Counters active for this long in seconds
pumptime = 10 * 10^-6; %10 us pumping
cushion = 20 * 10^-9; %20 ns pad after turning on or off laser
% trigwidth = 50 * 10^-9; %Reference trigger width in sec

%Setup for the night, no need for this again
% %% Connect to instruments and data file
% %% Sig Generator setup
% sg = instrfind('Type', 'serial', 'Port', PortMap('SG com'), 'Tag', '');
% 
% if isempty(sg)
%     sg = serial(PortMap('SG com'));
% else
%     fclose(sg);
%     sg = sg(1);
% end
% set(sg,'BaudRate',115200);
% fopen(sg);
% fprintf(sg, ['AMPR ',num2str(MWpower)]); %set power in dBm
% fprintf(sg, ['FREQ ', MWfreq, 'e9']);
% fclose(sg);
% 
% %% FPGA setup                
% Cyc10FunctionPool('RabiConfigC10', steps, dt, counttime, pumptime, cushion) %%Check reliability...

%% File settings
droot = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\';
ddir = [droot, '5-28_SimpleRabi_Distance4_MidG_-5dBm_150um_\'];
if ~exist(ddir, 'dir')
  mkdir(ddir);
end

%% DAQ setup
LoadNIDAQmx

Vx = 0.34;
Vy = 2.06;
%Fix Galvos in base position
DAQmxFunctionPool('WriteVoltage',PortMap('Galvo x'), Vx);
DAQmxFunctionPool('WriteVoltage',PortMap('Galvo y'), Vy);
samples = steps * 2;
% passes = idivide(int32(1) , counttime * 1000); %1 milisecond of integration time per point per run
passes = 1;

N = 2 * 2 * steps * passes;
[~, hCounter] = SetCounter(N, '/Dev1/PFI4'); 
DAQmx_Val_Rising = 10280;

%% Run multiple times and write data, they'll start on different times for now...
    % Add DAQ hardware trigger later
    % May need fourier for processing

dat = zeros(steps, 3);
ts = int32((0:steps - 1) * dt * 10^9); %Print in ns 
odds = (0:((steps * passes * 2) - 1))*2 + 1;
evens = odds + 1;
dat(:,1) = ts;
for run_idx = 0:1000000
    dfilename = ['Rabi_',num2str(MWpower),'dBm_', ...
        num2str(counttime * 10 ^9),'ns_run', num2str(run_idx),'.csv'];
    dfile = fopen([ddir, dfilename],'w+'); 
    fprintf(dfile, '%s\r\n', ['MW on time (s), Data Cts/',num2str(counttime), ', Ref cts']);

    for i = 1:20000
        status = DAQmxCfgDigEdgeStartTrig (hCounter, '/Dev1/PFI3', DAQmx_Val_Rising);
        arr = DAQmxFunctionPool('ReadCounter', hCounter, N, -1);
        DAQmxStopTask(hCounter);
%         while (arr == zeros(N))
%             [~, hCounter] = resetDAQ_Rabi(Vx, Vy, N)
%             status = DAQmxCfgDigEdgeStartTrig (hCounter, '/Dev1/PFI3', DAQmx_Val_Rising);
%             arr = DAQmxFunctionPool('ReadCounter', hCounter, N, -1);
%             DAQmxStopTask(hCounter);
%             mean(arr)
%         end
        arr = (arr(evens) - arr(odds));
        arr = reshape(arr, 2, steps * passes);
        cts = reshape(arr(1,:), steps, passes);
        ref = reshape(arr(2,:), steps, passes);
        dat(:,2) = dat(:,2) + sum(cts,2);
        dat(:,3) = dat(:,3) + sum(ref,2);
    end
	fprintf(dfile, '%d,%d,%d\r\n', transpose(dat));
    fclose(dfile);
    if (mod(run_idx,100) == 1)
       run_idx 
    end
end

DAQmxClearTask(hCounter)


function [status, task] = SetCounter(N, sampclk)
    DAQmx_Val_Volts= 10348; % measure volts
    DAQmx_Val_Rising = 10280; % Rising
    DAQmx_Val_FiniteSamps = 10178; % Finite Samples
    DAQmx_Val_CountUp = 10128; % Count Up
    DAQmx_Val_CountDown = 10124; % Count Down
    DAQmx_Val_GroupByChannel = 0; % Group per channel
    DAQmx_Val_ContSamps =10123; % Continuous Samples

    [ status, ~, task ] = DAQmxCreateTask([]);
    DAQmxErr(status);

    status = DAQmxCreateCICountEdgesChan(task,PortMap('Ctr in'),'',... %
        DAQmx_Val_Rising , 0, DAQmx_Val_CountUp);
    DAQmxErr(status);

    status = DAQmxCfgSampClkTiming(task,sampclk,1.0,... 
        DAQmx_Val_Rising,DAQmx_Val_FiniteSamps ,N);
    DAQmxErr(status);
end


function [stat, hCounter] = resetDAQ_Rabi(vx, vy, N)
    DAQmxResetDevice('Dev1')
    %Fix Galvos in base position
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo x'), vx);
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo y'), vy);
    [stat, hCounter] = SetCounter(N, '/Dev1/PFI4'); 
    DAQmx_Val_Rising = 10280;
end