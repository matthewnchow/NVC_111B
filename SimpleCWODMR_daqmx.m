%% Connect to instruments and data file
LoadNIDAQmx

Vx = 0.34;
Vy = 2.06;
%Fix Galvos in base position
DAQmxFunctionPool('WriteVoltage',PortMap('Galvo x'), Vx);
DAQmxFunctionPool('WriteVoltage',PortMap('Galvo y'), Vy);

samps = 21;
dt = 0.01; %10 ms bins, total integration time of 1 sec per pt for 5 runs, 20 samps
sampfreq = 1/dt; 
[~, hPulse] = DigInfPulse(sampfreq, 0.5);
[~, hCounter] = DAQmxFunctionPool('SetCounter', samps); 

SRS384 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 27, 'Tag', '');
if isempty(SRS384)
    SRS384 = gpib('NI', 0, 27);
else
    fclose(SRS384);
    SRS384 = SRS384(1);
end
fopen(SRS384);

droot = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\';
ddir = [droot, '5-27_cwODMR_Distance4_MidG_fullspec_150um_PowerBroad\'];
if ~exist(ddir, 'dir')
  mkdir(ddir);
end

%% Sweep frequency settings
freq_i = 2.78;
freq_f = 2.94;
N = 160;
M = 5;
delta_freq = (freq_f - freq_i)/N;

xs = linspace(freq_i + delta_freq, freq_f, N);

param = 'Power';
p_unit = 'dBm';
for run_idx = 1:60
    for prm_idx = -10:-2
        param_idx = prm_idx * 5;
        fprintf(SRS384, ['AMPR ',num2str(param_idx)]); %set power in dBm

        dfilename = [param, num2str(param_idx), '_', num2str(run_idx),'.csv'];
        dfile = fopen([ddir, dfilename],'w+');

        ys = zeros(N,1);
        fig = figure(abs(param_idx));
        line_handle = plot(xs, ys, '.');
        title([param, param_idx, p_unit, 'run',num2str(run_idx)]);
        
        fprintf(dfile, '%s\r\n', ['GHz,Cts/',num2str(dt),',std']);
        for j = 1:M
            for i = 1:N
                ghz = num2str(freq_i + delta_freq * i);
                str = ['FREQ ', ghz, 'e9'];
                fprintf(SRS384, str);
                DAQmxStartTask(hPulse);
                DAQmxStartTask(hCounter);
                arr = diff(DAQmxFunctionPool('ReadCounter', hCounter, samps, (samps + 1) * dt));
                DAQmxStopTask(hCounter);
                DAQmxStopTask(hPulse);
                u = mean(arr);
                if (u == 0)
                    [hPulse, hCounter] = resetDAQ_ODMR(Vx, Vy, sampfreq, samps)
                    DAQmxStartTask(hPulse);
                    DAQmxStartTask(hCounter);
                    arr = diff(DAQmxFunctionPool('ReadCounter', hCounter, samps, (samps + 1) * dt));
                    DAQmxStopTask(hCounter);
                    DAQmxStopTask(hPulse);
                    u = mean(arr)
                end
                ys(i) = ys(i) + u / M;
                fprintf(dfile, '%s\r\n', [ghz, ',', num2str(u),',',num2str(std(arr))]);
            end
            if (j == 1) %mod(j,5) == 0 || 
                j
                if isempty(line_handle)
                    line_handle = plot(xs, ys);
                end
                set(line_handle, 'Ydata', ys);
            end
        end
        set(line_handle, 'Ydata', ys);
%         param_idx
        fclose(dfile);
    end
    run_idx
end

'done'
%% Close and clean up
DAQmxClearTask(hCounter);
DAQmxClearTask(hPulse);
fclose(SRS384);
delete(SRS384)

function [hPulse, hCounter] = resetDAQ_ODMR(vx, vy, sampfreq, samps)

    DAQmxResetDevice('Dev1')
    %Fix Galvos in base position
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo x'), vx);
    DAQmxFunctionPool('WriteVoltage',PortMap('Galvo y'), vy);
    
    [~, hPulse] = DigInfPulse(sampfreq, 0.5);
    [~, hCounter] = DAQmxFunctionPool('SetCounter', samps); 
end