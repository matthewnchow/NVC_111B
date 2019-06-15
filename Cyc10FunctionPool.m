% Control script for custom FPGA Pulse Generator. 
% 8 channel output, (only 4 supported with line drivers at the moment)
% Programmable initial state
% Up to 32 'Edges' that invert a channel's output
% Hardware: Cyclone 10 LP Evaluation Kit with peripherals for:
%                   line driving (IL711)
%                   RS232 to UART level translation (MAX232)
%               Pulse clock: 100 MHz -> can be faster if replace IL711's
%                   Times in Matlab in seconds, converted here to clock
%                   cycles
% Written by Matt Chow 5/26/19. matthewnchow@berkeley.edu
%

function varargout = Cyc10FunctionPool(varargin)
    switch varargin{1}
        case 'ObjHand' % global C10 ?? Handle this as any other instrument,
                        % will want init and stuff later
            [varargout{1}] = ObjHand();
%% Primitive functions of FPGA
        case 'SetPeriod'
            [varargout{1}] = SetPeriod(varargin{2});
        case 'SetOuterPeriod'
            [varargout{1}] = SetOuterPeriod(varargin{2});
        case 'SetState0'
            [varargout{1}] = SetState0(varargin{2});
        case 'SetEdge'
            [varargout{1}] = SetEdge(varargin{2}, varargin{3}, ...
                    varargin{4}, varargin{5}, varargin{6});            
        case 'ClearAll'
            [varargout{1}] = ClearAll();
        case 'SoftTrig'
            [varargout{1}] = SoftTrig();
        case 'SetDPeriod'
            [varargout{1}] = SetDPeriod(varargin{2});            
%         case 'Print' %Extensions
%         case 'WriteSPI'

%% Functions built of primitives 
        case 'CWChannels'
            [varargout{1}] =  CWChannels(varargin{2});
        case 'RabiConfigC10'
            [varargout{1}] = RabiConfigC10(varargin{2}, varargin{3},...
                varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8});
        case 'PumpTime'
            [varargout{1}] = PumpTime(varargin{2}, varargin{3},... % steps, dt, AOMdelay, counttime, depump time, longpump
                varargin{4}, varargin{5}, varargin{6}, varargin{7});    
        case 'CWESR'
            [varargout{1}] = CWESR(varargin{2}, varargin{3}, varargin{4});
      
        case 'LoadSequence'
            [varargout{1}] = LoadSequence(varargin{2});
        otherwise
            [varargout{1}] = 'I do not get it, Cyc10';
    end
end

function status = CWChannels(onchannels)

    status = ClearAll();
    status = status + SetPeriod(0.001);
    status = status + SetState0(onchannels);

end

function [status] = LoadSequence(seq)
    status = ClearAll();
    status = status + SetState0(seq.state0);

    edID = 0;
    % loop through all the channels and program each edge
    period = 0;
    for i=1:numel(seq.ch)
        chID = seq.ch(i).ID;
        for j=1:length(seq.ch(i).x)
            x = seq.ch(i).x(j) * 10^-9; %Convert to seconds
            dx = seq.ch(i).dx(j) * 10^-9;
            period = max(period, x + max(seq.N * dx, 0));
            %SetEdge(edID, enable, chID, x, dx)
            status = status + SetEdge(edID, 1, chID, x, dx);
            edID = edID + 1;
            if status < 0
                return
            end
        end
    end
    cushion = 500 * 10^-9;
    status = status + SetPeriod(period + cushion);
    status = status + SetOuterPeriod(seq.N);
    global gmSEQ
    gmSEQ.period = period + cushion;
    % Add status = status + SetDPeriod(dperiod)?
end

%% LoadYaoSequence Abandoned by MChow 6/6/19
function [status] = LoadYaoSequence(SEQ)
    global gmSEQ;
    
    status = ClearAll();
    status = status + SetState0(0);
    % Given a SEQ structure, this function computes the delay times for all PB
    % channels and sends the instructions to Cyclone 10 FPGA board
    edID = 0;
    % loop through all the channels and program each edge
    period = 0;
    for i=1:numel(SEQ.CHN)
        chID = SEQ.CHN(i).PBN;
        risefallidx = 0;
        for j=1:SEQ.CHN(i).NRise
            x = (SEQ.CHN(i).T(j) + SEQ.CHN(i).Delays(risefallidx + 1)) * 10^-9; %Delays is 2 element array... rise and fall?
            dx = SEQ.CHN(i).DT(j) / SEQ.NSweepParam * 10^-9;
%             steps = 
            period = max(period, x + max(SEQ.NSweepParam * dx, 0));
            %SetEdge(edID, enable, chID, x, dx)
            status = status + SetEdge(edID, 1, chID, x, dx);
            edID = edID + 1;
            risefallidx = ~risefallidx;
        end
    end
    pad = 500 * 10^-9;
    status = status + SetPeriod(period + pad);
    status = status + SetOuterPeriod(SEQ.Repeat);
    % Add status = status + SetDPeriod(dperiod)?
end

function ohand = ObjHand()
    ohand = instrfind('Type', 'serial', 'Port', PortMap('C10'), 'Tag', '');
    if isempty(ohand)
        ohand = serial(PortMap('C10'));
        ohand.BaudRate = 9600;
        ohand.Timeout = 2.0;
        set(ohand, 'DataTerminalReady', 'off');
        set(ohand, 'RequestToSend', 'off');
        set(ohand, 'DataBits', 8);
        ohand.StopBits = 1;
    else
        fclose(ohand);
        ohand = ohand(1);
    end
    fopen(ohand);
end

function stat = SetPeriod(period)
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('period'), 'uint8')
        fwrite(C10, int32(sec2clk(period)), 'int32');
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('period')
            stat = 0;
            return
        end
    end
    stat = -1;
end

function stat = SetDPeriod(dperiod)
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('dperiod'), 'uint8')
        fwrite(C10, int32(sec2clk(dperiod)), 'int32');
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('period')
            stat = 0;
            return
        end
    end
    stat = -1;
end

function stat = SetOuterPeriod(period)
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('outer_period'), 'uint8')
        fwrite(C10, int32(period), 'int32');
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('outer_period')
            stat = 0;
            return
        end
    end
    stat = -1;
end

function stat = SetState0(state)
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('state0'), 'uint8')
        fwrite(C10, uint8(state), 'uint8');
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10,1);
        if ping == C10Const('state0')
            stat = 0;
            return
        end
    end
    stat = -1;
end

function stat = SetEdge(edID, enable, chID, x, dx)
    C10 = ObjHand();
    lastbyte = uint8(enable * 8 + chID);
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('ed'), 'uint8')
        fwrite(C10, uint8(edID), 'uint8'); % Edge ID 
        fwrite(C10, int32(sec2clk(x)), 'int32'); % x
        fwrite(C10, int32(sec2clk(dx)), 'int32'); % dx
        fwrite(C10, uint8(lastbyte), 'uint8'); % Channel 0, enabled on
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('ed')
            stat = 0;
            return 
        end
    end
    stat = -1;
end

function stat = ClearAll()
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('clr'), 'uint8')
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('clr')
            stat = 0;
            return
        end
    end
    stat = -1;    
end 

function stat = SoftTrig()
    C10 = ObjHand();
    for i = 1:C10Const('write_attempts')
        fwrite(C10, C10Const('soft_trig'), 'uint8')
        fwrite(C10, C10Const('CR'));
        fwrite(C10, C10Const('NL'));
        ping = fread(C10, 1);
        if ping == C10Const('soft_trig')
            stat = 0;
            return
        end
    end
    stat = -1;    
end 

%% Sequences


%%Simple Rabi seq
function stat = CWESR(steps, counttime, pumptime)
    % Note that counttime must be >= 5us
    % trig in own channel
    % Chan 0 = laser
    % Chan 1 = MW switch
    % Chan 2 = Gates for counters 1 & 2 (signal and reference)

    %% Rabi Parameters
    counterpw = 30 * 10^-9; %30 ns pulse time for counter edges
    trigwait = 5*10^-6;
    minperiod = trigwait + 2 * counttime + pumptime + counterpw;
    
    %% Set up FPGA
    stat = Cyc10FunctionPool('ClearAll');
    stat = stat + Cyc10FunctionPool('SetPeriod', minperiod);
    stat = stat + Cyc10FunctionPool('SetOuterPeriod', steps);

% Note, DAQ counter uses only the rising edges, so 2 rising edges for start
% and stop must be used
%     Cyc10FunctionPool('SetEdge', edID, enable, chID, x, dx)

    state0 = uint8(3); %0000_0000 start laser and MW on
    stat = stat + Cyc10FunctionPool('SetState0', state0);
    edges = {{1, 1, trigwait + counttime + counterpw, 0}, ... %Turn MW off
             {1, 2, trigwait, 0},...  %Turn on data counter
             {1, 2, trigwait + counterpw, 0},...%Falling edge for ^^
             {1, 2, trigwait + counttime, 0},...  %Turn off data counter
             {1, 2, trigwait + counttime + counterpw, 0},... % Falling edge for ^^
             {1, 2, trigwait + counttime + pumptime, 0},...  %Turn on ref counter
             {1, 2, trigwait + counttime + pumptime + counterpw, 0},... % Falling edge for ^^
             {1, 2, trigwait + 2 * counttime + pumptime, 0},...  %Turn off ref counter
             {1, 2, trigwait + 2 * counttime + pumptime + counterpw, 0} % Falling edge for ^^
             };
   
    for i = 1:(length(edges))
        ed = edges{i};
        stat = stat + Cyc10FunctionPool('SetEdge', i-1,ed{1}, ed{2}, ed{3}, ed{4});
    end
end    

%%Simple Rabi seq
function stat = DarkVaryConfigC10(steps, dt, counttime, pumptime, cushion, AOMdelay)
    % Note that pumptime must be greater than counttime
    % Chan 0 = trig
    % Chan 1 = laser
    % Chan 2 = MW
    % Chan 3 = Gates for counters 1 & 2 (signal and reference)

    %% Rabi Parameters
    counterpw = 30 * 10^-9; %30 ns pulse time for counter edges
    AOMrisefall = 500*10^-9; %Measured as 400, 500 just to be safe
%     trigwait = 0.75*10^-3;
    trigwait = 2*10^-6; %Testing HW trigger
    minperiod = trigwait + cushion + steps * dt + pumptime + AOMrisefall...
        + counttime + counterpw + cushion * 2 + AOMdelay; % This is really slow, but it has to be for slowTrig

    %% Set up FPGA
    stat = Cyc10FunctionPool('ClearAll');
    stat = stat + Cyc10FunctionPool('SetPeriod', minperiod);
    stat = stat + Cyc10FunctionPool('SetOuterPeriod', steps);

    state0 = uint8(1); %0000_0000 start laser on
    stat = stat + Cyc10FunctionPool('SetState0', state0);
    edges = {{1, 0, trigwait, 0}, ... % Turn Laser off
             {1, 1, trigwait + cushion / 2, 0}, ... %Turn MW on
             {1, 1, trigwait + cushion / 2, dt}, ... %Turn MW off, grows by dt each time
             {1, 0, trigwait + cushion + steps * dt, 0}, ... % Turn Laser on
             {1, 2, trigwait + cushion + steps * dt + AOMdelay, 0},...  %Turn on data counter
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + counterpw, 0},...%Falling edge for ^^
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + counttime, 0},...  %Turn off data counter
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + counttime + counterpw, 0},... % Falling edge for ^^
             {1, 0, trigwait + cushion + steps * dt + pumptime, 0}, ... % Turn Laser off
             {1, 0, trigwait + cushion + steps * dt + pumptime + AOMrisefall, 0}, ... % Turn Laser back on             
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + pumptime + AOMrisefall, 0},...  %Turn on ref counter
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + pumptime + AOMrisefall + counterpw, 0},... % Falling edge for ^^
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + pumptime + AOMrisefall + counttime, 0},...  %Turn off ref counter
             {1, 2, trigwait + cushion + steps * dt + AOMdelay + pumptime + AOMrisefall + counttime + counterpw, 0} % Falling edge for ^^
             };
    
    for i = 1:(length(edges))
        ed = edges{i};
        stat = stat + Cyc10FunctionPool('SetEdge', i-1,ed{1}, ed{2}, ed{3}, ed{4}); %i-1+badedcount
    end
end    



%%Simple Rabi seq
function stat = RabiConfigC10(steps, dt, counttime, pumptime, AOMrisefall, AOMdelay, MWdelay)
    % Note that pumptime must be greater than counttime
    % Chan 0 = trig
    % Chan 1 = laser
    % Chan 2 = MW
    % Chan 3 = Gates for counters 1 & 2 (signal and reference)

    %% Rabi Parameters
    counterpw = 30 * 10^-9; %30 ns pulse time for counter edges
%     AOMrisefall = 800*10^-9; %Measured as 400, 800 just to be safe
%     trigwait = 0.75*10^-3;
%     MWdelay = 100 *10^-9;
    pad = 1 * 10^-6;
    trigwait = 2*10^-6; %Testing HW trigger
    minperiod = trigwait + AOMrisefall + steps * dt + pumptime + AOMrisefall...
        + counttime + counterpw + pad;

    %% Set up FPGA
    stat = Cyc10FunctionPool('ClearAll');
    stat = stat + Cyc10FunctionPool('SetPeriod', minperiod);
    stat = stat + Cyc10FunctionPool('SetOuterPeriod', steps);

% Note, DAQ counter uses only the rising edges, so 2 rising edges for start
% and stop must be used
%     Cyc10FunctionPool('SetEdge', edID, enable, chID, x, dx)
%     Should do this from file later

    state0 = uint8(1); %0000_0000 start laser on
    stat = stat + Cyc10FunctionPool('SetState0', state0);
    edges = {{1, 0, trigwait, 0}, ... % Turn Laser off
             {1, 1, trigwait + AOMrisefall, 0}, ... %Turn MW on
             {1, 1, trigwait + AOMrisefall, dt}, ... %Turn MW off, grows by dt each time
             {1, 0, trigwait + AOMrisefall + steps * dt + MWdelay, 0}, ... % Turn Laser on
             {1, 2, trigwait + AOMrisefall + steps * dt + MWdelay + AOMdelay, 0},...  %Turn on data counter
             {1, 2, trigwait + AOMrisefall + steps * dt + MWdelay + AOMdelay + counterpw, 0},...%Falling edge for ^^
             {1, 2, trigwait + AOMrisefall + steps * dt + MWdelay + AOMdelay + counttime, 0},...  %Turn off data counter
             {1, 2, trigwait + AOMrisefall + steps * dt + MWdelay + AOMdelay + counttime + counterpw, 0},... % Falling edge for ^^
             {1, 0, trigwait + AOMrisefall + steps * dt + pumptime, 0}, ... % Turn Laser off
             {1, 0, trigwait + AOMrisefall + steps * dt + pumptime + AOMrisefall, 0}, ... % Turn Laser back on             
             {1, 2, trigwait + AOMrisefall + steps * dt + AOMdelay + pumptime + AOMrisefall, 0},...  %Turn on ref counter
             {1, 2, trigwait + AOMrisefall + steps * dt + AOMdelay + pumptime + AOMrisefall + counterpw, 0},... % Falling edge for ^^
             {1, 2, trigwait + AOMrisefall + steps * dt + AOMdelay + pumptime + AOMrisefall + counttime, 0},...  %Turn off ref counter
             {1, 2, trigwait + AOMrisefall + steps * dt + AOMdelay + pumptime + AOMrisefall + counttime + counterpw, 0} % Falling edge for ^^
             };
    
    for i = 1:(length(edges))
        ed = edges{i};
        stat = stat + Cyc10FunctionPool('SetEdge', i-1,ed{1}, ed{2}, ed{3}, ed{4}); %i-1+badedcount
    end
end    


%% Internals
function cycles = sec2clk(sec)
    clk_rate = 100 * 10 ^ 6;
    cycles = int32(sec * clk_rate);
end

function const = C10Const(what)
    switch what
        case 'state0'
            const = uint8(1);
        case 'period'
            const = uint8(2);
        case 'outer_period'
            const = uint8(3);
        case 'ed'
            const = uint8(4);
        case 'clr'
            const = uint8(5);
        case 'print'
            const = uint8(6);
        case 'soft_trig'
            const = uint8(7);
        case 'dperiod'
            const = uint8(8); %Not functional yet            
        case 'CR'
            const = 13;
        case 'NL'
            const = 10;
        case 'write_attempts'
            const = 3;
    end
end
