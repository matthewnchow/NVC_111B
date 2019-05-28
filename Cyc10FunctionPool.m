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
            [varargout{1}] = ObjHand(varargin{2});
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
%         case 'Print' %Extensions
%         case 'WriteSPI'

        case 'RabiConfigC10'
            [varargout{1}] = RabiConfigC10(varargin{2}, varargin{3},...
                varargin{4}, varargin{5}, varargin{6});
    end
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
function stat = RabiConfigC10(steps, dt, counttime, pumptime, cushion)
    % Note that pumptime must be greater than counttime
    % Chan 0 = trig
    % Chan 1 = laser
    % Chan 2 = MW
    % Chan 3 = Gates for counters 1 & 2 (signal and reference)

    %% Rabi Parameters
    state0 = uint8(0); %0000_0000 start all off
    counterpw = 30 * 10^-9; %30 ns pulse time for counter edges 
%     steps = 100; %Data points to take for Rabi cycle
    
    %% Set up FPGA
    stat = Cyc10FunctionPool('ClearAll');
    stat = stat + Cyc10FunctionPool('SetPeriod', pumptime + steps * dt...
            + counttime + cushion * 20);
    stat = stat + Cyc10FunctionPool('SetOuterPeriod', steps);
    stat = stat + Cyc10FunctionPool('SetState0', state0);

% Note, DAQ counter uses only the rising edges, so 2 rising edges for start
% and stop must be used
    % Cyc10FunctionPool('SetEdge', edID, enable, chID, x, dx)
    % Should do this from file later
    edges = {{1, 1, cushion, 0}, ... %Turn MW on
             {1, 1, cushion, dt}, ... %Turn MW off, grows by dt each time
             {1, 0, cushion + steps * dt, 0}, ... % Turn Laser on
             {1, 2, cushion + steps * dt + cushion, 0},...  %Turn on data counter
             {1, 2, cushion + steps * dt + cushion + counterpw, 0},...%Falling edge for ^^
             {1, 2, cushion + steps * dt + cushion + counttime, 0},...  %Turn off data counter
             {1, 2, cushion + steps * dt + cushion + counttime + counterpw, 0},... % Falling edge for ^^
             {1, 2, cushion + steps * dt + cushion + pumptime - counttime, 0},...  %Turn on ref counter
             {1, 2, cushion + steps * dt + cushion + pumptime - counttime + counterpw, 0},... % Falling edge for ^^
             {1, 2, cushion + steps * dt + cushion + pumptime, 0},...  %Turn off ref counter
             {1, 2, cushion + steps * dt + cushion + pumptime + counterpw, 0} % Falling edge for ^^
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
        case 'CR'
            const = 13;
        case 'NL'
            const = 10;
        case 'write_attempts'
            const = 3;
    end
end
