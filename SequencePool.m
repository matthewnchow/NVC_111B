function varargout = SequencePool(varargin)
global gmSEQ
if isfield(gmSEQ,'CHN')&& ~isequal(varargin{1},'PBDictionary')
    gmSEQ=rmfield(gmSEQ,'CHN');
end
if isfield(gmSEQ,'bLiO') && ~isequal(varargin{1},'PBDictionary')
    gmSEQ=rmfield(gmSEQ,'bLiO');
end
switch varargin{1}
    case 'PopulateSeq'
        varargout{1} = PopulateSeq();
        return
    case 'Rabi'
        varargout{1} = Rabi();
        return
    case 'PumpTime'
        varargout{1} = PumpTime();
        return
    case 'AOM Delay'
        varargout{1} = AOMDelay();
        return
    case 'CW ESR'
        varargout{1} = CWESR();
        return
    case 'Ramsey'
        varargout{1} = Ramsey();
        return
    case 'Hahn'
        varargout{1} = Hahn();
        return       
    case 'TPrimeEcho'
        varargout{1} = TPrimeEcho();
        return              
%% Unused or need to be rewritten for C10 pulse blaster        
    case 'Echo 4 counters'
        Echo4();
    case 'ESR'
        ESR();
    case 'T1'
        T1();
    case 'ODMR'
        ODMR();
    case 'CPMG-N'
        CPMGN();
    case 'T1JC'
        T1JC();
    case 'PBDictionary'
        varargout{1}=PBDictionary(varargin{2});
        return
    case 'CtrDur'
        CtrDur();
    case 'Counter calibration'
        CtrCal();
    case 'Pulsed ESR'
        PESR();
    case 'Rabi Dark'
        RabiDark();
    case 'ODMR Dark'
        ODMRDark();
    case 'Select Sequence'
        varargout{1} = {};
        return
end
varargout{1} = {};

function StrL = PopulateSeq
StrL{1} = 'Select Sequence';
StrL{numel(StrL)+1}='Rabi';
StrL{numel(StrL)+1}='PumpTime';
StrL{numel(StrL)+1}='AOM Delay';
StrL{numel(StrL)+1}='CW ESR';
StrL{numel(StrL)+1}='Ramsey';
StrL{numel(StrL)+1}='Hahn';
StrL{numel(StrL)+1}='TPrimeEcho';


%Should add these:
% StrL{numel(StrL)+1}='Counter calibration';
% StrL{numel(StrL)+1}='Pulsed ESR';
% Variable Echo
% StrL{numel(StrL)+1}='T1';

%Prob don't need these ~MChow
% StrL{numel(StrL)+1}='CPMG-T'; 
% StrL{numel(StrL)+1}='CPMG-N';
% StrL{numel(StrL)+1}='Echo 4 counters';
% StrL{numel(StrL)+1}='T1JC';
% StrL{numel(StrL)+1}='CtrDur';
% StrL{numel(StrL)+1}='Rabi Dark';
% StrL{numel(StrL)+1}='ODMR Dark';
% StrL{numel(StrL)+1}='ODMR';
% StrL{numel(StrL)+1}='ESR';


%% gmSEQ structure for C10:
% state0 = state0
% N = outerperiod
% dperiod = dperiod
% ch = channels (analagous to Yao CHN)
%     ch.ID
%     ch.x = time value of edge
%     ch.dx = difference in time value between sequences
%     ch.delay = delay for edges on this channel


function seq = Rabi()
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='0';
% gSG.bModSrc='External'; %Could do external pulse modulation of SRS384

%%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
seq.To = gmSEQ.To;
seq.From = gmSEQ.From;
seq.CtrGateDur = gmSEQ.CtrGateDur;
seq.N = gmSEQ.N;
%%%%%%%%

AOMrisefall = 500;
trigwait = 1000; %Make sure there is time for the DAQ to always get the first StartStop trigger first
counterpw = 50; %Counter start stop pulse width
stepsize = (seq.To - seq.From)/(seq.N - 1);
pumptime = gmSEQ.PumpT;
AOMcounterdelay = gmSEQ.AOMCounterDelay;
MWdelay = gmSEQ.MWSwitchDelay;

%%%%% Fixed sequence length %%%%%%
seq.state0 = 0; %% start with all off
%%%%% MW turn on first, turn off at stepsize * number of steps %%%
seq.ch(1).ID = PBDictionary('MWswitch');
seq.ch(1).x = [(trigwait) (seq.From + trigwait)];
seq.ch(1).dx = [0 stepsize];
%%%%% Turn on laser after MW is done %%%%% MW delay hard coded in 6/6/19
seq.ch(numel(seq.ch) + 1).ID = PBDictionary('AOM');
seq.ch(numel(seq.ch)).x  = [(seq.To + trigwait + MWdelay)...
    (seq.To + trigwait + pumptime) (seq.To + trigwait + pumptime...
     + AOMrisefall * 2)];
seq.ch(numel(seq.ch)).dx = [0 0 0];
%%%%% Counter edges 
seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
seq.ch(numel(seq.ch)).x  = [...
    (seq.To + trigwait + MWdelay + AOMcounterdelay)...
    (seq.To + trigwait + MWdelay + AOMcounterdelay + counterpw) ...
    (seq.To + trigwait + MWdelay + AOMcounterdelay + seq.CtrGateDur)...
    (seq.To + trigwait + MWdelay + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
    (seq.To + trigwait + pumptime + AOMrisefall * 2 + AOMcounterdelay)...
    (seq.To + trigwait + pumptime + AOMrisefall * 2 + AOMcounterdelay + counterpw)...
    (seq.To + trigwait + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur)...
    (seq.To + trigwait + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw)...
    ];
seq.ch(numel(seq.ch)).dx = [0, 0, 0, 0, 0, 0, 0, 0];

%%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
%%%%% sequence -- Note that this could alternately be compensated for with
%%%%% pad in LoadSequence (in Cyc10FunctionPool)
seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
seq.ch(numel(seq.ch)).x = [2*trigwait + seq.To + trigwait + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw];
seq.ch(numel(seq.ch)).dx = [0];

%%%%%%%% Global variables give me a headache part 2...~MChow
gmSEQ.ch = seq.ch;
%%%%%%%%
% ApplyNoDelays(); %FIXME, make the delays not hardcoded like they are here
%%%%% For variable length, program C10 for dperiod %%

function seq = PumpTime() %Meant to take a pretty long MW on time (from pi box) ~MChow
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=1;
    gSG.bMod='0';

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.To = gmSEQ.To;
    seq.From = gmSEQ.From;
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    seq.N = gmSEQ.N;
    seq.pi = gmSEQ.pi;
    %%%%%%%%

    AOMrisefall = 500; %This is pretty generous, it's more like 100ns
    trigwait = 1000; %Make sure there is time for the DAQ to always get the first StartStop trigger first
    counterpw = 50; %Counter start stop pulse width
    stepsize = (seq.To - seq.From)/(seq.N - 1);
    AOMcounterdelay = gmSEQ.AOMCounterDelay; %Negative since fluorescence (AOM on signal) appears to arrive first at the DAQ
    MWdelay = gmSEQ.MWSwitchDelay * 2; % times 2 for a bit extra cushion

    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 2^PBDictionary('MWswitch'); %% start with MW on
    %%%%% MW turn on first, depump state %%%
    seq.ch(1).ID = PBDictionary('MWswitch');
    seq.ch(1).x = [...
        (seq.pi)...
        (seq.pi + MWdelay + seq.To + AOMrisefall * 2 + seq.CtrGateDur * 2)...
        (seq.pi + MWdelay + seq.To + AOMrisefall * 2 + seq.CtrGateDur * 2 + seq.pi)...
        ];
    seq.ch(1).dx = [0 0 0];
    %%%%% Turn on laser after MW is done %%%%% MW delay hard coded in 6/6/19
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('AOM');
    seq.ch(numel(seq.ch)).x  = [...
        (seq.pi + MWdelay + seq.From)...
        (seq.pi + MWdelay + seq.From)...
        (seq.pi + MWdelay + seq.To + AOMrisefall * 2)...
        (seq.pi + MWdelay + seq.To + AOMrisefall * 2 + seq.CtrGateDur * 2)...
        (seq.pi + MWdelay + seq.To + AOMrisefall * 2 + seq.CtrGateDur * 2 ...
            + AOMrisefall + seq.pi + MWdelay)...        
        ];
    seq.ch(numel(seq.ch)).dx = [0 stepsize 0 0 0];
    %%%%% Counter edges %note that the reference is counted first here
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    %% counter edges for each counting laser edge = [...
%         (laseron + AOMcounterdelay)...
%         (laseron + AOMcounterdelay + counterpw) ...
%         (laseron + AOMcounterdelay + seq.CtrGateDur)...
%         (laseron + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
%          ];
    seq.ch(numel(seq.ch)).x  = [...
        (seq.ch(numel(seq.ch) - 1).x(3) + AOMcounterdelay)...
        (seq.ch(numel(seq.ch) - 1).x(3) + AOMcounterdelay + counterpw) ...
        (seq.ch(numel(seq.ch) - 1).x(3) + AOMcounterdelay + seq.CtrGateDur)...
        (seq.ch(numel(seq.ch) - 1).x(3) + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
        (seq.ch(numel(seq.ch) - 1).x(5) + AOMcounterdelay)...
        (seq.ch(numel(seq.ch) - 1).x(5) + AOMcounterdelay + counterpw) ...
        (seq.ch(numel(seq.ch) - 1).x(5) + AOMcounterdelay + seq.CtrGateDur)...
        (seq.ch(numel(seq.ch) - 1).x(5) + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
        ];
    seq.ch(numel(seq.ch)).dx = [0, 0, 0, 0, 0, 0, 0, 0];

    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [trigwait + 2 * seq.pi + 2 * MWdelay + seq.To + AOMrisefall * 4 + 4 * seq.CtrGateDur + counterpw];
    seq.ch(numel(seq.ch)).dx = [0];

    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;

    
function seq = AOMDelay()
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=1;
    gSG.bMod='0';

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.To = gmSEQ.To;
    seq.From = gmSEQ.From;
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    seq.N = gmSEQ.N;
    %%%%%%%%

    trigwait = 1000; %Important since the samp clk will be delayed quite a bit relative to the arm start trig
    per = max(10000, trigwait + seq.To + seq.CtrGateDur); %Make sure there is time for the DAQ to always get the first StartStop trigger first
    counterpw = 50; %Counter start stop pulse width
    stepsize = (seq.To - seq.From)/(seq.N - 1);

    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 0; %% start with all off
    %%%%% Sweep laser turn on time through sequence
    seq.ch(1).ID = PBDictionary('AOM');
    seq.ch(numel(seq.ch)).x  = [seq.From];
    seq.ch(numel(seq.ch)).dx = [stepsize];
    %%%%% Counter edges %note that the reference is counted first here
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    seq.ch(numel(seq.ch)).x  = [...
        ((seq.To + seq.From)/2 - seq.CtrGateDur/2)...
        ((seq.To + seq.From)/2 - seq.CtrGateDur/2 + counterpw) ...
        ((seq.To + seq.From)/2 + seq.CtrGateDur/2)...
        ((seq.To + seq.From)/2 + seq.CtrGateDur/2 + counterpw) ...
        ];
    seq.ch(numel(seq.ch)).dx = [0, 0, 0, 0];
    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [per];
    seq.ch(numel(seq.ch)).dx = [0];

    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;

function seq = CWESR()
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=0; %This will be the parameter that is swept each time
    gSG.bMod='0'; %Switching on and off with a MiniCircuits switch atm

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    CtrMWdelay = gmSEQ.MWSwitchDelay; 
                %This isn't quite right...fine for now
                %Should really have a delay time for each thing relative to
                %   the FPGA
    pumptime = gmSEQ.PumpT;
    %%%%%%%%

    trigwait = 1000;
    counterpw = 50; %Counter start stop pulse width

    per = 2 * (trigwait + seq.CtrGateDur + 2 * CtrMWdelay + counterpw) + pumptime;

    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 2^PBDictionary('MWswitch') + 2^PBDictionary('AOM'); %% start with laser and MW on, 
                            %leave laser on forever, switch MW off at per/2
    seq.N = 1; %Sets the outer period to 1
    
    %%%%% Switch MW off after done counting
	seq.ch(1).ID = PBDictionary('MWswitch');
    seq.ch(numel(seq.ch)).x = [trigwait + seq.CtrGateDur + 2 * CtrMWdelay];
    seq.ch(numel(seq.ch)).dx = [0];                            
    %%%%% Counter edges, first while MW on, second while MW is off
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    seq.ch(numel(seq.ch)).x  = [...
        (trigwait + CtrMWdelay)...
        (trigwait + CtrMWdelay + counterpw) ...
        (trigwait + CtrMWdelay + seq.CtrGateDur)...
        (trigwait + CtrMWdelay + seq.CtrGateDur + counterpw) ...
        (per/2 + pumptime / 2 + CtrMWdelay)...
        (per/2 + pumptime / 2 + CtrMWdelay + counterpw) ...
        (per/2 + pumptime / 2 + CtrMWdelay + seq.CtrGateDur)...
        (per/2 + pumptime / 2 + CtrMWdelay + seq.CtrGateDur + counterpw) ...
        ];
    seq.ch(numel(seq.ch)).dx = [0, 0, 0, 0, 0, 0, 0, 0];
    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [per];
    seq.ch(numel(seq.ch)).dx = [0];
    
    %Sweep over frequency
    seq.SweepParam = ((0:gmSEQ.N-1)*((gmSEQ.To-gmSEQ.From)/gmSEQ.N)) + gmSEQ.From;
    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;            
    gmSEQ.SweepParam = seq.SweepParam;
    gmSEQ.NSweepParam = gmSEQ.N;

function seq = Ramsey() %pi/2 -T - pi/2      
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=1;
    gSG.bMod='0';
    % gSG.bModSrc='External'; %Could do external pulse modulation of SRS384

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.To = gmSEQ.To; % Thing being swept is wait time
    seq.From = gmSEQ.From;
    seq.N = gmSEQ.N;
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    seq.pi = gmSEQ.pi;
    %%%%%%%%

    AOMrisefall = 600;
    trigwait = 1000; %Make sure there is time for the DAQ to always get the first StartStop trigger first
    counterpw = 50; %Counter start stop pulse width
    stepsize = (seq.To - seq.From)/(seq.N - 1);
    pumptime = gmSEQ.PumpT;
    AOMcounterdelay = gmSEQ.AOMCounterDelay;
    MWdelay = gmSEQ.MWSwitchDelay;
    
    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 0; %% start with all off
    %%%%% MW, two pi/2 pulses, with variable wait time between %%%
    seq.ch(1).ID = PBDictionary('MWswitch');
    seq.ch(1).x = [trigwait,...
                   trigwait + seq.pi / 2,...
                   trigwait + seq.pi / 2 + seq.From, ...
                   trigwait + seq.pi + seq.From, ...
                   ];
    seq.ch(1).dx = [0, 0, stepsize, stepsize];
    %%%%% Turn on laser after MW is done %%%%% MW delay hard coded in 6/6/19
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('AOM');
    seq.ch(numel(seq.ch)).x  = [trigwait + seq.pi + seq.To + MWdelay,...
        trigwait + seq.pi + seq.To + MWdelay + pumptime,...
        trigwait + seq.pi + seq.To + MWdelay + pumptime + AOMrisefall * 2 ... %Leave laser on till period ends
        ];
    seq.ch(numel(seq.ch)).dx = [0 0 0];
    %%%%% Counter edges 
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    seq.ch(numel(seq.ch)).x  = [...
        (trigwait + seq.pi + seq.To + MWdelay + AOMcounterdelay)...
        (trigwait + seq.pi + seq.To + MWdelay + AOMcounterdelay + counterpw) ...
        (trigwait + seq.pi + seq.To + MWdelay + AOMcounterdelay + seq.CtrGateDur)...
        (trigwait + seq.pi + seq.To + MWdelay + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
        (trigwait + seq.pi + seq.To + MWdelay + pumptime + AOMrisefall * 2 + AOMcounterdelay)...
        (trigwait + seq.pi + seq.To + MWdelay + pumptime + AOMrisefall * 2 + AOMcounterdelay + counterpw)...
        (trigwait + seq.pi + seq.To + MWdelay + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur)...
        (trigwait + seq.pi + seq.To + MWdelay + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw)...
        ];
    seq.ch(numel(seq.ch)).dx = [0, 0, 0, 0, 0, 0, 0, 0];

    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [4 * trigwait + seq.pi + seq.To + pumptime...
        + AOMrisefall * 2 + abs(AOMcounterdelay) + seq.CtrGateDur + counterpw];
    seq.ch(numel(seq.ch)).dx = [0];
    
    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;

function seq = Hahn() %pi/2 -T - pi - T'=T - pi/2   %Adding in the dPeriod function to C10 would be good for this
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=1;
    gSG.bMod='0';
    % gSG.bModSrc='External'; %Could do external pulse modulation of SRS384

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.To = gmSEQ.To; % Thing being swept is wait time
    seq.From = gmSEQ.From;
    seq.N = gmSEQ.N;
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    seq.pi = gmSEQ.pi;
    %%%%%%%%

    AOMrisefall = 600;
    trigwait = 1000; %Make sure there is time for the DAQ to always get the first StartStop trigger first
    counterpw = 50; %Counter start stop pulse width
    stepsize = (seq.To - seq.From)/(seq.N - 1);
    pumptime = gmSEQ.PumpT;
    AOMcounterdelay = gmSEQ.AOMCounterDelay;
    MWdelay = gmSEQ.MWSwitchDelay;
    
    lastMW = trigwait + seq.pi / 2 + seq.From + seq.pi + seq.From + seq.pi/2 + MWdelay;
    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 0; %% start with all off
    %%%%% MW, two pi/2 pulses, with variable wait time between %%%
    seq.ch(1).ID = PBDictionary('MWswitch');
    seq.ch(1).x = [trigwait,...
                   trigwait + seq.pi / 2,...
                   trigwait + seq.pi / 2 + seq.From, ...
                   trigwait + seq.pi / 2 + seq.From + seq.pi, ...
                   trigwait + seq.pi / 2 + seq.From + seq.pi + seq.From, ...
                   trigwait + seq.pi / 2 + seq.From + seq.pi + seq.From + seq.pi/2 ...
                   ];
    seq.ch(1).dx = stepsize * [0, 0, 1, 1, 2, 2];
    
    %%%%% Turn on laser after MW is done %%%%% MW delay hard coded in 6/6/19
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('AOM');
    seq.ch(numel(seq.ch)).x  = [lastMW + MWdelay,...
        lastMW + pumptime,...
        lastMW + pumptime + AOMrisefall * 2 ... %Leave laser on till period ends
        ];
    seq.ch(numel(seq.ch)).dx = stepsize * [2, 2, 2];
    %%%%% Counter edges 
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    seq.ch(numel(seq.ch)).x  = [...
        (lastMW + AOMcounterdelay)...
        (lastMW + AOMcounterdelay + counterpw) ...
        (lastMW + AOMcounterdelay + seq.CtrGateDur)...
        (lastMW + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + counterpw)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw)...
        ];
    seq.ch(numel(seq.ch)).dx = stepsize * [2, 2, 2, 2, 2, 2, 2, 2];

    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [4 * trigwait + lastMW + MWdelay + pumptime ...
        + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw + counterpw];
    seq.ch(numel(seq.ch)).dx = stepsize * [2];
    
    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;    
        
    
function seq = TPrimeEcho() %pi/2 - T - pi - T' - pi/2   %Adding in the dPeriod function to C10 would be good for this
    global gmSEQ gSG
    gSG.bfixedPow=1;
    gSG.bfixedFreq=1;
    gSG.bMod='0';
    % gSG.bModSrc='External'; %Could do external pulse modulation of SRS384

    %%%%%%%% Global variables give me a headache...this bit is unnecessary. ~MChow
    seq.To = gmSEQ.To; % Thing being swept is wait time
    seq.From = gmSEQ.From;
    seq.N = gmSEQ.N;
    seq.CtrGateDur = gmSEQ.CtrGateDur;
    seq.pi = gmSEQ.pi;
    %%%%%%%%

    AOMrisefall = 600;
    trigwait = 1000; %Make sure there is time for the DAQ to always get the first StartStop trigger first
    counterpw = 50; %Counter start stop pulse width
    stepsize = (seq.To - seq.From)/(seq.N - 1);
    pumptime = gmSEQ.PumpT;
    AOMcounterdelay = gmSEQ.AOMCounterDelay;
    MWdelay = gmSEQ.MWSwitchDelay;
    T = (seq.To + seq.From)/2;
    lastMW = trigwait + seq.pi / 2 + T + seq.pi + seq.From + seq.pi/2 + MWdelay;
    %%%%% Fixed sequence length %%%%%%
    seq.state0 = 0; %% start with all off
    %%%%% MW, two pi/2 pulses, with variable wait time between %%%
    seq.ch(1).ID = PBDictionary('MWswitch');
    seq.ch(1).x = [trigwait,...
                   trigwait + seq.pi / 2,...
                   trigwait + seq.pi / 2 + T, ...
                   trigwait + seq.pi / 2 + T + seq.pi, ...
                   trigwait + seq.pi / 2 + T + seq.pi + seq.From, ...
                   trigwait + seq.pi / 2 + T + seq.pi + seq.From + seq.pi/2 ...
                   ];
    seq.ch(1).dx = [0, 0, 0, 0, stepsize, stepsize];
    
    %%%%% Turn on laser after MW is done %%%%% MW delay hard coded in 6/6/19
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('AOM');
    seq.ch(numel(seq.ch)).x  = [lastMW + MWdelay,...
        lastMW + pumptime,...
        lastMW + pumptime + AOMrisefall * 2 ... %Leave laser on till period ends
        ];
    seq.ch(numel(seq.ch)).dx = stepsize * [1, 1, 1];
    %%%%% Counter edges 
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('ctr0');
    seq.ch(numel(seq.ch)).x  = [...
        (lastMW + AOMcounterdelay)...
        (lastMW + AOMcounterdelay + counterpw) ...
        (lastMW + AOMcounterdelay + seq.CtrGateDur)...
        (lastMW + AOMcounterdelay + seq.CtrGateDur + counterpw) ...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + counterpw)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur)...
        (lastMW + pumptime + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw)...
        ];
    seq.ch(numel(seq.ch)).dx = stepsize * [1, 1, 1, 1, 1, 1, 1, 1];

    %%%%% Dummy ensures that there is a pad of ~trigwait at the end of the
    %%%%% sequence -- Note that this could alternately be compensated for with
    %%%%% pad in LoadSequence (in Cyc10FunctionPool)
    seq.ch(numel(seq.ch) + 1).ID = PBDictionary('dummy');
    seq.ch(numel(seq.ch)).x = [4 * trigwait + lastMW + MWdelay + pumptime ...
        + AOMrisefall * 2 + AOMcounterdelay + seq.CtrGateDur + counterpw + counterpw];
    seq.ch(numel(seq.ch)).dx = [stepsize];
    
    %%%%%%%% Global variables give me a headache part 2...~MChow
    gmSEQ.ch = seq.ch;    
    
    
%%%%%%%%%%%%%%%%%%%%%%%Below is inconsistent with the rest of the program
function RabiDark
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';
delay=1050;
%%%%% Fixed sequence length %%%%%%
gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(1).NRise=3;
gmSEQ.CHN(1).T=[gmSEQ.To+delay+gmSEQ.pi*2 gmSEQ.To+delay+gmSEQ.pi*2+gmSEQ.readout-(gmSEQ.CtrGateDur+1000)+1000 gmSEQ.To+delay+gmSEQ.pi*2+gmSEQ.readout-(gmSEQ.CtrGateDur+1000)+1000+gmSEQ.CtrGateDur+2000];
if strcmp(gmSEQ.meas,'SPCM')
    gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
elseif strcmp(gmSEQ.meas,'APD')
    gmSEQ.CHN(1).DT=[1000 1000 1000];
end
gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[delay delay+gmSEQ.pi+gmSEQ.m];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.pi gmSEQ.pi];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('BNC645');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[delay+gmSEQ.pi];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.m];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('AOM');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[gmSEQ.CHN(1).T(1) gmSEQ.CHN(1).T(2)];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.readout-(gmSEQ.CtrGateDur+1000) gmSEQ.CtrGateDur+1000];
gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=1000;
ApplyDelays();


function CtrDur
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';

mwlength=1100+gmSEQ.pi;

%%%%% Fixed sequence length %%%%%%
gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=2;
gmSEQ.CHN(1).T=[100 100+gmSEQ.readout+mwlength];
gmSEQ.CHN(1).DT=[gmSEQ.readout gmSEQ.readout];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[gmSEQ.CHN(1).T(1)-100 gmSEQ.CHN(1).T(2)-100];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.m gmSEQ.m];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[gmSEQ.CHN(1).T(2)-gmSEQ.pi-100];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.pi];

ApplyDelays();


function T1JC
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';

mwlength=gmSEQ.m +gmSEQ.pi;

gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=3;
gmSEQ.CHN(1).T=[mwlength mwlength+gmSEQ.readout+mwlength mwlength+gmSEQ.readout+mwlength+gmSEQ.readout+200];
gmSEQ.CHN(1).DT=[gmSEQ.readout gmSEQ.readout gmSEQ.CtrGateDur+200];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=3;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=gmSEQ.CHN(1).T;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[gmSEQ.CHN(1).T(2)-gmSEQ.pi-50];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.pi];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[0];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[10];

ApplyDelays();


function CtrCal
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';

gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(1).NRise=2;
gmSEQ.CHN(1).T=[2000+gmSEQ.m 2000+gmSEQ.m+gmSEQ.readout+1000];
gmSEQ.CHN(1).DT=[12 12];

gmSEQ.CHN(2).PBN=PBDictionary('AOM');
gmSEQ.CHN(2).NRise=2;
gmSEQ.CHN(2).T=[2000 2000+gmSEQ.readout+1000];
gmSEQ.CHN(2).DT=[gmSEQ.readout gmSEQ.readout];

gmSEQ.CHN(3).PBN=PBDictionary('dummy');
gmSEQ.CHN(3).NRise=1;
gmSEQ.CHN(3).T=0;
gmSEQ.CHN(3).DT=12;

gmSEQ.CHN(4).PBN=PBDictionary('I');
gmSEQ.CHN(4).NRise=1;
gmSEQ.CHN(4).T=2000-50-gmSEQ.pi+500;
gmSEQ.CHN(4).DT=gmSEQ.pi;
ApplyNoDelays;

% function Ramsey
% global gmSEQ gSG
% gSG.bfixedPow=1;
% gSG.bfixedFreq=1;
% gSG.bMod='IQ';
% gSG.bModSrc='External';

%%%%% Fixed sequence length no ctr2%%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=1;
% gmSEQ.CHN(1).T=gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi;
% gmSEQ.CHN(1).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).PBN=PBDictionary('ctr1');
% gmSEQ.CHN(2).NRise=1;
% gmSEQ.CHN(2).T=gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout-gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(3).PBN=PBDictionary('I');
% gmSEQ.CHN(3).NRise=2;
% gmSEQ.CHN(3).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m];
% gmSEQ.CHN(3).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2];
% gmSEQ.CHN(4).PBN=PBDictionary('AOM');
% gmSEQ.CHN(4).NRise=2;
% gmSEQ.CHN(4).T=[0 gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi];
% gmSEQ.CHN(4).DT=[gmSEQ.readout gmSEQ.readout];
% ApplyDelays();

% %%%% Variable sequence length ctr 2%%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=1;
% gmSEQ.CHN(1).T=gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi;
% gmSEQ.CHN(1).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).PBN=PBDictionary('ctr1');
% gmSEQ.CHN(2).NRise=1;
% gmSEQ.CHN(2).T=gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi+gmSEQ.readout+1000;
% gmSEQ.CHN(2).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(3).PBN=PBDictionary('I');
% gmSEQ.CHN(3).NRise=3;
% gmSEQ.CHN(3).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+1000-50-gmSEQ.pi];
% gmSEQ.CHN(3).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 gmSEQ.pi];
% gmSEQ.CHN(4).PBN=PBDictionary('AOM');
% gmSEQ.CHN(4).NRise=4;
% gmSEQ.CHN(4).T=[0 gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi+gmSEQ.readout+1000 gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+1000];
% gmSEQ.CHN(4).DT=[gmSEQ.readout gmSEQ.readout 1000 1000];
% gmSEQ.CHN(5).PBN=PBDictionary('ctr2');
% gmSEQ.CHN(5).NRise=1;
% gmSEQ.CHN(5).T=gmSEQ.readout+gmSEQ.m+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+gmSEQ.pi+1000;
% gmSEQ.CHN(5).DT=gmSEQ.CtrGateDur;
% ApplyDelays();


%%%% Fixed sequence length ctr 2%%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=1;
% gmSEQ.CHN(1).T=gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi;
% gmSEQ.CHN(1).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).PBN=PBDictionary('ctr1');
% gmSEQ.CHN(2).NRise=1;
% gmSEQ.CHN(2).T=gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout+1000;
% gmSEQ.CHN(2).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(3).PBN=PBDictionary('I');
% gmSEQ.CHN(3).NRise=3;
% gmSEQ.CHN(3).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+1000-50-gmSEQ.pi];
% gmSEQ.CHN(3).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 gmSEQ.pi];
% gmSEQ.CHN(4).PBN=PBDictionary('AOM');
% gmSEQ.CHN(4).NRise=4;
% gmSEQ.CHN(4).T=[0 gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout+1000 gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+1000];
% gmSEQ.CHN(4).DT=[gmSEQ.readout gmSEQ.readout 1000 1000];
% gmSEQ.CHN(5).PBN=PBDictionary('ctr2');
% gmSEQ.CHN(5).NRise=1;
% gmSEQ.CHN(5).T=gmSEQ.readout+gmSEQ.To+1100+gmSEQ.pi+gmSEQ.readout+1000+1000+1000;
% gmSEQ.CHN(5).DT=gmSEQ.CtrGateDur;
% ApplyDelays();

% %%%% Fixed sequence length, IDENTICAL DUTY CYCLE %%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=3;
% gmSEQ.CHN(1).T=[gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+gmSEQ.pi+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100];
% gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
% gmSEQ.CHN(2).PBN=PBDictionary('AOM');
% gmSEQ.CHN(2).NRise=6;
% gmSEQ.CHN(2).T=[0 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100];
% gmSEQ.CHN(2).DT=[gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000];
% gmSEQ.CHN(3).PBN=PBDictionary('I');
% gmSEQ.CHN(3).NRise=3;
% gmSEQ.CHN(3).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100-gmSEQ.pi-50];
% gmSEQ.CHN(3).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 gmSEQ.pi];
% ApplyDelays();

% %%%% Variable sequence length, IDENTICAL DUTY CYCLE %%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=3;
% gmSEQ.CHN(1).T=[gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+gmSEQ.pi+1100 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100];
% gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
% gmSEQ.CHN(2).PBN=PBDictionary('AOM');
% gmSEQ.CHN(2).NRise=6;
% gmSEQ.CHN(2).T=[0 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100];
% gmSEQ.CHN(2).DT=[gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000];
% gmSEQ.CHN(3).PBN=PBDictionary('I');
% gmSEQ.CHN(3).NRise=3;
% gmSEQ.CHN(3).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.m+(gmSEQ.pi)+1100-gmSEQ.pi-50];
% gmSEQ.CHN(3).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 gmSEQ.pi];
% ApplyDelays();


% %%%% Fixed sequence length, IDENTICAL DUTY CYCLE WITH EXTRA PI %%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=3;
% gmSEQ.CHN(1).T=[gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+1100];
% gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
% gmSEQ.CHN(2).PBN=PBDictionary('I');
% gmSEQ.CHN(2).NRise=4;
% gmSEQ.CHN(2).T=[gmSEQ.readout+1050 gmSEQ.readout+1050+(gmSEQ.pi)/2+gmSEQ.m gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+1050 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+1050];
% gmSEQ.CHN(2).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 gmSEQ.pi gmSEQ.pi];
% gmSEQ.CHN(3).PBN=PBDictionary('AOM');
% gmSEQ.CHN(3).NRise=8;
% gmSEQ.CHN(3).T=[0 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000 gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+(gmSEQ.pi)+1100+1000+1000+gmSEQ.readout+gmSEQ.To+1100];
% gmSEQ.CHN(3).DT=[gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000];
% ApplyDelays();

function Echo4
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';

precess=gmSEQ.m+1100+2*(gmSEQ.pi);
precessTo=gmSEQ.To+1100+2*(gmSEQ.pi);
readout1=gmSEQ.readout-(gmSEQ.CtrGateDur+500);
readout2=gmSEQ.CtrGateDur+500;
readoutTot=readout1+readout2+200;
% 
% gmSEQ.CHN(1).PBN=PBDictionary('AOM');
% gmSEQ.CHN(1).NRise=3;
% gmSEQ.CHN(1).T=[precess+300 precess+readout1+200+300 precess+readoutTot+1100+gmSEQ.pi+300];
% gmSEQ.CHN(1).DT=[readout1 readout2 gmSEQ.readout];

gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=5;
gmSEQ.CHN(1).T(1)=0;
gmSEQ.CHN(1).T(2)=readout1+200;
gmSEQ.CHN(1).T(3)=readoutTot+precess;
gmSEQ.CHN(1).T(4)=gmSEQ.CHN(1).T(3)+gmSEQ.readout+precess+gmSEQ.pi;
gmSEQ.CHN(1).T(5)=gmSEQ.CHN(1).T(4)+gmSEQ.readout+1100+gmSEQ.pi;
gmSEQ.CHN(1).DT=[readout1 readout2 gmSEQ.readout gmSEQ.readout readout2];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=4;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=gmSEQ.CHN(1).T(2:end);
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=5;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)=gmSEQ.CHN(1).T(2)+readout2+1050;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(2)=gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)+gmSEQ.pi*3/2+gmSEQ.m;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(3)=gmSEQ.CHN(1).T(3)+gmSEQ.readout+1050;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(4)=gmSEQ.CHN(numel(gmSEQ.CHN)).T(3)+gmSEQ.pi*3/2+gmSEQ.m;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(5)=gmSEQ.CHN(1).T(5)-50-gmSEQ.pi;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[(gmSEQ.pi)/2 (gmSEQ.pi)/2 (gmSEQ.pi)/2 3/2*(gmSEQ.pi) gmSEQ.pi];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('Q');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)=gmSEQ.CHN(numel(gmSEQ.CHN)-1).T(1)+gmSEQ.pi/2+gmSEQ.m/2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(2)=gmSEQ.CHN(numel(gmSEQ.CHN)-1).T(3)+gmSEQ.pi/2+gmSEQ.m/2;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.pi gmSEQ.pi];
if gmSEQ.bFixDuty
    gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
    gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
    gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
    gmSEQ.CHN(numel(gmSEQ.CHN)).DT=readoutTot+precessTo+gmSEQ.readout+precessTo+gmSEQ.pi+gmSEQ.readout+1100+gmSEQ.pi+readout2;
end

ApplyDelays();

function ESR
global gSG gmSEQ
gSG.bfixedPow=1;
gSG.bfixedFreq=0;
gSG.bMod='Sweep';
gSG.bModSrc='External';
gSG.sweepRate=10;
gmSEQ.bLiO=1;
gmSEQ.ctrN=1;
% dummy sequence for DrawSequence
gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=1;
gmSEQ.CHN(1).T=0;
gmSEQ.CHN(1).DT=1;

function pbn = PBDictionary(type)
switch type
    case 'AOM'
        pbn=0;
    case 'MWswitch'
        pbn=1;
    case 'ctr0'
        pbn=2;
    case 'dummy' %this is used to maintain a given duty cycle
        pbn=3;
    case 'APDGate'
        pbn=4;
    case 'I'
        pbn=5;
    case 'Q'
        pbn=6;
    case 'ctr1'
        pbn=7;
end


function ApplyDelays
global gmSEQ;

    aom_delay=800; %%Rise time is approx 400 ns, start of rise at ~ -100ns 6/6/19 MChow. See AOM data
    
    if strcmp(gmSEQ.meas,'APD')
        detector_delay=-880;
    else
        detector_delay=0;
    end
    
for i=1:numel(gmSEQ.CHN)
    if gmSEQ.CHN(i).PBN==PBDictionary('AOM')
        gmSEQ.CHN(i).Delays=ones(1,2)*aom_delay;
    elseif gmSEQ.CHN(i).PBN==PBDictionary('ctr0')
        gmSEQ.CHN(i).Delays=ones(1,2)*detector_delay;
    else
        gmSEQ.CHN(i).Delays=zeros(1,2);
    end
end

function ApplyNoDelays
global gmSEQ
for i=1:numel(gmSEQ.CHN)
    gmSEQ.CHN(i).Delays=zeros(1,max(2,gmSEQ.CHN(i).NRise));
end

function T1
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';

%%%%% Fixed sequence length %%%%%%
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=1;
% gmSEQ.CHN(1).T=gmSEQ.readout+gmSEQ.m;
% gmSEQ.CHN(1).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).PBN=PBDictionary('ctr1');
% gmSEQ.CHN(2).NRise=1;
% gmSEQ.CHN(2).T=gmSEQ.readout+gmSEQ.m+gmSEQ.readout-gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(3).PBN=PBDictionary('AOM');
% gmSEQ.CHN(3).NRise=2;
% gmSEQ.CHN(3).T=[0 gmSEQ.readout+gmSEQ.m];
% gmSEQ.CHN(3).DT=[gmSEQ.readout gmSEQ.readout];
% gmSEQ.CHN(4).PBN=PBDictionary('dummy');
% gmSEQ.CHN(4).NRise=1;
% gmSEQ.CHN(4).T=(gmSEQ.readout)*2+gmSEQ.To-20;
% gmSEQ.CHN(4).DT=20;
% ApplyDelays();

%%%%% Variable sequence length%%%%%%
gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(1).NRise=3;
gmSEQ.CHN(1).T=[gmSEQ.readout+gmSEQ.m gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi+1000+1000+gmSEQ.readout+1000+gmSEQ.pi];
gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
gmSEQ.CHN(2).PBN=PBDictionary('AOM');
% gmSEQ.CHN(2).NRise=6;
% gmSEQ.CHN(2).T=[0 gmSEQ.readout+gmSEQ.m gmSEQ.readout+gmSEQ.m+1000+1000 gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+gmSEQ.m gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+gmSEQ.m+1000+1000 gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+gmSEQ.m];
% gmSEQ.CHN(2).DT=[gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000];
gmSEQ.CHN(2).NRise=6;
gmSEQ.CHN(2).T=[0 gmSEQ.readout+gmSEQ.m gmSEQ.readout+gmSEQ.m+1000+1000 gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi+1000+1000 gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi+1000+1000+gmSEQ.readout+1000+gmSEQ.pi];
gmSEQ.CHN(2).DT=[gmSEQ.readout 1000 gmSEQ.readout 1000 gmSEQ.readout 1000];
gmSEQ.CHN(3).PBN=PBDictionary('I');
gmSEQ.CHN(3).NRise=1;
gmSEQ.CHN(3).T=gmSEQ.readout+gmSEQ.m+1000+1000+gmSEQ.readout+1000+gmSEQ.pi+1000+1000+gmSEQ.readout+1000-50;
gmSEQ.CHN(3).DT=gmSEQ.pi;
ApplyDelays();

function ODMR
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=0;
gSG.bMod='IQ';
gSG.bModSrc='External';

delay=1050;%1050 or 100000

%%%%% Fixed sequence length %%%%%%

gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=2;
gmSEQ.CHN(1).T(1)=gmSEQ.pi+delay+50;
gmSEQ.CHN(1).T(2)=gmSEQ.pi+delay+50+gmSEQ.readout-(gmSEQ.CtrGateDur+1000)+gmSEQ.CtrGateDur+1000;
gmSEQ.CHN(1).DT=[gmSEQ.readout-(gmSEQ.CtrGateDur+1000) gmSEQ.CtrGateDur+1000];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=3;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)=gmSEQ.pi+delay+50;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(2)=gmSEQ.pi+delay+50+gmSEQ.readout-(gmSEQ.CtrGateDur+1000);
gmSEQ.CHN(numel(gmSEQ.CHN)).T(3)=gmSEQ.CHN(1).T(2);
if strcmp(gmSEQ.meas,'SPCM')
    gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
elseif strcmp(gmSEQ.meas,'APD')
    gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[1000 1000 1000];
end

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=delay;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=gmSEQ.pi;



gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=100;
ApplyDelays();

function ODMRDark
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=0;
gSG.bMod='IQ';
gSG.bModSrc='External';

delay=1050;
delay2=200;
%%%%% Fixed sequence length %%%%%%

gmSEQ.CHN(1).PBN=PBDictionary('AOM');
gmSEQ.CHN(1).NRise=2;
gmSEQ.CHN(1).T(1)=gmSEQ.pi+gmSEQ.pidark+gmSEQ.pi+delay+50+delay2+delay2;
gmSEQ.CHN(1).T(2)=gmSEQ.CHN(1).T(1)+gmSEQ.readout-(gmSEQ.CtrGateDur+1000)+gmSEQ.CtrGateDur+1000;
gmSEQ.CHN(1).DT=[gmSEQ.readout-(gmSEQ.CtrGateDur+1000) gmSEQ.CtrGateDur+1000];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=3;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)=gmSEQ.CHN(1).T(1);
gmSEQ.CHN(numel(gmSEQ.CHN)).T(2)=gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)+gmSEQ.readout-(gmSEQ.CtrGateDur+500);
gmSEQ.CHN(numel(gmSEQ.CHN)).T(3)=gmSEQ.CHN(1).T(2);
if strcmp(gmSEQ.meas,'SPCM')
    gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
elseif strcmp(gmSEQ.meas,'APD')
    gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[1000 1000 1000];
end

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=2;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)=delay;
gmSEQ.CHN(numel(gmSEQ.CHN)).T(2)=gmSEQ.CHN(numel(gmSEQ.CHN)).T(1)+gmSEQ.pi+gmSEQ.pidark+delay2+delay2;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=[gmSEQ.pi gmSEQ.pi];

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('BNC645');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=delay+gmSEQ.pi+delay2;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=gmSEQ.pidark;

gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=10000;
ApplyDelays();
function PESR
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=0;
gSG.bMod='IQ';
gSG.bModSrc='External';

%%%%% Fixed sequence length %%%%%%
gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(1).NRise=1;
gmSEQ.CHN(1).T=[0];
if strcmp(gmSEQ.meas,'SPCM')
    gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur];
elseif strcmp(gmSEQ.meas,'APD')
    gmSEQ.CHN(1).DT=[1000 1000];
end
gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('I');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=gmSEQ.CtrGateDur;
gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('AOM');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=[0];
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=gmSEQ.CtrGateDur;
gmSEQ.CHN(numel(gmSEQ.CHN)+1).PBN=PBDictionary('dummy');
gmSEQ.CHN(numel(gmSEQ.CHN)).NRise=1;
gmSEQ.CHN(numel(gmSEQ.CHN)).T=0;
gmSEQ.CHN(numel(gmSEQ.CHN)).DT=100000;
ApplyDelays();

function CPMGN %not really CPMGN, just testing
global gmSEQ gSG
gSG.bfixedPow=1;
gSG.bfixedFreq=1;
gSG.bMod='IQ';
gSG.bModSrc='External';
% gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
% gmSEQ.CHN(1).NRise=1;
% gmSEQ.CHN(1).T=10000;
% gmSEQ.CHN(1).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(2).PBN=PBDictionary('ctr1');
% gmSEQ.CHN(2).NRise=1;
% gmSEQ.CHN(2).T=10000;
% gmSEQ.CHN(2).DT=gmSEQ.CtrGateDur;
% gmSEQ.CHN(3).PBN=PBDictionary('AOM');
% gmSEQ.CHN(3).NRise=1;
% gmSEQ.CHN(3).T=[0];
% gmSEQ.CHN(3).DT=[(gmSEQ.readout)*2];
% gmSEQ.CHN(4).PBN=PBDictionary('ctr2');
% gmSEQ.CHN(4).NRise=1;
% gmSEQ.CHN(4).T=10000;
% gmSEQ.CHN(4).DT=gmSEQ.CtrGateDur;
% ApplyDelays();
gmSEQ.CHN(1).PBN=PBDictionary('ctr0');
gmSEQ.CHN(1).NRise=3;
gmSEQ.CHN(1).T=[10000 13000 16000];
gmSEQ.CHN(1).DT=[gmSEQ.CtrGateDur gmSEQ.CtrGateDur gmSEQ.CtrGateDur];
gmSEQ.CHN(2).PBN=PBDictionary('AOM');
gmSEQ.CHN(2).NRise=1;
gmSEQ.CHN(2).T=[0];
gmSEQ.CHN(2).DT=[(gmSEQ.readout)*2];
ApplyDelays();

function ApplyAPDGate()
global gmSEQ
if strcmp(gmSEQ.readout,'APD')
    NCHN=numel(gmSEQ.CHN)+1;
    gmSEQ.CHN(NCHN).PBN=PBDictionary('APDGate');
    CHNctr0=0;
    for i=1:numel(gmSEQ.CHN) %find which CHN is ctr0
        if gmSEQ.CHN(i).PBN==PBDictionary('ctr0')
            CHNctr0=i;
            break
        end
    end
    gmSEQ.CHN(NCHN).T=gmSEQ.CHN(CHNctr0).T;
    gmSEQ.CHN(NCHN).DT=gmSEQ.CHN(CHNctr0).DT;
end

