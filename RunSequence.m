function RunSequence(hObject, eventdata, handles)
global gmSEQ gSG tmax hCPS

audioread('WindowsNotify.wav');
BackupFile = PortMap('BackupFile');

InitializeData(handles);
% gmSEQ.integrate=1;
gmSEQ.integrate=str2double(get(handles.Passes,'String'));
try
disp(strcat('Commencing ',{' '},string(gmSEQ.name), ' sequence...'))
drawnow;

if gSG.bfixedPow && gSG.bfixedFreq
    gmSEQ.refCounts=Track('Init');
    SignalGeneratorFunctionPool('SetMod');
    SignalGeneratorFunctionPool('WritePow');
    SignalGeneratorFunctionPool('WriteFreq');
    gSG.bOn=1; SignalGeneratorFunctionPool('RFOnOff');
    SequencePool(string(gmSEQ.name));

    Nstartstop =  gmSEQ.ctrN * 2 * gmSEQ.N * gmSEQ.Passes;
    evens = (1:Nstartstop/2)* 2;
    odds = evens - 1;
    runtime = Nstartstop / 2 * gmSEQ.period / gmSEQ.ctrN;
    disp(['Runtime per sequence = ', num2str(runtime)]);
    [status, hCounter] = DAQmxFunctionPool('SetStartStopCounter', Nstartstop);
    DAQmxErr(status);
    for i=1:gmSEQ.Average
        gmSEQ.iAverage=i;
        handles.biAverage.String=num2str(gmSEQ.iAverage);
%         j=1; %Goes with AOM warmup
%         iwarmup=1;
        DAQmxErr(DAQmxWaitUntilTaskDone(hCounter, 10 + runtime * 1.2));
        vec = DAQmxFunctionPool('ReadCounter', hCounter, Nstartstop, 10 + runtime * 1.2); %Might not want to have infinite timeout
        vec = vec(evens) - vec(odds); %Get difference in counts between start and stop pulses
        sigDatum = transpose(sum(reshape(vec, gmSEQ.ctrN, gmSEQ.N, gmSEQ.Passes), 3)./gmSEQ.Passes);
%         [sigDatum] = ProcessData(vec);
        %% changed this part to do the whole array at once
        if i==1
            gmSEQ.signal=sigDatum;
        else
            gmSEQ.signal=((gmSEQ.signal*(i-1))+sigDatum)/i; %Simple averaging
        end
        % save a backup of the data here in case matlab crashes
        TemporarySave(BackupFile);
        ExperimentFunctionPool('Plot',hObject, eventdata, handles,0);
        drawnow;
        if ~gmSEQ.bGo
            break
        end
        if gmSEQ.bTrack
            Track('Run');
        end
        %Took this part out since AOM warm up is not too critical here
%         if (gmSEQ.bWarmUpAOM && iwarmup==2)||~gmSEQ.bWarmUpAOM||i~=1||j~=1
%             j=j+1;
%         elseif gmSEQ.bWarmUpAOM
%             iwarmup=iwarmup+1;
%         end
        if ~gmSEQ.bGo || ~gmSEQ.bGoAfterAvg
            break
        end
    end
    DAQmxClearTask(hCounter);

elseif gSG.bfixedPow && ~gSG.bfixedFreq %Chow: CWESR %%Yao: ODMR or ODMR Dark
    gmSEQ.refCounts=Track('Init');
    gSG.bOn=1; SignalGeneratorFunctionPool('RFOnOff');
    SignalGeneratorFunctionPool('WritePow');
    
    SignalGeneratorFunctionPool('SetMod');
    SequencePool(string(gmSEQ.name));
    gmSEQ.SweepParam=gmSEQ.SweepParam*1e9; %This is GHz
    
    Nstartstop = gmSEQ.ctrN * 2 * gmSEQ.Passes;
    evens = (1:Nstartstop/2)* 2;
    odds = evens - 1;
    runtime = gmSEQ.Passes * gmSEQ.period;
    disp(['Runtime per point = ', num2str(runtime)]);
    [status, hCounter] = DAQmxFunctionPool('SetStartStopCounter', Nstartstop);
    DAQmxErr(status);    

    gmSEQ.refCounts=Track('Run');
    
    for i=1:gmSEQ.Average
        gmSEQ.iAverage=i;
        handles.biAverage.String=num2str(gmSEQ.iAverage);
        j=1;
%         iwarmup=1;
        while j<=gmSEQ.NSweepParam
            gmSEQ.iSweep=j;
            gSG.Freq=gmSEQ.SweepParam(j);
            SignalGeneratorFunctionPool('WriteFreq');
            DAQmxErr(DAQmxWaitUntilTaskDone(hCounter, 10 + runtime * 1.2));
            vec = DAQmxFunctionPool('ReadCounter', hCounter, Nstartstop, 10 + runtime * 1.2);
            vec = vec(evens) - vec(odds);
            vec = reshape(vec, gmSEQ.ctrN, gmSEQ.Passes);
            [sigDatum] = mean(vec, 2);
            if i==1
                gmSEQ.signal(j, 1:gmSEQ.ctrN)=sigDatum(1:gmSEQ.ctrN);
            else
                gmSEQ.signal(j, 1:gmSEQ.ctrN)=(gmSEQ.signal(j, 1:gmSEQ.ctrN)*(i-1)+transpose(sigDatum(1:gmSEQ.ctrN)))/i;
            end
            % save a backup of the data here in case matlab crashes
			TemporarySave(BackupFile);
            drawnow;
            if ~gmSEQ.bGo
                break
            end
            j = j + 1;
%             if (gmSEQ.bWarmUpAOM && iwarmup==2)||~gmSEQ.bWarmUpAOM||i~=1||j~=1
%                 j=j+1;
%             elseif gmSEQ.bWarmUpAOM
%                 iwarmup=iwarmup+1;
%             end
            if gmSEQ.bTrack
                Track('Run');
            end
        end
        ExperimentFunctionPool('Plot',hObject, eventdata, handles,0);
        if ~gmSEQ.bGo || ~gmSEQ.bGoAfterAvg
            break
        end
    end
    DAQmxErr(DAQmxClearTask(hCounter));
end

gSG.bOn=0; SignalGeneratorFunctionPool('RFOnOff');
gmSEQ.bGo=0;
disp('Experiment completed!')
%sound(y,Fs);
catch ME
    KillAllTasks;
    gSG.bOn=0; SignalGeneratorFunctionPool('RFOnOff');
    rethrow(ME);
end
	
function TemporarySave(BackupFile)
global gSG gmSEQ

% convert relevant globals to a bigger structure
BackupExp.gmSEQ=gmSEQ;
BackupExp.gSG=gSG;
%save the data in matlab binary format
save(BackupFile,'BackupExp');
    
function [vec, NN, Freq] = MakeSweepVector
global gSG gmSEQ

if strcmp(gmSEQ.meas,'SPCM')
    vecA=-1:(2/(gmSEQ.NSweepParam-1)):1;
    vec=[vecA fliplr(vecA)];
    vec=repmat(vec,1,ceil(gSG.sweepRate*gmSEQ.misc/2));
    vec=[vec(1) vec];
    % Frequency of the timing is set such that analog output sweeps through
    % the range [-1 1] (one direction only) once every 1/sweepRate.
    Freq=gmSEQ.NSweepParam*gSG.sweepRate;
else % APD
    repeats=floor(500000/(gmSEQ.NSweepParam*gSG.sweepRate)); % 500 kHz the max analog input freq.
    vecA=-1:(2/(gmSEQ.NSweepParam-1)):1;
    vecB=repmat(vecA,repeats,1);
    vecB=reshape(vecB,[1,length(vecA)*repeats]);
    vec=[vecB fliplr(vecB)];
    vec=repmat(vec,1,ceil(gSG.sweepRate*gmSEQ.misc/2));
    Freq=repeats*gmSEQ.NSweepParam*gSG.sweepRate;
end
NN=length(vec);

function [sigDatum] = ProcessData(varargin)
global gmSEQ gSG
if strcmp(gmSEQ.meas,'SPCM')
    AA=diff(varargin{1});
else
    AA=varargin{1};
end

if isfield(gmSEQ,'bLiO') % ESR
    if strcmp(gmSEQ.meas,'APD')
        NN=length(AA);
        repeats=floor(500000/(gmSEQ.NSweepParam*gSG.sweepRate));
        AA=reshape(AA,[repeats NN/repeats]);
        AA=mean(AA);
    end
    sigDatum=zeros(1,gmSEQ.NSweepParam);
    samps=gmSEQ.misc*gSG.sweepRate;
    AA=reshape(AA,[gmSEQ.NSweepParam samps]);
    AA(:,2:2:(samps))=flipud(AA(:,2:2:(samps)));
    for i=1:gmSEQ.NSweepParam
        sigDatum(i)=sum(AA(i,:))/(samps);
    end
else
    sigDatum=NaN(1:gmSEQ.ctrN);
    for i=1:(gmSEQ.ctrN)
        sigDatum(i)=sum(AA(i:gmSEQ.ctrN:end));
    end
end

function refCounts = Track(what)
global gmSEQ gSG gScan
if gmSEQ.bTrack || strcmp(what,'CheckBadSignal')
     % get the handle of Gui1
     h = findobj('Tag','ImageNVCGUI');

     % if exists (not empty)
     if ~isempty(h)
        % get handles and other user-defined data associated to Gui1
        handles_ImageNVC = guidata(h);
     end

%     Leave laser at whatever duty cycle for tracking
%     PBFunctionPool('PBON',2^SequencePool('PBDictionary','AOM'));
    currentCounts = ImageFunctionPool('RunCPSOnce',0, 0, handles_ImageNVC);
%     if ~isfield(gmSEQ,'bLiO')
%         ExperimentFunctionPool('PBOFF',0, 0, handles_ImageNVC);
%     end
else
    refCounts=0;
    return
end
switch what
    case {'Init','CheckBadSignal'}
        refCounts=currentCounts;
        gmSEQ.trackCoords=[];
        return
    case 'Run'
        if currentCounts<.9*gmSEQ.refCounts
            for i=1:2
                currentCounts = ImageFunctionPool('NewTrackFast',0, 0, handles_ImageNVC);
                drawnow;
                if gmSEQ.bGo==0
                    refCounts=currentCounts;
                    return
                end
            end
            if currentCounts<.6*gmSEQ.refCounts
                KillAllTasks;
                gSG.bOn=0; SignalGeneratorFunctionPool('RFOnOff');
                gmSEQ.bGo=0;
                error('Tracking failed! Aborting...')
            end
            refCounts=currentCounts;
        else
            refCounts=gmSEQ.refCounts;
            gmSEQ.trackCoords=[gmSEQ.trackCoords; gScan.FixVx gScan.FixVy gScan.FixVz];
        end
end
