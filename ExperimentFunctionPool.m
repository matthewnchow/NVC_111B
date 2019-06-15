function varargout=ExperimentFunctionPool(what,hObject, eventdata, handles, misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Written by Jeronimo Maze, July 2007 %%%%%%%%%%%%%%%%%%
%%%%%%%%%% Harvard University, Cambridge, USA  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch what
    case 'PBON'
        PBON(hObject, eventdata, handles);
    case 'PBOFF'
        Cyc10FunctionPool('CWChannels', 0);
    case 'Initialize'
        Initialize(hObject, eventdata, handles);
    case 'LoadSEQ'
        LoadSEQ(hObject, eventdata, handles,misc);
    case 'Run'
        LoadUserInputs(hObject,eventdata,handles);
        RunFirstPoint(hObject, eventdata, handles);
    case 'SaveData'
        SaveIgorText(handles); %In separate file
    case 'Plot'
        PlotRaw(hObject,eventdata,handles,misc);
        PlotData(hObject,eventdata,handles,misc);
    case 'Normalize'
        varargout{1}=Normalize();
    case 'PlotRaw'
        PlotRaw(hObject,eventdata,handles,misc);
    case 'PlotData'
        PlotData(hObject,eventdata,handles,misc);
    case 'RunFirstPoint'
        seq = LoadSEQ(hObject, eventdata, handles,handles.axes1);
        status = PrepC10(seq);
        if (status ~= 0)
            disp('FPGA not loaded properly!');
        else 
            disp('FPGA done loading. Verify on scope')
        end
    case 'LoadUserInputs'
        LoadUserInputs(hObject,eventdata,handles);
    otherwise
        disp('No Matches found in Pool Function');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PBON(hObject, eventdata, handles)

LioPB.PBN(1) = 0;
LioPB.OnOff(1) = (get(handles.LioPB0,'Value'));

LioPB.PBN(2) = 1;
LioPB.OnOff(2) = (get(handles.LioPB1,'Value'));

LioPB.PBN(3) = 2;
LioPB.OnOff(3) = (get(handles.LioPB2,'Value'));

LioPB.PBN(4) = 3;
LioPB.OnOff(4) = (get(handles.LioPB3,'Value'));

LioPB.PBN(5) = 4;
LioPB.OnOff(5) = (get(handles.LioPB4,'Value'));

LioPB.PBN(6) = 5;
LioPB.OnOff(6) = (get(handles.LioPB5,'Value'));

LioPB.PBN(7) = 6;
LioPB.OnOff(7) = (get(handles.LioPB6,'Value'));

LioPB.PBN(8) = 7;
LioPB.OnOff(8) = (get(handles.LioPB7,'Value'));

%Binary number
OutPuts = 0;
for ipbn = 1:8
    OutPuts = OutPuts + LioPB.OnOff(ipbn)*2^(LioPB.PBN(ipbn));
end

Cyc10FunctionPool('CWChannels', uint8(OutPuts))
%Jero 2008-07-10, changed by M Chow 6/6/19

function Initialize(hObject, eventdata, handles)
a = instrfind;
delete(a);
global gmSEQ 
StrL = SequencePool('PopulateSeq');
set(handles.sequence,'String',StrL);
clear StrL;

set(handles.axes1,'FontSize',8);
set(handles.axes2,'FontSize',8);
set(handles.axes3,'FontSize',8);

gmSEQ.bGo=0;

% load NI-DAQ MX DLL
LoadNIDAQmx;
% load SRS SG384 DLL
SignalGeneratorFunctionPool('Init',PortMap('SG com'));
gmSEQ.meas=PortMap('meas');

function varargout = LoadSEQ(hObject, eventdata, handles,ax)
global gmSEQ
LoadUserInputs(hObject,eventdata,handles);
[varargout{1}] = SequencePool(string(gmSEQ.name));
DrawSequence(gmSEQ, hObject, eventdata, ax);

function LoadUserInputs(hObject,eventdata,handles)
global gmSEQ gSG
str=get(handles.sequence, 'String');
val=get(handles.sequence, 'Value');
gmSEQ.name= str(val);
gmSEQ.From= str2double(get(handles.FROM1, 'String'));
gmSEQ.To= str2double(get(handles.TO1, 'String'));
gmSEQ.N= str2double(get(handles.SweepNPoints, 'String'));
gmSEQ.pi= str2double(get(handles.pi, 'String'));
gmSEQ.pidark=str2double(get(handles.pidark, 'String'));
gmSEQ.readout= str2double(get(handles.readout, 'String'));
gmSEQ.integrate= str2double(get(handles.integrate, 'String')); %Changed this from misc... ~MChow
%%%%% MChow added these buttons
gmSEQ.Passes= str2double(get(handles.Passes, 'String'));
gmSEQ.CtrGateDur=str2double(get(handles.CtrGateDur,'String'));
gmSEQ.AOMCounterDelay=str2double(get(handles.AOMCounterDelay,'String'));
gmSEQ.PumpT=str2double(get(handles.PumpT,'String'));
gmSEQ.MWSwitchDelay=str2double(get(handles.MWSwitchDelay,'String'));
%%%%%
gmSEQ.bWarmUpAOM=get(handles.bWarmUpAOM,'Value');
gmSEQ.bTrack=get(handles.bTrack,'Value');
gmSEQ.bFixDuty=get(handles.bFixDuty,'Value');
gSG.Pow = str2double(get(handles.fixPow, 'String'));
gSG.VDark = str2double(get(handles.fixVDark, 'String'));
gSG.Freq = str2double(get(handles.fixFreq, 'String'))*1e9;
gSG.FreqDark = str2double(get(handles.fixFreqDark, 'String'))*1e6;
% 
% function PlotExtRaw(hObject,eventdata,handles)
% global gmSEQ ScaleT ScaleStr
% figure;
% 
% 
% % set(cla,'FontSize',8);        %todo: figure out how to make this work...
% plot(single(gmSEQ.SweepParam).*ScaleT,single(gmSEQ.reference),'-b','LineStyle','--')
% 
% hold('on')
% 
% plot(single(gmSEQ.SweepParam).*ScaleT,single(gmSEQ.signal),'-r') 
% if gmSEQ.ctrN==3
%     hold('on')
%     plot(single(gmSEQ.SweepParam)*ScaleT,single(gmSEQ.reference2),'-m','LineStyle','--')
% end
% ylabel('Fluorescence counts');
%  xlabel(ScaleStr);
% hold('off')

% 
% 
% function PlotExt(hObject,eventdata,handles)
% global gmSEQ ScaleT ScaleStr
% figure;
% if gmSEQ.ctrN==3
%     data=(gmSEQ.signal(~isnan(gmSEQ.signal))-gmSEQ.reference2(~isnan(gmSEQ.reference2)))./(gmSEQ.reference(~isnan(gmSEQ.reference))-gmSEQ.reference2(~isnan(gmSEQ.reference2)));
% else
%     data=gmSEQ.signal(~isnan(gmSEQ.signal))./gmSEQ.reference(~isnan(gmSEQ.reference));
% end
% 
% plot(gmSEQ.SweepParam(1:length(data)).*ScaleT,data,'-g')
% ylabel('Fluorescence contrast');
% xlabel(ScaleStr);
% xlim([gmSEQ.SweepParam(1)*ScaleT gmSEQ.SweepParam(gmSEQ.NSweepParam)*ScaleT]);
 
function SaveData(handles)
global gSaveData gmSEQ  gSG
%global gScan
now = clock;
date = [num2str(now(1)),'-',num2str(now(2)),'-',num2str(round(now(3)))];
fullPath=fullfile(PortMap('Data'),date,'\');
if ~exist(fullPath,'dir')
    mkdir(fullPath);
end
gSaveData.path = fullPath;


gSaveData.file = ['_' date '.txt'];
name=regexprep(gmSEQ.name,'\W',''); % rewrite the sequence name without spaces/weird characters
%File name and prompt
B=fullfile(gSaveData.path, strcat(name, gSaveData.file));
file = strcat(name, gSaveData.file);

%Prevent overwriting
mfile = strrep(B,'.txt','*');
mfilename = strrep(gSaveData.file,'.txt','');

A = ls(char(mfile));
ImgN = 0;
for f = 1:size(A,1)
    sImgN = sscanf(A(f,:),strcat(name, string(mfilename), '_%d.txt'));
    if ~isempty(sImgN)
        if sImgN > ImgN
            ImgN = sImgN;
        end
    end
end
ImgN = ImgN + 1;
file = strrep(file,'.txt',sprintf('_%03d.txt',ImgN));
final= fullfile(gSaveData.path, file);
%Save File as Data
fnSEQ=fieldnames(gmSEQ);
fnSG=fieldnames(gSG);

fid = fopen(string(final),'wt');
fprintf(fid,'Sweep vector\n');
fprintf(fid,'%d\t',gmSEQ.SweepParam);
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid,'Signal vector\n');
fprintf(fid,'%d\t',gmSEQ.signal);
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid,'Reference vector\n');
fprintf(fid,'%d\t',gmSEQ.reference);
fprintf(fid, '\n');
fprintf(fid, '\n');

%fprintf(fid, 'Galvo Position: Vx = %.4f Vy = %.4f Vz = %.4f\n',[gScan.FixVx gScan.FixVy gScan.FixVz]);

fprintf(fid,'SEQUENCE PARAMETERS');
fprintf(fid, '\n'); 
fprintf(fid, '\n'); 

for i=1:length(fnSEQ)
    st=fnSEQ(i);
    if ~strcmp(st,'CHN') &&~strcmp(st,'SweepParam')&&~strcmp(st,'reference')&&~strcmp(st,'signal')
        fprintf(fid,string(st));
        fprintf(fid, '\n');      
        fprintf(fid,string(gmSEQ.(char(st))));
        fprintf(fid, '\n');        
        fprintf(fid, '\n'); 
    end
end

fprintf(fid,'SIGNAL GENERATOR PARAMETERS');
fprintf(fid, '\n'); 
fprintf(fid, '\n'); 

for i=1:length(fnSG)
    st=fnSG(i);
    if ~strcmp(st,'serial')&&~strcmp(st,'qErr')
        fprintf(fid,string(st));
        fprintf(fid, '\n');      
        fprintf(fid,string(gSG.(char(st))));
        fprintf(fid, '\n');        
        fprintf(fid, '\n'); 
    end
end

fclose(fid);

handles.textFileName.String = file;

function status = PrepC10(seq)
    status = Cyc10FunctionPool('LoadSequence', seq);
    attemptsleft = 2;
    while (status < 0 && attemptsleft > 0)
        disp(['Load FPGA failed! Attempts left:', num2str(attemptsleft)]);
        attemptsleft = attemptsleft - 1;
        status = Cyc10FunctionPool('LoadSequence', seq);
    end

function RunFirstPoint(hObject, eventdata, handles)
global gSG gmSEQ
gmSEQ.bGo=1;
% SignalGeneratorFunctionPool('SetIQ'); %We don't have IQ ports at the
% moment
if isfield(gmSEQ,'bLiO')   
%     PBFunctionPool('PBON',2^SequencePool('PBDictionary','AOM'));
    Cyc10FunctionPool('CWChannels',2^SequencePool('PBDictionary','AOM'));
    SignalGeneratorFunctionPool('WritePow');
    gSG.bOn=1; SignalGeneratorFunctionPool('RFOnOff');
    gSG.Freq=gmSEQ.m;
    SignalGeneratorFunctionPool('WriteFreq');
    while gmSEQ.bGo
        pause(.1);
        drawnow;
    end
else
    %if gSG.bfixedPow && gSG.bfixedFreq
    SignalGeneratorFunctionPool('WritePow');
    SignalGeneratorFunctionPool('WriteFreq');
    gSG.bOn=1; SignalGeneratorFunctionPool('RFOnOff');
%     for k=1:numel(gmSEQ.CHN)
%         gmSEQ.CHN(k).T=gmSEQ.CHN(k).T/1e9;
%         gmSEQ.CHN(k).DT=gmSEQ.CHN(k).DT/1e9;
%         gmSEQ.CHN(k).Delays=gmSEQ.CHN(k).Delays/1e9;
%     end
%     gmSEQ.Repeat=1e4;
    gmSEQ.Repeat=str2double(get(handles.Repeat,'String'));
    gmSEQ.Passes=str2double(get(handles.Passes,'String'));
    % Draw out the sequence once first
%     for i=1:50
%         pause(.1); drawnow;
%         if ~gmSEQ.bGo
%             break
%         end
%     end
    while gmSEQ.bGo
        RunSequence(hObject, eventdata, handles);
    end
end
    gSG.bOn=0; SignalGeneratorFunctionPool('RFOnOff');
    ExperimentFunctionPool('PBOFF',hObject, eventdata, handles);
%SaveData;

function PlotRaw(hObject,eventdata,handles,ext)
global gmSEQ
if ext
    figure
    ax=gca;
else
    ax=handles.axes2;
end
[ScaleT, ScaleStr, Color] = ScalePlot;

% if isfield(gmSEQ,'bLiO')
%     plot(ax, gmSEQ.SweepParam*ScaleT,gmSEQ.signal,Color{1})
% else
    plot(ax, gmSEQ.SweepParam*ScaleT,gmSEQ.signal(:,1),Color{1},'LineStyle','--')
% end

hold(ax, 'on')
for i=2:gmSEQ.ctrN
    plot(ax, gmSEQ.SweepParam*ScaleT,gmSEQ.signal(:,i),Color{i},'LineStyle','--')
end

 set(ax,'FontSize',8);
 ylabel(ax, 'Fluorescence counts');
 xlabel(ax, ScaleStr);
 xlim(ax, [gmSEQ.SweepParam(1)*ScaleT gmSEQ.SweepParam(gmSEQ.NSweepParam)*ScaleT]);
hold(ax, 'off')

function PlotData(hObject,eventdata,handles,ext)
    global gmSEQ
    if ext
        figure
        ax=gca;
    else
        ax=handles.axes3;
    end
    [ScaleT, ScaleStr, ~] = ScalePlot;
if gmSEQ.ctrN~=1 %~isfield(gmSEQ,'bLiO')&& %Taken out by MChow
    data=Normalize;
    plot(ax, gmSEQ.SweepParam(1:length(data)).*ScaleT,data,'.k')
    set(ax,'FontSize',8);
    ylabel(ax, 'Fluorescence contrast');
    xlabel(ax, ScaleStr);
    xlim(ax, [gmSEQ.SweepParam(1)*ScaleT gmSEQ.SweepParam(gmSEQ.NSweepParam)*ScaleT]);
end

function data=Normalize
    global gmSEQ
    
    %Won't have NaN's since sweeping through outer period in first run
    %~MChow
%             %Remove NaNs from signal so that they don't interfere with arithemetic.
%     if(gmSEQ.iAverage==1)
%         sig=NaN(gmSEQ.iSweep,gmSEQ.ctrN);
%     else
%         sig=NaN(gmSEQ.NSweepParam,gmSEQ.ctrN);
%     end
    
    for i=1:gmSEQ.ctrN
        sig0=gmSEQ.signal(:, i);
        sig(:,i)=sig0(~isnan(sig0));
    end
    
    if gmSEQ.ctrN==3
       if strcmp(gmSEQ.name,'T1JC')
            data=(sig(:,1)-sig(:,2))./sig(:,3);
            gmSEQ.normalize='(sig(:,1)-sig(:,2))./sig(:,3);';
       elseif or(strcmp(gmSEQ.name,'ODMR'),strcmp(gmSEQ.name,'ODMR Dark'))
           data=(sig(:,1)-sig(:,2))./(sig(:,3)-sig(:,2));
           gmSEQ.normalize='(sig(:,1)-sig(:,2))./(sig(:,3)-sig(:,2));';
       else
            data=(sig(:,1)-sig(:,3))./(sig(:,2)-sig(:,3));
            gmSEQ.normalize='(sig(:,1)-sig(:,3))./(sig(:,2)-sig(:,3));';
        end
    elseif gmSEQ.ctrN==4 && strcmp(gmSEQ.name,'Echo 4 counters')
        data=(sig(:,2)-sig(:,3))./(sig(:,1)-sig(:,4));
        gmSEQ.normalize='(sig(:,2)-sig(:,3))./(sig(:,1)-sig(:,4));';
    else
        data=sig(:,1)./sig(:,2);
        gmSEQ.normalize='sig(:,1)./sig(:,2);';
    end
    
function [ScaleT, ScaleStr, Color] = ScalePlot
global gmSEQ
if strcmp(gmSEQ.name,'ODMR') ||strcmp(gmSEQ.name,'ESR')||strcmp(gmSEQ.name,'CW ESR')||strcmp(gmSEQ.name,'Pulsed ESR')
    ScaleT=1e-9;
    ScaleStr='Frequency (GHz)';
elseif gmSEQ.To<1e3
    ScaleT=1;
    ScaleStr='ns';
elseif gmSEQ.To<1e6
    ScaleT=1e-3;
    ScaleStr='\mus';
else
    ScaleT=1e-6;
    ScaleStr='ms';
end
Color = {'r','b','k','m','c','y','g','r','b','k','m','c','y'};
