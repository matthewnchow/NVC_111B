function SaveIgorText(handles)
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
%%%%%%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,comment('SEQUENCE PARAMETERS'));

for i=1:length(fnSEQ)
    st=fnSEQ(i);
    if ~strcmp(st,'CHN')&& ~strcmp(st,'ch') &&~strcmp(st,'SweepParam')&&~strcmp(st,'reference')&&~strcmp(st,'signal')&&~strcmp(st,'reference2')&&~strcmp(st,'trackCoords')
        fprintf(fid,comment(string(st)));
        fprintf(fid,comment(string(gmSEQ.(char(st)))));
    end
end

fprintf(fid,comment('SIGNAL GENERATOR PARAMETERS'));

for i=1:length(fnSG)
    st=fnSG(i);
    if ~strcmp(st,'serial')&&~strcmp(st,'qErr')
        fprintf(fid,comment(string(st)));
        fprintf(fid,comment(string(gSG.(char(st)))));
    end
end
%%%%%%%%%%%%%%%%%%%Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove NaNs from signal so that they don't interfere with arithemetic.
    if(gmSEQ.iAverage==1)
        sig=NaN(gmSEQ.iSweep,gmSEQ.ctrN);
    else
        sig=NaN(gmSEQ.NSweepParam,gmSEQ.ctrN);
    end
    
    if gmSEQ.ctrN==1
        sig=gmSEQ.signal(~isnan(gmSEQ.signal));
        sweep=gmSEQ.SweepParam(1:length(sig));
        organized=[sweep; sig.'];
    else
        for i=1:gmSEQ.ctrN
            sig0=gmSEQ.signal(:, i);
            sig(:,i)=sig0(~isnan(sig0));
        end
        sweep=gmSEQ.SweepParam(1:length(sig(:,1)));
        organized=[sweep; sig.'];
    end
    

% if gmSEQ.ctrN==3
%     fprintf(fid,'IGOR\nWAVES/D/O sweep,sig,ref,ref2\nBEGIN\n');
% elseif gmSEQ.ctrN==2
%     fprintf(fid,'IGOR\nWAVES/D/O sweep,sig,ref\nBEGIN\n');
% elseif gmSEQ.ctrN==1
%     fprintf(fid,'IGOR\nWAVES/D/O sweep,sig\nBEGIN\n');
% end
temp='';
for i=1:gmSEQ.ctrN
    temp=[temp '%d '];
end
temp=[temp '%d\n'];
fprintf(fid,temp,organized);
fprintf(fid, '\n');
% if gmSEQ.ctrN==3
%     fprintf(fid, 'END\nX Duplicate/O sig, data; data=(sig-ref2)/(ref-ref2)\nX Display data vs sweep; ModifyGraph width=200, height=100\nX ModifyGraph rgb(data)=(0,65535,0)\nX Display sig vs sweep; AppendToGraph ref vs sweep; ModifyGraph width=200, height=100\nX ModifyGraph rgb(ref)=(0,0,65535);\n');
%     fprintf(fid, 'X AppendToGraph ref2 vs sweep;\n');
% elseif gmSEQ.ctrN==2
%     fprintf(fid, 'END\nX Duplicate/O sig, data; data=sig/ref\nX Display data vs sweep; ModifyGraph width=200, height=100\nX ModifyGraph rgb(data)=(0,65535,0)\nX Display sig vs sweep; AppendToGraph ref vs sweep; ModifyGraph width=200, height=100\nX ModifyGraph rgb(ref)=(0,0,65535);\n');
% else
%     
%     fprintf(fid, 'END\nX Display sig vs sweep; ModifyGraph width=200, height=100\nX ModifyGraph rgb(ref)=(0,0,65535);\n');
% end

fclose(fid);

handles.textFileName.String = file;

function outStr = comment(inStr)
outStr=strcat('X// ',inStr,'\n');