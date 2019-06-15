function InitializeData(handles)

global gmSEQ

gmSEQ.Repeat=str2double(get(handles.Repeat,'String'));
gmSEQ.Passes=str2double(get(handles.Passes,'String'));

if get(handles.bAverage,'Value')
    gmSEQ.Average=str2double(get(handles.Average,'String'));
else
    gmSEQ.Average=1;
end

gmSEQ.SweepParam=gmSEQ.From:(gmSEQ.To-gmSEQ.From)/(gmSEQ.N-1):gmSEQ.To;
gmSEQ.NSweepParam=length(gmSEQ.SweepParam);


gmSEQ.bGo=1;
gmSEQ.bGoAfterAvg=1;
for i=1:numel(gmSEQ.ch)
    if gmSEQ.ch(i).ID==SequencePool('PBDictionary','ctr0')
        gmSEQ.ctrN=int32(length(gmSEQ.ch(i).x)/4); %number of startstop counters per sequence
    end
end
gmSEQ.signal=NaN(gmSEQ.NSweepParam, gmSEQ.ctrN);