function EEG = createDummyEvents(EEG,labels,duration,removeOld)

if removeOld
    EEG.event = {};
end


num = floor(EEG.xmax/duration)-1;      
pnts = 1;
ep = length(EEG.event);

% create label list
labelList = repmat(1:length(labels),1,floor(num/length(labels)));
tmp = randperm(length(labelList));
labelList = labelList(tmp);
%keyboard;
% create dummy epochs
for i = 1:length(labelList)
    EEG.event(ep+i).type = labels{labelList(i)};
    EEG.event(ep+i).latency = pnts;
    EEG.event(ep+i).urevent = ep;
    pnts = pnts + (EEG.srate * duration);
end



