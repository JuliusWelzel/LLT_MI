% function that adds a further stream to an eeg dataset stored in
% xdf-format

% file: the filepath/filename of the .xdf-file
%
% eegStream: the name of eeg stream within LSL
%
% otherstream: the name of the other stream within LSL that shall be added
% to the eeg structure
%
% interpolate: 1= do interpolate; 0 = do not interpolate
%
% SR: number to which the dataset shall be up- or downsampled (optional)
%       --> caution: to some unknown reason this destroys the bcilab online
%        classifier. Do not use this funtionality within a bci paradigm.





function [EEG] = addStream(file,eegStream,otherStream,interpolate,SR)
 

    if strcmp(eegStream,otherStream)
       disp('Achtung: beide Streams können nicht gleich sein!'); 
    end
    
    tmp = load_xdf(file);
    % find the right streams in the tmp cell
    for i=1:length(tmp)
        if strcmp(tmp{i}.info.name,otherStream)
            idxOther = i;
        end
        if strcmp(tmp{i}.info.name,eegStream)
            idxEEG = i;
        end
    end
    
    
    EEG = eeg_load_xdf(file,'streamname',eegStream);
    
    % add lsl timeline to the eeg dataset
    EEG.lslTimes = tmp{idxEEG}.time_stamps;
    
    % if requested, fill in all datapoints of the other stream at the
    % correct positions
    if ~isempty(otherStream)
        EEG.myStream = nan(1,length(EEG.data));
        for i=1:length(tmp{idxOther}.time_stamps)
            val = tmp{idxOther}.time_stamps(i);
            tmp2 = abs(EEG.lslTimes-val);
            [idx idx] = min(tmp2); %index of closest value
            EEG.myStream(idx)=tmp{idxOther}.time_series(1,i);   
        end
    end
    
    % interpolate, if requested:
    if interpolate && ~isempty(otherStream)
        EEG.myStream = inpaint_nans(double(EEG.myStream));
    end
    clear tmp
        
    % change sampling rate, if requested
    if SR && ~isempty(otherStream)
       EEG.chanlocs(end+1).labels = 'tmp';
       EEG.data(end+1,:) = EEG.myStream;
       EEG.nbchan = size(EEG.data,1);
       EEG = pop_resample(EEG,SR);
       EEG.myStream = EEG.data(end,:);
       EEG.data(end,:)=[];
       EEG.nbchan = size(EEG.data,1);
       %keyboard;
       EEG.chanlocs = EEG.chanlocs(1:end-1);
       %keyboard;
    end
    

