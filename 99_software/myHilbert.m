% this auxilary function calculates the absolute values and hilbert
% transform of a signal. It can either take continuous data 
% (channels x samples) or epoched data (channels x trials x samples)
function dataOut = myHilbert(dataIn)

    % check whether the data is continuous or epoched
    epoched = 0; if length(size(dataIn))== 3; epoched = 1; end;
           
    %% process continuous data    
    if epoched == 0
        disp('Calculate abs and hilbert for continuous data'); dataOut =[];
        for c=1:size(dataIn,1)
            tmpChannel = dataIn(c,:)';
            tmpChannel = abs(hilbert(tmpChannel')');
            dataOut(c,:) = tmpChannel';
        end           
    end
    
    %% process epoched data 
    if epoched == 1
        disp('Calculate abs and hilbert for epoched data'); dataOut =[];
        for c=1:size(dataIn,1)
            tmpChannel = squeeze(dataIn(c,:,:))';
            tmpChannel = abs(hilbert(tmpChannel')');
            dataOut(c,:,:) = tmpChannel';
        end       
        
    end
    
   
    
    