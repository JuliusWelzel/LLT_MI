% Calculates the power spectral density of an epoched EEG datasets.
%
% Input:
% data: 3D array (channels x samples x trials) of class1
% srate: integer specifying the sampling rate of the recorded signal
% channels (optional): integer array for specifying only a subset of channels.
% trials (optional): integer for only specifying a subset of trials.
%
%
% Output: 
% psd: 3D array containing the PSD values (channels, trials, frequency bin)

function [psd,f] = calcPSD(data,srate,channels,trials)

if isempty(channels)
    channels = 1:size(data,1);
end
if isempty(trials)
    trials = size(data,3);
end
    
psd = [];
%psd2 = [];

for c = 1:length(channels)
    for t = 1:trials
        [psd(c,t,:),f] = pwelch(data(channels(c),1:end,t),[],[],[],srate);
    end
end



%for c = 1:59
%    for t=1:EEGm.trials
%        [psd2(c, t,:),f] = pwelch(data2.data(c,1:end,t),[],[],[],srate);
%    end
%end

