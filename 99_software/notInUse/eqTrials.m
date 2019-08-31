% Reduces the number of trials of all specified eeg data structures
% to the number of trials of the smallest dataset

function varargout  = eqTrials(varargin)

trials = [];

for i=1:length(varargin)
    varargin{i}.trials;
    trials = [trials, varargin{i}.trials];
end

trials = min(trials);


for i=1:length(varargin)
    varargin{i}.data = varargin{i}.data(:,:,1:trials);
    varargin{i}.trials = trials;
end

varargout = varargin;


%% Equalize number of trials for both groups (needed for some of the remaining steps)
%if EEGnm.trials < EEGm.trials
%    EEGm.data = EEGm.data(:,:,1:EEGnm.trials)
%    EEGm.trials = EEGnm.trials;
%elseif EEGm.trials< EEGnm.trials
%    EEGnm.data = EEGnm.data(:,:,1:EEGm.trials)
%    EEGnm.trials = EEGm.trials;
%end
