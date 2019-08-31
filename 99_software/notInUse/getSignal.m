% Produces and outputs a specified signal as well as its corresponding
% trigger values
% Inputs:
%   fillLeft: double specifying the window width of the
%   resting state at the beginning of a new trial (in seconds)
%
%   opened: double specifying the window width of the
%   open state
%
%   flexion: double specifying the window width of the flexion (closing) state
%
%   fillBefClose: double specifying the window width of the resting state
%   before the closed imaging phase.
%
%   closed: double specifying the window width of the
%   closed state
%
%   extension: double specifying the window width of the extension (opening) state
%
%   fillRight: double specifying the window width of the
%   resting state at the end of a trial (might be slightly larger when fillLeft is jittered)
%
%   trials: integer specifying the number of trials to be created
%
%   jitter: double specifying the jitter before LED onset in seconds
%
%   jitterBeforeRest: double specifying the jitter before the resting
%   periods

%
% Outputs:
%   signal: 1 x n array constituting the whole signal
%
%   trigger: 1 x n array specifying the trigger positions of the different
%   mental state phases. Coding:
%       0 : nothing
%       1 : onset of Trial (fillLeft)
%       2 : onset of open state
%       3 : onset of flexion state
%       4 : onset of resting state before closed state
%       5 : onset of closed state
%       6 : onset of extension state
%       7 : onset of resting state at trial end

%
%   trialLength: integer indicating the number of samples of one trial

% examples:
% [signal, trigger] = getSignal(3.5,2,0.8,2,0.8,2,50);
% [signal,trigger] = getSignal([3 4.2],2,0.8,2,0.8,2,50);

function [signal, trigger, trialLength] = getSignal(fillLeft,opened,flexion,fillBefClosed,closed,extension,fillRight,trials,jitter,jitterBeforeRest)

Fs = 100;


% generate time-vecs
t.fillLeft = 0 : 1/Fs: fillLeft-1/Fs;
t.fillRight = 0 : 1/Fs: fillRight-1/Fs;
t.fillBefClosed = 0 : 1/Fs:fillBefClosed-1/Fs;
n =Fs*flexion; t.flexion = linspace(0,1,n);
%t.flexion	= 0 : 1/Fs: flexion-1/Fs;
t.opened= 0 : 1/Fs: opened-1/Fs;
n =Fs*extension; t.extension = linspace(0,1,n);
t.closed= 0 : 1/Fs: closed-1/Fs;

% generate signals
s.fillLeft = zeros(1,length(t.fillLeft)); %nan(1,length(t.fillLeft));
s.fillRight = zeros(1,length(t.fillRight)); %nan(1,length(t.fillRight));
s.fillBefClosed = ones(1,length(t.fillBefClosed))+1; %nan(1,length(t.fillBefClosed));
s.opened= zeros(1,length(t.opened));
s.flexion	= (cos(pi*t.flexion)*-1)+1;
s.closed = ones(1,length(t.closed))+1;
s.extension	= (cos(pi*t.extension))+1;



if isempty(jitter)
     
    % compute signal for one trial
    signal = [s.fillLeft s.opened s.flexion s.fillBefClosed s.closed s.extension s.fillRight];
    %figure; plot(signal);
    trialLength = length(signal);

    % compute trigger for one trial;
    trigger = zeros(1,length(signal));
    trigger(1,1) = 1;
    trigger(1,length(s.fillLeft)) = 2;
    trigger(1,length(s.fillLeft)+length(s.opened))=3;
    trigger(1,length(s.fillLeft)+length(s.opened)+length(s.flexion))=4;
    trigger(1,length(s.fillLeft)+length(s.opened)+length(s.flexion)+length(s.fillBefClosed))=5;
    trigger(1,length(s.fillLeft)+length(s.opened)+length(s.flexion)+length(s.fillBefClosed)+length(s.closed))=6;
    trigger(1,length(s.fillLeft)+length(s.opened)+length(s.flexion)+length(s.fillBefClosed)+length(s.closed)+length(s.extension))=7;

    % computer overall signal and overall trigger
    signal = repmat(signal,1,trials);
    trigger = repmat(trigger,1,trials);

else
    signal = [];
    trigger = [];
    jitter = jitter/(1/Fs);
    jitterBeforeRest = jitterBeforeRest/(1/Fs);
    
    
    for i = 1:trials
       tmp.jitter = randi(jitter);
       tmp.opened = [s.opened zeros(1,tmp.jitter)];
       tmp.closed = [s.closed ones(1,jitter-tmp.jitter)+1];
       
       tmp.jitterBeforeRest = randi(jitterBeforeRest);
       tmp.fillLeft = [s.fillLeft zeros(1,tmp.jitterBeforeRest)];
       tmp.fillBefClosed = [s.fillBefClosed ones(1,jitter-tmp.jitterBeforeRest)+1];
       
       % compute signal for current trial
       tmp.signal = [];
       tmp.signal = [tmp.fillLeft tmp.opened s.flexion tmp.fillBefClosed tmp.closed s.extension s.fillRight];
       trialLength = length(tmp.signal);
       
       
       % compute trigger for one trial;
       tmp.trigger = [];
       tmp.trigger = zeros(1,length(tmp.signal));
       tmp.trigger(1,1) = 1;
       tmp.trigger(1,length(tmp.fillLeft)) = 2;
       tmp.trigger(1,length(tmp.fillLeft)+length(tmp.opened))=3;
       tmp.trigger(1,length(tmp.fillLeft)+length(tmp.opened)+length(s.flexion))=4;
       tmp.trigger(1,length(tmp.fillLeft)+length(tmp.opened)+length(s.flexion)+length(tmp.fillBefClosed))=5;
       tmp.trigger(1,length(tmp.fillLeft)+length(tmp.opened)+length(s.flexion)+length(tmp.fillBefClosed)+length(tmp.closed))=6;
       tmp.trigger(1,length(tmp.fillLeft)+length(tmp.opened)+length(s.flexion)+length(tmp.fillBefClosed)+length(tmp.closed)+length(s.extension))=7;

       % concatente new trial to overall signal and trigger array
       signal = [signal,tmp.signal];
       trigger = [trigger,tmp.trigger];
       
       
    end    
        
    
end