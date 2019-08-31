function play_ep( OBJ, epoch, factor )
% INPUT:
%   OBJ: Object containing audio epoched data and Speech onset data
%   epoch: number of epoch to look at
%   factor of downsampling audio

% AUTHOR: Julius Welzel

audio_epoch = OBJ.audio_ep;%{epoch};
x = audio_epoch{epoch}; 
fs = 48000*(1/factor);
time=(0:length(x)-1)/fs;
figure(1)
set(gcf, 'Units','normalized','Position',[0 0 1 1]); % set figure size in cm
plot(time, x)
title(['Audio epoch ' num2str(epoch) OBJ.stim{epoch,1}]);
vline (OBJ.RT_SO_ms(epoch)/(1000*1/factor));
xlabel('Time (s)')
grid on
end_time = length(x)/fs;
h=line([0,0],[min(x) max(x)],'color','r','marker', 'o', 'linewidth', 2);
sound(x, fs) % starts playing the sound
tic % Starts Matlab timer
t=toc; % Gets the time since the timer started
while t<end_time
   set(h, 'xdata', t*[1,1]) % Moves the line to the time indicated by t
   drawnow % necessary to get figure updated
   t=toc; % get the current time for the next update
end

end

