
function [onset] =  rms_cepst_onset_detection( audio_vec, Fs )



%
% rms_cepst_onset_detection( audio_file, Fs ) computes the speech onset and offset
% of a spoken word in the audiovector 'audio_vec'
%
% Inputs:
%
%       audio_vec   :  vector containing audio samples
%       Fs          :  sampling rate
%
% Outputs:
%
%       onset       :  speech onset in samples
%       offset      :  speech offset in samples
%       to_check    :  flag is set to true, if there has been a problem in
%                      finding the offset
% Requirements:
%
%       voicebox speech processing toolbox (can be downloaded here: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%


downsample_factor = 30; % downsample audio by factor 30;
minthreshold = 20;  %threshold for the 'findchangepts' function
window_size_env = 90; %window size for the envelope function
window_size_ceps = 15; % window size for the cepstrum

audio = struct('samples', '', 'x_changes','', 'cepband_one', '', 'cepband_one_smooth', '', 'chosen_mark', '', 'umschlag', '', 'cepband_two_smooth', '', 'window_factor', '');

sound.samples = audio_vec;

%% highpass filter audio signal to rid signal of baseline fluctuations
d = designfilt('highpassiir','FilterOrder',8, ...
    'PassbandFrequency',35,'PassbandRipple',0.2, ...
    'SampleRate',44100);

sound.samples = filter(d,sound.samples);


%% run the algorithm twice; once for onset, once for offset
for onset_offset = 1%:2
    
    if onset_offset == 3
        
        sound.samples = flip(sound.samples); % invert the audio signal when looking for offset
        
        %% check if there is any audio at the end of the file and if so, remove it
        
        file_ok = false;
        start_sample = 1;
        window_len = 500;
        window_len_pause = 3000; % approx 70ms
        checked = 1;
        to_check = false;
        
        while file_ok == false
            
            if start_sample >length(sound.samples)*0.4
                
                
                start_sample = 1;
                to_check = true;
                file_ok = true;
                
            end
            
            if mean(abs(sound.samples(start_sample:(start_sample+window_len))))<0.0004
                
                
                for check_for_pause = 1:3
                    
                    if mean(abs(sound.samples(start_sample+window_len_pause*(check_for_pause-1):(start_sample+window_len_pause*check_for_pause))))<0.0005
                        
                        
                        checked = checked +1;
                        
                        if checked == 5  %% If there has been audio detected 
                            
                            if start_sample<length(sound.samples)*0.4
                                
                                file_ok = true;
                                
                                sound.samples = sound.samples(start_sample:end);
                            
                            else
                                
                                start_sample = 1;
                                file_ok = true;
                                to_check = true;
                                
                            end
                            
                        end
                        
                    else
                        
                       start_sample = start_sample+window_len; 
                        
                    end
                    
                end
                
                
            else
                
                start_sample = start_sample+window_len;
                
            end
            
            
            
            
        end
        
    end
    
    
    audio_vec_down = downsample(sound.samples,downsample_factor); %Downsample for speed improvement
    
    %% compute envelope and apply lowpass filter
    
    sound.umschlag = envelope(audio_vec_down,300);
    
    [a_a,b_a] = butter(10,0.999);
    
    sound.umschlag = filter(b_a,10,sound.umschlag);
    
    %% find abrupt changes in the rms of the downsampled audio vector, which will be rejected or confirmed in the loop below
    
    sound.x_changes = findchangepts(audio_vec_down,'Statistic','rms','MinThreshold',minthreshold, 'MinDistance',150);
    
    %% compute cepstrum, extract fundamental frequency an apply low pass filter to it
    
    [c, tc] = melcepst(sound.samples, Fs, 'p', 10);
    cepst_spacing = diff(tc);
    cepst_spacing = cepst_spacing(1);
    sound.cepband_one = c(:,1);
    
    [a_b,b_b] = butter(10,0.80);
    sound.cepband_one_smooth = filter(b_b,10,sound.cepband_one);
    max_unfilt = max(sound.cepband_one);
    sound.cepband_one_smooth = (sound.cepband_one_smooth/max(sound.cepband_one_smooth))*max_unfilt;
    
    %% extract second frequency band and apply low pass filter
    
    [a_c,b_c] = butter(10,0.80);
    sound.cepband_two_smooth = filter(b_c,10,c(:,2));
    max_unfilt = max(c(:,2));
    sound.cepband_two_smooth = (sound.cepband_two_smooth/max(sound.cepband_two_smooth))*max_unfilt;
    
    
    %% define parameters for the following loop
    sound.chosen_mark = 1;
    check_marks = true;
    marker_i = 0;
    
    
    %% check the onset markers in  'sound.x_changes'
    
    
    while check_marks == true && marker_i < length(sound.x_changes)
        
        keep_searching = true;
        marker_i = marker_i +1;
        
        if round(sound.x_changes(marker_i)*downsample_factor/cepst_spacing)+window_size_ceps  < length(sound.cepband_one_smooth)
            
            %% define windows
            current_window_ceps_one = abs(sound.cepband_one_smooth(round(sound.x_changes(marker_i)*downsample_factor/cepst_spacing): round(sound.x_changes(marker_i)*downsample_factor/cepst_spacing)+window_size_ceps)-2);
            current_window_ceps_two = abs(sound.cepband_two_smooth(round(sound.x_changes(marker_i)*downsample_factor/cepst_spacing): round(sound.x_changes(marker_i)*downsample_factor/cepst_spacing)+window_size_ceps)-6);
            current_window_env = (sound.umschlag(sound.x_changes(marker_i): sound.x_changes(marker_i)+window_size_env));

        else 
            
            check_marks = false;
            keep_searching = false;
            to_check = true;
            
        end
            
            
        %% conditions in which the onset marker is not rejected
        
        if length(find(current_window_env >  0.4)) > 70 || length(find(abs(current_window_ceps_one)>10))>=4 || length(find(abs(current_window_ceps_two)>20))>=2
            
            if mean(current_window_env)>0.25
                
                keep_searching = false;
                
            end
            
        end
        
        
        %% if the marker has been rejected, check the next one
        if keep_searching == true
            
            if  sound.chosen_mark <= length(sound.x_changes)-1
                
                sound.chosen_mark = sound.chosen_mark +1;
                
            else
                
                sound.chosen_mark = length(sound.x_changes);
                
            end
            
        elseif keep_searching == false
            
            check_marks = false;
            sound.chosen_mark = marker_i;
                     
        end
        
    end
    
    %% If we are looking for the offset, the audiovector is reversed
    if onset_offset == 2
        
       try 
           
           offset = ((length(audio_vec_down)+start_sample/downsample_factor) - sound.x_changes(sound.chosen_mark))*downsample_factor;
        
       catch 
           
           to_check = true;
           offset = 0;
           
       end
       
    else
        
        onset = sound.x_changes(sound.chosen_mark)*downsample_factor;
        
    end
    
end


end

% end of function