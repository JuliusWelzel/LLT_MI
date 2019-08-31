 function [RT_BP_ms, RT_SO_ms, audio_ep, error_ep] = RT_LLT(OBJ, MAIN);

%          RT_LLT - This function extract button press (BP) and
%          speech onset (SO) reaction times (RT) of a dataset
%          following a LLT paradigm

% 
%         [RT_BP_ms, RT_SO_ms, audio_ep, error_ep] = RT_all(OBJ, MAIN);
% 
%         INPUT
%             - OBJ: data structure containing latencies of trigger
%             - MAIN: path to store figures
% 
%         OUTPUT
%             - RT_BP_ms: Button pressed RT to a trigger
%             - RT_SO_ms: RT to a trigger found with speech onset 
%             - audio_ep: time series of audio data of single epoch
%             - error_ep: time series of audio data of single epoch where SO failed
% 

% 
%         Author: Julius Welzel
%         email: julius.welzel@uol.de
%         Date: 03.02.18

    PATHOUT = [MAIN '02_data\03_RTs\find_RT_SO\'];
    %Check if PATHOUT-folder is already there; if not: create
    if ~isdir(PATHOUT)
        mkdir(PATHOUT)
    end
    
    %% get all triggers
    all_trig = {OBJ.event(:).type};
    an = {'1', '2', '3','empty'};
    not_pic = ismember(all_trig,an);
    pic_i = find(~ismember(all_trig,an));

    % get BP RT in ms
    for e = 1:size(pic_i,2);
        RT_BP_ms(e) = (OBJ.event(pic_i(e)+1).latency-OBJ.event(pic_i(e)).latency)*(1000/OBJ.srate);
    end

    % load original audio data
    audio = load_xdf([MAIN '02_data\00_raw\LLT_' OBJ.ID '.xdf']);
    for s = 1:size(audio,2)
        if strcmp(audio{1,s}.info.name,'AudioCaptureWin') 
           ab = s;
        end
    end

    audio = audio{ab};

    %% epoch audio data to hand to onset function

    %merge stereo to mono
    audio.time_series = mean(audio.time_series);
    audio_lat_vec = 1:size(audio.time_series,2);
    audio_lat_vec = audio_lat_vec*(1000/str2num(audio.info.nominal_srate));

    error_ep = {};

    for e = 1:size(pic_i,2);
        try
        %find latency of marker in audio
        [~,pic_lat] = min(abs(audio_lat_vec-(OBJ.event(pic_i(e)).latency)*(1000/OBJ.srate)));
        [~,bp_lat] = min(abs(audio_lat_vec-(OBJ.event(pic_i(e)+1).latency)*(1000/OBJ.srate)));
        % find SO in single epochs
        audio_ep {e} = audio.time_series(pic_lat:bp_lat);
        [RT_SO] =  rms_cepst_onset_detection_beta_21_11_17( audio_ep{e}, str2num(audio.info.nominal_srate) );
        RT_SO_ms(e) = RT_SO *(1000/str2num(audio.info.nominal_srate));

        catch % store epoch when SO did not work
            error_ep{end+1} = e;
            RT_SO_ms(e) = length(audio_ep{e})*(1000/str2num(audio.info.nominal_srate));

        end

    end


    %% figure to plot time differences between BP & SO for all trials
    error_ep = [error_ep{:}];

    % declare epoch to be checked manually if:
    % diff BP & SO is outlier or bigger than 600 ms
    % SO is smaller than 300 ms or bigger than 6000 ms
    diff_RT = RT_BP_ms-RT_SO_ms;
    crit_ep = isoutlier(diff_RT,'median')| diff_RT>600 | RT_SO_ms<300 |RT_SO_ms>6000; 
    crit_ep(error_ep) = 1;

    figure;
    hold on
        for i = 1:length(diff_RT)
            h=bar(i,diff_RT(i));
            if crit_ep(i) ==0
                set(h,'FaceColor','k');
            else
                set(h,'FaceColor','r');
            end
        end
    hold off
    xlim ([1 length(audio_ep)]);
    xlabel('trials');ylabel('ms');
    title(['MEAN: ' num2str(round(mean(diff_RT))) ' ms']);
    hline(mean(diff_RT));

    save_fig(gcf,[MAIN '02_data\03_RTs\find_RT_SO\'],[OBJ.ID '_SO_overview'],'figtype','.fig')
    close(gcf);

    %% plot time outlier epochs time series
    ce_idx = find(crit_ep);
    for c = 1:length(ce_idx)
        figure
        plot(audio_ep{ce_idx(c)});
        title(['Epoch ' num2str(ce_idx(c))]);
        xlabel('Time [ms]');ylabel('Amplitude [muV]');
        vline(RT_SO_ms(ce_idx(c))/(1000/48000));
        save_fig(gcf,[MAIN '02_data\03_RTs\find_RT_SO\'],[OBJ.ID '_epoch_' num2str(ce_idx(c))],'figtype','.fig');
    end
end