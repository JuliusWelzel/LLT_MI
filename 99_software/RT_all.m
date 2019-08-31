 function [RT_BP_ms, RT_SO_ms, audio_epoch, error_ep] = RT_all(OBJ, audio);
    
    MAIN = ['Z:\projects\juw_LLT\'];

    all_trig = {OBJ.event(:).type};
    an = {'1', '2', '3','empty'};
    not_pic = ismember(all_trig,an);
    pic_i = find(~ismember(all_trig,an));

    for e = 1:size(pic_i,2);
        RT_BP_ms(e) = (OBJ.event(pic_i(e)+1).latency-OBJ.event(pic_i(e)).latency)*(1000/OBJ.srate);
    end

    % cut audio data to hand to onset function

    %merge stereo to mono
    audio.time_series = mean(audio.time_series);
    audio_lat_vec = 1:size(audio.time_series,2);
    audio_lat_vec = audio_lat_vec*(1000/str2num(audio.info.nominal_srate));


    for e = 1:size(pic_i,2);
        try
        %find latency of marker in audio
            [~,pic_lat] = min(abs(audio_lat_vec-(OBJ.event(pic_i(e)).latency)*(1000/OBJ.srate)));
            [~,bp_lat] = min(abs(audio_lat_vec-(OBJ.event(pic_i(e)+1).latency)*(1000/OBJ.srate)));
            audio_epoch {e} = audio.time_series(pic_lat:bp_lat);
            [RT_SO] =  rms_cepst_onset_detection_beta_21_11_17( audio_epoch{e}, str2double(audio.info.nominal_srate) );
            RT_SO_ms(e) = RT_SO *(1000/str2num(audio.info.nominal_srate));

        catch
            RT_SO_ms(e) = length(audio_epoch{e})*(1000/str2num(audio.info.nominal_srate));

        end
    end % all pics



    diff_RT = RT_BP_ms-RT_SO_ms;
    crit_ep = isoutlier(diff_RT,'median');%diff_RT>1000 |diff_RT>mean(diff_RT);

    figure; % subject overview of all epochs
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
    xlim ([1 length(audio_epoch)]);
    xlabel('trials');ylabel('ms');
    title(['mean: ' num2str(round(mean(diff_RT)))]);
    hline(mean(diff_RT));

    save_fig(gcf,[MAIN '02_data\03_RTs\find_RT_SO\'],[OBJ.ID '_SO_overview'],'figtype','.fig')
    close(gcf);

    for c = 1:length(crit_ep) % save single epcoh where difference RT_BP and RT_SO is BIIIIG
        if crit_ep(c) ==1
            figure;
            plot(audio_epoch{c});
            title(['Epoch ' num2str(c)]);
            xlabel('trials');ylabel('ms');
            vline(RT_SO_ms(c)/(1000/48000));
            save_fig(gcf,[MAIN '02_data\03_RTs\find_RT_SO\'],[OBJ.ID '_epoch_' num2str(crit_ep(c))],'figtype','.fig');
            close(gcf);
        end       
    end
 end %end fun