% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];

PATHOUT_WAVELETS = [MAIN '02_data\04_wavelets\'];
PATHOUT_plots = [PATHOUT_WAVELETS 'TF_plots\'];


% Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_plots)
end

% assign necessary parameters for following analysis
cfg.ep_length   = [-1 2];
cfg.ep_prune    = 3; % pruning parameters for rejection
cfg.freqs       = [1:40]; % freqs OI
cfg.srate       = 500; 
cfg.times       = cfg.ep_length(1)*1000:1000/cfg.srate:cfg.ep_length(2)*1000-1; % time vec for epoch
cfg.cycles      = [5 8]; %cycle range from 5 to 8, Debener ea., 2005
cfg.bl          = [-500 -200];
cfg.bl_erp      = [-600 -400];



%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*finICA_clean.set']));
SUBJ = extractBefore({list.name},'_');


%%
for sub = 1:length(SUBJ)
    %% check if dataset already exists
    if exist([SUBJ{sub} '_ep_filt.set'],'file') == 2 % update for current evaluation
        disp(['WAVELETS for ' SUBJ{sub} ' have already been calculated. Continue with next dataset.']);
        continue;
    end
        
    %% load ICA cleaned sets
        
    EEG     = pop_loadset('filename',[SUBJ{sub} '_finICA_clean.set'],'filepath',PATHIN_eeg);
    EEG     = eeg_checkset(EEG, 'eventconsistency' );
    EEG.ID  = SUBJ{sub};
    

    % wavelets according to Debener ea., 2005
    EEG             = eeg_checkset( EEG );
    EEG             = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG             = pop_reref( EEG, []);
    nms_trig        = unique({EEG.event.type});
    idx_trig_OI     = contains(nms_trig,'E_');
    EEG             = pop_epoch( EEG, nms_trig(idx_trig_OI),cfg.ep_length, 'epochinfo', 'yes');
    EEG             = pop_rmbase( EEG, cfg.bl_erp ,[]);
    
    % Remove remaining non-stereotyped artifacts 
    EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg.ep_prune,cfg.ep_prune,0,0);
    EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg.ep_prune,cfg.ep_prune,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = pop_rejepoch(EEG,find(EEG.reject.rejglobal),0);
    EEG = eeg_checkset( EEG );


   for ch = 1:size(EEG.data,1)
        tmpWAVELET = nan(numel(cfg.freqs),length(cfg.times),size(EEG.data,3)); % freqs x samples x trials
        % and each trial
        for t=1:size(EEG.data,3)
            tmpCh = EEG.data(ch,:,t);
            [tmpWAVELET(:,:,t), fres, tres, foi ] = jt_st_wave(tmpCh,cfg.srate,cfg.freqs,cfg.cycles(1),cfg.cycles(2),cfg.times);
        end
        
        % convert to dB
        base = squeeze(mean(tmpWAVELET,3));     % across trials --> freqs x times
        base = mean(base(:,TFP.BL(1)+125:TFP.BL(2)),2);             % across BL samples --> freqs
        base = repmat(base,1,size(tmpWAVELET,2));       % freqs x times
        base = repmat(base,1,1,size(EEG.data,3)); % freqs x times x trials
    
        % baseline correct and transform into percentage
        wv_ab = 100*(((tmpWAVELET)-base)./base);
        WAVELET_ep(:,:,:,ch) = wv_ab;      %Frequencies x TP x trials x channels

   end

    end
    end
end

originalEEG = EEG;
    %identify bad channels
    meanstdAllPs  = mean(std(EEG.data(:,:),0,2));
    stdstdAllPs  = std(std(EEG.data(:,:),0,2));
    [EEG] = trimOutlier(EEG, (meanstdAllPs-(3*stdstdAllPs)), (meanstdAllPs+(3*stdstdAllPs)), Inf, 0);
    
    %re-ref the data
    EEG = pop_reref(EEG,[17 18]);
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

%     % get individual RTs for BP & SO & stimuli parameters
%     [EEG.RT_BP_ms, EEG.RT_SO_ms, EEG.audio_ep, EEG.error_ep] = RT_LLT(EEG,MAIN); % see RT_all for further detail
%     [EEG] = LLT_stim(EEG);
%     
%     
%     e_type = {EEG.event.type};
%     e_idx = contains({EEG.event.type},'E_');
%     ep_e = e_type(e_idx);
%     EEG = pop_epoch( EEG,ep_e, [-0.2  0.7], 'newname', 'filt eps', 'epochinfo', 'yes');
%     EEG = pop_rmbase( EEG, [-150  0]);
%     figure; pop_plottopo(EEG, [1:24] , ['P300 ALL STIMULI  // ' EEG.ID], 0, 'ydir',1);
%     save_fig(gcf,PATHOUT_ERPplot,[EEG.ID '_chanERP']);
    
    % save set 
    EEG.setname = [SUBJ{sub},'_clean_ERP_SO'];
    EEG = pop_saveset( EEG, 'filename',[SUBJ{sub},'_ep_filt.set'],'filepath',PATHOUT_WAVELETS);    

   
    %% Correct trials with big difference BP & SO
if correcttrial
    % check epoch with containing artifacts
    % find bad audio epochs
    figs = dir(fullfile([MAIN '02_data\03_RTs\find_RT_SO\' EEG.ID  '*epoch*.fig']));
    EEG.speech_flag = ones(1,length(EEG.audio_ep));
    EEG.error_ep = logical(zeros(1,length(EEG.audio_ep)));

    name_figs = {figs.name};

    for f = 1:length(figs)

        onset = 0;
        while onset < 0.1 && onset > -1
            under = strfind(name_figs{f}, '_');
            dot = strfind(name_figs{f}, '.');
            fig_num = str2num(figs(f).name(under(end)+1:dot-1));

            ep = fig_num;
            down_s = 0.7; %audio downsampling factor
            play_ep(EEG,ep,down_s);

            % enter onset per hand
            onset = input('Enter onset in seconds (0 if repeat, 1 if correct, -1 if no speech): ');

        end
        
        if onset == 1 % SO found correct onset
            EEG.RT_SO_ms(ep) = EEG.RT_SO_ms(ep);
            EEG.speech_flag(ep) = 1;
        elseif onset == 0 % repeat play of epoch
            onset = 0;
            EEG.speech_flag(ep) = 0;
        elseif onset == -1; % no valuable speech
            onset = -1;
            EEG.speech_flag(ep) = 0;
        else % store RT by hand
            EEG.RT_SO_ms(ep) = onset*(1000/down_s);
            EEG.speech_flag(ep) = 1;
        end

    end

    % declare epoch to be checked manually if:
    % diff BP & SO is outlier or bigger than 600 ms
    % SO is smaller than 300 ms or bigger than 6000 ms
    diff_RT = EEG.RT_BP_ms-EEG.RT_SO_ms;
    crit_ep = isoutlier(diff_RT,'median')| diff_RT>600 | EEG.RT_SO_ms<300 |EEG.RT_SO_ms>6000; 
    crit_ep(EEG.error_ep) = 1;

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
    xlim ([1 length(EEG.audio_ep)]);
    xlabel('trials');ylabel('ms');
    title(['mean: ' num2str(round(mean(diff_RT)))]);
    hline(mean(diff_RT));
    save_fig(gcf,[MAIN '02_data\03_RTs\find_RT_SO\'],[EEG.ID '_SO_overview_new']);
    close(gcf);

    
    % save object with updated times
    EEG.setname = [SUBJ{sub},'_clean_ERP_SO'];
    EEG = pop_saveset( EEG, 'filename',[SUBJ{sub},'_ep_filt.set'],'filepath',PATHOUT_WAVELETS);    

    %% save clean RTs
    RT_ALL(sub).ID = EEG.ID;
    RT_ALL(sub).BP_ms = EEG.RT_BP_ms;
    RT_ALL(sub).SO_ms = EEG.RT_SO_ms;
    RT_ALL(sub).con_ind = EEG.con_ind;
    RT_ALL(sub).stim_ind = EEG.stim_ind;
    RT_ALL(sub).acc = EEG.acc;
    RT_ALL(sub).speech_flag = EEG.speech_flag;
    RT_ALL(sub).stim = EEG.stim;

end


    catch
        disp(['Error with ' SUBJ{sub}])
        docError(end+1) = SUBJ(sub);
        close
        break;

    end
   
%         % continue?
%     weiter = input('Enter 1 if next participant, 0 if abort: ');
%     if weiter == 0
%         break
%     end  
    
    end
    

    
    save([PATHOUT_RTs 'RT_ALL.mat'],'RT_ALL');
end
