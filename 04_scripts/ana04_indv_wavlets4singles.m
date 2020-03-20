% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];
PATHIN_rts = [MAIN '02_data\03_RTs\'];

PATHOUT_WAVELETS = [MAIN '02_data\04_wavelets\'];
PATHOUT_plots = [PATHOUT_WAVELETS 'beta_plots\'];


% Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_plots)
end

% assign necessary parameters for following analysis in EEGLAB *yay*
cfg_el.ep_length   = [-1.2 8];
cfg_el.ep_prune    = 3; % pruning parameters for rejection
cfg_el.freqs       = 1:1:35; % freqs OI
cfg_el.srate       = 500; 
cfg_el.times       = cfg_el.ep_length(1)*1000:1000/cfg_el.srate:cfg_el.ep_length(2)*1000-1; % time vec for epoch
cfg_el.cycles      = [5 8]; %cycle range from 5 to 8, Debener ea., 2005
cfg_el.bl_wvlt     = [-500 -200];
cfg_el.bl_erp      = [-600 -400];
cfg_el.postBETAt   = [0 1000];

save([PATHOUT_WAVELETS 'cfg.mat'],'cfg_el');


% create idices for analysis
idx_bl      = dsearchn(cfg_el.times',cfg_el.bl_wvlt');
idx_alpha   = cfg_el.freqs <= 12 & cfg_el.freqs >= 8;
idx_beta    = cfg_el.freqs <= 25 & cfg_el.freqs >= 15;

chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);


%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*finICA_clean.set']));
SUBJ = extractBefore({list.name},'_');

load([PATHIN_rts 'RT_ALL.mat']);
SUBJ = SUBJ(contains(SUBJ,{RT_ALL.ID})); % only condiser subjects with valid RT data

%%
doc_error = {};

for sub = 1:length(SUBJ)
    
   
    % check if dataset already exists
    if contains(SUBJ(sub),[dat_wvlt_all.id])% update for current evaluation
        disp(['WAVELETS for ' SUBJ{sub} ' have already been calculated. Continue with next dataset.']);
        continue;
    end
        
    %% load ICA cleaned sets
        
    EEG     = pop_loadset('filename',[SUBJ{sub} '_finICA_clean.set'],'filepath',PATHIN_eeg);
    EEG     = eeg_checkset(EEG, 'eventconsistency' );
    EEG.ID  = SUBJ{sub};
     
    EEG             = eeg_checkset( EEG );
    EEG             = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG             = pop_reref( EEG, []); % reref CAR
    nms_trig        = unique({EEG.event.type});
    idx_trig_OI     = contains(nms_trig,'E_');
    EEG             = pop_epoch( EEG, nms_trig(idx_trig_OI),cfg_el.ep_length, 'epochinfo', 'yes');
    EEG             = pop_rmbase( EEG, cfg_el.bl_erp ,[]);
    
    % Remove remaining non-stereotyped artifacts 
    EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg_el.ep_prune,cfg_el.ep_prune,0,0);
    EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg_el.ep_prune,cfg_el.ep_prune,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = eeg_checkset( EEG );
    
    EEG.idx_ep_prune    = EEG.reject.rejglobal;
    EEG.idx_art         = squeeze(max(EEG.data,[],2))>100 | squeeze(min(EEG.data,[],2))<-100;
    
    % epoching to correct response
    try
        idx_SUB_inRTall = strcmp({RT_ALL.ID},SUBJ(sub));
        EEG.idx_cor     = RT_ALL(idx_SUB_inRTall).acc(11:end);
        EEG.SO_ms       = RT_ALL(idx_SUB_inRTall).SO_ms(11:end);

    catch 
        doc_error{end+1} = [SUBJ{sub} ' // no sufficient RTs from audio'];
    end

    
    
    %% Convert to fieldtrip struct for TF anaylsis
    
    dat_raw = eeglab2fieldtrip(EEG, 'preprocessing'); 
    
    % MRLT-WVLT analysis according to Debener ea., 2005
    % create powspec size( 96/192 trials, 24 channel, 35 freqs, 4600 sample [-1.2-8s]))
    cfg             = [];
    cfg.channel     = 'EEG';
    cfg.method      = 'wavelet';
    cfg.output      = 'pow';
    cfg.keeptrials  = 'yes';
    cfg.foi         = 1:1:35; %frequency range
    cfg.toi         = dat_raw.time{1};
    cfg.width       = linspace(cfg_el.cycles(1),cfg_el.cycles(2),numel(cfg.foi));
    cfg.pad         = 'nextpow2';
    dat_wvlt        = ft_freqanalysis(cfg, dat_raw);
    
    %normalize trial data 
    tmp_dat     = dat_wvlt.powspctrm(:,:,:,:); %trials x 1channel x freqs x time
    tmp_bl      = squeeze(nanmean(tmp_dat(:,:,:,idx_bl(1):idx_bl(2)),4));
    tmp_dat_dB  = real(10*log10(tmp_dat./tmp_bl)); %transform to dB // dB = 10*log10 (signal/baseline)
    
    %update for new var
    dat_wvlt_bl = dat_wvlt;
    dat_wvlt_bl.powspctrm = tmp_dat_dB;
    
    cfg = [];
    cfg.showlabels   = 'yes';
%     cfg.baseline     = cfg_el.bl_wvlt/1000;
%     cfg.baselinetype = 'db';
    ft_multiplotTFR(cfg, dat_wvlt_bl) 
    save_fig(gcf,PATHOUT_plots,[SUBJ{sub} '_wvlts']);
    

    % add additional information
    % stimuli
    unique_ep = unique([EEG.event.epoch]);
    
    for st = 1:length(unique_ep)
        idx_ep_stim = contains({EEG.event.type},'E_') & [EEG.event.epoch] == unique_ep(st);
        stim_pic = EEG.event(find(idx_ep_stim,1)).type; % only use first important stimuli
        stim_split = strsplit(stim_pic,{'_','.'});
        stim(st,:) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=depth rotation, 4 = lateral rotation
    end
    
    
    dat_wvlt_bl.id = SUBJ{sub};
    dat_wvlt_bl.events = stim;
    clear stim
    
    dat_wvlt_bl.idx_ep_prune       = EEG.idx_ep_prune;
    dat_wvlt_bl.idx_art            = EEG.idx_art;
    dat_wvlt_bl.idx_cor            = EEG.idx_cor;
    dat_wvlt_bl.SO_ms              = EEG.SO_ms;

    save([PATHOUT_WAVELETS 'wvlts_' SUBJ{sub} '.mat'],'dat_wvlt', '-v7.3');
    
    %% extract post stim beta ERS & beta bursts
    
    idx_chan = find(strcmp(dat_wvlt_bl.label,'CPz')); 
    idx_med  = contains(dat_wvlt_bl.events(:,1),'lh') & contains(dat_wvlt_bl.events(:,4),{'60','120'}) |...
               contains(dat_wvlt_bl.events(:,1),'rh') & contains(dat_wvlt_bl.events(:,4),{'240','300'});

    idx_lat  = contains(dat_wvlt_bl.events(:,1),'lh') & contains(dat_wvlt_bl.events(:,4),{'240','300'}) |...
               contains(dat_wvlt_bl.events(:,1),'rh') & contains(dat_wvlt_bl.events(:,4),{'60','120'});
                     
    idx_cor  = dat_wvlt_bl.idx_cor;
    

    % extract beta rebound after SO + 1s if avaliable, otherwise NaN
    for t = 1:length(dat_wvlt_bl.events)
        [v idx_SO_ms] = min(abs(dat_wvlt_bl.time-dat_wvlt_bl.SO_ms(t)/1000));
        if idx_SO_ms+500 <= 4600 %maximum possible epoich length is 4600 samples
            dat_beta_st(t,:) = squeeze(mean(dat_wvlt_bl.powspctrm(t,idx_chan,idx_beta,idx_SO_ms:idx_SO_ms+500)))';   
        else
            dat_beta_st(t,:) = nan(1,501);
            dat_beta_st(t,1:length(squeeze(mean(dat_wvlt_bl.powspctrm(t,idx_chan,idx_beta,idx_SO_ms:end)))')) = ...
                squeeze(mean(dat_wvlt_bl.powspctrm(t,idx_chan,idx_beta,idx_SO_ms:end)))';  
        end
    end

    
    % absolute beta rebound from dat_beta_st [SO+1s] value per trial splited for correct and
    % incorrect trials
    dat_betaReb_cor = nanmean(dat_beta_st(~dat_wvlt_bl.idx_art(idx_chan,:) & dat_wvlt_bl.idx_cor,:),2);
    dat_betaReb_err = nanmean(dat_beta_st(~dat_wvlt_bl.idx_art(idx_chan,:) & ~dat_wvlt_bl.idx_cor,:),2);
    clear dat_beta_st
    
    
    % single trial beta timecourse
    [sort_SO idx_ep_length] = sort(dat_wvlt_bl.SO_ms(~dat_wvlt_bl.idx_art(idx_chan,:) & dat_wvlt_bl.idx_cor));
    dat_st_beta             = squeeze(mean(dat_wvlt_bl.powspctrm(~dat_wvlt_bl.idx_art(idx_chan,:) & dat_wvlt_bl.idx_cor,idx_chan,idx_beta,:),3));
    dat_st_beta             = dat_st_beta(idx_ep_length,:);
    dat_beta_topo           = squeeze(nanmean(nanmean(nanmean(dat_wvlt_bl.powspctrm(~dat_wvlt_bl.idx_art(idx_chan,:) & dat_wvlt_bl.idx_cor,:,idx_beta,idx_SO_ms:idx_SO_ms+500),1),3),4));
    

    %%
    close all
    figure
    subplot(1,3,1)
    singleBoxplot({dat_betaReb_cor',dat_betaReb_err'})
        xticks ([1,2]);
        xticklabels ({'Correct','Error'});
        tune_BP([c_lh;c_rh])
        [h p] = ttest2(dat_betaReb_cor,dat_betaReb_err);
        sigstar([1,2],p)
        ylabel '\beta-Power [dB]'

    subplot(1,3,2)
    topoplot(dat_beta_topo,chanlocs,...
    'maplimits',[-6 3],'emarker2',{idx_chan,'.','w',15})
        title '\beta-ERS for correct trials'
    
    subplot(1,3,3)
    imagesc(dat_wvlt_bl.time,1:size(dat_st_beta,1),dat_st_beta,[-20 20])
    hold on
    plot((sort_SO + abs(cfg_el.times(1)))/1000,1:size(dat_st_beta,1),'k','LineWidth',1)
        xlim ([-0.5 7])
        ylabel 'Trials sorted by RT'
        xlabel 'Time [s]'
        cb = colorbar;
        ylabel(cb,'\beta-Power [dB]')
        v = vline(0,{'k:', 'LineWidth', 1});

    save_fig(gcf,PATHOUT_plots,[SUBJ{sub} '_betaERS'],'FigSize',[0 0 30 15]);   
    
    %% store mean for all trials for alpha & beta for all participants
        
    dat_wvlt_all(sub).id                  = SUBJ(sub);
    dat_wvlt_all(sub).dat_alpha(:,:,:)     = squeeze(mean(dat_wvlt_bl.powspctrm(:,:,idx_alpha,:),3));
    dat_wvlt_all(sub).dat_beta(:,:,:)      = squeeze(mean(dat_wvlt_bl.powspctrm(:,:,idx_beta,:),3));
    
    dat_wvlt_all(sub).events        = dat_wvlt_bl.events;
    dat_wvlt_all(sub).idx_art       = dat_wvlt_bl.idx_ep_prune';
    dat_wvlt_all(sub).cor           = dat_wvlt_bl.idx_cor';

    dat_wvlt_all(sub).times             = dat_wvlt_bl.time;
    dat_wvlt_all(sub).freqs             = dat_wvlt_bl.freq;

    clear dat_wvlt_bl

    
end

save([PATHOUT_WAVELETS 'wvlts_all.mat'],'dat_wvlt_all', '-v7.3');
save([PATHOUT_WAVELETS 'error_all.mat'],'doc_error');












































