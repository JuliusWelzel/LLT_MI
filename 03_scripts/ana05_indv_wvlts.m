% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];
PATHIN_rts = [MAIN '02_data\03_RTs\'];

PATHOUT_WAVELETS = [MAIN '02_data\05_wavelets\'];
PATHOUT_plots = [MAIN '04_plots\05_tf_plots\'];


% Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_plots)
end

% assign necessary parameters for following analysis in EEGLAB *yay*
load([PATHIN_eeg 'cfg.mat']);
cfg.wvlt.ep_length   = [-1.5 7];
cfg.wvlt.ep_prune    = 3; % pruning parameters for rejection
cfg.wvlt.freqs       = 1:1:35; % freqs OI
cfg.wvlt.srate       = 250;
cfg.wvlt.times       = cfg.wvlt.ep_length(1)*1000:1000/cfg.wvlt.srate:cfg.wvlt.ep_length(2)*1000-1; % time vec for epoch
cfg.wvlt.cycles      = [5 8]; %cycle range from 5 to 8, Debener ea., 2005
cfg.wvlt.bl_wvlt     = [-500 -200];
cfg.wvlt.bl_erp      = [-600 -400];

% create idices for analysis
idx_bl      = dsearchn(cfg.wvlt.times',cfg.wvlt.bl_wvlt');
idx_alpha   = cfg.wvlt.freqs <= 12 & cfg.wvlt.freqs >= 8;
idx_beta    = cfg.wvlt.freqs <= 25 & cfg.wvlt.freqs >= 15;

chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);


%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*finICA_clean.set']));
SUBJ = extractBefore({list.name},'_');

load([PATHIN_rts 'RT_ALL.mat']);
SUBJ = SUBJ(contains(SUBJ,{RT_ALL.ID})); % only condiser subjects with valid RT data

%%
doc_error = {};

for sub = 1:length(SUBJ)
    %% load ICA cleaned sets

    EEG     = pop_loadset('filename',[SUBJ{sub} '_finICA_clean.set'],'filepath',PATHIN_eeg);
    EEG     = eeg_checkset(EEG, 'eventconsistency' );
    EEG.ID  = SUBJ{sub};

    EEG             = eeg_checkset( EEG );
    EEG             = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG             = pop_resample(EEG,cfg.wvlt.srate);
    EEG             = pop_reref( EEG, []); % reref CAR
    nms_trig        = unique({EEG.event.type});
    idx_trig_OI     = contains(nms_trig,'E_');
    EEG             = pop_epoch( EEG, nms_trig(idx_trig_OI),cfg.wvlt.ep_length, 'epochinfo', 'yes');

    % Remove remaining non-stereotyped artifacts
    EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg.wvlt.ep_prune,cfg.wvlt.ep_prune,0,0);
    EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg.wvlt.ep_prune,cfg.wvlt.ep_prune,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = eeg_checkset( EEG );

    EEG.idx_ep_prune    = EEG.reject.rejglobal;
    EEG.idx_art         = squeeze(max(EEG.data,[],2))>100 | squeeze(min(EEG.data,[],2))<-100;

    idx_SUB_inRTall = strcmp({RT_ALL.ID},SUBJ(sub));
    EEG.idx_cor     = RT_ALL(idx_SUB_inRTall).acc(11:end);
    EEG.SO_ms       = RT_ALL(idx_SUB_inRTall).SO_ms(11:end);

    %% Convert to fieldtrip struct for TF anaylsis

    dat_wvlt = eeglab2fieldtrip(EEG, 'preprocessing');


    % MRLT-WVLT analysis according to Debener ea., 2005
    % create powspec size( 96/192 trials, 24 channel, 35 freqs, 4600 sample [-1.2-8s]))
    cfg_ft              = [];
    cfg_ft.channel      = 'EEG';
    cfg_ft.method       = 'wavelet';
    cfg_ft.output       = 'pow';
    cfg_ft.keeptrials   = 'yes';
    cfg_ft.foi          = cfg.wvlt.freqs; %frequency range
    cfg_ft.toi          = dat_wvlt.time{1};
    cfg_ft.width        = linspace(cfg.wvlt.cycles(1),cfg.wvlt.cycles(2),numel(cfg_ft.foi));
    cfg_ft.pad          = 'nextpow2';
    dat_wvlt            = ft_freqanalysis(cfg_ft, dat_wvlt);

    %normalize trial data
    tmp_bl      = squeeze(nanmean(dat_wvlt.powspctrm(:,:,:,idx_bl(1):idx_bl(2)),4));
    dat_wvlt.powspctrm  = real(10*log10(dat_wvlt.powspctrm./tmp_bl)); %transform to dB // dB = 10*log10 (signal/baseline)

    %update for new var
    cfg_ft = [];
    cfg_ft.showlabels   = 'yes';
    ft_multiplotTFR(cfg_ft, dat_wvlt)
    axis xy
    save_fig(gcf,PATHOUT_plots,[SUBJ{sub} '_wvlts']);

    % add additional information
    % stimuli
    unique_ep = unique([EEG.event.epoch]);

    for st = 1:length(unique_ep)
        idx_ep_stim = contains({EEG.event.type},'E_') & [EEG.event.epoch] == unique_ep(st);
        stim_pic = EEG.event(find(idx_ep_stim,1)).type; % only use first important stimuli
        stim_split = strsplit(stim_pic,{'_','.'});
        ep_stim(:,st) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=depth rotation, 4 = lateral rotation
    end


    dat_wvlt.id = SUBJ{sub};
    dat_wvlt.events = ep_stim;
    clear stim

    dat_wvlt.idx_ep_prune       = EEG.idx_ep_prune;
    dat_wvlt.idx_art            = EEG.idx_art;
    dat_wvlt.idx_cor            = EEG.idx_cor;
    dat_wvlt.SO_ms              = EEG.SO_ms;

        % find additional indices
    % idx lat/med
    dat_wvlt.idx_lat      = contains(ep_stim(1,:),'l')' & contains(ep_stim(4,:),{'240','300'})'|...
                            contains(ep_stim(1,:),'r')' & contains(ep_stim(4,:),{'60','120'})';
    dat_wvlt.idx_med      = contains(ep_stim(1,:),'l')' & contains(ep_stim(4,:),{'60','120'})'|...
                            contains(ep_stim(1,:),'r')' & contains(ep_stim(4,:),{'240','300'})';

    % idx hand/foot
    dat_wvlt.idx_hand      = contains(ep_stim(1,:),'h')';
    dat_wvlt.idx_foot      = contains(ep_stim(1,:),'f')';

    %angle
    dat_wvlt.idx_s_60      = strcmp(ep_stim(4,:),'60')'| strcmp(ep_stim(4,:),'300')';
    dat_wvlt.idx_s_120     = strcmp(ep_stim(4,:),'120')'|  strcmp(ep_stim(4,:),'240')';

    % idx in depth rotation
    dat_wvlt.idx_id_0      = contains(ep_stim(3,:),'0')';
    dat_wvlt.idx_id_60     = contains(ep_stim(3,:),'60')';
    dat_wvlt.idx_id_300    = contains(ep_stim(3,:),'300')';
    save([PATHOUT_WAVELETS 'wvlts_' SUBJ{sub} '.mat'],'dat_wvlt', '-v7.3');

    clear EEG

end
