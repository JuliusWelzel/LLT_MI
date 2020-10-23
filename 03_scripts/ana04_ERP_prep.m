% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];
PATHIN_rts = [MAIN '02_data\03_RTs\'];

PATHOUT_ERP = [MAIN '02_data\04_ERPs\'];


%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_ERP)
    mkdir(PATHOUT_ERP)
end

%creating variables for ICA
cfg.ERP.LP = 20;
cfg.ERP.HP = 1;
cfg.ERP.PRUNE = 3;
cfg.ERP.resam = 100;
cfg.ERP.ep_time = [-1 3];
cfg.ERP.BL = [-500 -100];

load([PATHIN_rts 'RT_ALL.mat']);

%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*_finICA_clean.set']));
SUBJ = extractBefore({list.name},'_');

%% process all subjects for ERPs

for sub = 1:length(SUBJ)

        % load xdf files        
        EEG = pop_loadset('filename',[SUBJ{sub} '_finICA_clean.set'],'filepath',PATHIN_eeg);
        EEG = eeg_checkset(EEG, 'eventconsistency' );
        EEG.ID = SUBJ{sub};

        %low-pass filter
        EEG = pop_eegfiltnew(EEG, [],cfg.ERP.LP,[],0,[],0);

        %resample to 100 Hz
        EEG = pop_resample(EEG, cfg.ERP.resam);

        %high-pass filter
        EEG = pop_eegfiltnew(EEG, [],cfg.ERP.HP,[],1,[],0);

        %re-ref the data
        EEG = pop_reref(EEG,{'TP9','TP10'});
        EEG = pop_interp(EEG,EEG.chanlocs, 'spherical');


        e_type = {EEG.event.type};
        e_idx = contains({EEG.event.type},'E_');
        ep_e = e_type(e_idx);
        EEG = pop_epoch( EEG,ep_e,cfg.ERP.ep_time, 'newname', 'filt eps', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, cfg.ERP.BL);

        % Remove remaining non-stereotyped artifacts 
        EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg.ERP.PRUNE,cfg.ERP.PRUNE,0,0);
        EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg.ERP.PRUNE,cfg.ERP.PRUNE,0,0);
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);

        EEG.idx_ep_prune    = EEG.reject.rejglobal;
        EEG.idx_art         = squeeze(max(EEG.data,[],2))>100 | squeeze(min(EEG.data,[],2))<-100;

        % extract information from RT_ALL for SO onset and behaviour
        idx_SUB_inRTall = strcmp({RT_ALL.ID},SUBJ(sub));
        if sum(idx_SUB_inRTall) == 0;continue;end
        EEG.idx_cor     = RT_ALL(idx_SUB_inRTall).acc(11:end);
        EEG.SO_ms       = RT_ALL(idx_SUB_inRTall).SO_ms(11:end);
        
        % save set 
        EEG.setname = [SUBJ{sub} '_clean_ERP_SO'];
        EEG = pop_saveset( EEG, 'filename',[SUBJ{sub} '_ep_filt.set'],'filepath',PATHOUT_ERP);    
        
        % store part mean ERP
        ERP_all(sub).ID = EEG.ID;
        ERP_all(sub).pics = ep_e;
        ERP_all(sub).mERP = EEG.data;
        
        %store ep info
        ERP_all(sub).idx_prune  = EEG.idx_ep_prune;
        ERP_all(sub).idx_art    = EEG.idx_art;
        ERP_all(sub).idx_cor    = EEG.idx_cor;
        ERP_all(sub).RT_SO_s    = EEG.SO_ms;

end % SUB LOOP

cfg.ERP.chanlocs = EEG.chanlocs;
idx_badpart = cellfun(@isempty,{ERP_all.idx_art});
ERP_all(idx_badpart) = [];
save([PATHOUT_ERP 'ERPall.mat'],'ERP_all');
save ([PATHOUT_ERP 'cfg.mat'],'cfg');


