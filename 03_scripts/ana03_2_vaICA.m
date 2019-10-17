% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

global MAIN
PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];

PATHOUT_ERP = [MAIN '02_data\04_ERPs\'];
PATHOUT_ERPplot = [PATHOUT_ERP 'ERP_plots\'];


%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_ERP)
    mkdir(PATHOUT_ERP)
    mkdir(PATHOUT_ERPplot)
end

%creating variables for ICA
load 'cfg.mat'
cfg.ERP.LP = 20;
cfg.ERP.HP = 0.1;
cfg.ERP.PRUNE = 3;
cfg.ERP.resam = 500;
cfg.ERP.ep_time = [-0.7 1.3];
cfg.ERP.BL = [-500 -100];

iclab_nms = {'finICA'};

% document potential problems with a subject
docError = {};   

%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*_brain_clean.set']));
SUBJ = extractBefore({list.name},'_');

%%
for ii = 1:size(iclab_nms,2)

    for sub =1:length(SUBJ)
    %% check if dataset already exists
        if exist([SUBJ{sub} '_ep_filt.set'],'file') == 2 % update for current evaluation
            disp(['RTs for ' SUBJ{sub} ' have already been calculated. Continue with next dataset.']);
        else 
        
    try
    % load xdf files        
    EEG = pop_loadset('filename',[SUBJ{sub} '_' iclab_nms{ii} '_clean.set'],'filepath',PATHIN_eeg);
    EEG = eeg_checkset(EEG, 'eventconsistency' );
    EEG.ID = SUBJ{sub};
    
    %low-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ERP.LP,[],0,[],0);

    %resample to 500 Hz
    EEG = pop_resample(EEG, cfg.ERP.resam);
    
    %high-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ERP.HP,[],1,[],0);
            
    %re-ref the data
    EEG = pop_reref(EEG,[]);
    EEG = pop_interp(EEG,EEG.chanlocs, 'spherical');
    
    
    e_type = {EEG.event.type};
    e_idx = contains({EEG.event.type},'E_');
    ep_e = e_type(e_idx);
    EEG = pop_epoch( EEG,ep_e,cfg.ERP.ep_time, 'newname', 'filt eps', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, cfg.ERP.BL);
    
    % save set 
    EEG.setname = [SUBJ{sub} '_clean_ERP_SO'];
    EEG = pop_saveset( EEG, 'filename',[SUBJ{sub} '_' iclab_nms{ii} '_ep_filt.set'],'filepath',PATHOUT_ERP);    
   
    % store part mean ERP
    ERP_all(sub).ID = EEG.ID;
    ERP_all(sub).pics = ep_e;
    ERP_all(sub).mERP = EEG.data;
    
    catch
        disp(['Error with ' SUBJ{sub}])
        docError(end+1) = SUBJ(sub);
        close
        break;

    end % error catch 
    end % if already processed
    end % SUB LOOP
    
    save([PATHOUT_ERP 'ERPall_' iclab_nms{ii} '.mat'],'ERP_all');
    save ([PATHOUT_ERP 'cfg.mat'],'cfg');
    
end % for different condition ICA

