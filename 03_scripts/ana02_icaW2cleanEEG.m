% ICAw to clean data
% Identify Bad Components on 1-40 Hz filtered data and
% backproject them on rawdata


PATHIN_eeg = [MAIN '02_data\00_raw\'];
PATHIN_ica = [MAIN '02_data\01_ICAw\'];
PATHOUT_eeg = [MAIN '02_data\02_cleanEEG\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_eeg)
    mkdir(PATHOUT_eeg);
end

% set possible bad ICs
iclab_cfg = {[2:7],[2,3,4,5],[3,5]};
iclab_nms = {'Brain','mdrt','cnsvtv'};

% document potential problems with a subject
docError = {};      
if exist([PATHIN_ica,'cfg.mat'],'file') == 2;load([PATHIN_ica,'cfg.mat']),end;

%% Gather all available datasets
list = dir(fullfile([PATHIN_ica '*Ica*.set']));
SUBJ = extractBefore({list.name},'_');

%% Automatically reject ICs (EegLAB, ICLabel)

for sub =1:length(SUBJ)
    if exist([PATHOUT_eeg,SUBJ{sub},'_',iclab_nms{1}, '_clean.set'],'file') == 2
        disp(['ICA for ' SUBJ{sub} ' has already been cleaned. Continue with next dataset.']);
        pause(0.2);
    else
    try
    for ii = 1:size(iclab_cfg,2)
        
    % load dataset with ICA weights
    EEG = pop_loadset('filename',[SUBJ{sub} '_ICAw.set'],'filepath',PATHIN_ica);
    EEG = eeg_checkset(EEG,'ica' );
    EEG.icaact = eeg_getica(EEG);
    
    %% ICLabel automated IC detection
    EEG = iclabel(EEG);
    
    %% enter bad components
    
    [v IC_idx] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);    
    badcomps = find(~ismember(IC_idx,iclab_cfg{ii}));
    
    %% store bad components & ICA relevant info in tmp
    TMP.icaweights = EEG.icaweights;
    TMP.icawinv = EEG.icawinv;
    TMP.icasphere = EEG.icasphere;
    TMP.badcomps = badcomps;
    TMP.icachansind = EEG.icachansind;
    TMP.etc = EEG.etc;
    
    %% Load "original" data set
    EEG = pop_loadxdf([PATHIN_eeg 'LLT_' SUBJ{sub} '.xdf'], 'streamtype', 'EEG', 'exclude_markerstreams', {});

    %add chanlocs
    EEG = pop_chanedit(EEG, 'lookup',[MAIN '99_software\mobile24_proper_labels.elp']);    
    
    %% overwrite tmp ICA info in original data
    EEG.icaweights = TMP.icaweights;
    EEG.icawinv = TMP.icawinv;
    EEG.icasphere = TMP.icasphere;
    EEG.badcomps = TMP.badcomps;
    EEG.icachansind = TMP.icachansind;
    EEG.etc = TMP.etc;
            
    %% delete bad components
    EEG = pop_subcomp(EEG, [EEG.badcomps], 0);
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset( EEG, 'filename',[SUBJ{sub} '_' iclab_nms{ii} '_clean.set'],'filepath',PATHOUT_eeg);    
   
    end
   catch % write in Error doc
        disp(['Error with ' SUBJ{sub}])
        docError(end+1) = SUBJ(sub);
        close
        continue;
    end
    end
end

save([PATHOUT_eeg 'docError'], 'docError');
