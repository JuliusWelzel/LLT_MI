% ICAw to clean data
% Identify Bad Components on 1-40 Hz filtered data and
% backproject them on rawdata


PATHIN_eeg = [MAIN '02_data\00_raw\'];
PATHIN_ica = [MAIN '02_data\01_ICAw\'];
PATHOUT_eeg = [MAIN '02_data\02_cleanEEG\'];
PATHOUT_ica = [PATHIN_ica 'indv_topoplot_and_BC\'];
PATHOUT_topop = [PATHOUT_ica 'indv_topoplot_and_BC\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_eeg)
    mkdir(PATHOUT_eeg);
end

if ~isdir(PATHOUT_ica)
    mkdir(PATHOUT_ica);
end

if ~isdir(PATHOUT_topop)
    mkdir(PATHOUT_topop);
end

% document potential problems with a subject
docError = {};      


%% Gather all available datasets
list = dir(fullfile([PATHIN_ica '*Ica*.set']));
SUBJ = extractBefore({list.name},'_');

%% Automatically reject ICs (EegLAB, ICLabel)

for sub =1:3%length(SUBJ)
    if exist([PATHOUT_eeg,SUBJ{sub},'_clean.set'],'file') == 9
        disp(['ICA for ' SUBJ{sub} ' has already been cleaned. Continue with next dataset.']);
        pause(0.2);
    else
    try
    
    %% load dataset with ICA weights
    EEG = pop_loadset('filename',[SUBJ{sub} '_ICAw.set'],'filepath',PATHIN_ica);
    EEG = eeg_checkset( EEG );

    %% ICLabel automated IC detection
    EEG = iclabel(EEG);
    
    %% enter bad components
    
    [v IC_idx] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    badcomps = find(IC_idx ~=1)';
    
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
    
    %% plot ICs with ICLabel
    IC_nms = EEG.etc.ic_classification.ICLabel.classes;
    IC_num = 1:length(EEG.icaweights);
    IC_cls = IC_idx;
    
    tempt = [''];
    for i = 1:length(IC_nms)
        tempt = [tempt IC_nms{i} ': ' num2str(IC_num(IC_cls ==i)) ', '];
    end  
    
    pop_viewprops( EEG, 0, [1:24], {'freqrange', [1 40]}, {}, 1, 'ICLabel' )
    save_fig(gcf,PATHOUT_ica,['ICs_' SUBJ{sub}])
        
        
    %% delete bad components
    EEG = pop_subcomp(EEG, [EEG.badcomps], 0);
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset( EEG, 'filename',[SUBJ{sub},'_clean.set'],'filepath',PATHOUT_eeg);    
    
   catch % write in Error doc
        disp(['Error with ' SUBJ{sub}])
        docError(end+1) = SUBJ(sub);
        close
        continue;
    end
    end
end

save([PATHOUT_eeg 'docError'], 'docError');
