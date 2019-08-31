%% nic_LLT eeg preprocessing for ICA
% all avaliable LLT datasets
% author: mareike daeglau, 06.10.17
% edited by Julius Welzel, 23.10.17 for nic_LLT data, only LLT data 


% # of datapoints needed for ICA: 3*N(umber of channels)^2
% for 64 channel cap: 12.288

%Define MAINPATH, concatenate with In- and Outputfolder

PATHIN = [MAIN '02_data\00_raw\'];
PATHOUT = [MAIN '02_data\01_ICAw\'];
PATHOUTeye = [PATHOUT 'EyeCatchComponents\'];
PATHOUTrejB = [PATHOUT 'delEEGdata_forICA\'];


%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT)
    mkdir(PATHOUT);
end

if ~isdir(PATHOUTeye)
    mkdir(PATHOUTeye);
end

if ~isdir(PATHOUTrejB)
    mkdir(PATHOUTrejB);
end

%creating variables for ICA
cfg.ICA.LP = 40;
cfg.ICA.HP = 1;
cfg.ICA.PRUNE = 3;

% document potential problems with a subject
docError = {};      


%% Gather all available datasets
% datasets are storde as LLT_SUBJXX.xdf
list = dir(fullfile([PATHIN '*LLT*SUBJ*.xdf']));
SUBJ = extractBetween({list.name},'_','.');


%% loop over datasets


for sub = 1:length(SUBJ)
    
    try % write in doc_error when failing
    
    %% check if dataset already exists
    if exist([SUBJ{sub} '_ICAw.set'],'file') == 9
        disp(['ICA for ' SUBJ{sub} ' has already been cleaned. Continue with next dataset.']);
    else 

%% load xdf files
    EEG = pop_loadxdf([PATHIN 'LLT_' SUBJ{sub} '.xdf'], 'streamtype', 'EEG', 'exclude_markerstreams', {});
    EEG = eeg_checkset(EEG, 'eventconsistency' );

    %add chanlocs
    EEG = pop_chanedit(EEG, 'lookup',[MAIN '99_software\mobile24_proper_labels.elp']);

    %low-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ICA.LP,[],0,[],0);

    %resample to 250 Hz
    EEG = pop_resample(EEG, 250);
    
    %high-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ICA.HP,[],1,[],0);
    
    
%% Find the start & end marker of the to be rejected data --> just ICA
    
%  find all presentation markers which are NOT pictures and store them

    trig = {EEG.event.type};
    not_pic = ismember(trig,{'1','2','3','empty'}); %all triggers which are not a pic
    lat_pic = [EEG.event(~not_pic).latency]; % latency of pictures
    type_pic = {EEG.event(~not_pic).type}; %type of stimulus
    
    %find training & experimental pictures
    T_pic_all = strfind(trig,'T');
    E_pic_all = strfind(trig,'E');
    
    % check for number of trials trianing & experiment
    if sum([E_pic_all{:}]) == (96) % check if number of experimental triggers is okay
        num_blocks(1).type = 'T';
        num_blocks(1).n_trig = sum([T_pic_all{:}]);
        num_blocks(2).type = 'E';
        num_blocks(2).n_trig = sum([E_pic_all{:}]);
    else
        docError(end+1) = {'Numbers of triggers not correct'}
    end
    EEG.num_blocks = num_blocks;
    
    %number of blocks to be deleted
    num_blocks_del = floor((sum([E_pic_all{:}])/32)+2);
    if sub == 30| sub == 25
        num_blocks_del = num_blocks_del +1;
    elseif sub == 39
            num_blocks_del = num_blocks_del -1;
    end

    % define start & end position of blocks to be deleted in samples (from
    % EEG.event.latency)
    % define the first bloch to be deleted (before training block)
    EEG.EndRejTrig(1,:) = [1 lat_pic(1)];

    % define other blocks to be deleted (training block 1st to 10th Picture, after that a complete block consists of 32 pictures)
    indx_s = sum([T_pic_all{:}]);
    indx_end = sum([T_pic_all{:}])+1;
    if sub == 39
            indx_s = 32;
            indx_end = 32+1;
    end
        
    for nb = 2:max(num_blocks_del)-1
            EEG.EndRejTrig(nb,:) = [lat_pic(indx_s) lat_pic(indx_end)];
            indx_s = indx_s+32;
            indx_end = indx_end+32;
    end
        
    % define the block to be deleted after the last block of pictures
    EEG.EndRejTrig(max(num_blocks_del),:) = [EEG.event(end).latency size(EEG.data,2)];

    % Reject the Interblock Data
    % save original dataset in new Dataset nEEG
    oEEG = EEG;
    % reject before definded parts of the data
    EEG = eeg_eegrej(EEG,EEG.EndRejTrig);

    % plot the original mean of all EEG channels & the deleted windows
    % Content of plot 
    plot_con = mean(oEEG.data(:,:),1);
    % transform position of marker to sec for plotting
    EEG.EndRejTrig = EEG.EndRejTrig./EEG.srate ;


    figure 
    % create timevector
    timevec = (1/EEG.srate:(1/EEG.srate):size(oEEG.data,2)/EEG.srate);
    plot(timevec,plot_con)
    hold on
    % plot rectangles around the deleted blocks
    for q = 1:length(EEG.EndRejTrig);
    rectangle ('Position',[EEG.EndRejTrig(q,1) min(plot_con) [EEG.EndRejTrig(q,2)-EEG.EndRejTrig(q,1)] [max(plot_con)-min(plot_con)]],'FaceColor', [1 0 0 0.5]);
    hold on
    end
    xlabel ('[sec]')
    ylabel ('\muV')
    %Set titel which includes the time differences
    if sub == 39 
        num_blocks_del = num_blocks_del +1;
    end
    title({[num2str(round((length(oEEG.data)-length(EEG.data))/EEG.srate)) ' sec of data were deleted in // ' '( -' round(num2str((100-size(EEG.data)/size(oEEG.data)*100))) '%) '],'LLT included ' num2str(num_blocks_del-2) 'blocks'}) 

    %save sec of deleted data & number of LLT blocks for further analysis
    SUBJinfo(sub).numsec_del = num2str(round((length(oEEG.data)-length(EEG.data))/EEG.srate));
    SUBJinfo(sub).num_LLT_bl = num2str(num_blocks_del-2);
    SUBJinfo(sub).num_pics = size(type_pic,2);
    SUBJinfo(sub).num_trig = length(EEG.event);

    % Save figure
    print([PATHOUTrejB, SUBJ{sub} ,'_deletedDATA.png'],'-dpng');  % raster image
    close
    %     
    %% ICA
    run_ICA = 1;
    
    if run_ICA 
    
    
    % making dummy epochs for ICA
    EEG = eeg_regepochs(EEG,'recurrence',1);
    
    N_EPOCHS_BEFORE = size(EEG.data,3);
    
    % pruning
    EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg.ICA.PRUNE,cfg.ICA.PRUNE,0,0);
    EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg.ICA.PRUNE,cfg.ICA.PRUNE,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = pop_rejepoch(EEG, EEG.reject.rejglobal ,0);

    N_EPOCHS_AFTER = size(EEG.data,3);
    
   %% ICA    
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'pca',EEG.nbchan);

    % saving variables
    cfg.ICA.N_EPOCHS_BEFORE(sub) = N_EPOCHS_BEFORE;
    cfg.ICA.N_EPOCHS_AFTER(sub) = N_EPOCHS_AFTER;

    
    EEG.setname = [SUBJ{sub},'_ICAw'];
    EEG = pop_saveset(EEG, 'filename',[SUBJ{sub},'_ICAw','.set'],'filepath',PATHOUT);
    end %if run ICA
    end %if subj ICAw
    catch % write in Error doc
        disp(['Error with ' SUBJ{sub}])
        docError(end+1) = SUBJ(sub);
        close
        continue;

    end
    clear trig;
    clear not_pic;
    clear lat_pic;
    clear type_pic;
    clear T_pic;
    clear T_pic_all:
    clear E_pic;
    clear E_pic_all;
    clear EEG.EndRejTrig;
    
end
save([PATHOUT 'docError'], 'docError');

