%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Prep of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich
% Supervisor: Cornelia Kranczioch (University of Oldenburg)

PATHIN_raw      = fullfile(MAIN,'02_data','00_raw');
PATHIN_RTs      = fullfile(MAIN,'02_data','03_RTs');
PATHOUT_data    = fullfile(MAIN,'02_data','03_RTs');

load(fullfile(PATHIN_RTs,'RT_ALL.mat'));

%% Gather all available datasets
% datasets are storde as LLT_SUBJXX.xdf
list = dir(fullfile(PATHIN_raw,'*LLT.xdf'));
SUBJ = extractBefore({list.name},'_');

%% loop over datasets

c = length(RT_ALL) + 1;

for sub = 1:length(SUBJ)
    
    % load xdf files
    EEG = pop_loadxdf(fullfile(PATHIN_raw,[SUBJ{sub} '_LLT.xdf']), 'streamtype', 'EEG');
    EEG = eeg_checkset(EEG, 'eventconsistency' );
      
    %% get all triggers
    all_trig = {EEG.event(:).type};
    an = {'1', '2', '3','4','empty'};
    not_pic = ismember(all_trig,an);
    pic_i = find(~ismember(all_trig,an));

    % get BP RT in ms
    for e = 1:size(pic_i,2);
        RT_ALL(c).BP_ms(e)  = (EEG.event(pic_i(e)+1).latency-EEG.event(pic_i(e)).latency)*(1000/EEG.srate);
        tmp  = strsplit(EEG.event(pic_i(e)).type, {'_','.'});
    end
    
    %% get all trigger information
    EEG     = LLT_stim( EEG );   
    
    % transfer to RT_ALL
    RT_ALL(c).ID        = SUBJ{sub};
    RT_ALL(c).con_ind   = EEG.con_ind;
    RT_ALL(c).stim_ind  = EEG.con_ind;
    RT_ALL(c).acc       = EEG.acc;
    RT_ALL(c).stim      = EEG.stim;
    c = c+1;
   
    
end

%% Transfer to table

%% Prep data for JASP long table format
SUBJ = str2double(extractAfter({RT_ALL.ID},'SUBJ'));

% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;

c = 1;
for s = 1:numel(RT_ALL) % Loop over all subjects
    
    for e = 10:size(RT_ALL(s).BP_ms,2) % Loop over all epochs pers subject
        
        % ident vars
        rts(c).ID       = RT_ALL(s).ID;
        rts(c).group    = LLTgroupname(RT_ALL(s).ID);
        
        rts(c).BP_ms    = RT_ALL(s).BP_ms(e);
        rts(c).con_idx  = RT_ALL(s).con_ind(e);
        rts(c).stim_idx = RT_ALL(s).stim_ind(e);
        rts(c).acc      = RT_ALL(s).acc(e);
        c = c+1;
    end
    
end
        
rts = struct2table(rts);
writetable(rts,fullfile(PATHOUT_data,'compare_nic_juw_llt.csv'));