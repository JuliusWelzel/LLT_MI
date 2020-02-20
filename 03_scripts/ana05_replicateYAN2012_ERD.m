%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_RTs  = [MAIN '02_data\03_RTs\'];
PATHIN_WVLT = [MAIN '02_data\04_wavelets\'];

PATHOUT_YAN = [MAIN '02_data\05_yan2012\'];
PATHOUT_plots = [MAIN '02_data\05_yan2012\plots\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_YAN)
    mkdir(PATHOUT_YAN);
    mkdir(PATHOUT_plots);
end

%% Gather all available datacfg_elts with clean EEG data
list = dir(fullfile([PATHIN_WVLT 'wvlts_*']));
list(strcmp({list.name},'wvlts_all.mat')) = []; 
nms_SUBJ = {list.name};
SUBJ = str2double(extractBetween({list.name},'wvlts_SUBJ','.mat'));

% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;

% groups vars for plotting
g_nms   = {'stroke','old','young'}; 

for i = 1:length(g_nms)
    groups(i).names = g_nms{i};
end

groups(1).idx = idx_stroke;
groups(2).idx = idx_old;
groups(3).idx = idx_young;


%% replicate behavioural analysis

% load RTs 
load([PATHIN_RTs 'RT_ALL.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                CONTINUE WITH EXTRACTING TIMED ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(SUBJ)
    
    load([PATHIN_WVLT nms_SUBJ{s}]);
    
    % handiness
    idx_lh      = contains(ERP_all(s).pics,'lh');
    idx_rh      = contains(ERP_all(s).pics,'rh');
    idx_sub_in_RTall = strcmp({RT_ALL.ID},{ERP_all(s).ID});
    idx_cor     = logical(RT_ALL(idx_sub_in_RTall).acc(11:end));
    
    % create 3D matrix [subs x channel x time_period] -> 56 x 24 x 4(BL,BEG,MDL,END)
    dat_erp_lh(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_lh & idx_cor),3),2)); 
    dat_erp_lh(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_lh & idx_cor),3),2)); 
    dat_erp_lh(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_lh & idx_cor),3),2)); 
    dat_erp_lh(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_lh & idx_cor),3),2)); 
    
    dat_erp_rh(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_rh & idx_cor),3),2)); 
    dat_erp_rh(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_rh & idx_cor),3),2)); 
    dat_erp_rh(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_rh & idx_cor),3),2)); 
    dat_erp_rh(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_rh & idx_cor),3),2)); 



    
end





