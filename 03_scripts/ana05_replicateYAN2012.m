%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_ERP  = [MAIN '02_data\03_ERPs\'];
PATHIN_RTs  = [MAIN '02_data\03_RTs\'];
PATHIN_WVLT = [MAIN '02_data\04_wavelets'];

PATHOUT_YAN = [MAIN '02_data\05_yan2012\'];
PATHOUT_plots = [MAIN '02_data\05_yan2012\plots\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_YAN)
    mkdir(PATHOUT_YAN);
    mkdir(PATHOUT_plots);
end

%% Gather all available datacfg_elts with clean EEG data
list = dir(fullfile([PATHIN_ERP '*_ep_filt.set']));
SUBJ = str2double(extractAfter({RT_ALL.ID},'SUBJ'));

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

%% get data behavioural

for s = 1:numel(RT_ALL)
    
    % handiness
    idx_lh      = strcmp(RT_ALL(s).stim(11:end,1),'lh');
    idx_rh      = strcmp(RT_ALL(s).stim(11:end,1),'rh');
    idx_cor     = logical(RT_ALL(s).acc(11:end)');
   
    % create 3D matrix [subs x hand x mod] -> 56 x 2(LH,RH) x 2(RT,ACC)
    dat_b_hand(s,1,1) = nanmean(RT_ALL(s).SO_ms(idx_lh & idx_cor)); 
    dat_b_hand(s,2,1) = nanmean(RT_ALL(s).SO_ms(idx_rh & idx_cor));
    dat_b_hand(s,1,2) = nanmean(RT_ALL(s).acc(idx_lh));
    dat_b_hand(s,2,2) = nanmean(RT_ALL(s).acc(idx_rh));
    
    
    %angle
    idx_60      = strcmp(RT_ALL(s).stim(11:end,4),'60') | strcmp(RT_ALL(s).stim(11:end,4),'300');    
    idx_120     = strcmp(RT_ALL(s).stim(11:end,4),'120') | strcmp(RT_ALL(s).stim(11:end,4),'240');
    
    % create 3D matrix [subs x ang x mod] -> 56 x 2(60,120) x 2(RT,ACC)
    dat_b_ang(s,1,1) = nanmean(RT_ALL(s).SO_ms(idx_60 & idx_cor));
    dat_b_ang(s,2,1) = nanmean(RT_ALL(s).SO_ms(idx_120 & idx_cor));
    dat_b_ang(s,1,2) = nanmean(RT_ALL(s).acc(idx_60));
    dat_b_ang(s,2,2) = nanmean(RT_ALL(s).acc(idx_120));

end


%% %%%%%%%%%%%%%%%%%% plot behaviroual results handiness %%%%%%%%%%%%%%%%%%

%prepare data
RT_LH_s = dat_b_hand(idx_stroke,1,1);
RT_LH_o = dat_b_hand(idx_old,1,1);
RT_RH_s = dat_b_hand(idx_stroke,2,1);
RT_RH_o = dat_b_hand(idx_old,2,1);

ACC_LH_s = dat_b_hand(idx_stroke,1,2);
ACC_LH_o = dat_b_hand(idx_old,1,2);
ACC_RH_s = dat_b_hand(idx_stroke,2,2);
ACC_RH_o = dat_b_hand(idx_old,2,2);


% RT compare old stroke
[m_LH_o se_LH_o] = mean_SEM(RT_LH_o');
[m_LH_s se_LH_s] = mean_SEM(RT_LH_s');
[m_RH_o se_RH_o] = mean_SEM(RT_RH_o');
[m_RH_s se_RH_s] = mean_SEM(RT_RH_s');

close all
figure
subplot(2,2,1)
peb = errorbar([m_LH_o m_LH_s;m_RH_o m_RH_s],[se_LH_o se_LH_s;se_RH_o se_RH_s],'-s','MarkerSize',10);
peb(1).Color = c_lh;
peb(2).Color = c_rh;
peb(1).MarkerFaceColor = c_lh;
peb(2).MarkerFaceColor = c_rh;

xlim ([0.5 2.5])
ylim ([1000 2500])
ylabel 'RT [ms]'
xticks ([1 2])
xticklabels ({'LH','RH'})
legend({'Control','Stroke'})
box off

% ACC compare old stroke
[m_LH_o se_LH_o] = mean_SEM(ACC_LH_o');
[m_LH_s se_LH_s] = mean_SEM(ACC_LH_s');
[m_RH_o se_RH_o] = mean_SEM(ACC_RH_o');
[m_RH_s se_RH_s] = mean_SEM(ACC_RH_s');

subplot(2,2,2)
peb = errorbar([m_LH_o m_LH_s;m_RH_o m_RH_s],[se_LH_o se_LH_s;se_RH_o se_RH_s],'-s','MarkerSize',10);
peb(1).Color = c_lh;
peb(2).Color = c_rh;
peb(1).MarkerFaceColor = c_lh;
peb(2).MarkerFaceColor = c_rh;

xlim ([0.5 2.5])
ylim ([0.7 1])
ylabel 'ACC [%]'
xticks ([1 2])
xticklabels ({'LH','RH'})
legend({'Control','Stroke'})
box off


%%%%%%%%%%%%%%% plot behaviroual results angle %%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare data
RT_LH_s = dat_b_ang(idx_stroke,1,1);
RT_LH_o = dat_b_ang(idx_old,1,1);
RT_RH_s = dat_b_ang(idx_stroke,2,1);
RT_RH_o = dat_b_ang(idx_old,2,1);

ACC_LH_s = dat_b_ang(idx_stroke,1,2);
ACC_LH_o = dat_b_ang(idx_old,1,2);
ACC_RH_s = dat_b_ang(idx_stroke,2,2);
ACC_RH_o = dat_b_ang(idx_old,2,2);

% RT compare old stroke
[m_LH_o se_LH_o] = mean_SEM(RT_LH_o');
[m_LH_s se_LH_s] = mean_SEM(RT_LH_s');
[m_RH_o se_RH_o] = mean_SEM(RT_RH_o');
[m_RH_s se_RH_s] = mean_SEM(RT_RH_s');

subplot(2,2,3)
peb = errorbar([m_LH_o m_LH_s;m_RH_o m_RH_s],[se_LH_o se_LH_s;se_RH_o se_RH_s],'-s','MarkerSize',10);
peb(1).Color = c_lh;
peb(2).Color = c_rh;
peb(1).MarkerFaceColor = c_lh;
peb(2).MarkerFaceColor = c_rh;

xlim ([0.5 2.5])
ylim ([1000 2500])
ylabel 'RT [ms]'
xticks ([1 2])
xticklabels ({'+-60�','+-120�'})
legend({'Control','Stroke'})
box off

% ACC compare old stroke
[m_LH_o se_LH_o] = mean_SEM(ACC_LH_o');
[m_LH_s se_LH_s] = mean_SEM(ACC_LH_s');
[m_RH_o se_RH_o] = mean_SEM(ACC_RH_o');
[m_RH_s se_RH_s] = mean_SEM(ACC_RH_s');

subplot(2,2,4)
peb = errorbar([m_LH_o m_LH_s;m_RH_o m_RH_s],[se_LH_o se_LH_s;se_RH_o se_RH_s],'-s','MarkerSize',10);
peb(1).Color = c_lh;
peb(2).Color = c_rh;
peb(1).MarkerFaceColor = c_lh;
peb(2).MarkerFaceColor = c_rh;

xlim ([0.5 2.5])
ylim ([0.7 1])
ylabel 'ACC [%]'
xticks ([1 2])
xticklabels ({'+-60�','+-120�'})
legend({'Control','Stroke'})
box off
sgtitle 'Figure 2. Behavior performance.'
save_fig(gcf,PATHOUT_plots,'rep_yan12_fig2','FigSize',[0 0 24 20]);



%% replicate ERP analysis

%load ERPs
load ([PATHIN_ERP 'ERPall_finICA.mat']);
idx_bad_part = ~contains({ERP_all.ID},{RT_ALL.ID});
ERP_all(idx_bad_part) = []; %remove bad participants
load ([PATHIN_ERP 'cfg.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
SUBJ = str2double(extractAfter({ERP_all.ID},'SUBJ'));

% define groups
idx_old     = SUBJ < 70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;

%% get data ERP
ep_time = linspace(cfg.ERP.ep_time(1),cfg.ERP.ep_time(2),200)*1000;
idx_bl = ep_time >= -500 & ep_time <= -100;
idx_beg = ep_time >= 0 & ep_time <= 300;
idx_mdl = ep_time >= 300 & ep_time <= 800;
idx_end = ep_time >= 800 & ep_time <= 1200;


for s = 1:numel(ERP_all)
    
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

    dat_erp(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_cor),3),2)); 
    dat_erp(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_cor),3),2)); 
    dat_erp(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_cor),3),2)); 
    dat_erp(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_cor),3),2)); 


end

%% Plot ERP results

nms_time_period = {'P200 [0:300ms]','P300 [300:800ms]'};

% replicate Fig 4
dat_erp_f4= dat_erp(:,:,[2,3]);

close all
figure
i = 1;
for p = [1 3]

    subplot(2,2,p)
    topoplot(mean(dat_erp_f4(idx_old,:,i)),chanlocs,'maplimits',[-6 6]);
    title (['O: ' nms_time_period{i}])
    
    subplot(2,2,1+p)
    topoplot(mean(dat_erp_f4(idx_stroke,:,i)),chanlocs,'maplimits',[-6 6]);
    title (['S: ' nms_time_period{i}])    
    
    i = i+1;
end

sgtitle 'Figure 4. Brain mappings of ERP'
save_fig(gcf,PATHOUT_plots,'rep_yan12_fig4','FigSize',[0 0 25 20]);


%% replicate ERD results


