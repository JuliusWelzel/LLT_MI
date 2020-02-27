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
% load RTs 
load([PATHIN_RTs 'RT_ALL.mat']);

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


for s = 1:numel(RT_ALL)
    
    nms_f = fieldnames(RT_ALL);
    for f = 2:numel(nms_f)-1
        RT_ALL(s).(nms_f{f})(1:10) = [];
    end
    RT_ALL(s).stim(1:10,:) = [];
        
    
    % handiness
    idx_lh      = strcmp(RT_ALL(s).stim(:,1),'lh');
    idx_rh      = strcmp(RT_ALL(s).stim(:,1),'rh');
    idx_cor     = logical(RT_ALL(s).acc(:));
   
    % create 3D matrix [subs x hand x mod] -> 56 x 2(LH,RH) x 2(RT,ACC)
    dat_b_hand(s,1,1) = nanmean(RT_ALL(s).SO_ms(idx_lh & idx_cor)); 
    dat_b_hand(s,2,1) = nanmean(RT_ALL(s).SO_ms(idx_rh & idx_cor));
    dat_b_hand(s,1,2) = nanmean(RT_ALL(s).acc(idx_lh));
    dat_b_hand(s,2,2) = nanmean(RT_ALL(s).acc(idx_rh));
    
    
    %angle
    idx_60      = strcmp(RT_ALL(s).stim(:,4),'60') & contains(RT_ALL(s).stim(:,1),'h') |...
                    strcmp(RT_ALL(s).stim(:,4),'300') & contains(RT_ALL(s).stim(:,1),'h');   
    idx_120     = strcmp(RT_ALL(s).stim(:,4),'120') & contains(RT_ALL(s).stim(:,1),'h') | ...
                    strcmp(RT_ALL(s).stim(:,4),'240') & contains(RT_ALL(s).stim(:,1),'h');
    
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
xticklabels ({'+-60°','+-120°'})
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
xticklabels ({'+-60°','+-120°'})
legend({'Control','Stroke'})
box off
sgtitle 'Figure 2. Behavior performance.'
% save_fig(gcf,PATHOUT_plots,'rep_yan12_fig2','FigSize',[0 0 24 20]);



%% replicate ERP analysis

%load ERPs
load ([PATHIN_ERP 'ERPall_finICA.mat']); % from ana03_preERP // filter [0.1-20Hz] // 100 srate // CAR // BL [-500 -100]
idx_bad_part = ~contains({ERP_all.ID},{RT_ALL.ID});
ERP_all(idx_bad_part) = []; %remove bad participants
load ([PATHIN_ERP 'cfg.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
idx_chanref = contains({chanlocs.labels},{'TP9','TP10'});
chanlocs(idx_chanref) =  [];
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
    idx_cor     = logical(RT_ALL(idx_sub_in_RTall).acc(:))';
   
    % create 3D matrix [subs x channel x time_period] -> 56 x 24 x 4(BL,BEG,MDL,END)
    dat_tp_lh(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_lh & idx_cor),3),2)); 
    dat_tp_lh(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_lh & idx_cor),3),2)); 
    dat_tp_lh(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_lh & idx_cor),3),2)); 
    dat_tp_lh(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_lh & idx_cor),3),2)); 
    
    dat_tp_rh(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_rh & idx_cor),3),2)); 
    dat_tp_rh(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_rh & idx_cor),3),2)); 
    dat_tp_rh(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_rh & idx_cor),3),2)); 
    dat_tp_rh(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_rh & idx_cor),3),2)); 

    dat_tp(s,:,1) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_bl,idx_cor),3),2)); 
    dat_tp(s,:,2) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_beg,idx_cor),3),2)); 
    dat_tp(s,:,3) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_mdl,idx_cor),3),2)); 
    dat_tp(s,:,4) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_end,idx_cor),3),2)); 


end

%% Plot ERP results
load('cfg.mat');

nms_time_period = {'P200 [0:300ms]','P300 [300:800ms]'};

% replicate Fig 4
dat_erp_f4= dat_tp(:,:,[2,3]);

close all
figure
i = 1;
for p = [1 3]

    subplot(2,2,p)
    topoplot(mean(dat_erp_f4(idx_old,:,i)),chanlocs,'maplimits',[-10 10]);
    title (['O: ' nms_time_period{i}])
    
    subplot(2,2,1+p)
    topoplot(mean(dat_erp_f4(idx_stroke,:,i)),chanlocs,'maplimits',[-10 10]);
    title (['S: ' nms_time_period{i}])    
    
    i = i+1;
end

sgtitle 'Figure 4. Brain mappings of ERP'
save_fig(gcf,PATHOUT_plots,'rep_yan12_fig4','FigSize',[0 0 25 20]);


%% replicate yan12 figure 3 eith P4 & Pz
% from ana03_preERP // filter [0.1-20Hz] // 100 srate // CAR // BL [-500 -100]

idx_Pz = strcmp({chanlocs.labels},'Pz');
idx_P4 = strcmp({chanlocs.labels},'P4');


for s = 1:numel(ERP_all)
    
    % get stimuli information
    
    stim = splitLLTstim(ERP_all(s).pics); % 4 positions - 1. Laterality [LH, RH, LF, RF] // 2. View [back,palm] // ...
                                                       % 3. in depth rotation // 4. lateral rotation 
    
    idx_lh      = contains(stim(1,:),'lh');
    idx_rh      = contains(stim(1,:),'rh');
    idx_lf      = contains(stim(1,:),'lf');
    idx_rf      = contains(stim(1,:),'rf');
 
    idx_0       = strcmp(stim(4,:),'0');
    idx_60      = strcmp(stim(4,:),'60') | strcmp(stim(4,:),'300');
    idx_120     = strcmp(stim(4,:),'120')| strcmp(stim(4,:),'240');
    idx_180     = strcmp(stim(4,:),'180');
    
    % correct stimuli
    idx_sub_in_RTall    = strcmp({RT_ALL.ID},{ERP_all(s).ID});
    idx_cor             = logical(RT_ALL(idx_sub_in_RTall).acc(:))';
    idx_art             = squeeze(sum(ERP_all(s).mERP(idx_Pz,:,:) > 100,2))';
   
    % create 2D matrix [subs x samples] -> 56 x 200
    % laterality
    dat_erp_lh(s,:) = mean(ERP_all(s).mERP(idx_Pz,:,idx_lh & ~idx_art & idx_cor),3); 
    dat_erp_rh(s,:) = mean(ERP_all(s).mERP(idx_Pz,:,idx_rh & ~idx_art & idx_cor),3); 
    dat_erp_lf(s,:) = mean(ERP_all(s).mERP(idx_Pz,:,idx_lf & ~idx_art & idx_cor),3); 
    dat_erp_rf(s,:) = mean(ERP_all(s).mERP(idx_Pz,:,idx_rf & ~idx_art & idx_cor),3); 
     
    %rotation only hands
    dat_erp_0(s,:)      = mean(ERP_all(s).mERP(idx_Pz,:,idx_0 & ~idx_art & idx_cor),3); 
    dat_erp_60(s,:)     = mean(ERP_all(s).mERP(idx_Pz,:,idx_60 & ~idx_art & idx_cor),3); 
    dat_erp_120(s,:)    = mean(ERP_all(s).mERP(idx_Pz,:,idx_120 & ~idx_art & idx_cor),3); 
    dat_erp_180(s,:)    = mean(ERP_all(s).mERP(idx_Pz,:,idx_180 & ~idx_art & idx_cor),3); 

   
end


%% plot results
close all
figure

%degree
subplot(3,2,1)
plot(ep_time,mean(dat_erp_0(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_60(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_120(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_180(idx_stroke,:)))

title 'Stroke patient'
erp_leg({'0°','\pm60°','\pm120°','180°'})


subplot(3,2,3)
plot(ep_time,mean(dat_erp_0(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_60(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_120(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_180(idx_old,:)))

title 'Old subjects'
erp_leg({'0°','\pm60°','\pm120°','180°'})

subplot(3,2,5)
plot(ep_time,mean(dat_erp_0(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_60(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_120(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_180(idx_young,:)))

title ' Young subjects'
erp_leg({'0°','\pm60°','\pm120°','180°'})


%Hands
subplot(3,2,2)
plot(ep_time,mean(dat_erp_lh(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_rh(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_lf(idx_stroke,:)))
hold on
plot(ep_time,mean(dat_erp_rf(idx_stroke,:)))

title 'Stroke patient'
erp_leg({'LH','RH','LF','RF'})

subplot(3,2,4)
plot(ep_time,mean(dat_erp_lh(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_rh(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_lf(idx_old,:)))
hold on
plot(ep_time,mean(dat_erp_rf(idx_old,:)))

title 'Old subjects'
erp_leg({'LH','RH','LF','RF'})

subplot(3,2,6)
plot(ep_time,mean(dat_erp_lh(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_rh(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_lf(idx_young,:)))
hold on
plot(ep_time,mean(dat_erp_rf(idx_young,:)))

title ' Young subjects'
erp_leg({'LH','RH','LF','RF'})


sgtitle 'Pz'


save_fig(gcf,PATHOUT_plots,'rep_yan12_fig4_Pz','FigSize',[0 0 25 20]);




%% HAND ERP all groups

c_stroke    = c_rh;
c_old       = c_lh;
c_young     = [149,165,166]/255;

% hand
% stroke
mERP_f_s        = mean([dat_erp_lh(idx_stroke,:);dat_erp_lh(idx_stroke,:)]);
std_mERP_f_s    = std([dat_erp_lh(idx_stroke,:);dat_erp_lh(idx_stroke,:)])/sqrt(sum((idx_stroke)));
%old
mERP_h_o        = mean([dat_erp_lh(idx_old,:);dat_erp_lh(idx_old,:)]);
std_mERP_h_o    = std([dat_erp_lh(idx_old,:);dat_erp_lh(idx_old,:)])/sqrt(sum((idx_old)));
% young
mERP_h_y        = mean([dat_erp_lh(idx_young,:);dat_erp_lh(idx_young,:)]);
std_mERP_h_y    = std([dat_erp_lh(idx_young,:);dat_erp_lh(idx_young,:)])/sqrt(sum((idx_young)));


close all
figure
bl = boundedline(   ep_time, mERP_f_s, std_mERP_f_s,...
                    ep_time, mERP_h_o, std_mERP_h_o,...
                    ep_time, mERP_h_y, std_mERP_h_y,...
    'alpha','cmap',[c_stroke;c_old;c_young]); % alpha makes bounds transparent

legend ([bl(1) bl(2) bl(3)],{'stroke','old','young'})
legend('boxoff')
xlabel 'time [ms]'
ylabel 'amplitude [\muV]'
xlim ([-300 1200])
ylim ([-20 20])
vline(0,'--k');
bl(1).LineWidth = 1;
bl(2).LineWidth = 1;
bl(3).LineWidth = 1;

% Add all our previous improvements:
ax = gca();
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
title 'Hand ERP at Pz'

save_fig(gcf,PATHOUT_plots,['ERP_hand']);


%% HAND ERP all groups

c_stroke    = c_rh;
c_old       = c_lh;
c_young     = [149,165,166]/255;

% hand
% stroke
mERP_f_s        = mean([dat_erp_lf(idx_stroke,:);dat_erp_lf(idx_stroke,:)]);
std_mERP_f_s    = std([dat_erp_lf(idx_stroke,:);dat_erp_lf(idx_stroke,:)])/sqrt(sum((idx_stroke)));
%old
mERP_f_o        = mean([dat_erp_lf(idx_old,:);dat_erp_lf(idx_old,:)]);
std_mERP_f_o    = std([dat_erp_lf(idx_old,:);dat_erp_lf(idx_old,:)])/sqrt(sum((idx_old)));
% young
mERP_f_y        = mean([dat_erp_lf(idx_young,:);dat_erp_lf(idx_young,:)]);
std_mERP_f_y    = std([dat_erp_lf(idx_young,:);dat_erp_lf(idx_young,:)])/sqrt(sum((idx_young)));


close all
figure
bl = boundedline(   ep_time, mERP_f_s, std_mERP_f_s,...
                    ep_time, mERP_f_o, std_mERP_f_o,...
                    ep_time, mERP_f_y, std_mERP_f_y,...
    'alpha','cmap',[c_stroke;c_old;c_young]); % alpha makes bounds transparent

legend ([bl(1) bl(2) bl(3)],{'stroke','old','young'})
legend('boxoff')
xlabel 'time [ms]'
ylabel 'amplitude [\muV]'
xlim ([-300 1200])
ylim ([-20 20])
vline(0,'--k');
bl(1).LineWidth = 1;
bl(2).LineWidth = 1;
bl(3).LineWidth = 1;

% Add all our previous improvements:
ax = gca();
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
title 'Feet ERP at Pz'

save_fig(gcf,PATHOUT_plots,['ERP_feet']);































