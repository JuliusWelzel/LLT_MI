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


%% replicate figure 6 

load([PATHIN_WVLT 'cfg.mat'])

% get time epcoh data and indices for time period
ep_time     = cfg_el.times;
idx_bl      = ep_time >= -500 & ep_time <= -100;
idx_beg     = ep_time >= 0 & ep_time <= 300;
idx_mdl     = ep_time >= 300 & ep_time <= 800;
idx_end     = ep_time >= 800 & ep_time <= 2500;

% indeices for freq band
idx_aplha   = cfg_el.freqs >= 8 & cfg_el.freqs <= 12;
idx_beta    = cfg_el.freqs >= 15 & cfg_el.freqs <= 25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                CONTINUE WITH EXTRACTING TIMED ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load RTs 
load([PATHIN_RTs 'RT_ALL.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);

% for ERD analysis
dat_lh_alpha    = []; %sub x chan x t_period [56 x 24 x 4 (BL,BEG,MID,END)]
dat_rh_alpha    = [];
dat_lh_beta     = [];
dat_rh_beta     = [];

% for RT vs ERD analysis
dat_RTERD_alpha     = []; %sub x t_period x con [56 x 4 (BL,BEG,MID,END) x 2 (RT, ERD)]
dat_RTERD_beta      = [];


for s = 1:numel(SUBJ)
    
    % laod wvlt data
    load([PATHIN_WVLT nms_SUBJ{s}]);
    display (['Working on Subject ' num2str(SUBJ(s))])
    
    %normalize trial data 
    tmp_dat     = dat_wvlt.powspctrm;  % -> size powspctrm [96 x 24 x 35 x 2000]
    tmp_bl      = squeeze(mean(tmp_dat(:,:,:,idx_bl),4));
    tmp_dat_dB  = 10*log10(tmp_dat./tmp_bl); %transform to dB // dB = 10*log10 (signal/baseline)


    % grab relevant RT data form RT_ALL
    % experimental stimuli as in ana04_indv_wavelets
    
    idx_sub_in_RTall    = strcmp({RT_ALL.ID},dat_wvlt.id);
    if sum(idx_sub_in_RTall) == 0; continue; end; % skip this participant if not found in RT_ALL
    idx_cor             = logical(RT_ALL(idx_sub_in_RTall).acc(11:end));
    
    % inidces laterality
    idx_lat      = contains(dat_wvlt.events(:,1),'lh') & contains(dat_wvlt.events(:,4),{'240','300'}) |...
                   contains(dat_wvlt.events(:,1),'rh') & contains(dat_wvlt.events(:,4),{'60','120'});
    idx_med      = contains(dat_wvlt.events(:,1),'lh') & contains(dat_wvlt.events(:,4),{'60','120'}) |...
                   contains(dat_wvlt.events(:,1),'rh') & contains(dat_wvlt.events(:,4),{'240','300'});
    % inidces handiness
    idx_lh      = contains(dat_wvlt.events(:,1),'lh');
    idx_rh      = contains(dat_wvlt.events(:,1),'rh');
        
    
    %% HANDINESS
    % create 3D matrix -> sub x chan x t_period [56 x 24 x 4 (BL,BEG,MID,END)]
    % alpha_lh all time periods
    dat_lh_alpha(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_aplha,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_lh_alpha(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_aplha,idx_beg),4),3)));
    dat_lh_alpha(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_aplha,idx_mdl),4),3))); 
    dat_lh_alpha(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_aplha,idx_end),4),3))); 
   
    % alpha_rh all time periods
    dat_rh_alpha(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_rh',:,idx_aplha,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_rh_alpha(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_rh',:,idx_aplha,idx_beg),4),3)));
    dat_rh_alpha(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_rh',:,idx_aplha,idx_mdl),4),3))); 
    dat_rh_alpha(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_rh',:,idx_aplha,idx_end),4),3))); 

    % beta_lh all time periods
    dat_lh_beta(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_lh_beta(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_beg),4),3)));
    dat_lh_beta(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_mdl),4),3))); 
    dat_lh_beta(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_end),4),3))); 

    % alpha_rh all time periods
    dat_rh_beta(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_rh_beta(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_beg),4),3)));
    dat_rh_beta(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_mdl),4),3))); 
    dat_rh_beta(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lh',:,idx_beta,idx_end),4),3))); 

    
    %% LATERALITY
    % create 3D matrix -> sub x chan x t_period [56 x 24 x 4 (BL,BEG,MID,END)]
    % alpha_lat all time periods
    dat_lat_alpha(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_aplha,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_lat_alpha(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_aplha,idx_beg),4),3)));
    dat_lat_alpha(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_aplha,idx_mdl),4),3))); 
    dat_lat_alpha(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_aplha,idx_end),4),3))); 
   
    % alpha_med all time periods
    dat_med_alpha(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_aplha,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_med_alpha(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_aplha,idx_beg),4),3)));
    dat_med_alpha(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_aplha,idx_mdl),4),3))); 
    dat_med_alpha(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_aplha,idx_end),4),3))); 

    % beta_lat all time periods
    dat_lat_beta(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_beta,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_lat_beta(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_beta,idx_beg),4),3)));
    dat_lat_beta(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_beta,idx_mdl),4),3))); 
    dat_lat_beta(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_lat',:,idx_beta,idx_end),4),3))); 

    % alpha_med all time periods
    dat_med_beta(s,:,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_beta,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_med_beta(s,:,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_beta,idx_beg),4),3)));
    dat_med_beta(s,:,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_beta,idx_mdl),4),3))); 
    dat_med_beta(s,:,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor & idx_med',:,idx_beta,idx_end),4),3))); 

    
    %% RTs at CPz
    % create 3D matrix -> sub x t_period x con [56 x 4 (BL,BEG,MID,END) x 2 (RT, ERD)]
    % for RT vs ERD analysis
    idx_chan                = find(strcmp(dat_wvlt.label,'CPz'));
    tmp_RT                  = RT_ALL(idx_sub_in_RTall).SO_ms(11:end);
    dat_RTERD.RT(s)         = mean(tmp_RT(idx_cor)); % time period & RT
    
    %alpha ERD
    dat_RTERD.alpha(s,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_aplha,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.alpha(s,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_aplha,idx_beg),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.alpha(s,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_aplha,idx_mdl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.alpha(s,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_aplha,idx_end),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    
    % beta ERD
    dat_RTERD.beta(s,1)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_beta,idx_bl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.beta(s,2)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_beta,idx_beg),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.beta(s,3)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_beta,idx_mdl),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]
    dat_RTERD.beta(s,4)    = squeeze(mean(mean(mean(tmp_dat_dB(idx_cor,idx_chan,idx_beta,idx_end),4),3))); %-> size powspctrm [96 x 24 x 35 x 2000]

end


%% save data

save([PATHOUT_YAN 'ERD_hand.mat'],'dat_lh_alpha','dat_rh_alpha','dat_lh_beta','dat_rh_beta');
save([PATHOUT_YAN 'ERD_latmed.mat'],'dat_lat_alpha','dat_med_alpha','dat_lat_beta','dat_med_beta');
save([PATHOUT_YAN 'ERD_RT.mat'],'dat_RTERD');

%% Plot results // replicate figure 6

nms_time_period = {'BL [-500:-200ms]','BEG [0:300ms]','MDL [300:800ms]','END [800:2500ms]'};
map_limits = [-8 8];


% LH aplha
close all
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_lh_alpha(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB % function changed
    
    subplot(3,4,p+4)
    topoplot(mean(dat_lh_alpha(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_lh_alpha(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB    

    i = i+1;
end

sgtitle 'Figure 6. Time-course \alpha-ERD LH'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_alpha_LH','FigSize',[0 0 30 20]);

%RH alpha
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_rh_alpha(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_rh_alpha(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_rh_alpha(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \alpha-ERD RH'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_alpha_RH','FigSize',[0 0 30 20]);

% LH beta
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_lh_beta(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_lh_beta(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_lh_beta(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \beta-ERD LH'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_beta_LH','FigSize',[0 0 30 20]);

%RH alpha
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_rh_beta(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_rh_beta(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_rh_beta(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \beta-ERD RH'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_beta_RH','FigSize',[0 0 30 20]);


%% plot time course for lateral/ medial
map_limits = [-8 8];

% Lat alpha
close all
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_lat_alpha(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB
    
    subplot(3,4,p+4)
    topoplot(mean(dat_lat_alpha(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_lat_alpha(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB    

    i = i+1;
end

sgtitle 'Figure 6. Time-course \alpha-ERD lateral'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_alpha_lateral','FigSize',[0 0 30 20]);

%medial alpha
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_med_alpha(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_med_alpha(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_med_alpha(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \alpha-ERD medial'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_alpha_medial','FigSize',[0 0 30 20]);

% lateral beta
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_lat_beta(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_lat_beta(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_lat_beta(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \beta-ERD lateral'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_beta_lateral','FigSize',[0 0 30 20]);

%RH alpha
figure
i = 1;
for p = 1:numel(nms_time_period)

    subplot(3,4,p)
    topoplot(mean(dat_med_beta(idx_stroke,:,i)),chanlocs,'maplimits',map_limits);
    title (['S: ' nms_time_period{i}])
    put_CB

    subplot(3,4,p+4)
    topoplot(mean(dat_med_beta(idx_old,:,i),1),chanlocs,'maplimits',map_limits);
    title (['O: ' nms_time_period{i}])    
    put_CB

    subplot(3,4,p+8)
    topoplot(mean(dat_med_beta(idx_young,:,i),1),chanlocs,'maplimits',map_limits);
    title (['Y: ' nms_time_period{i}])    
    put_CB

    i = i+1;
end

sgtitle 'Figure 6. Time-course \beta-ERD medial'
save_fig(gcf,PATHOUT_plots,'ERD_timecourse_beta_medial','FigSize',[0 0 30 20]);




%% CONNY ERD STATS
% Könntest du für die TF-Maps Differenz-Plots machen: medial-lateral jeweils 
% für jede Gruppe und Alt-Jung bzw. Alt-Stroke jeweils für medial und lateral getrennt? 
% Mann könnte die einzelnen Elektroden-Positionen auch gleich noch pasierend auf 
% p-Werten markieren (dafür müsste man jeweils einen T-Test für unabhängige oder 
% abhängige Stichproben pro Kanal rechnen und dann könnte man zB jeden Kanal ja 
% nach p-Wert farblich markieren oder in der Größe variiren. Man kann p-WErte natürlich auch mappen). 



































