%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data
% Data: nicolor.c_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_ERP  = [MAIN '02_data\04_ERPs\'];
PATHIN_RTs  = [MAIN '02_data\03_RTs\'];

PATHOUT_plots   = [MAIN '04_plots\04_ERPs\'];
PATHOUT_data    = [MAIN '02_data\04_ERPs\JASP\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_plots);mkdir(PATHOUT_plots);end
if ~isdir(PATHOUT_data);mkdir(PATHOUT_data);end

%% Gather all available data with clean EEG & RTs
% load RTs

% load all ERP data
load ([PATHIN_ERP 'ERPall.mat']);
load ([PATHIN_ERP 'cfg.mat']);
cfg.ep_time = linspace(cfg.ERP.ep_time(1),cfg.ERP.ep_time(2),diff(cfg.ERP.ep_time)*cfg.ERP.resam)*1000;


SUBJ = str2double(extractAfter({ERP_all.ID},'SUBJ'));

% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;
idx_group   = SUBJ;
idx_group(idx_stroke)   = 0;
idx_group(idx_old)      = 1;
idx_group(idx_young)    = 2;

% groups vars for plotting
nms_group   = {'stroke','old','young'};

%% Explore ERP dist
%add relevant indices
idx_time_P2 = dsearchn(cfg.ep_time',[180]'):dsearchn(cfg.ep_time',[260]');
idx_time_P3 = dsearchn(cfg.ep_time',[300]'):dsearchn(cfg.ep_time',[600]');
idx_ROI = ismember({cfg.ERP.chanlocs.labels},{'P3','Pz','P4'});


for s = 1:length(ERP_all)

    clear ep_stim
    clear idx_art

    ep_stim = splitLLTstim(ERP_all(s).pics);
    for t = 1:length(ep_stim)
        idx_time_bl_start   = dsearchn(cfg.ep_time',cfg.ERP.BL(1));
        idx_ep_end          = dsearchn(cfg.ep_time',1500);
        if idx_ep_end == 1;idx_ep_end = dsearchn(cfg.ep_time',cfg.ep_time(end));end
        idx_art(:,t) = squeeze(max(ERP_all(s).mERP(:,idx_time_bl_start:idx_ep_end,t),[],2))>100 | squeeze(min(ERP_all(s).mERP(:,idx_time_bl_start:idx_ep_end,t),[],2))<-100;
    end

    ERP_all(s).idx_art = idx_art;
    idx_art = mean(ERP_all(s).idx_art(idx_ROI,:));
    idx_cor = ERP_all(s).idx_cor' & ~idx_art';

    % find indices
    % idx lat/med
    ERP_all(s).idx_lat      = contains(ep_stim(1,:),'l')' & contains(ep_stim(4,:),{'240','300'})'|...
                            contains(ep_stim(1,:),'r')' & contains(ep_stim(4,:),{'60','120'})';
    ERP_all(s).idx_med      = contains(ep_stim(1,:),'l')' & contains(ep_stim(4,:),{'60','120'})'|...
                            contains(ep_stim(1,:),'r')' & contains(ep_stim(4,:),{'240','300'})';

    % idx hand/foot
    ERP_all(s).idx_hand      = contains(ep_stim(1,:),'h')';
    ERP_all(s).idx_foot      = contains(ep_stim(1,:),'f')';

    %angle
    ERP_all(s).idx_s_60      = strcmp(ep_stim(4,:),'60')'| strcmp(ep_stim(4,:),'300')';
    ERP_all(s).idx_s_120     = strcmp(ep_stim(4,:),'120')'|  strcmp(ep_stim(4,:),'240')';

    % idx in depth rotation
    ERP_all(s).idx_id_0      = contains(ep_stim(3,:),'0')';
    ERP_all(s).idx_id_60     = contains(ep_stim(3,:),'60')';
    ERP_all(s).idx_id_300    = contains(ep_stim(3,:),'300')';


    % group assignment, only consider artifact free epochs and ID
    ERP_all(s).group        = ones(1,size(ep_stim,2))*idx_group(s);
    dat_trial_hERP(s).ID    = cell(1,sum(idx_cor & ERP_all(s).idx_hand));    
    dat_trial_hERP(s).ID(:) = {ERP_all(s).ID};
    dat_trial_fERP(s).ID    = cell(1,sum(idx_cor & ERP_all(s).idx_foot));    
    dat_trial_fERP(s).ID(:) = {ERP_all(s).ID};

    %single trial vs single sub
    dat_trial_hERP(s).P3 = squeeze(nanmean(nanmean(ERP_all(s).mERP(idx_ROI,idx_time_P3,idx_cor' & ERP_all(s).idx_hand'))))';
    dat_trial_fERP(s).P3 = squeeze(nanmean(nanmean(ERP_all(s).mERP(idx_ROI,idx_time_P3,idx_cor' & ERP_all(s).idx_foot'))))';

    % add conditions for ANOVA
    dat_trial_hERP(s).group = ERP_all(s).group(idx_cor & ERP_all(s).idx_hand);
    dat_trial_fERP(s).group = ERP_all(s).group(idx_cor & ERP_all(s).idx_foot);

    % tmp condn idx
    tmp_con(ERP_all(s).idx_lat) = string('lat');
    tmp_con(ERP_all(s).idx_med) = string('med');
    dat_trial_hERP(s).condn = tmp_con(idx_cor & ERP_all(s).idx_hand);
    dat_trial_fERP(s).condn = tmp_con(idx_cor & ERP_all(s).idx_foot);

    tmp_ang(ERP_all(s).idx_s_60) = string('60°');
    tmp_ang(ERP_all(s).idx_s_120) = string('120°');
    dat_trial_hERP(s).angle = tmp_ang(idx_cor & ERP_all(s).idx_hand);
    dat_trial_fERP(s).angle = tmp_ang(idx_cor & ERP_all(s).idx_foot);

    %single sub average ERP data hand
    dat_trial_hERP_lat(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_hand' & ERP_all(s).idx_lat')),3));
    dat_trial_hERP_med(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_hand' & ERP_all(s).idx_med')),3));

    dat_trial_hERP_s60(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_hand' & ERP_all(s).idx_s_60')),3));
    dat_trial_hERP_s120(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_hand' & ERP_all(s).idx_s_120')),3));

    %single sub average ERP data feet
    dat_trial_fERP_lat(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_foot' & ERP_all(s).idx_lat')),3));
    dat_trial_fERP_med(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_foot' & ERP_all(s).idx_med')),3));

    dat_trial_fERP_s60(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_foot' & ERP_all(s).idx_s_60')),3));
    dat_trial_fERP_s120(s,:) = squeeze(mean(mean(ERP_all(s).mERP(idx_ROI,:,idx_cor' & ERP_all(s).idx_foot' & ERP_all(s).idx_s_120')),3));

end


%% Grand average

cfg.color   = color;
channel_oi  = {'P3','Pz','P4'};
grnd_avg_plot(ERP_all,channel_oi,cfg)
save_fig(gcf,PATHOUT_plots,'ERP_overview')

%% Hand ERPs

close all
figure
%%%%%%%%%%%%%%%%%%%%%% plot lateral/medial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1)
plot(cfg.ep_time,mean(dat_trial_hERP_lat(idx_old,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_med(idx_old,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,3)
plot(cfg.ep_time,mean(dat_trial_hERP_lat(idx_young,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_med(idx_young,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,5)
plot(cfg.ep_time,mean(dat_trial_hERP_lat(idx_stroke,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_med(idx_stroke,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper;
    title 'STROKE'

%%%%%%%%%%%%%%%%%%%%%% plot 60/120  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)
plot(cfg.ep_time,mean(dat_trial_hERP_s60(idx_old,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_s120(idx_old,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,4)
plot(cfg.ep_time,mean(dat_trial_hERP_s60(idx_young,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_s120(idx_young,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,6)
plot(cfg.ep_time,mean(dat_trial_hERP_s60(idx_stroke,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_hERP_s120(idx_stroke,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    title 'STROKE'
    legend_ERP_paper

sgtitle (['Hand ERP [300-600 ms] at P3,Pz,P4'])
save_fig(gcf,PATHOUT_plots,'grp_hERPs_ANOVA')


%% Feet ERPs

close all
figure
%%%%%%%%%%%%%%%%%%%%%% plot lateral/medial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1)
plot(cfg.ep_time,mean(dat_trial_fERP_lat(idx_old,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_med(idx_old,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,3)
plot(cfg.ep_time,mean(dat_trial_fERP_lat(idx_young,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_med(idx_young,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,5)
plot(cfg.ep_time,mean(dat_trial_fERP_lat(idx_stroke,:)),'Color',color.c_lat)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_med(idx_stroke,:)),'Color',color.c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper;
    title 'STROKE'

%%%%%%%%%%%%%%%%%%%%%% plot 60/120  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)
plot(cfg.ep_time,mean(dat_trial_fERP_s60(idx_old,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_s120(idx_old,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,4)
plot(cfg.ep_time,mean(dat_trial_fERP_s60(idx_young,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_s120(idx_young,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,6)
plot(cfg.ep_time,mean(dat_trial_fERP_s60(idx_stroke,:)),'Color',color.c_lh)
hold on
plot(cfg.ep_time,mean(dat_trial_fERP_s120(idx_stroke,:)),'Color',color.c_rh)
    legend ({'60°','120°'})
    title 'STROKE'
    legend_ERP_paper

sgtitle (['Feet ERP [300-600 ms] at P3,Pz,P4'])
save_fig(gcf,PATHOUT_plots,'grp_fERPs_ANOVA')


%% Prep data for JASP

c = 1;
for s = 1:numel(ERP_all) % Loop over all subjects
    

    for e = 1:size(ERP_all(s).mERP,3) % Loop over all epochs pers subject
        
        % ident vars
        ERP(c).ID       = ERP_all(s).ID;
        ERP(c).group    = LLTgroupname(ERP_all(s).ID);
        
        % epoch vars
        clear ep_stim % delete from previous iteration
        ep_stim             = splitLLTstim(ERP_all(s).pics(e)); % Get all stim infos per participant
        [ERP(c).extremity ERP(c).condition ERP(c).angle]  = LLTepochinfo(ep_stim);
        ERP(c).response     = ERP_all(s).idx_cor(e);
        
        % check for artifact
        idx_time_bl_start   = dsearchn(cfg.ep_time',cfg.ERP.BL(1));
        idx_ep_end          = dsearchn(cfg.ep_time',1500);
        if idx_ep_end == 1;idx_ep_end = dsearchn(cfg.ep_time',cfg.ep_time(end));end
        ERP(c).artefact = max(squeeze(max(ERP_all(s).mERP(idx_ROI,idx_time_bl_start:idx_ep_end,e),[],2))>100 | squeeze(min(ERP_all(s).mERP(idx_ROI,idx_time_bl_start:idx_ep_end,e),[],2))<-100);
        
        ERP(c).P2 = squeeze(nanmean(nanmean(ERP_all(s).mERP(idx_ROI,idx_time_P2,e))))';
        ERP(c).P3 = squeeze(nanmean(nanmean(ERP_all(s).mERP(idx_ROI,idx_time_P3,e))))';

        c = c+1;
    end
    
end
        
ERP = struct2table(ERP);
writetable(ERP,[PATHOUT_data 'ERPs.csv']);

ERPcon   = unstack(ERP,{'P2','P3'},'condition','AggregationFunction',@mean);
ERPang   = unstack(ERP,{'P2','P3'},'angle','AggregationFunction',@mean);
ERPext   = unstack(ERP,{'P2','P3'},'extremity','AggregationFunction',@mean);

writetable(ERPcon,[PATHOUT_data 'ERP_condition.csv']);
writetable(ERPang,[PATHOUT_data 'ERP_angle.csv']);
writetable(ERPext,[PATHOUT_data 'ERP_extremity.csv']);



%% Plot number of trials which are rejected

for s = 1:size(ERP_all,2)
    var_trl_art(s) = mean(sum(ERP_all(s).idx_art(idx_ROI,:)) > 0) * 100;
end


dat_art_bp = {var_trl_art(idx_stroke),var_trl_art(idx_old),var_trl_art(idx_young)};
close all
singleBoxplot(dat_art_bp)
tune_BP([color.c_stroke;color.c_old;color.c_young])
    xticks ([1 2 3])
    xticklabels (nms_group)
    ylabel 'Number of artifactual epochs [%]'

save_fig(gcf,PATHOUT_plots,'grp_art_eps')
