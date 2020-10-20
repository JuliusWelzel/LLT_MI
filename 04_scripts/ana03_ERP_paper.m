%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_ERP  = [MAIN '02_data\ana03_ERPs\'];
PATHIN_RTs  = [MAIN '02_data\03_RTs\'];

PATHOUT_plots = [MAIN '03_plots\ana03_ERPs\ANOVA\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_plots)
    mkdir(PATHOUT_plots);
end

%% Gather all available data with clean EEG & RTs
% load RTs

% load all ERP data
load ([PATHIN_ERP 'ERPall.mat']);
load ([PATHIN_ERP 'cfg.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
ep_time = linspace(cfg.ERP.ep_time(1),cfg.ERP.ep_time(2),diff(cfg.ERP.ep_time)*cfg.ERP.resam)*1000;


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
idx_time_P3 = dsearchn(ep_time',[450]'):dsearchn(ep_time',[600]');
idx_ROI = ismember({chanlocs.labels},{'P3','Pz','P4'});


for s = 1:length(ERP_all)

    clear ep_stim
    clear idx_art

    ep_stim = splitLLTstim(ERP_all(s).pics);
    for t = 1:length(ep_stim)
        idx_time_bl_start   = dsearchn(ep_time',cfg.ERP.BL(1));
        idx_ep_end          = dsearchn(ep_time',1500);
        if idx_ep_end == 1;idx_ep_end = dsearchn(ep_time',ep_time(end));end
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


    % group assignment, only consider artifact free epochs
    ERP_all(s).group = ones(1,size(ep_stim,2))*idx_group(s);

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

    tmp_ang(ERP_all(s).idx_s_60) = string('60�');
    tmp_ang(ERP_all(s).idx_s_120) = string('120�');
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


%% Hand ERPs

close all
figure
%%%%%%%%%%%%%%%%%%%%%% plot lateral/medial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1)
plot(ep_time,mean(dat_trial_hERP_lat(idx_old,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_hERP_med(idx_old,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,3)
plot(ep_time,mean(dat_trial_hERP_lat(idx_young,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_hERP_med(idx_young,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,5)
plot(ep_time,mean(dat_trial_hERP_lat(idx_stroke,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_hERP_med(idx_stroke,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper;
    title 'STROKE'

%%%%%%%%%%%%%%%%%%%%%% plot 60/120  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)
plot(ep_time,mean(dat_trial_hERP_s60(idx_old,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_hERP_s120(idx_old,:)),'Color',c_rh)
    legend ({'60�','120�'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,4)
plot(ep_time,mean(dat_trial_hERP_s60(idx_young,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_hERP_s120(idx_young,:)),'Color',c_rh)
    legend ({'60�','120�'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,6)
plot(ep_time,mean(dat_trial_hERP_s60(idx_stroke,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_hERP_s120(idx_stroke,:)),'Color',c_rh)
    legend ({'60�','120�'})
    title 'STROKE'
    legend_ERP_paper

sgtitle (['Hand ERP [450-600 ms] at P3,Pz,P4'])
save_fig(gcf,PATHOUT_plots,'grp_hERPs_ANOVA')


% ANOVA peak & mean (group [3] x con [2] x ang [2])

% Hand
dat_amp_hand_P3 = [dat_trial_hERP.P3];

idx_all_stroke = [dat_trial_hERP.group] == 0;
idx_all_old = [dat_trial_hERP.group] == 1;
idx_all_young = [dat_trial_hERP.group] == 2;
nms_all_grp(idx_all_stroke)     = string('stroke');
nms_all_grp(idx_all_old)        = string('old');
nms_all_grp(idx_all_young)      = string('young');

nms_all_con = [dat_trial_hERP.condn];
nms_all_ang = [dat_trial_hERP.angle];


[~,tbl,stats] = anovan(dat_amp_hand_P3,{nms_all_grp,nms_all_con,nms_all_ang},'model','interaction','varnames',{'Group','Condition','Angle'})
multcompare(stats,'Dimension',[1 2 3])
title 'Multicompare hand ERPs'
xlabel 'mean Amplitude [\muV]'
save_fig(gcf,PATHOUT_plots,'mltcmpr_hERPs_ANOVA');
writecell(tbl,[PATHOUT_plots 'tbl_hERPs_ANOVA.xls']);

tmp_tab = struct2table(dat_trial_hERP);
writetable(tmp_tab,[PATHOUT_plots 'st_hERP.xls']);


%% Feet ERPs

close all
figure
%%%%%%%%%%%%%%%%%%%%%% plot lateral/medial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1)
plot(ep_time,mean(dat_trial_fERP_lat(idx_old,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_fERP_med(idx_old,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,3)
plot(ep_time,mean(dat_trial_fERP_lat(idx_young,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_fERP_med(idx_young,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,5)
plot(ep_time,mean(dat_trial_fERP_lat(idx_stroke,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_trial_fERP_med(idx_stroke,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend_ERP_paper;
    title 'STROKE'

%%%%%%%%%%%%%%%%%%%%%% plot 60/120  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)
plot(ep_time,mean(dat_trial_fERP_s60(idx_old,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_fERP_s120(idx_old,:)),'Color',c_rh)
    legend ({'60�','120�'})
    legend_ERP_paper
    title 'OLD'

subplot(3,2,4)
plot(ep_time,mean(dat_trial_fERP_s60(idx_young,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_fERP_s120(idx_young,:)),'Color',c_rh)
    legend ({'60�','120�'})
    legend_ERP_paper
    title 'YOUNG'

subplot(3,2,6)
plot(ep_time,mean(dat_trial_fERP_s60(idx_stroke,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_trial_fERP_s120(idx_stroke,:)),'Color',c_rh)
    legend ({'60�','120�'})
    title 'STROKE'
    legend_ERP_paper

sgtitle (['Feet ERP [450-600 ms] at P3,Pz,P4'])
save_fig(gcf,PATHOUT_plots,'grp_fERPs_ANOVA')

%% Prep data for JASP
% ANOVA mean ERP

% Feet
clear nms_all_grp
clear nms_all_con
clear nms_all_ang

dat_amp_feet_P3 = [dat_trial_fERP.P3];

clear nms_all_grp
idx_all_stroke  = [dat_trial_hERP.group] == 0;
idx_all_old     = [dat_trial_hERP.group] == 1;
idx_all_young   = [dat_trial_hERP.group] == 2;
nms_all_grp(idx_all_stroke)     = string('stroke');
nms_all_grp(idx_all_old)        = string('old');
nms_all_grp(idx_all_young)      = string('young');

nms_all_con = [dat_trial_fERP.condn];
nms_all_ang = [dat_trial_fERP.angle];


tab_hand    = table(dat_amp_hand_P3',nms_all_id',nms_all_grp',nms_all_con',nms_all_ang','VariableNames',{'P3','ID','Group','Condition','Angle'});
tab_hand_con   = unstack(tab_hand,'P3','Condition','AggregationFunction',@mean)
writetable(tab_hand_con,[PATHOUT_data 'hERP_P3_condition.csv']);

tab_hand    = table(dat_amp_hand_P3',nms_all_id',nms_all_grp',nms_all_con',nms_all_ang','VariableNames',{'P3','ID','Group','Condition','Angle'});
tab_hand_ang   = unstack(tab_hand,'P3','Angle','AggregationFunction',@mean)
writetable(tab_hand_ang,[PATHOUT_data 'hERP_P3_angle.csv']);

% Foot
dat_amp_foot_P3 = [dat_trial_fERP.P3];

clear nms_all_grp
idx_all_stroke  = [dat_trial_fERP.group] == 0;
idx_all_old     = [dat_trial_fERP.group] == 1;
idx_all_young   = [dat_trial_fERP.group] == 2;
nms_all_grp(idx_all_stroke)     = string('stroke');
nms_all_grp(idx_all_old)        = string('old');
nms_all_grp(idx_all_young)      = string('young');

nms_all_con = [dat_trial_fERP.condn];
nms_all_ang = [dat_trial_fERP.angle];
nms_all_id  = [dat_trial_fERP.ID];

tab_foot    = table(dat_amp_foot_P3',nms_all_id',nms_all_grp',nms_all_con',nms_all_ang','VariableNames',{'P3','ID','Group','Condition','Angle'});
tab_foot_con   = unstack(tab_foot,'P3','Condition','AggregationFunction',@mean)
writetable(tab_foot_con,[PATHOUT_data 'fERP_P3_condition.csv']);

tab_foot    = table(dat_amp_foot_P3',nms_all_id',nms_all_grp',nms_all_con',nms_all_ang','VariableNames',{'P3','ID','Group','Condition','Angle'});
tab_foot_ang   = unstack(tab_foot,'P3','Angle','AggregationFunction',@mean)
writetable(tab_foot_ang,[PATHOUT_data 'fERP_P3_angle.csv']);


%% Plot number of trials which are rejected

for s = 1:size(ERP_all,2)
    var_trl_art(s) = mean(sum(ERP_all(s).idx_art(idx_ROI,:)) > 0) * 100;
end


dat_art_bp = {var_trl_art(idx_stroke),var_trl_art(idx_old),var_trl_art(idx_young)};
close all
singleBoxplot(dat_art_bp)
tune_BP([c_stroke;c_old;c_young])
    xticks ([1 2 3])
    xticklabels (nms_group)
    ylabel 'Number of artifactual epochs [%]'

save_fig(gcf,PATHOUT_plots,'grp_art_eps')
