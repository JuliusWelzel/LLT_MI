%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_WVLT = [MAIN '02_data\ana04_wavelets\'];

PATHOUT_dat_wvlt = [MAIN '02_data\ana05_terhorst2013\'];
PATHOUT_plots = [MAIN '03_plots\ana05_terhorst2013\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_dat_wvlt)
    mkdir(PATHOUT_dat_wvlt);
    mkdir(PATHOUT_plots);
end

%% Gather all available datacfg_wvltts with clean EEG data
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

%% Extract single trials

% ...   [ We expected to replicate the finding of a stronger central mu ERD for
%       medial stimuli during mental rotation as compared to lateral
%       stimuli in younger participants and extend it to older participants] ...

load([PATHIN_WVLT 'cfg.mat'])

% get time epcoh data and indices for time period
ep_time         = cfg_wvlt.times;
idx_bl          = ep_time >= -500 & ep_time <= -200;
idx_time_erd    = ep_time >= 300 & ep_time <= 800;

% indeices for freq band
idx_alpha   = cfg_wvlt.freqs >= 8 & cfg_wvlt.freqs <= 12;
idx_beta    = cfg_wvlt.freqs >= 15 & cfg_wvlt.freqs <= 25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                CONTINUE WITH EXTRACTING TIMED ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load Chanlocs & wavlets all
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);

%% run loop for all participants

for s = 1:numel(SUBJ)
    
    % load wvlt data
    display (['Working on Subject ' num2str(SUBJ(s))])
    load([PATHIN_WVLT nms_SUBJ{s}]);
    
    % experimental stimuli as in ana04_indv_wavelets
    idx_cor = logical(dat_wvlt.idx_cor)';
    idx_art = dat_wvlt.idx_ep_prune';
    idx_plm = strcmp(dat_wvlt.events(:,2),'palm');
    idx_good = idx_cor & ~idx_art & idx_plm;
    
    % indeices for freq band
    idx_alpha   = dat_wvlt.freq >= 8 & dat_wvlt.freq <= 12;

    % inidces laterality
    idx_lat      = contains(dat_wvlt.events(:,1),'lh') & contains(dat_wvlt.events(:,4),{'240','300'}) |...
                   contains(dat_wvlt.events(:,1),'rh') & contains(dat_wvlt.events(:,4),{'60','120'});
    idx_med      = contains(dat_wvlt.events(:,1),'lh') & contains(dat_wvlt.events(:,4),{'60','120'}) |...
                   contains(dat_wvlt.events(:,1),'rh') & contains(dat_wvlt.events(:,4),{'240','300'});
     
    %% create 3D matrix -> sub x chan x time for alpha/mu

    % for ERD analysis
    dat_lat(s,:,:)   = squeeze(mean(mean(dat_wvlt.powspctrm(idx_lat & idx_good,:,idx_alpha,:),3))); % -> 96 tr x 24 chan x 35 freq x 2000 time
    dat_med(s,:,:)   = squeeze(mean(mean(dat_wvlt.powspctrm(idx_med & idx_good,:,idx_alpha,:),3))); 
    
    
        %% create 2D matrix for corr & err beta stuff -> sub x trial
%     
%     idx_ROI_lft = ismember({chanlocs.labels},{'C4','CP2','CP6','P4'});
%     idx_ROI_rgt = ismember({chanlocs.labels},{'C3','CP1','CP5','P3'});
%     idx_lft = contains(dat_wvlt.events(:,1),'lh')';
%     idx_rgt = contains(dat_wvlt.events(:,1),'rh')';
% 
%     idx_err = find(~idx_cor & ~idx_art);
% 
%     for t = 1:length(idx_err)
%         
%         idx_time_ERD = ep_time >= 300 & ep_time <= dat_wvlt.SO_ms(idx_err(t));
%         idx_time_ERS = ep_time >= dat_wvlt.SO_ms(idx_err(t)) & ep_time <= dat_wvlt.SO_ms(idx_err(t))+1000
% 
%         if idx_lft(idx_err(t)) == 1
%             dat_mu_ERD_err(s,c) = dat_wvlt.powspctrm(idx_err(t),idx_ROI_lft,idx_alpha,idx_time_ERD) % -> 96 tr x 24 chan x 35 freq x 2000 time
%         elseif idx_rgt(idx_err(t)) == 1
%             dat_mu_ERD_err(s,c) = dat_wvlt.powspctrm(idx_err(t),idx_ROI_rgt,idx_alpha,idx_time_ERD) % -> 96 tr x 24 chan x 35 freq x 2000 time
% 
%     var_RT_err = dat_wvlt.SO_ms(~idx_cor & ~idx_art);
%     var_RT_cor = dat_wvlt.SO_ms(idx_cor & ~idx_art);
%     var_RT_cor = unique(var_RT_cor(dsearchn(dat_wvlt.SO_ms(idx_cor & ~idx_art)',var_RT_err')));
%     
%     idx_cor = ismember(dat_wvlt.SO_ms,var_RT_cor);
%     idx_err = ismember(dat_wvlt.SO_ms,var_RT_err);
%     
%     
%     for t = 1:length(idx_art)
%         if idx_art(t) == 1; continue;end;
%         
%         idx_time_ERD = ep_time >= 300 & ep_time <= dat_wvlt.SO_ms(t);
%         idx_time_ERS = ep_time >= dat_wvlt.SO_ms(t) & ep_time <= dat_wvlt.SO_ms(t)+1000
%         if idx_cor(t) == 1
%             dat_mu_ERD_cor(s,c) = dat_wvlt.powspctrm(t,idx_ROI_lft,idx_alpha,idx_time_ERD) % -> 96 tr x 24 chan x 35 freq x 2000 time
%         elseif idx_err(t) == 1
%             dat_mu_ERD_err(s,:)
%         
%     end

end


save([PATHOUT_dat_wvlt 'ERD_latmed.mat'],'dat_lat','dat_med');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%               Mu ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([PATHOUT_dat_wvlt 'ERD_latmed.mat']);

% time course
idx_ROI_mu = ismember({chanlocs.labels},{'P4','P3'})
idx_ROI_alpha = ismember({chanlocs.labels},{'O1','O2'})
map_limits      = [-8 8];

close all
figure

%%%%%%%%%% YOUNG %%%%%%%%%%%%%%%%%%%%%%%
subplot(3,4,1)
plot(ep_time,squeeze(mean(mean(dat_med(idx_young,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_young,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013 % first medial than lateral in legend
    title 'mu_{young}'

subplot(3,4,2)
plot(ep_time,squeeze(mean(mean(dat_med(idx_young,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_young,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013
    title '\alpha_{young}'
    
% topoplots
subplot(3,4,3)
topoplot(squeeze(mean(mean(dat_med(idx_young,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'medial_{young}'
    set(findobj(gca,'type','patch'),'facecolor',c_med)
subplot(3,4,4)
topoplot(squeeze(mean(mean(dat_lat(idx_young,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'lateral_{young}'
    set(findobj(gca,'type','patch'),'facecolor',c_lat)


        
%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%
subplot(3,4,5)
plot(ep_time,squeeze(mean(mean(dat_med(idx_old,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_old,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013 % first medial than lateral in legend
    title 'mu_{old}'

subplot(3,4,6)
plot(ep_time,squeeze(mean(mean(dat_med(idx_old,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_old,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013
    title '\alpha_{old}'
 
% topoplots
subplot(3,4,7)
topoplot(squeeze(mean(mean(dat_med(idx_old,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'medial_{old}'
    set(findobj(gca,'type','patch'),'facecolor',c_med)
subplot(3,4,8)
topoplot(squeeze(mean(mean(dat_lat(idx_old,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'lateral_{old}'
    set(findobj(gca,'type','patch'),'facecolor',c_lat)


%%%%%%%%%% STROKE %%%%%%%%%%%%%%%%%%%%%%%    
subplot(3,4,9)
plot(ep_time,squeeze(mean(mean(dat_med(idx_stroke,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_stroke,idx_ROI_mu,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013 % first medial than lateral in legend
    title 'mu_{stroke}'

subplot(3,4,10)
plot(ep_time,squeeze(mean(mean(dat_med(idx_stroke,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time,squeeze(mean(mean(dat_lat(idx_stroke,idx_ROI_alpha,:),2))),'LineWidth',2,'Color',c_lat)
    legend_ERD_terhorst2013
    title '\alpha_{stroke}'

%topoplots
subplot(3,4,11)
topoplot(squeeze(mean(mean(dat_med(idx_stroke,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'medial_{stroke}'
    set(findobj(gca,'type','patch'),'facecolor',c_med)
subplot(3,4,12)
topoplot(squeeze(mean(mean(dat_lat(idx_stroke,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{find(idx_ROI_mu | idx_ROI_alpha),'.','w',15})
    title 'lateral_{stroke}'
    set(findobj(gca,'type','patch'),'facecolor',c_lat)


save_fig(gcf,PATHOUT_plots,'avgERD_timecourse_latmed','FigSize',[0 0 45 30]);


%% statistic

dat_ERD_stat = [squeeze(mean(mean(dat_lat(:,idx_ROI_mu,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_med(:,idx_ROI_mu,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_lat(:,idx_ROI_alpha,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_med(:,idx_ROI_alpha,idx_time_erd),2),3))];
            
nms_group(idx_stroke)   = string('stroke');
nms_group(idx_old)      = string('old');
nms_group(idx_young)    = string('young');
nms_group = [nms_group nms_group nms_group nms_group];

nms_lat(1:size(dat_lat,1))   = string('lat');
nms_med(1:size(dat_med,1))   = string('med');
nms_con = [nms_lat nms_med nms_lat nms_med];

nms_mu(1:size(dat_lat,1))   = string('mu');
nms_alpha(1:size(dat_lat,1))   = string('alpha');
nms_loc = [nms_mu nms_mu nms_alpha nms_alpha];

[~,tbl,stats] = anovan(dat_ERD_stat,{nms_group,nms_con,nms_loc})
multcompare(stats,'Dimension',[1 2 3])

%% Prep for JASP
dat_JASP = [squeeze(mean(mean(dat_lat(:,idx_ROI_mu,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_med(:,idx_ROI_mu,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_lat(:,idx_ROI_alpha,idx_time_erd),2),3));...
                squeeze(mean(mean(dat_med(:,idx_ROI_alpha,idx_time_erd),2),3))];
            
nms_group(idx_stroke)   = string('stroke');
nms_group(idx_old)      = string('old');
nms_group(idx_young)    = string('young');
nms_group = [nms_group nms_group nms_group nms_group];

nms_lat(1:size(dat_lat,1))   = string('lat');
nms_med(1:size(dat_med,1))   = string('med');
nms_con = [nms_lat nms_med nms_lat nms_med];

nms_mu(1:size(dat_lat,1))   = string('mu');
nms_alpha(1:size(dat_lat,1))   = string('alpha');
nms_loc = [nms_mu nms_mu nms_alpha nms_alpha];

nms_sub = extractBetween(nms_SUBJ,'wvlts_','.mat');
nms_sub = repmat(nms_sub,1,4)

dat_JASP = table(dat_ERD_stat,nms_group',nms_con',nms_loc',nms_sub','VariableNames',{'ERD','Group','Condition','ROI','ID'});
dat_JASP_ = unstack(dat_JASP,{'ERD'},'Condition');
writetable(dat_JASP_,[PATHOUT_dat_wvlt 'dat_JASP.csv'])

%% Plots results for statistic
c_stroke    = [149,165,166]/255;
c_old       = c_lh;
c_young     = c_rh;

dat_m_lat_y = mean(mean(dat_lat(idx_young,idx_chan_mu,idx_time_erd),3),2);
dat_m_med_y = mean(mean(dat_med(idx_young,idx_chan_mu,idx_time_erd),3),2);

dat_m_lat_o = mean(mean(dat_lat(idx_old,idx_chan_mu,idx_time_erd),3),2);
dat_m_med_o = mean(mean(dat_med(idx_old,idx_chan_mu,idx_time_erd),3),2);


%%

close all
figure

subplot(1,2,1)
singleBoxplot([dat_m_lat_y,dat_m_med_y])
tune_BP([c_lat;c_med])
    ylabel 'Power [dB]'
    ylim ([-8 0])
    xticks ([1 2])
    xticklabels ({'Lateral','Medial'})
    [h p] = ttest2(dat_m_lat_y,dat_m_med_y);
    sigstar({{'Lateral','Medial'}},p)
    title (['young // p_{lateral-medial} = ' num2str(round(p,3))])

subplot(1,2,2)
singleBoxplot([dat_m_lat_o,dat_m_med_o])
    tune_BP([c_lat;c_med])
    ylabel 'Power [dB]'
    ylim ([-8 0])
    xticks ([1 2])
    xticklabels ({'Lateral','Medial'})
    [h p] = ttest2(dat_m_lat_o,dat_m_med_o);
    sigstar({{'Lateral','Medial'}},p)
    title (['old // p_{lateral-medial} = ' num2str(round(p,3))])
    
save_fig(gcf,PATHOUT_plots,'BP_timecourse_mu','FigSize',[0 0 20 20]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%               Alpha ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% time course
idx_chan_alpha     = [20 21 22 23 24];
map_limits      = [-8 8];

close all
figure

subplot(2,4,[1:3])
plot(ep_time, squeeze(mean(mean(dat_med(idx_young,idx_chan_alpha,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time, squeeze(mean(mean(dat_lat(idx_young,idx_chan_alpha,:),2))),'LineWidth',2,'Color',c_lat)
    legend ({'medial','lateral'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'Power [dB]'
    xlim ([-500 2000])
    ylim ([-6 3])
    vline(0,'k');
    vline([300 1500],'--k');
    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'young'
        

subplot(2,4,[5:7])
plot(ep_time, squeeze(mean(mean(dat_med(idx_old,idx_chan_alpha,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time, squeeze(mean(mean(dat_lat(idx_old,idx_chan_alpha,:),2))),'LineWidth',2,'Color',c_lat)
    legend ({'medial','lateral'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'Power [dB]'
    xlim ([-500 2000])
    ylim ([-6 3])
    vline(0,'k');
    vline([300 1500],'--k');
    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'old'



% topoplots
subplot(2,4,[4])
topoplot(squeeze(mean(mean(dat_med(idx_young,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{idx_chan_alpha,'.','w',15})
subplot(2,4,[8])
topoplot(squeeze(mean(mean(dat_lat(idx_old,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{idx_chan_alpha,'.','w',15})


save_fig(gcf,PATHOUT_plots,'avgERD_timecourse_alpha','FigSize',[0 0 30 20]);

%% Plots results for statistic
c_stroke    = [149,165,166]/255;
c_old       = c_lh;
c_young     = c_rh;

dat_m_lat_y = mean(mean(dat_lat(idx_young,idx_chan_alpha,idx_time_erd),3),2);
dat_m_med_y = mean(mean(dat_med(idx_young,idx_chan_alpha,idx_time_erd),3),2);

dat_m_lat_o = mean(mean(dat_lat(idx_old,idx_chan_alpha,idx_time_erd),3),2);
dat_m_med_o = mean(mean(dat_med(idx_old,idx_chan_alpha,idx_time_erd),3),2);


%%

close all
figure

subplot(1,2,1)
singleBoxplot([dat_m_lat_y,dat_m_med_y])
tune_BP([c_lat;c_med])
    ylabel 'Power [dB]'
    ylim ([-8 0])
    xticks ([1 2])
    xticklabels ({'Lateral','Medial'})
    [h p] = ttest2(dat_m_lat_y,dat_m_med_y);
    sigstar({{'Lateral','Medial'}},p)
    title (['young // p_{lateral-medial} = ' num2str(round(p,3))])

subplot(1,2,2)
singleBoxplot([dat_m_lat_o,dat_m_med_o])
    tune_BP([c_lat;c_med])
    ylabel 'Power [dB]'
    ylim ([-8 0])
    xticks ([1 2])
    xticklabels ({'Lateral','Medial'})
    [h p] = ttest2(dat_m_lat_o,dat_m_med_o);
    sigstar({{'Lateral','Medial'}},p)
    title (['old // p_{lateral-medial} = ' num2str(round(p,3))])
    
save_fig(gcf,PATHOUT_plots,'BP_timecourse_alpha','FigSize',[0 0 20 20]);




%%

























