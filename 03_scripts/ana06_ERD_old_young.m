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
PATHIN_WVLT = [MAIN '02_data\05_wvlts\'];

PATHOUT_ERD = [MAIN '02_data\05_ERD\'];
PATHOUT_plots = [MAIN '04_plots\05_ERD\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_ERD)
    mkdir(PATHOUT_ERD);
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


%% Hypothesis 1

% ...   [ We expected to replicate the finding of a stronger central mu ERD for
%       medial stimuli during mental rotation as compared to lateral
%       stimuli in younger participants and extend it to older participants] ...

load([PATHIN_WVLT 'cfg.mat'])
% assign necessary parameters for following analysis in EEGLAB *yay*
cfg_el.ep_length   = [-1 3];
cfg_el.ep_prune    = 3; % pruning parameters for rejection
cfg_el.freqs       = [1:40]; % freqs OI
cfg_el.srate       = 500; 
cfg_el.times       = cfg_el.ep_length(1)*1000:1000/cfg_el.srate:cfg_el.ep_length(2)*1000-1; % time vec for epoch
cfg_el.cycles      = [5 8]; %cycle range from 5 to 8, Debener ea., 2005
cfg_el.bl          = [-500 -200];
cfg_el.bl_erp      = [-600 -400];

% get time epcoh data and indices for time period
ep_time         = cfg_el.times;
idx_bl          = ep_time >= -500 & ep_time <= -200;
idx_time_erd    = ep_time >= 300 & ep_time <= 1500;

% indeices for freq band
idx_alpha   = cfg_el.freqs >= 8 & cfg_el.freqs <= 12;
idx_beta    = cfg_el.freqs >= 15 & cfg_el.freqs <= 25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                CONTINUE WITH EXTRACTING TIMED ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load RTs 
load([PATHIN_RTs 'RT_ALL.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);

dat_lat = [];
dat_med = [];

for s = 48%1:numel(SUBJ)
    
    % laod wvlt data
    load([PATHIN_WVLT nms_SUBJ{s}]);
    display (['Working on Subject ' num2str(SUBJ(s))])
    
    % indeices for freq band
    idx_alpha   = dat_wvlt.freq >= 8 & dat_wvlt.freq <= 12;

    %normalize trial data 
    tmp_dat     = dat_wvlt.powspctrm;  % -> size powspctrm [96 x 24 x 35 x 2000]
    tmp_bl      = squeeze(nanmean(tmp_dat(:,:,:,idx_bl),4));
    tmp_dat_dB  = 10*log10(tmp_dat./tmp_bl); %transform to dB // dB = 10*log10 (signal/baseline)
    plot(ep_time,squeeze(tmp_dat(:,13,8,:)),'LineWidth',0.3)
    hold on
    plot(ep_time,squeeze(mean(tmp_dat_dB(:,13,8,:))),'-k','LineWidth',2)
    vline([-500 -200],'--k');

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
     
    %% complexity
    % create 2D matrix -> sub x chan

    % for ERD analysis
    dat_lat(s,:,:)   = squeeze(mean(mean(tmp_dat_dB(idx_lat & idx_cor',:,idx_alpha,:),3))); % -> 96 tr x 24 chan x 35 freq x 2000 time
    dat_med(s,:,:)   = squeeze(mean(mean(tmp_dat_dB(idx_med & idx_cor',:,idx_alpha,:),3))); 
    

end


%% save data

idx_bad_sub = nanmean(nanmean(dat_med(:,idx_chan_mu,:),3),2) == 0;
SUBJ(idx_bad_sub) = [];
% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;

dat_lat(idx_bad_sub,:,:) = [];
dat_med(idx_bad_sub,:,:) = [];


save([PATHOUT_ERD 'ERD_latmed.mat'],'dat_lat','dat_med');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%               Mu ERD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% time course
idx_chan_mu     = [9 16 14 21 10 15 17 22];
map_limits      = [-8 8];

close all
figure

subplot(2,4,[1:3])
plot(ep_time, squeeze(mean(mean(dat_med(idx_young,idx_chan_mu,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time, squeeze(mean(mean(dat_lat(idx_young,idx_chan_mu,:),2))),'LineWidth',2,'Color',c_lat)
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
plot(ep_time, squeeze(mean(mean(dat_med(idx_old,idx_chan_mu,:),2))),'LineWidth',2,'Color',c_med)
hold on
plot(ep_time, squeeze(mean(mean(dat_lat(idx_old,idx_chan_mu,:),2))),'LineWidth',2,'Color',c_lat)
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
    'maplimits',[-6 3],'emarker2',{idx_chan_mu,'.','w',15})
subplot(2,4,[8])
topoplot(squeeze(mean(mean(dat_lat(idx_old,:,idx_time_erd),3))),chanlocs,...
    'maplimits',[-6 3],'emarker2',{idx_chan_mu,'.','w',15})


save_fig(gcf,PATHOUT_plots,'avgERD_timecourse_mu','FigSize',[0 0 30 20]);

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
































