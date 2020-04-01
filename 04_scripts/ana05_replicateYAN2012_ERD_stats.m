%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Stats of ERD 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERD statistical analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_RTs  = [MAIN '02_data\03_RTs\'];
PATHIN_WVLT = [MAIN '02_data\ana04_wavelets\'];
PATHOUT_WVLT_plots = [MAIN '03_plots\ana04_wavelets\TF_plots\'];

PATHOUT_YAN = [MAIN '02_data\ana05_yan2012\'];
PATHOUT_plots = [MAIN '03_plots\ana05_yan2012\'];

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

% timing information
load([PATHIN_WVLT 'cfg.mat'])

% get time epcoh data and indices for time period
ep_time     = cfg_wvlt.times;
idx_bl      = ep_time >= -500 & ep_time <= -100;
idx_beg     = ep_time >= 0 & ep_time <= 300;
idx_mdl     = ep_time >= 300 & ep_time <= 800;
idx_end     = ep_time >= 800 & ep_time <= 2500;

nms_time_period = {'BL [-500:-200ms]','BEG [0:300ms]','MDL [300:800ms]','END [800:2500ms]'};
map_limits = [-8 8];


%% CONNY ERD STATS
% K�nntest du f�r die TF-Maps Differenz-Plots machen: medial-lateral jeweils 
% f�r jede Gruppe und Alt-Jung bzw. Alt-Stroke jeweils f�r medial und lateral getrennt? 
% Mann k�nnte die einzelnen Elektroden-Positionen auch gleich noch basierend auf 
% p-Werten markieren (daf�r m�sste man jeweils einen T-Test f�r unabh�ngige oder 
% abh�ngige Stichproben pro Kanal rechnen und dann k�nnte man zB jeden Kanal ja 
% nach p-Wert farblich markieren oder in der Gr��e variiren. Man kann p-Werte nat�rlich auch mappen). 

chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
load([PATHOUT_YAN 'ERD_latmed.mat']); %subj x chan x  t_period [56 x 24 x 4 (BL,BEG,MID,END)]

dat_lat_alpha(dat_lat_alpha == 0) = NaN;
dat_med_alpha(dat_med_alpha == 0) = NaN;
%%

%---------------------- Per group Med vs. Lat ----------------------

%Alpha
vals_map_lim = [-1 1];

close all
figure

for t = 1:numel(nms_time_period)
    
    %prep data
    dat_diff_old = nanmean(dat_lat_alpha(idx_old,:,t)-dat_med_alpha(idx_old,:,t));
    [h_old p_old] = ttest(dat_lat_alpha(idx_old,:,t),dat_med_alpha(idx_old,:,t));
    dat_diff_young = nanmean(dat_lat_alpha(idx_young,:,t)-dat_med_alpha(idx_young,:,t));
    [h_young p_young] = ttest(dat_lat_alpha(idx_young,:,t),dat_med_alpha(idx_young,:,t));
    dat_diff_stroke = nanmean(dat_lat_alpha(idx_stroke,:,t)-dat_med_alpha(idx_stroke,:,t));
    [h_stroke p_stroke] = ttest(dat_lat_alpha(idx_stroke,:,t),dat_med_alpha(idx_stroke,:,t));
    
        
    subplot(3,4,0+t)
    topoplot(dat_diff_old,chanlocs,'electrodes','off','emarker2',{find(h_old),'.','k'},'maplimits',vals_map_lim)
        title (['O: ' nms_time_period{t}])    
    subplot(3,4,4+t)
    topoplot(dat_diff_young,chanlocs,'electrodes','off','emarker2',{find(h_young),'.','k'},'maplimits',vals_map_lim)
        title (['Y: ' nms_time_period{t}])    
    subplot(3,4,8+t)
    topoplot(dat_diff_stroke,chanlocs,'electrodes','off','emarker2',{find(h_stroke),'.','k'},'maplimits',vals_map_lim)
        title (['S: ' nms_time_period{t}])    
    
end

sgtitle 'ERD [8-12 Hz] differences medial-lateral'
save_fig(gcf,PATHOUT_WVLT_plots,'ERD_dif_grp_alpha','FigSize',[0 0 30 25]);

% beta
vals_map_lim = [-1 1];

close all
figure

for t = 1:numel(nms_time_period)
    
    %prep data
    dat_diff_old = nanmean(dat_lat_beta(idx_old,:,t)-dat_med_beta(idx_old,:,t));
    [h_old p_old] = ttest(dat_lat_beta(idx_old,:,t),dat_med_beta(idx_old,:,t));
    dat_diff_young = nanmean(dat_lat_beta(idx_young,:,t)-dat_med_beta(idx_young,:,t));
    [h_young p_young] = ttest(dat_lat_beta(idx_young,:,t),dat_med_beta(idx_young,:,t));
    dat_diff_stroke = nanmean(dat_lat_beta(idx_stroke,:,t)-dat_med_beta(idx_stroke,:,t));
    [h_stroke p_stroke] = ttest(dat_lat_beta(idx_stroke,:,t),dat_med_beta(idx_stroke,:,t));
    
        
    subplot(3,4,0+t)
    topoplot(dat_diff_old,chanlocs,'electrodes','off','emarker2',{find(h_old),'.','k'},'maplimits',vals_map_lim)
        title (['O: ' nms_time_period{t}])    
    subplot(3,4,4+t)
    topoplot(dat_diff_young,chanlocs,'electrodes','off','emarker2',{find(h_young),'.','k'},'maplimits',vals_map_lim)
        title (['Y: ' nms_time_period{t}])    
    subplot(3,4,8+t)
    topoplot(dat_diff_stroke,chanlocs,'electrodes','off','emarker2',{find(h_stroke),'.','k'},'maplimits',vals_map_lim)
        title (['S: ' nms_time_period{t}])    
    
end

sgtitle 'ERD [15-25 Hz] differences medial-lateral'
save_fig(gcf,PATHOUT_WVLT_plots,'ERD_dif_grp_beta','FigSize',[0 0 30 25]);


%% ---------------------- Per Con group compare ----------------------

%Alpha
vals_map_lim = [-3 3];

close all
figure

for t = 1:numel(nms_time_period)
    
    %prep data
    dat_diff_lat_oy = nanmean(dat_lat_alpha(idx_old,:,t))-nanmean(dat_lat_alpha(idx_young,:,t));
    [h_lat_oy p_lat_oy] = ttest2(dat_lat_alpha(idx_old,:,t),dat_lat_alpha(idx_young,:,t));
    dat_diff_med_oy = nanmean(dat_med_alpha(idx_old,:,t))-nanmean(dat_med_alpha(idx_young,:,t));
    [h_med_oy p_med_oy] = ttest2(dat_med_alpha(idx_old,:,t),dat_med_alpha(idx_young,:,t));

    dat_diff_lat_os = nanmean(dat_lat_alpha(idx_old,:,t))-nanmean(dat_lat_alpha(idx_stroke,:,t));
    [h_lat_os p_lat_os] = ttest2(dat_lat_alpha(idx_old,:,t),dat_lat_alpha(idx_stroke,:,t));
    dat_diff_med_os = nanmean(dat_med_alpha(idx_old,:,t))-nanmean(dat_med_alpha(idx_stroke,:,t));
    [h_med_os p_med_os] = ttest2(dat_med_alpha(idx_old,:,t),dat_med_alpha(idx_stroke,:,t));

        
    subplot(4,4,0+t)
    topoplot(dat_diff_lat_oy,chanlocs,'electrodes','off','emarker2',{find(h_lat_oy),'.','k'},'maplimits',vals_map_lim)
        title (['LAT_{o-y}: ' nms_time_period{t}])    
    subplot(4,4,4+t)
    topoplot(dat_diff_lat_os,chanlocs,'electrodes','off','emarker2',{find(h_lat_os),'.','k'},'maplimits',vals_map_lim)
        title (['LAT_{o-s}: ' nms_time_period{t}])    
    subplot(4,4,8+t)
    topoplot(dat_diff_med_oy,chanlocs,'electrodes','off','emarker2',{find(h_med_oy),'.','k'},'maplimits',vals_map_lim)
        title (['MED_{o-y}: ' nms_time_period{t}])  
    subplot(4,4,12+t)
    topoplot(dat_diff_med_os,chanlocs,'electrodes','off','emarker2',{find(h_med_os),'.','k'},'maplimits',vals_map_lim)
        title (['MED_{o-s}: ' nms_time_period{t}])    

    
end

sgtitle 'ERD [8-12 Hz] differences group'
save_fig(gcf,PATHOUT_WVLT_plots,'ERD_dif_con_alpha','FigSize',[0 0 30 35]);

%beta

close all
figure

for t = 1:numel(nms_time_period)
    
    %prep data
    dat_diff_lat_oy = nanmean(dat_lat_beta(idx_old,:,t))-nanmean(dat_lat_beta(idx_young,:,t));
    [h_lat_oy p_lat_oy] = ttest2(dat_lat_beta(idx_old,:,t),dat_lat_beta(idx_young,:,t));
    dat_diff_med_oy = nanmean(dat_med_beta(idx_old,:,t))-nanmean(dat_med_beta(idx_young,:,t));
    [h_med_oy p_med_oy] = ttest2(dat_med_beta(idx_old,:,t),dat_med_beta(idx_young,:,t));

    dat_diff_lat_os = nanmean(dat_lat_beta(idx_old,:,t))-nanmean(dat_lat_beta(idx_stroke,:,t));
    [h_lat_os p_lat_os] = ttest2(dat_lat_beta(idx_old,:,t),dat_lat_beta(idx_stroke,:,t));
    dat_diff_med_os = nanmean(dat_med_beta(idx_old,:,t))-nanmean(dat_med_beta(idx_stroke,:,t));
    [h_med_os p_med_os] = ttest2(dat_med_beta(idx_old,:,t),dat_med_beta(idx_stroke,:,t));

        
    subplot(4,4,0+t)
    topoplot(dat_diff_lat_oy,chanlocs,'electrodes','off','emarker2',{find(h_lat_oy),'.','k'},'maplimits',vals_map_lim)
        title (['LAT_{o-y}: ' nms_time_period{t}])    
    subplot(4,4,4+t)
    topoplot(dat_diff_lat_os,chanlocs,'electrodes','off','emarker2',{find(h_lat_os),'.','k'},'maplimits',vals_map_lim)
        title (['LAT_{o-s}: ' nms_time_period{t}])    
    subplot(4,4,8+t)
    topoplot(dat_diff_med_oy,chanlocs,'electrodes','off','emarker2',{find(h_med_oy),'.','k'},'maplimits',vals_map_lim)
        title (['MED_{o-y}: ' nms_time_period{t}])  
    subplot(4,4,12+t)
    topoplot(dat_diff_med_os,chanlocs,'electrodes','off','emarker2',{find(h_med_os),'.','k'},'maplimits',vals_map_lim)
        title (['MED_{o-s}: ' nms_time_period{t}])    

    
end

sgtitle 'ERD [15-25 Hz] differences group'
save_fig(gcf,PATHOUT_WVLT_plots,'ERD_dif_con_beta','FigSize',[0 0 30 35]);














%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   CSPs for med vs. lat
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([PATHIN_WVLT 'cfg.mat'])

% get time epcoh data and indices for time period
ep_time         = cfg_el.times;
idx_bl          = ep_time >= -500 & ep_time <= -200;
idx_time_erd    = ep_time >= 300 & ep_time <= 1500;

% indeices for freq band
idx_alpha   = cfg_el.freqs >= 8 & cfg_el.freqs <= 12;
idx_beta    = cfg_el.freqs >= 15 & cfg_el.freqs <= 25;

% load RTs 
load([PATHIN_RTs 'RT_ALL.mat']);

for s = 48%1:numel(SUBJ)
    
    % laod wvlt data
    display (['Working on Subject ' num2str(SUBJ(s))])
    load([PATHIN_WVLT nms_SUBJ{s}]);
    
    % indeices for freq band
    idx_alpha   = dat_wvlt.freq >= 8 & dat_wvlt.freq <= 12;

    %normalize trial data 
    tmp_dat     = dat_wvlt.powspctrm;  % -> size powspctrm [96 x 24 x 35 x 2000]
    tmp_bl      = squeeze(nanmean(tmp_dat(:,:,:,idx_bl),4));
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
     
    %% CSPs
    
    % for CSP analysis
    dat_lat  = squeeze(mean(tmp_dat_dB(idx_lat & idx_cor',:,idx_alpha,:),3)); % -> 96 tr x 24 chan x 35 freq x 2000 time
    dat_med  = squeeze(mean(tmp_dat_dB(idx_med & idx_cor',:,idx_alpha,:),3)); 
 
    covar1 = cov(X1');   % S1~[C x C]
    covar2 = cov_shrink(X2');   % S2~[C x C]
    covar1(~isfinite(covar1)) = 0; covar2(~isfinite(covar2)) = 0;

    [V,tmp] = eig(covar1,covar1+covar2);          % V = eigenvectoren; D = eigenwerte

    tmpFilters = V(:,[1:nrPat end-nrPat+1:end]);  % Select the last n patterns and the first n patterns
    P = inv(V);
    tmpPatterns = P([1:nrPat end-nrPat+1:end],:);

    % If necessary, remove imaginary part
    tmpFilters = real(tmpFilters); tmpPatterns = real(tmpPatterns);
        
        
        
    % create 2D matrix -> sub x chan

    % for ERD analysis
    dat_lat(s,:,:)   = squeeze(mean(mean(tmp_dat_dB(idx_lat & idx_cor',:,idx_alpha,:),3))); % -> 96 tr x 24 chan x 35 freq x 2000 time
    dat_med(s,:,:)   = squeeze(mean(mean(tmp_dat_dB(idx_med & idx_cor',:,idx_alpha,:),3))); 
    

end






























