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


%% CONNY ERD STATS
% Könntest du für die TF-Maps Differenz-Plots machen: medial-lateral jeweils 
% für jede Gruppe und Alt-Jung bzw. Alt-Stroke jeweils für medial und lateral getrennt? 
% Mann könnte die einzelnen Elektroden-Positionen auch gleich noch basierend auf 
% p-Werten markieren (dafür müsste man jeweils einen T-Test für unabhängige oder 
% abhängige Stichproben pro Kanal rechnen und dann könnte man zB jeden Kanal ja 
% nach p-Wert farblich markieren oder in der Größe variiren. Man kann p-Werte natürlich auch mappen). 

chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
load([PATHOUT_YAN 'ERD_latmed.mat']); %subj x chan x  t_period [56 x 24 x 4 (BL,BEG,MID,END)]

dat_lat_alpha(dat_lat_alpha == 0) = NaN;
dat_med_alpha(dat_med_alpha == 0) = NaN;
%%
close all
figure

%%%%%%%%%%%%%%%%%%%%%% Lateral %%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
[h p ci tstat] = ttest2(dat_lat_alpha(idx_old,:,3),dat_lat_alpha(idx_young,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'old vs. young (lateral)'

subplot(2,3,2)
[h p ci tstat] = ttest2(dat_lat_alpha(idx_old,:,3),dat_lat_alpha(idx_stroke,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'old vs. stroke (lateral)'

subplot(2,3,3)
[h p ci tstat] = ttest2(dat_lat_alpha(idx_stroke,:,3),dat_lat_alpha(idx_young,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'stroke vs. young (lateral)'

%%%%%%%%%%%%%%%%%%%%%% medial %%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4)
[h p ci tstat] = ttest2(dat_med_alpha(idx_old,:,3),dat_med_alpha(idx_young,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'old vs. young (medial)'

subplot(2,3,5)
[h p ci tstat] = ttest2(dat_med_alpha(idx_old,:,3),dat_med_alpha(idx_stroke,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'old vs. stroke (medial)'

subplot(2,3,6)
[h p ci tstat] = ttest2(dat_med_alpha(idx_stroke,:,3),dat_med_alpha(idx_young,:,3));
topoplot(tstat.tstat,chanlocs)
put_CB('t-value')
title 'stroke vs. young (medial)'

sgtitle 't-values ERD [8-12 Hz] at 300-800 ms'
save_fig(gcf,PATHOUT_plots,'ERD_t_vals','FigSize',[0 0 30 20]);


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






























