% Get RTs of LLT data
% Use speech onset function to evaluate RTs of every trial

PATHIN_eeg = [MAIN '02_data\02_cleanEEG\'];

PATHOUT_WAVELETS = [MAIN '02_data\04_wavelets\'];
PATHOUT_plots = [PATHOUT_WAVELETS 'TF_plots\'];


% Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_WAVELETS)
    mkdir(PATHOUT_plots)
end

% assign necessary parameters for following analysis in EEGLAB *yay*
cfg_el.ep_length   = [-1 3];
cfg_el.ep_prune    = 3; % pruning parameters for rejection
cfg_el.freqs       = [1:40]; % freqs OI
cfg_el.srate       = 500; 
cfg_el.times       = cfg_el.ep_length(1)*1000:1000/cfg_el.srate:cfg_el.ep_length(2)*1000-1; % time vec for epoch
cfg_el.cycles      = [5 8]; %cycle range from 5 to 8, Debener ea., 2005
cfg_el.bl          = [-500 -200];
cfg_el.bl_erp      = [-600 -400];

save([PATHOUT_WAVELETS 'cfg.mat'],'cfg_el');

%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_eeg '*finICA_clean.set']));
SUBJ = extractBefore({list.name},'_');

%%
doc_error = {};

for sub = 1:length(SUBJ)
    
    try
    % check if dataset already exists
    if exist([PATHOUT_WAVELETS 'wvlts_' SUBJ{sub} '.mat'],'file') == 2 % update for current evaluation
        disp(['WAVELETS for ' SUBJ{sub} ' have already been calculated. Continue with next dataset.']);
        continue;
    end
        
    %% load ICA cleaned sets
        
    EEG     = pop_loadset('filename',[SUBJ{sub} '_finICA_clean.set'],'filepath',PATHIN_eeg);
    EEG     = eeg_checkset(EEG, 'eventconsistency' );
    EEG.ID  = SUBJ{sub};
    

    % epoching 
    EEG             = eeg_checkset( EEG );
    EEG             = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG             = pop_reref( EEG, []); % reref CAR
    nms_trig        = unique({EEG.event.type});
    idx_trig_OI     = contains(nms_trig,'E_');
    EEG             = pop_epoch( EEG, nms_trig(idx_trig_OI),cfg_el.ep_length, 'epochinfo', 'yes');
    EEG             = pop_rmbase( EEG, cfg_el.bl_erp ,[]);
    
    % Remove remaining non-stereotyped artifacts 
    prune = false;
    if prune
        EEG = pop_jointprob(EEG,1, 1:EEG.nbchan ,cfg_el.ep_prune,cfg_el.ep_prune,0,0);
        EEG = pop_rejkurt(EEG,1, 1:EEG.nbchan ,cfg_el.ep_prune,cfg_el.ep_prune,0,0);
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch(EEG,find(EEG.reject.rejglobal),0);
    end
    EEG = eeg_checkset( EEG );

    
    %% Convert to fieldtrip struct for TF anaylsis
    
    dat_wvlt = eeglab2fieldtrip(EEG, 'preprocessing'); 
    
    % MRLT-WVLT analysis according to Debener ea., 2005
    cfg             = [];
    cfg.channel     = 'EEG';
    cfg.method      = 'wavelet';
    cfg.output      = 'pow';
    cfg.keeptrials  = 'yes';
    cfg.foi         = 1:1:35; %frequency range
    cfg.toi         = dat_wvlt.time{1};
    cfg.width       = linspace(5,8,numel(cfg.foi));
    cfg.pad         ='nextpow2'; %quicker computatio time
    dat_wvlt        = ft_freqanalysis(cfg, dat_wvlt);

    cfg = [];
    cfg.baseline     = [-0.5 -0.2];
    cfg.baselinetype = 'db';
    cfg.showlabels   = 'yes';
    figure
    ft_multiplotTFR(cfg, dat_wvlt)    
    save_fig(gcf,PATHOUT_plots,[SUBJ{sub} '_wvlts']);
    
    % add additional information
    % stimuli
    unique_ep = unique([EEG.event.epoch]);
    
    for st = 1:length(unique_ep)
        idx_ep_stim = contains({EEG.event.type},'E_') & [EEG.event.epoch] == unique_ep(st);
        stim_pic = EEG.event(find(idx_ep_stim,1)).type; % only use first important stimuli
        stim_split = strsplit(stim_pic,{'_','.'});
        stim(st,:) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=depth rotation, 4 = lateral rotation
    end
    dat_wvlt.id = SUBJ{sub};
    dat_wvlt.events = stim;
    clear stim
    
    save([PATHOUT_WAVELETS 'wvlts_' SUBJ{sub} '.mat'],'dat_wvlt', '-v7.3');
    
    %% store for all participants
    
    idx_chan = find(strcmp(dat_wvlt.label,'CPz')); 
    
    dat_wvlt_all.id(sub)        = SUBJ(sub);
    dat_wvlt_all.data(sub,:,:)  = squeeze(mean(dat_wvlt.powspctrm(:,idx_chan,:,:),1));
    dat_wvlt_all.times          = dat_wvlt.time;
    dat_wvlt_all.freqs          = dat_wvlt.freq;

    catch
        doc_error{end+1} = SUBJ(sub);
    end %try subj
    
end

save([PATHOUT_WAVELETS 'wvlts_all.mat'],'dat_wvlt_all');
save([PATHOUT_WAVELETS 'error_all.mat'],'doc_error');

%% Plot single subject CPz results
if ~exist('dat_wvlt_all');load([PATHOUT_WAVELETS 'wvlts_all.mat']);end;

psp     = numSubplots(size(dat_wvlt_all.data,1)+1);
wvlt_bl = [-0.5 -0.2];


figure
for s = 1:size(dat_wvlt_all.data,1)
    
    subplot(psp(1),psp(2),s)

    %normalize trial data 
    tmp_dat     = squeeze(dat_wvlt_all.data(s,:,:)); %trials x 1channel x freqs x time
    idx_bl      = dsearchn(dat_wvlt_all.times',wvlt_bl');
    tmp_bl      = squeeze(mean(tmp_dat(:,idx_bl(1):idx_bl(2)),2));
    tmp_dat_dB  = real(10*log10(tmp_dat./tmp_bl)); %transform to dB // dB = 10*log10 (signal/baseline)
    
    %store for all subs
    dat_wvlt_all.data_dB(s,:,:)     = tmp_dat_dB;

    % plot data
    imagesc(dat_wvlt_all.times,dat_wvlt_all.freqs,tmp_dat_dB,[-10 10])
    axis xy
    xlim ([-0.55 2.5])
    ylim ([6 35])
    vline(0,{'k:', 'LineWidth', 0.5}); % epoch start
    title ([dat_wvlt_all.id(s)])

end

% reate legend
subplot(psp(1),psp(2),[s+2:psp(1)*psp(2)])
caxis ([-10 10])
cb = colorbar;
cb.Limits = [-10 10];
ylabel(cb,'dB')
xlabel 'time [s]'
ylabel 'frequency [Hz]' 
xlim ([-0.55 2.5])
ylim ([6 35])
v = vline(0,{'k:', 'LineWidth', 0.5});
legend ([v],'Stimulus')
legend boxoff

% overall title
sgtitle 'CPz'

% save plot
save_fig(gcf,PATHOUT_plots,'single_sub_wvlt_CPz')


%% group for CPz all conditions

% extract RT data
vec_RTs = [];
for s = 1:numel(RT_ALL)
    
    vec_RTs = [vec_RTs RT_ALL(s).SO_ms(11:end)];

end

%%
close all
figure
% plot wvlt data
subplot(3,1,[1:2])
imagesc(dat_wvlt_all.times,dat_wvlt_all.freqs,squeeze(mean(dat_wvlt_all.data_dB,1)),[-8 8])
axis xy
xlim ([-0.55 2.5])
ylim ([6 35])
vline(0,{'k:', 'LineWidth', 2}); % epoch start

r = patch([-.5 -.5 -.2 -.2],[6 35 35 6],[.4 .4 .4]);
r.FaceAlpha = .5;
r.EdgeAlpha = 0;

% cb = colorbar;
% ylabel(cb,'dB')
% xlabel 'time [s]'
ylabel 'frequency [Hz]' 
l = legend ('BL')
l.Box = 'off'
% overall title
title 'Grand average all trials all participants [CPz]'


subplot(3,1,3)
%prepare single data
jitterAmount = 0.3;
vals_jitter = 2*(rand(size(vec_RTs))-0.5)*jitterAmount;

plot(vec_RTs/1000,vals_jitter+1,'.','Color',[.8 .8 .8],'LineWidth',0.7,'MarkerSize',5)
hold on
boxplot(vec_RTs/1000,1,'Orientation','horizontal','OutlierSize',0.00001,'Symbol','k.','Widths',0.25)

tune_BP(c_rh)
ax = gca;
ax.YAxis.Visible = 'off'; % remove y-axis
hold on
box off
xlim ([-0.55 2.5])
xlabel 'time [s]'
title 'Reaction times'


save_fig(gcf,PATHOUT_plots,'GrAvg_wvlt_CPz')












































