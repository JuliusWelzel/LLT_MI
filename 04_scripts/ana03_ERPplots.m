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

PATHOUT_plots = [MAIN '03_plots\ana03_ERPs\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_plots)
    mkdir(PATHOUT_plots);
end

%% Gather all available data with clean EEG & RTs
% load RTs 

list = dir(fullfile([PATHIN_ERP '*_ep_filt.set']));
SUBJ = str2double(extractBetween({list.name},'SUBJ','_ep'));

% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;

% groups vars for plotting
g_nms   = {'stroke','old','young'}; 

% load all ERP data
load ([PATHIN_ERP 'ERPall.mat']);
load ([PATHIN_ERP 'cfg.mat']);
chanlocs = readlocs([MAIN '99_software\mobile24_proper_labels.elp']);
ep_time = linspace(cfg.ERP.ep_time(1),cfg.ERP.ep_time(2),diff(cfg.ERP.ep_time)*cfg.ERP.resam)*1000;

%% replicate behavioural analysis
idx_Pz = strcmp({chanlocs.labels},'Pz');

for s = 1:length(ERP_all)
    
    clear ep_stim
    idx_best = ~ERP_all(s).idx_art(idx_Pz,:)' & ERP_all(s).idx_cor';
    ep_stim = splitLLTstim(ERP_all(s).pics)';
    
    % find indices
    % idx lat/med
    idx_lat      = contains(ep_stim(:,1),'lh') & contains(ep_stim(:,4),{'240','300'}) & idx_best |...
                   contains(ep_stim(:,1),'rh') & contains(ep_stim(:,4),{'60','120'}) & idx_best;
    idx_med      = contains(ep_stim(:,1),'lh') & contains(ep_stim(:,4),{'60','120'}) & idx_best |...
                   contains(ep_stim(:,1),'rh') & contains(ep_stim(:,4),{'240','300'}) & idx_best;

    % idx hand/foot
    idx_hand      = contains(ep_stim(:,1),'h') & idx_best;
    idx_foot      = contains(ep_stim(:,1),'f') & idx_best;

    % idx in depth rotation
    idx_id_0      = contains(ep_stim(:,3),'0') & idx_best; 
    idx_id_60     = contains(ep_stim(:,3),'60') & idx_best; 
    idx_id_300    = contains(ep_stim(:,3),'300') & idx_best; 
    
    % store single subject ERPs
    dat_lat(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_lat),3));
    dat_med(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_med),3));
    
    dat_hand(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_hand),3));
    dat_foot(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_foot),3));

    dat_id_0(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_id_0),3));
    dat_id_60(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_id_60),3));
    dat_id_300(s,:) = squeeze(mean(ERP_all(s).mERP(idx_Pz,:,idx_id_300),3));


end %SUBJ loop


%% HAND ERP all groups

close all
figure
%%%%%%%%%%%%%%%%%%%%%% plot lateral/medial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,1)
plot(ep_time,mean(dat_lat(idx_old,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_med(idx_old,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'OLD'

subplot(3,3,4)
plot(ep_time,mean(dat_lat(idx_young,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_med(idx_young,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'YOUNG'

subplot(3,3,7)
plot(ep_time,mean(dat_lat(idx_stroke,:)),'Color',c_lat)
hold on
plot(ep_time,mean(dat_med(idx_stroke,:)),'Color',c_med)
    legend ({'lateral','medial'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'STROKE'

    
%%%%%%%%%%%%%%%%%%%%%% plot hand_feet  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,2)
plot(ep_time,mean(dat_hand(idx_old,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_foot(idx_old,:)),'Color',c_rh)
    legend ({'hand','foot'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'OLD'

subplot(3,3,5)
plot(ep_time,mean(dat_hand(idx_young,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_foot(idx_young,:)),'Color',c_rh)
    legend ({'hand','foot'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'YOUNG'

subplot(3,3,8)
plot(ep_time,mean(dat_hand(idx_stroke,:)),'Color',c_lh)
hold on
plot(ep_time,mean(dat_foot(idx_stroke,:)),'Color',c_rh)
    legend ({'hand','foot'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])

    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'STROKE'   
    
    
    
%%%%%%%%%%%%%%%%%%%%%% plot in depth rotation  %%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,3)
plot(ep_time,mean(dat_id_0(idx_old,:)),'Color',c_young)
hold on
plot(ep_time,mean(dat_id_60(idx_old,:)),'Color',c_old)
hold on
plot(ep_time,mean(dat_id_300(idx_old,:)),'Color',c_stroke)
    legend ({'0°','60°','300°'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])


    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'OLD'

subplot(3,3,6)
plot(ep_time,mean(dat_id_0(idx_young,:)),'Color',c_young)
hold on
plot(ep_time,mean(dat_id_60(idx_young,:)),'Color',c_old)
hold on
plot(ep_time,mean(dat_id_300(idx_young,:)),'Color',c_stroke)
    legend ({'0°','60°','300°'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])

    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'YOUNG'

subplot(3,3,9)
plot(ep_time,mean(dat_id_0(idx_stroke,:)),'Color',c_young)
hold on
plot(ep_time,mean(dat_id_60(idx_stroke,:)),'Color',c_old)
hold on
plot(ep_time,mean(dat_id_300(idx_stroke,:)),'Color',c_stroke)
    legend ({'0°','60°','300°'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-20 20])

    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    title 'STROKE' 
    
    
    
    
save_fig(gcf,PATHOUT_plots,['ERP_overview'],'FigSize',[0 0 30 20]);































