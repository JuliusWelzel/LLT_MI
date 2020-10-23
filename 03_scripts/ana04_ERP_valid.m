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

PATHOUT_plots = [MAIN '04_plots\04_ERPs\validity\'];

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

idx_time_P2 = dsearchn(ep_time',[180]'):dsearchn(ep_time',[280]');
idx_time_P3 = dsearchn(ep_time',[300]'):dsearchn(ep_time',[600]');

idx_Pz = strcmp({chanlocs.labels},'Pz');

%add relevant indices
for s = 1:length(ERP_all)
    
    clear ep_stim
    ep_stim = splitLLTstim(ERP_all(s).pics)';
    
    % find indices
    % idx lat/med
    ERP_all(s).idx_lat      = contains(ep_stim(:,1),'l')' & contains(ep_stim(:,4),{'240','300'})' |...
                            contains(ep_stim(:,1),'r')' & contains(ep_stim(:,4),{'60','120'})' ;
    ERP_all(s).idx_med      = contains(ep_stim(:,1),'l')' & contains(ep_stim(:,4),{'60','120'})' |...
                            contains(ep_stim(:,1),'r')' & contains(ep_stim(:,4),{'240','300'})' ;
                        
    ERP_all(s).nms_con(ERP_all(s).idx_lat)     = string('lateral');
    ERP_all(s).nms_con(ERP_all(s).idx_med)     = string('medial');
    ERP_all(s).nms_con(ERP_all(s).idx_art(idx_Pz,:)) = [];

                        
    % idx hand/foot
    ERP_all(s).idx_hand      = contains(ep_stim(:,1),'h');
    ERP_all(s).idx_foot      = contains(ep_stim(:,1),'f');

    % idx in depth rotation
    ERP_all(s).idx_id_0      = contains(ep_stim(:,3),'0'); 
    ERP_all(s).idx_id_60     = contains(ep_stim(:,3),'60'); 
    ERP_all(s).idx_id_300    = contains(ep_stim(:,3),'300'); 
    
    % group assignment, only consider artifact free epochs
    ERP_all(s).group = ones(1,sum(~ERP_all(s).idx_art(idx_Pz,:)))*idx_group(s);

    dat_ERP_peak(s).P2 = squeeze(max(ERP_all(s).mERP(idx_Pz,idx_time_P2,~ERP_all(s).idx_art(idx_Pz,:))))';
    dat_ERP_peak(s).P3 = squeeze(max(ERP_all(s).mERP(idx_Pz,idx_time_P3,~ERP_all(s).idx_art(idx_Pz,:))))';
    
    dat_ERP_mean(s).P2 = squeeze(mean(ERP_all(s).mERP(idx_Pz,idx_time_P2,~ERP_all(s).idx_art(idx_Pz,:)),2))';
    dat_ERP_mean(s).P3 = squeeze(mean(ERP_all(s).mERP(idx_Pz,idx_time_P3,~ERP_all(s).idx_art(idx_Pz,:)),2))';

         
end

%% Extract ERP info

%%% P2 %%%
close all
figure
scatterhist([dat_ERP_mean.P2],[dat_ERP_peak.P2],...
    'Group',[ERP_all.group],'Kernel','on','Color',[color.c_stroke;color.c_old;color.c_young],...
    'LineWidth',[2,2,2],'MarkerSize',[1 1 1]);
    
xlabel 'Mean Amplitude [\muV]'
ylabel 'Peak Amplitude [\muV]'
box off
sgtitle 'P2 [180-280ms] at Pz'

hLegend = findobj(gcf, 'Type', 'Legend');
%get text
hLegend.String = nms_group;
hLegend.Box = 'off'

save_fig(gcf,PATHOUT_plots,'cmpr_P2_Pz')

%%% P3 %%%
close all
figure
scatterhist([dat_ERP_mean.P3],[dat_ERP_peak.P3],...
    'Group',[ERP_all.group],'Kernel','on','Color',[color.c_stroke;color.c_old;color.c_young],...
    'LineWidth',[2,2,2],'MarkerSize',[1 1 1]);
    
xlabel 'Mean Amplitude [\muV]'
ylabel 'Peak Amplitude [\muV]'
box off
sgtitle 'P3 [300-600 ms] at Pz'

hLegend = findobj(gcf, 'Type', 'Legend');
%get text
hLegend.String = nms_group;
hLegend.Box = 'off'

save_fig(gcf,PATHOUT_plots,'cmpr_P3_Pz')


%% Boxplot P2 & P3

close all
figure
subplot(1,2,1)
singleBoxplot({[dat_ERP_peak(idx_old).P2],[dat_ERP_peak(idx_young).P2],[dat_ERP_peak(idx_stroke).P2]})
tune_BP([color.c_old;color.c_young;color.c_stroke])
    ylabel 'Amplitude [\muV]'
    xticks ([1 2 3])
    xticklabels ({'Old','Young','Stroke'})
    title 'Peak amplitude'
    box off
    
subplot(1,2,2)
singleBoxplot({[dat_ERP_mean(idx_old).P2],[dat_ERP_mean(idx_young).P2],[dat_ERP_mean(idx_stroke).P2]})
tune_BP([color.c_old;color.c_young;color.c_stroke])
    ylabel 'Amplitude [\muV]'
    xticks ([1 2 3])
    xticklabels ({'Old','Young','Stroke'})
    title 'Mean amplitude'
    box off

sgtitle 'P2 [180-280ms] at Pz'
save_fig(gcf,PATHOUT_plots,'cmpr_grp_P2_Pz')

close all
figure
subplot(1,2,1)
singleBoxplot({[dat_ERP_peak(idx_old).P3],[dat_ERP_peak(idx_young).P3],[dat_ERP_peak(idx_stroke).P3]})
tune_BP([color.c_old;color.c_young;color.c_stroke])
    ylabel 'Amplitude [\muV]'
    xticks ([1 2 3])
    xticklabels ({'Old','Young','Stroke'})
    title 'Peak amplitude'
    box off
    
subplot(1,2,2)
singleBoxplot({[dat_ERP_mean(idx_old).P3],[dat_ERP_mean(idx_young).P3],[dat_ERP_mean(idx_stroke).P3]})
tune_BP([color.c_old;color.c_young;color.c_stroke])
    ylabel 'Amplitude [\muV]'
    xticks ([1 2 3])
    xticklabels ({'Old','Young','Stroke'})
    title 'Mean amplitude'
    box off

sgtitle 'P3 [300-600 ms] at Pz'
save_fig(gcf,PATHOUT_plots,'cmpr_grp_P3_Pz')

%% ANOVA peak & mean (group [3] x con [2])
idx_all_stroke = [ERP_all.group] == 0;
idx_all_old = [ERP_all.group] == 1;
idx_all_young = [ERP_all.group] == 2;

nms_all_grp(idx_all_stroke)     = string('stroke');
nms_all_grp(idx_all_old)        = string('old');
nms_all_grp(idx_all_young)      = string('young');

nms_all_con = [ERP_all.nms_con];


% P2
dat_mean_P2 = [dat_ERP_mean.P2];
[~,~,stats] = anovan(dat_mean_P2,{nms_all_grp,nms_all_con},'model','interaction','varnames',{'Group','Condition'})
figure
multcompare(stats,'Dimension',[1 2])
    xlabel 'Amplitude [\muV]'

% P3
dat_mean_P3 = [dat_ERP_mean.P3];
[~,~,stats] = anovan(dat_mean_P3,{nms_all_grp,nms_all_con},'model','interaction','varnames',{'Group','Condition'})
figure
multcompare(stats,'Dimension',[1 2])













































