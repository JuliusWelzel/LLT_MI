%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_ERP = [MAIN '02_data\04_ERPs\'];
PATHOUT_ERPplot = [MAIN '02_data\04_ERPs\ERP_plots\'];
PATHOUT_ERP_ICA = [MAIN '02_data\04_ERPs_ICA\'];


%Check if PATHOUT-folder is already there; if not: create

if ~isdir(PATHOUT_ERP_ICA)
    mkdir(PATHOUT_ERP_ICA);
end

% IC cons
iclab_nms = {'finICA'};

%idx channel
idx_chan = 9;


% document potential problems with a subject
docError = {};   
load ([PATHIN_ERP 'cfg.mat'])

%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_ERP '*_ep_filt.set']));
SUBJ = str2double(extractBetween({list.name},'SUBJ','_'));

old = SUBJ<70 & SUBJ>=30;
young = SUBJ>= 70;
stroke = SUBJ <=30;

g_nms = {'stroke','old','young'};

for i = 1:length(g_nms)
groups(i).names = g_nms{i};
end

groups(1).idx = stroke;
groups(2).idx = old;
groups(3).idx = young;

%%
for i = 1:length(iclab_nms)
    load ([PATHIN_ERP 'ERPall_' iclab_nms{i} '.mat']);
for g = 1:length(groups)
    ss = 1;
    avgERP(g).group = groups(g).names;
    
for sub = find(groups(g).idx)
    
    % correction
    if length(ERP_all(sub).pics) <50;continue;end;% correct if more pics are empty
    
    if length(ERP_all(sub).pics) > size(ERP_all(sub).mERP,3) % correct if more pics than epochs
        ERP_all(sub).pics(size(ERP_all(sub).mERP,3):end) = [];
    elseif length(ERP_all(sub).pics) < size(ERP_all(sub).mERP,3) % correct if more pics than epochs
        ERP_all(sub).mERP(size(ERP_all(sub).pics,3):end) = [];
    end
    
    % correct artifact epochs
    art_ep = squeeze(max(ERP_all(sub).mERP(idx_chan,:,:))) > 100 | squeeze(min(ERP_all(sub).mERP(idx_chan,:,:))) < -100;
    ERP_all(sub).mERP(:,:,art_ep) = [];
    ERP_all(sub).pics(art_ep) = [];
    
    ERP_all(sub).art_ep(i) = sum(art_ep);
    

    % stimuli
    for st = 1:length(ERP_all(sub).pics)
        stim_pics = ERP_all(sub).pics;
        stim_split = strsplit(stim_pics{st},{'_','.'});
        stim(st,:) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=depth rotation, 4 = lateral rotation
    end

    %handiness
    idx_LH = contains(stim(:,1),'lh');
    idx_RH = contains(stim(:,1),'rh');
    avgERP(g).LH(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_LH),3);
    avgERP(g).RH(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_RH),3);

    % angle
    idx_0 = contains(stim(:,4),{'0'});
    idx_60 = contains(stim(:,4),{'60','300'});
    idx_120 = contains(stim(:,4),{'120','240'});
    idx_180 = contains(stim(:,4),{'180'});
    avgERP(g).a0(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_0),3);
    avgERP(g).a60(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_60),3);
    avgERP(g).a120(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_120),3);
    avgERP(g).a180(:,:,ss) = mean(ERP_all(sub).mERP(:,:,idx_180),3);

    
    
    ss = ss+1;
    clear stim
    clear art_ep
end

end % groups
ICA(i).avgERP = avgERP;
clear avgERP

end % ICA cons

%% compare ERPs
tvec = cfg.ERP.ep_time(1):1000/(1000*cfg.ERP.resam):cfg.ERP.ep_time(2)-1000/(1000*cfg.ERP.resam);
idx_chan = [6,3,7,9,8,10,21,20,22];
   
for g = 1:length(groups)
    figure 
    for c = 1:length(idx_chan)
        subplot(3,3,c)
        mERP_rh = double(mean(ICA(1).avgERP(g).RH(idx_chan(c),:,:),3));
        std_mERP_rh = double(std(ICA(1).avgERP(g).RH(idx_chan(c),:,:),[],3));
        
        mERP_lh = double(mean(ICA(1).avgERP(g).LH(idx_chan(c),:,:),3));
        std_mERP_lh = double(std(ICA(1).avgERP(g).LH(idx_chan(c),:,:),[],3));
        
        bl = boundedline(tvec, mERP_rh, std_mERP_rh,...
                    tvec, mERP_lh, std_mERP_lh,...
            'alpha','cmap',[c_rh;c_lh]); % alpha makes bounds transparent
        
        title (cfg.EEG.chanlocs(idx_chan(c)).labels)
        legend ([bl(1) bl(2)],{'RH','LH'})
        xlabel 'time [s]'
        ylabel 'amplitude [\muV]'
        xlim ([-0.3 1.0])
        ylim ([-4 4])
        vline(0,'k');
        vline([0.2 0.3],'--k');
        % Add all our previous improvements:
        ax = gca();
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        box off;
    end
    sgtitle ([groups(g).names])
    save_fig(gcf,PATHOUT_ERPplot,['ERP_overview_' groups(g).names]);
end
        






