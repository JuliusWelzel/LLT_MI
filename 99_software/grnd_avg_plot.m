function grnd_avg_plot(ERP_all,channel_oi,cfg)
%grnd_avg_plot Plot grand average ERPs in time and space for LLT-MI data
%INPUT:
%     - ERP_all -> MatLab struct with all ERPs
%     - channel -> channel to average over and highliht in topo
%     - time_wins -> time windows for topoplots
%     - cfg -> config file which contains all info
% 
% Author: Julius Welzel (j.welzel@neurologi.uni-kiel.de)

    idx_chan = ismember({cfg.ERP.chanlocs.labels},channel_oi);
    t_wins_oi   = [180, 260; 300, 600];

    idx_time_P2 = dsearchn(cfg.ep_time',[t_wins_oi(1,1)]'):dsearchn(cfg.ep_time',[t_wins_oi(1,2)]');
    idx_time_P3 = dsearchn(cfg.ep_time',[t_wins_oi(2,1)]'):dsearchn(cfg.ep_time',[t_wins_oi(2,2)]');
    
    
    
    for s = 1:numel(ERP_all)
        
        erp_time(s,:)   = squeeze(mean(mean(ERP_all(s).mERP(idx_chan,:,~any(ERP_all(s).idx_art(idx_chan,:),1)),3),1));
        erp_map_P2(s,:) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_time_P2,~any(ERP_all(s).idx_art(idx_chan,:),1)),3),2));
        erp_map_P3(s,:) = squeeze(mean(mean(ERP_all(s).mERP(:,idx_time_P3,~any(ERP_all(s).idx_art(idx_chan,:),1)),3),2));
        
        idx_group(s) = ERP_all(s).group(1);
        
    end
    
    %% Plotting per group

    figure
    subplot(2,2,[1 2])
    plot(cfg.ep_time,erp_time,'Color',[0.7 0.7 0.7],'LineWidth',0.003)
    hold on
    pm = plot(cfg.ep_time,mean(erp_time),'k','LineWidth',2)

    legend([pm],'Grand average')
    legend  boxoff
    legend_ERP_valid
    
    vline(t_wins_oi(1,:),{'Color',cfg.color.c_med,'LineStyle',':','LineWidth',2})
    vline(t_wins_oi(2,:),{'Color',cfg.color.c_lat,'LineStyle',':','LineWidth',2})

    subplot(2,2,3)
    topoplot(mean(erp_map_P2),cfg.ERP.chanlocs,...
    'maplimits',[-2 5],'emarker2',{find(idx_chan),'.','w',15})
        set(findobj(gca,'type','patch'),'facecolor',cfg.color.c_med)

    subplot(2,2,4)
    topoplot(mean(erp_map_P3),cfg.ERP.chanlocs,...
    'maplimits',[-5 10],'emarker2',{find(idx_chan),'.','w',15})
        set(findobj(gca,'type','patch'),'facecolor',cfg.color.c_lat)

end
