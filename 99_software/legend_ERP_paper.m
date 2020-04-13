function  legend_ERP_paper

    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-10 15])
    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    vline([450 600],'--k');
    
end

