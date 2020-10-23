function  legend_ERP_valid

    xlabel 'time [ms]'
    ylabel 'amplitude [\muV]'
    xlim ([-300 1200])
    ylim ([-10 15])
    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    
end
