function legend_ERD_terhorst2013

    legend ({'medial','lateral'})
    legend('boxoff')
    xlabel 'time [ms]'
    ylabel 'Power [dB]'
    xlim ([-500 1000])
    ylim ([-6 3])
    vline(0,'k');
    vline([300 800],'--k');
    % Add all our previous improvements:
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off;
    
    
end

