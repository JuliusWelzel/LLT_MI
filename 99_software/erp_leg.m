function erp_leg(nms_legend)


legend(nms_legend)
legend ('boxoff')
ylim ([-5 15])
xlim ([-200 1200])
xlabel 'time [ms]'
ylabel 'ERP [\muV]'
vline(0,'--k')

% center axis around 0
ax = gca();
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;


end

