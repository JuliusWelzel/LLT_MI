function erp_leg(nms_legend)


legend(nms_legend)
ylim ([-12 12])
xlim ([-200 1200])
xlabel 'time [ms]'
ylabel 'ERP [\muV]'
vline(0,'--k')


end

