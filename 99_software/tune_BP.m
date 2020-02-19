function tune_BP(color)


h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color,'FaceAlpha',.5);
 end

h = findobj(gca,'Tag','Median');
set(h,'Color',color)
set(h,'LineWidth',2)


end
