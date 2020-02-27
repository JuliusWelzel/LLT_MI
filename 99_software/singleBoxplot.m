function singleBoxplot(dat_tp)
%singleBoxplot - create boxplot (in groups even) with single subjects data
%   
%   Input: dat_tp <- data per group each colum
%           


%prepare single subjects data
jitterAmount = 0.2;
vals_jitter = 2*(rand(size(dat_tp))-0.5)*jitterAmount;

for g = 1:size(dat_tp,2)
    plot(vals_jitter(:,g)+g,dat_tp,'.','Color',[.8 .8 .8],'LineWidth',0.7,'MarkerSize',10)
    hold on
end
boxplot(dat_tp,'OutlierSize',0.00001,'Symbol','k.','Widths',0.25)



end

