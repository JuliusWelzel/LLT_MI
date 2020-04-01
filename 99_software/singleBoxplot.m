function singleBoxplot(dat_tp)
%singleBoxplot - create boxplot (in groups even) with single subjects data
%   
%   Input: dat_tp <- data per group each colum
%           

var_g = {};

for g = 1:size(dat_tp,2)
    %prepare single subjects data
    jitterAmount = 0.2;
    vals_jitter = 2*(rand(size(dat_tp{:,g}))-0.5)*jitterAmount;

    plot(vals_jitter+g,dat_tp{:,g},'.','Color',[.8 .8 .8],'LineWidth',0.7,'MarkerSize',10)
    hold on
    
    var_g(g) = {ones(size(dat_tp{:,g}))*g};

end

boxplot([dat_tp{:}],[var_g{:}],'OutlierSize',0.00001,'Symbol','k.','Widths',0.25)




end

