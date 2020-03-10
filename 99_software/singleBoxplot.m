function singleBoxplot(dat_tp)
%singleBoxplot - create boxplot (in groups even) with single subjects data
%   
%   Input: dat_tp <- data per group each colum
%           

if iscell(dat_tp)
    
        
    jitterAmount = 0.05;

    
    for g = 1:size(dat_tp,2)
        tmp_dat_bp(g).data = dat_tp{:,g};
        tmp_dat_bp(g).idx = g*ones(1,length(dat_tp{:,g}));
        tmp_dat_bp(g).vals_jitter = 2*(rand(size(dat_tp{:,g}))-0.5)*jitterAmount;
    end
    
    dat_tp = [tmp_dat_bp.data];
    idx_group = [tmp_dat_bp.idx];
    vals_jitter = [tmp_dat_bp.vals_jitter];
    
    
    plot(vals_jitter+idx_group,dat_tp,'.','Color',[.8 .8 .8],'LineWidth',0.7,'MarkerSize',10)
    hold on
    boxplot(dat_tp,idx_group,'OutlierSize',0.00001,'Symbol','k.','Widths',0.25)

else
    for g = 1:size(dat_tp,2)
        %prepare single subjects data
        jitterAmount = 0.2;
        vals_jitter = 2*(rand(size(dat_tp))-0.5)*jitterAmount;

        plot(vals_jitter(:,g)+g,dat_tp,'.','Color',[.8 .8 .8],'LineWidth',0.7,'MarkerSize',10)
        hold on
    end

    G = [ones(size(A1))  2*ones(size(A2)) 3*ones(size(A2)) 4*ones(size(A2))];
        X = [A1, A2, A3, A4];

        boxplot(X,G)

    boxplot(dat_tp,'OutlierSize',0.00001,'Symbol','k.','Widths',0.25)

end

end

